
dbNSFP_cleanup <- function(genie_maf_annotated){
  genie_maf_annotated <- setnames(genie_maf_annotated, "CADD_phred_hg19", "CADD_score")
  
  score_cols <- grep("_score", colnames(genie_maf_annotated), value=T)
  min_cols <-  c("SIFT_score", "FATHMM_score", "ESM1b_score")
  max_cols <- setdiff(score_cols, min_cols)
  pred_cols <- grep("_pred", colnames(genie_maf_annotated), value=T)
  
  genie_maf_annotated <- genie_maf_annotated[, (score_cols) := lapply(.SD, as.numeric), .SDcols = score_cols]
  
  genie_maf_annotated <- genie_maf_annotated[, (max_cols) := lapply(.SD, function(x) max(x, na.rm=T)),
                                             by = .(HGVSp_Short, Hugo_Symbol), .SDcols = max_cols]
  
  genie_maf_annotated <- genie_maf_annotated[, (min_cols) := lapply(.SD, function(x) min(x, na.rm=T)), 
                                             by = .(HGVSp_Short, Hugo_Symbol), .SDcols = min_cols]
  
  genie_maf_annotated <- genie_maf_annotated[, (score_cols) := lapply(.SD, function(x) ifelse(is.infinite(x), NA, x)),
                                             .SDcols = score_cols]
  
  genie_maf_annotated <- genie_maf_annotated[, (pred_cols) := lapply(.SD, function(x) ifelse(x==".", NA, x)), 
                                             .SDcols = pred_cols]
  
  genie_maf_annotated <- genie_maf_annotated[, (pred_cols) := lapply(.SD, function(x) ifelse(all(is.na(x)), NA, na.omit(x)[1])), 
                                             by = .(HGVSp_Short, Hugo_Symbol), .SDcols = pred_cols]
  
  genie_maf_annotated <- unique(genie_maf_annotated, by=c("Chromosome","Start_Position","Reference_Allele", "Tumor_Seq_Allele2","HGVSc_ANNOVAR", "HGVSp_Short", "Tumor_Sample_Barcode"))
  
  genie_maf_annotated <- genie_maf_annotated[, g := .GRP, by=c("Chromosome","Start_Position","Reference_Allele", "Tumor_Seq_Allele2","HGVSc_ANNOVAR", "HGVSp_Short")]
  genie_maf_annotated <- genie_maf_annotated[, g_count := .N, by="g"]
  
  return(genie_maf_annotated)
  
}


# Function to process dbNSFP label and standardize to a binary nomenclature system

pred_to_label <- function(maf){
  pred_cols <- grep("pred", colnames(maf), value = T)
  for (i in 1:length(pred_cols)){
    method <- sub("_pred", "_pathogenic", pred_cols[i])
    method_numeric <- sub("_pred", "_numeric", pred_cols[i])
    maf <- maf[, eval(method) := case.(get(pred_cols[i]) %in% c("T", "B", "L", "N"), "Non_pathogenic",
                                       get(pred_cols[i]) %in% c("D", "P", "H", "M"), "Pathogenic",
                                       default = "Unknown")]
    maf <- maf[, eval(method_numeric) := case.(get(pred_cols[i]) %in% c("T", "B", "L", "N"), 0,
                                               get(pred_cols[i]) %in% c("D", "P", "H", "M"), 1,
                                               default = -1)]
  }
  return(maf)
}

determine_cutpoint <- function(genie_maf_annotated, method_col, direction){
  
  genie_maf_annotated <- setnames(genie_maf_annotated, "ONCOGENIC", "oncogenic", skip_absent = T)
  
  cadd_subset <- unique(genie_maf_annotated[oncogenic %in% c("Oncogenic", "Likely Oncogenic", "Resistance"), ], by=c("HGVSp_Short", "Hugo_Symbol"))
  
  cadd_subset <- cadd_subset[, oncogenic_binary := 1][, colnames(cadd_subset) %in% c("oncogenic_binary", "HGVSp_Short", "Hugo_Symbol", method_col), with=F]
  
  cadd_subset <- cadd_subset
  
  cadd_subset_binary <- rbind(cadd_subset[!is.na(get(method_col)),], 
                              common_scores[, colnames(common_scores) == method_col, with=F][, oncogenic_binary := 0], fill=T)
  
  cadd_subset_binary <- cadd_subset_binary[, oncogenic_binary := as.character(oncogenic_binary)]
  cadd_subset_binary <- cadd_subset_binary[, eval(method_col) := as.numeric(get(method_col))]
  
  cadd_subset_binary <- cadd_subset_binary[!is.na(get(method_col)),]
  
  cadd_cutpoint <- cutpointr(
                             x=as.numeric(cadd_subset_binary[[method_col]]),
                             class=cadd_subset_binary$oncogenic_binary,
                             pos_class=1,
                             neg_class=0,
                             direction = direction,
                             metric = sum_sens_spec,
                             break_ties=median)
  
  pred_col <- sub("_score", "_pred", method_col)
  numeric_col <- sub("_score", "_numeric", method_col)
  
  if (direction == ">="){
    genie_maf_annotated <- genie_maf_annotated[, eval(pred_col) := ifelse(get(method_col) >= cadd_cutpoint$optimal_cutpoint, "D", "T")]
    genie_maf_annotated <- genie_maf_annotated[, eval(numeric_col) := ifelse(get(method_col) >= cadd_cutpoint$optimal_cutpoint, 1, 0)]
  } else {
    genie_maf_annotated <- genie_maf_annotated[, eval(pred_col) := ifelse(get(method_col) <= cadd_cutpoint$optimal_cutpoint, "D", "T")]
    genie_maf_annotated <- genie_maf_annotated[, eval(numeric_col) := ifelse(get(method_col) <= cadd_cutpoint$optimal_cutpoint, 1, 0)]
  }
 
  
  g1 <- plot_metric(cadd_cutpoint) +
    theme_classic()+
    labs(y="Sum of sensitivity and specificity", x="Prediction score", title=sub("_score", "", method_col), subtitle="")+
    geom_vline(xintercept = cadd_cutpoint$optimal_cutpoint, color="blue", linetype="dashed")+
    annotate("text", x = cadd_cutpoint$optimal_cutpoint - cadd_cutpoint$optimal_cutpoint*0.05, y = 1.03, 
             label = paste0("Optimal cutoff:\n",round(cadd_cutpoint$optimal_cutpoint, 2)), 
             color="blue", hjust=1, size=5)
  
  return(list(plot=g1, annotated_df=genie_maf_annotated))
  
}

# Function to label MAF file and standardize OncoKB/VEP annotation columns into 3 main columns:
# a "_score" column with raw prediction scores outputted by each method
# a "_pathogenic" column where a variant is labeled non-pathogenic, pathogenic or unknown, resulting from applying pre-determined cutoff on the prediction scores  
# a "_numeric" column where a variant is labeled 0 (non-pathogenic), 1 (pathogenic) or -1 (unknown) 
maf_process <- function(maf){
  maf <- maf[HGVSp_Short != "", ][, colnames(maf) %in% c("Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_Position", "HGVSc", 
                                                         "HGVSp_Short", "SIFT", "PolyPhen", "Hotspot", "indel_hotspot", 
                                                         "indel_hotspot_type", "oncogenic", "Highest_level", "PHRED", 
                                                         grep("numeric|pathogenic|score", colnames(maf), ignore.case=T, value=T)), with=F]
  
  
  maf <- maf[, PATIENT_ID := substr(Tumor_Sample_Barcode, 1, 9)]
  
  maf <- maf[, oncogenic_numeric := case.(oncogenic %in% c("Likely Oncogenic", "Oncogenic", "Resistance"), 1,
                                          oncogenic%in% c("Neutral", "Likely Neutral"), 0,
                                          default = -1)]
  
  label_cols <- grep("label|pathogenic", colnames(maf), value = T)
  numeric_cols <- grep("numeric", colnames(maf), value=T)
  score_cols <- grep("score", colnames(maf), value=T)
  
  maf <- maf[, (label_cols) := lapply(.SD, function(x) ifelse(is.na(x), "Unknown", x)), .SDcols = label_cols]
  maf <- maf[, (numeric_cols) := lapply(.SD, function(x) case.(x %in% c(2,10, 0), 0,
                                                               is.na(x), -1, 
                                                               x == 1, 1,
                                                               default = -1)), .SDcols = numeric_cols]
  
  maf <- maf[, (score_cols) := lapply(.SD, function(x) as.numeric(x)), .SDcols = score_cols]
  
  if (all(c("SIFT", "PolyPhen") %in% colnames(maf))){
    maf <- maf[, c("SIFT_label", "SIFT_score") := tstrsplit(SIFT, split="\\(")]
    maf <- maf[, SIFT_score := as.numeric(sub("\\)", "", SIFT_score))]
    maf <- maf[, c("PolyPhen_label", "PolyPhen_score") := tstrsplit(PolyPhen, split="\\(")]
    maf <- maf[, PolyPhen_score := as.numeric(sub("\\)", "", PolyPhen_score))]
    
    maf <- maf[, SIFT_numeric := ifelse(SIFT_label %in% c("deleterious_low_confidence", "deleterious"), 1,0)]
    maf <- maf[, PolyPhen_numeric := ifelse(PolyPhen_label %in% c("probably_damaging", "possibly_damaging"), 1, 0)]
    
    maf <- maf[, `:=`( PolyPhen_label_consolidated = str_to_title(gsub("possibly_|probably_", "", PolyPhen_label)),
                       SIFT_label_consolidated = str_to_title(gsub("_low_confidence", "", SIFT_label)))]
    
  }
  


  maf <- maf[, oncogenic_consolidated := case.(oncogenic == "Likely Oncogenic", gsub("Likely ", "", oncogenic),
                                               oncogenic == "Likely Neutral", "Neutral",
                                               oncogenic %in% c("Inconclusive", "Unknown"), "Unknown",
                                               is.na(oncogenic), "Unknown",
                                               default = oncogenic)]
  
  maf <- maf[, oncogenic_binary := ifelse(oncogenic %in% c("Likely Oncogenic", "Oncogenic", "Resistance"), "Pathogenic",
                                          "Non_pathogenic")]
  
  return(maf)
  
}


# Function to define reclassified/"rescued" mutation
# This adds another column "_Rescue" to every VEP where VUSs (variants labeled as "Inconclusive" or not annotated by OncoKB) can be "Rescued_Oncogenic" 
# if a VEP predicts that they are pathogenic, "Rescued_Benign" if a VEP predicts that they are non-pathogenic
# Variants with OncoKB annotations are labeled with their OncoKB annotations
define_rescue <- function(maf_process){
  method_cols <- grep("pathogenic", colnames(maf_process), value=T)
  for (i in 1:length(method_cols)){
    tmp_col <- method_cols[i]
    maf_process <- maf_process[, Rescue := case.(oncogenic_consolidated == "Unknown" & get(tmp_col) == "Pathogenic", "Rescued_Oncogenic",
                                                 oncogenic_consolidated == "Unknown" & get(tmp_col) == "Non_pathogenic", "Rescued_Benign",
                                                 oncogenic_consolidated %in% c("Oncogenic", "Resistance"), "OncoKB_Oncogenic",
                                                 oncogenic_consolidated == "Neutral", "OncoKB_Benign",
                                                 default = "Unknown")]
    
    maf_process <- setnames(maf_process, "Rescue", sub("_pathogenic", "_Rescue", tmp_col))
  }
  
  if (all(c("SIFT_label_consolidated", "PolyPhen_label_consolidated") %in% colnames(maf_process))){
    maf_process <- maf_process[ , SIFT_Rescue := case.(oncogenic_consolidated  %in% c("Oncogenic", "Resistance"),"OncoKB_Oncogenic",
                                                       oncogenic_consolidated == "Neutral", "OncoKB_Benign",
                                                       SIFT_label_consolidated == "Deleterious" & oncogenic_consolidated == "Unknown",
                                                       "Rescued_Oncogenic",
                                                       SIFT_label_consolidated == "Tolerated" & oncogenic_consolidated == "Unknown",
                                                       "Rescued_Benign",
                                                       default="Unknown")]
    
    maf_process <- maf_process[ , PolyPhen_Rescue := case.(oncogenic_consolidated  %in% c("Oncogenic", "Resistance"),"OncoKB_Oncogenic",
                                                           oncogenic_consolidated == "Neutral", "OncoKB_Benign",
                                                           PolyPhen_label_consolidated == "Damaging" & oncogenic_consolidated == "Unknown",
                                                           "Rescued_Oncogenic",
                                                           PolyPhen_label_consolidated == "Benign" & oncogenic_consolidated == "Unknown",
                                                           "Rescued_Benign",
                                                           default="Unknown")]
  }
  
 
  return(maf_process)
  
  
}

# Function to merge annotated MAF file with clinical file on a per-gene basis, fill in all columns for cases without mutation
# For GENIE only: subset patients sequenced by a panel with the gene included only
prepareGenomicDF <- function(maf, clinical, gene=NULL, panel_presence=NULL){
  
  # subset MAF by gene of interest
  if (!is.null(gene)){
    maf_subset <- maf[Hugo_Symbol %in% gene, ]
  } else {
    maf_subset <- maf
  }
  
  maf_subset <- maf_subset[Tumor_Sample_Barcode %in% clinical$SAMPLE_ID, ]
  maf_subset <- maf_subset[, oncogenic_numeric := case.(oncogenic_consolidated %in% c("Likely Oncogenic", "Oncogenic", "Resistance"), 1,
                                                        oncogenic_consolidated%in% c("Neutral", "Likely Neutral"), 0,
                                                        default = -1)]

  # merge in with clinical df
  pt_subset <- merge(clinical, maf_subset,
                     by.y=c("Tumor_Sample_Barcode", "PATIENT_ID"), by.x=c("SAMPLE_ID", "PATIENT_ID"),
                     all.x=T)
  
  if (!is.null(panel_presence)){
    pt_subset <- merge(pt_subset,
                       panel_presence[Hugo_Symbol == gene & !is.na(inPanel), ][, .(SEQ_ASSAY_ID)],
                       by="SEQ_ASSAY_ID")
  }

  # fill in the blanks for cases with no mutations for all columns 
  
  text_cols <- grep("consolidated|pathogenic|Rescue|oncogenic_binary", colnames(pt_subset), value=T)
  
  numeric_cols <- grep("numeric", colnames(pt_subset), value=T)
  
  pt_subset <- pt_subset[, (text_cols) := lapply(.SD, function(x) factor(ifelse(is.na(x), "No_mutation", x))), .SDcols = text_cols]
  
  pt_subset <- pt_subset[, (numeric_cols) := lapply(.SD, function(x) ifelse(is.na(x), -1, x)), .SDcols = numeric_cols]
  
  # convert to factors
  
  pt_subset <- pt_subset[, gene_of_interest := factor(ifelse(is.na(Hugo_Symbol), 0, 1), levels=c(0,1))]

  pt_subset <- pt_subset[, (text_cols) := lapply(.SD, function(x) relevel(x, ref="No_mutation")), .SDcols = text_cols]
  
  # select only mutation with highest OncoKB annotation in case there are multiple mutations
  if (!is.null(gene)){
    pt_subset <- pt_subset[pt_subset[, .I[which.max(oncogenic_numeric)], by=.(PATIENT_ID, Hugo_Symbol)]$V1]
  }
  
  # calculate correlation between VEP vs. oncogenic prediction
  corr_matrix <- pt_subset[, grep("numeric", colnames(pt_subset), value=T), with=F]
  corr_matrix <- setnames(corr_matrix, colnames(corr_matrix), gsub("_numeric", "", colnames(corr_matrix)))
  corr_matrix <- setnames(corr_matrix, c("AM", "MA", "Clinvar", "oncogenic"),
                          c("AlphaMissense", "MutationAssessor", "ClinVar", "OncoKB"),
                          skip_absent = T)
  
  p <- corrplot(cor(corr_matrix, method="spearman"), type="upper",tl.col = "black")
  
  return_list <- list(pt_subset, p)
  names(return_list) <- c("pt_subset", "corrplot")
  
  return(return_list)
  
}


# Function to do inverse probability of treatment weighting
# Return weights to be used for downstream Cox's proportional hazard models or Kaplan-Meier curve
calculateIPTWeights <- function(lhs, rhs, df, focal_level=NULL) {
  
  # Specify matching formula
  matching_formula <- as.formula(paste0(lhs, "~", paste0(rhs, collapse=" + ")))

  # Calculate weights using WeightIt package
  weighted <- weightit(matching_formula,
                         data = df, estimand = "ATE",
                         method = "ps")
  
  # Diagnostic table for matching
  match_summary <-  bal.tab(weighted, stats = c("m", "v"), thresholds = c(m = .25))
  
  # Love plot showing absolute mean differences before and after matching
  p <- love.plot(weighted, stars = "std",
                 abs=TRUE, line=T,
                 thresholds = c(m = .2),
                 limits=list(m = c(0, .5)),
                 sample.names = c("Unweighted", "PS Weighted"),
                 colors = c("grey", "black"))
  
  return(list(p,weighted))
}

# Function to run Cox's proportional hazard model and KM curve following IPTW for a single gene/outcome
# Takes a merged clinical and genomics dataframe, a formula to be passed to coxph() or Surv(), and an object output by calculateIPTWeights()
weightedSurvival <- function(df, survform, weighted) {
  
  cat("\nKM curve")
  
  # Fit KM curve
  weighted_km <- do.call(survfit, list(formula=as.formula(survform), data=df,
                                       weights = weighted$weights))
  
  # Plot
  s1 <- do.call(ggsurvplot, list(fit=weighted_km, risk.table = T, marks=T, conf.int=T,
                                 xlabs = "Time (months)", ystrataname = "", 
                                 data=df))
  
  cat("\nModel PFS with Cox PH")
  
  # Fit weighted Cox PH
  weighted_cox=coxph(as.formula(survform),
                     data=df,weights=weighted$weights)
  
  # Forest plot
  g1 <- ggforest(weighted_cox, data=df, main="")
  
  return_list <- list(s1, g1, weighted_km, weighted_cox)
  names(return_list) <- c("km_plot_weighted", "coxph_forest_plots", "km_object", "coxph_models")
  
  return(return_list)
  
}

# Function to run the whole pipeline for OS effect for a list of genes and a list of stratified variables
# Starting with an annotated MAF and a clinical dataframe, this function merges the two on a per-gene basis
# run IPTW and Cox PH, then return results as a list

weightedSurvivalGeneList <- function(maf, clinical, gene_to_run, 
                                     outcomes, matching_vars,
                                     survobjform, 
                                     panel_presence=NULL, #for GENIE only: allowing the subset of cases with a specific gene sequenced 
                                     count_threshold, # minimum number of patients in each strata to be included
                                     baseline="other") { # in case the baseline is different from "No_mutation" 
  
  # loop through gene list
  surv_list_snv <- list()
  for (g in 1:length(gene_to_run)){
    temp_gene <- gene_to_run[g]
    cat("\n", temp_gene)
    
    # merge clinical and genomic, fill in missing values for cases with no mutation
    temp_df <- prepareGenomicDF(maf = maf, clinical = clinical, gene=temp_gene, panel_presence=panel_presence)
    
    # loop through VEP/left-hand side variables
    surv_list <- list()
    
    for (i in 1:length(outcomes)){
      tmp_outcome <- outcomes[i]
      print(tmp_outcome)
      
      # This allows the specification of other comparisons instead of the usual mutation vs. no mutation comparison
      # If run default, leave baseline blank
      if (baseline == "oncogenic"){
        rundf <- temp_df[[1]][get(tmp_outcome) %in% c("OncoKB_Oncogenic", "Rescued_Oncogenic"), ]
      } else if (baseline == "benign") {
        rundf <- temp_df[[1]][get(tmp_outcome) %in% c("No_mutation", "Rescued_Benign"), ]
      } else {
        rundf <- temp_df[[1]]
      }
      
      # Remove strata with too few cases (below count_threshold)
      outcome_level <- as.character(rundf[, .N, by=tmp_outcome][N >= count_threshold, ][[1]])
      rundf <- rundf[rundf[[tmp_outcome]] %in% outcome_level, ]
      rundf[[tmp_outcome]] <- droplevels(rundf[[tmp_outcome]])
      
      # calculate weights
      tmp_weights <- tryCatch(calculateIPTWeights(lhs = tmp_outcome,
                                                  rhs = matching_vars,
                                                  df = rundf),
                              error=function(e) NULL)
      
      # set baseline
      if (!is.null(tmp_weights)){
        if (baseline == "oncogenic"){
          rundf[[tmp_outcome]] <- relevel(rundf[[tmp_outcome]], ref="OncoKB_Oncogenic") 
        } else if (baseline == "benign") {
          rundf[[tmp_outcome]] <- relevel(rundf[[tmp_outcome]], ref="No_mutation") 
        }
        
        # run single model CoxPH
        surv_list[[i]] <- weightedSurvival(df=rundf, 
                                           survform=paste0(survobjform, "~", tmp_outcome), 
                                           weighted=tmp_weights[[2]])
      } else { surv_list[[i]] <- NULL }
    }
    surv_list_snv[[g]] <- surv_list
  }
  
  # remove genes with blank coefficients due to having too few cases before returning
  gene_to_run_subset <- gene_to_run[sapply(surv_list_snv, function(x) length(x) > 0)]
  surv_list_snv <- surv_list_snv[sapply(surv_list_snv, function(x) length(x) > 0)]
  names(surv_list_snv) <-  gene_to_run_subset
  
  surv_list_snv <- lapply(surv_list_snv, function(x) x[unlist(sapply(x, function(y) !is.null(y)))])
  
  return(surv_list_snv)
}


# Function to format the result table from weightedSurvivalGeneList into a table with nice variable names for plotting
createPlotDFs <- function(surv_list, gene_to_run){
  coxph_coeffs <- lapply(surv_list, function(x) lapply(x, function(y) as.data.frame(tidy(y$coxph_models, conf.int = TRUE))))
  count_dfs <- lapply(surv_list, function(x) lapply(x, function(y) unique(setDT(y$km_plot_weighted[["table"]][["data"]])[, .(strata, strata_size)])))
  names(coxph_coeffs) <- names(count_dfs) <- gene_to_run
  
  temp_dfs <- list()
  for (i in 1:length(coxph_coeffs)){
    if (length(coxph_coeffs[[i]]) != 0){
      tmp_coefs <- rbindlist(coxph_coeffs[[i]])[, gene := names(coxph_coeffs)[i]]
      tmp_count <- rbindlist(count_dfs[[i]])[, term := gsub("=", "", strata)]
      temp_dfs[[i]] <- merge(tmp_coefs, tmp_count[, .(term, strata_size)], by="term")
    }
  }
  
  coefs_df <- rbindlist(temp_dfs)
  coefs_df <- coefs_df[, `:=`(hr = exp(estimate),
                              hr.conf.low=exp(conf.low),
                              hr.conf.high=exp(conf.high))]
  coefs_df <- coefs_df[, pval_adj := p.adjust(p.value, method = "fdr")]
  coefs_df <- coefs_df[, estimate_sig := ifelse(pval_adj <= 0.1, hr, NA)]
  #gene_with_coefs <- unique(coefs_df[!is.na(estimate_sig),]$gene)
  
  
  coefs_df <- coefs_df[, y_lab := gsub("_label_|consolidated|_|_of_interest1", " ", term)]
  coefs_df <- coefs_df[, y_lab := sub("  ", ": ", y_lab)]
  coefs_df <- coefs_df[, y_lab := sub("oncogenic: ", "OncoKB: ", y_lab)]
  coefs_df <- coefs_df[, y_lab := sub(" pathogenic", ": ", y_lab, ignore.case = F)]
  coefs_df <- coefs_df[, y_lab := ifelse(y_lab == "gene ", "Any alteration", y_lab)]
  coefs_df <- coefs_df[, ann_type := case.(grepl("Unknown", y_lab), "Unknown",
                                           grepl("Benign|Tolerated|Non pathogenic|Neutral", y_lab), "Benign",
                                           grepl("Oncogenic|Deleterious|Pathogenic|Damaging|Resistance", y_lab), "Damaging",
                                           default="Gene")]
  
  coefs_df <- coefs_df[, ann_type := factor(ann_type,
                                            levels=c("Gene", "Damaging", "Benign", "Unknown"))]

  coefs_rescued <- coefs_df[grepl("RescueRescued_Oncogenic|RescueRescued_Benign|oncogenic_binaryPathogenic|oncogenic_binaryNon_pathogenic|gene_of_interest1|any_VUS1", term), ]
  coefs_rescued <- coefs_rescued[, ann_type := case.(grepl("OncoKB|oncogenic_binary", term), "OncoKB",
                                                                   grepl("Clinvar", term), "ClinVar",
                                                                   grepl("SP", term), "SIFT/PolyPhen",
                                                                   grepl("MA", term), "MutationAssessor",
                                                                   grepl("AM", term), "AlphaMissense",
                                                                   grepl("gene_of_interest", term), "Gene",
                                                                   default=sub("_RescueRescued_Benign|_RescueRescued_Oncogenic", "", term))]
  
  #coefs_rescued <- coefs_rescued[, ann_type := relevel(as.factor(ann_type), ref="Gene")]
  
  coefs_rescued <- coefs_rescued[, y_lab := gsub(" RescueRescued", ": Rescued", y_lab)]
  coefs_rescued <- coefs_rescued[, y_lab := gsub(".*RescueOncoKB", "OncoKB", y_lab)]
  coefs_rescued <- coefs_rescued[, y_lab := ifelse(grepl("Unknown", y_lab), "Unknown", y_lab)]

  if ("Unknown" %in% unique(coefs_rescued$y_lab)){
    coefs_rescued <- coefs_rescued[, y_lab := relevel(factor(y_lab), ref="Unknown")]
    }
  if ("OncoKB Oncogenic" %in% unique(coefs_rescued$y_lab)){
    coefs_rescued <- coefs_rescued[, y_lab := relevel(factor(y_lab), ref="OncoKB Oncogenic")]
  }
  
  returnlist <- list(coefs_df, coefs_rescued)
  names(returnlist) <- c("all", "rescued")
  
  return(returnlist)
  
}
