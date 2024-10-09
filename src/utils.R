# Constants

ref_matrix <- matrix(c(1,1,2,2), nrow=2)

# Function to process dbNSFP label and standardize to a binary nomenclature system
pred_to_label <- function(maf){
  pred_cols <- grep("pred", colnames(maf), value = T)
  for (i in 1:length(pred_cols)){
    method <- sub("_pred", "_pathogenic", pred_cols[i])
    maf <- maf[, eval(method) := case.(get(pred_cols[i]) %in% c("T", "B", "L", "N"), "Non_pathogenic",
                                       get(pred_cols[i]) %in% c("D", "P", "H", "M"), "Pathogenic",
                                       default = "Unknown")]
  }
  return(maf)
}


# Function to label MAF file and standardize OncoKB/VEP annotation columns into 3 main columns:
# a "_score" column with raw prediction scores outputted by each method
# a "_pathogenic" column where a variant is labeled non-pathogenic, pathogenic or unknown, resulting from applying pre-determined cutoff on the prediction scores  
# a "_numeric" column where a variant is labeled 0 (non-pathogenic), 1 (pathogenic) or -1 (unknown) 
maf_process <- function(maf){
  maf <- maf[HGVSp_Short != "", ][, colnames(maf) %in% c("Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_Position", "HGVSc", 
                                                         "HGVSp_Short", "Hotspot", "indel_hotspot", 
                                                         "indel_hotspot_type", "oncogenic", "Highest_level", "PHRED", 
                                                         grep("numeric|pathogenic|score", colnames(maf), ignore.case=T, value=T)), with=F]
  
  
  maf <- maf[, PATIENT_ID := substr(Tumor_Sample_Barcode, 1, 9)]
  
  maf <- maf[, oncogenic_numeric := ifelse(oncogenic %in% c("Likely Oncogenic", "Oncogenic", "Resistance"), 1, 0)]
  
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


prepareGenomicDF <- function(maf, clinical, gene=NULL, panel_presence=NULL){
  
  # subset MAF by gene of interest
  if (!is.null(gene)){
    maf_subset <- maf[Hugo_Symbol %in% gene, ]
  } else {
    maf_subset <- maf
  }
  
  maf_subset <- maf_subset[Tumor_Sample_Barcode %in% clinical$SAMPLE_ID, ]
  maf_subset <- maf_subset[, oncogenic_numeric := ifelse(oncogenic_consolidated %in% c("Likely Oncogenic", "Oncogenic", "Resistance"), 1,
                                                       0)]

  # merge in with clinical df
  pt_subset <- merge(clinical, maf_subset,
                     by.y=c("Tumor_Sample_Barcode", "PATIENT_ID"), by.x=c("SAMPLE_ID", "PATIENT_ID"),
                     all.x=T)
  
  if (!is.null(panel_presence)){
    pt_subset <- merge(pt_subset,
                       panel_presence[Hugo_Symbol == gene & !is.na(inPanel), ][, .(SEQ_ASSAY_ID)],
                       by="SEQ_ASSAY_ID")
  }

  # fill in the blanks for all columns
  
  text_cols <- grep("consolidated|pathogenic|Rescue|oncogenic_binary", colnames(pt_subset), value=T)
  numeric_cols <- grep("numeric", colnames(pt_subset), value=T)
  
  pt_subset <- pt_subset[, (text_cols) := lapply(.SD, function(x) factor(ifelse(is.na(x), "No_mutation", x))), .SDcols = text_cols]
  
  pt_subset <- pt_subset[, (numeric_cols) := lapply(.SD, function(x) ifelse(is.na(x), -1, x)), .SDcols = numeric_cols]
  
  # convert to factors
  
  pt_subset <- pt_subset[, gene_of_interest := factor(ifelse(is.na(Hugo_Symbol), 0, 1), levels=c(0,1))]

  
  pt_subset <- pt_subset[, (text_cols) := lapply(.SD, function(x) relevel(x, ref="No_mutation")), .SDcols = text_cols]
  
  # select only mutation with highest oncokb annotation in case there are multiple mutations
  if (!is.null(gene)){
    pt_subset <- pt_subset[pt_subset[, .I[which.max(oncogenic_numeric)], by=.(PATIENT_ID, Hugo_Symbol)]$V1]
  }
  
  # calculate correlation between protein vs. oncogenic prediction
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


calculateIPTWeights <- function(lhs, rhs, df, focal_level=NULL) {
  
  matching_formula <- as.formula(paste0(lhs, "~", paste0(rhs, collapse=" + ")))
  
  # pre_matched <- bal.tab(matching_formula,
  #         data = clinical_pfs, estimand = "ATT", focal="LEVEL_1", thresholds = c(m = .2))
  
  if (!is.null(focal_level)) {
    weighted <- weightit(matching_formula,
                         data = df, estimand = "ATT",
                         focal=focal_level,
                         method = "ps")
    
  } else {
    weighted <- weightit(matching_formula,
                         data = df, estimand = "ATE",
                         method = "ps")
  }
  
  match_summary <-  bal.tab(weighted, stats = c("m", "v"), thresholds = c(m = .25))
  
  p <- love.plot(weighted, stars = "std",
                 abs=TRUE, line=T,
                 thresholds = c(m = .2),
                 limits=list(m = c(0, .5)),
                 sample.names = c("Unweighted", "PS Weighted (ATT)"),
                 colors = c("grey", "black"))
  
  #non_balanced <- setDT(match_summary$Balance.Across.Pairs)[M.Threshold == "Not Balanced, >0.2", ]
  
  # if (match_summary$Balanced.mean.diffs["Not Balanced, >0.25",] > 5){
  #    if (!is.null(focal_level)) {
  #     weighted <- weightit(matching_formula,
  #     data = df, estimand = "ATT",
  #     focal=focal_level,
  #     method = "gbm")
  # 
  # } else {
  #     weighted <- weightit(matching_formula,
  #     data = df, estimand = "ATE",
  #     method = "gbm")
  # }
  # 
  #   match_summary <- bal.tab(weighted, stats = c("m", "v"), thresholds = c(m = .2))
  # }
  
  # p <- love.plot(weighted, stars = "std",
  #           abs=TRUE, line=T,
  #           thresholds = c(m = .2),
  #           limits=list(m = c(0, .5)),
  #           sample.names = c("Unweighted", "PS Weighted (ATT)"),
  #           colors = c("grey", "black"))
  
  return(list(p,weighted))
}


weightedSurvival <- function(df, survform, weighted) {
  
  cat("\nKM curve")
  
  # unweighted_km <- survfit(as.formula(survform), data=df)
  # 
  # s1 <- ggsurvplot(unweighted_km, risk.table = T, marks=T,
  #                 xlabs = "Time (months)", ystrataname = "",
  #                 data=df)
  
  weighted_km <- do.call(survfit, list(formula=as.formula(survform), data=df,
                                       weights = weighted$weights))
  
  s2 <- do.call(ggsurvplot, list(fit=weighted_km, risk.table = T, marks=T, conf.int=T,
                                 xlabs = "Time (months)", ystrataname = "", 
                                 data=df))
  
  cat("\nModel PFS with Cox PH")
  
  print(survform)
  
  weighted_cox=coxph(as.formula(survform),
                     data=df,weights=weighted$weights)
  
  cox_noweight=coxph(as.formula(survform),
                     data=df)
  
  g1 <- ggforest(weighted_cox, data=df, main="")
  #g2 <- ggforest(cox_noweight, data=df, main="Hazard ratio (Unweighted)")
  
  #g_grid <- plot_grid(g2, g1, nrow=1)
  
  return_list <- list(s2, g1, weighted_km, weighted_cox)
  names(return_list) <- c("km_plot_weighted", "coxph_forest_plots", "km_object", "coxph_models")
  
  return(return_list)
  
}


weightedSurvivalGeneList <- function(maf, clinical, gene_to_run, 
                                     outcomes, matching_vars,
                                     survobjform, 
                                     panel_presence=NULL,
                                     count_threshold,
                                     baseline="other"){
  
  surv_list_snv <- list()
  for (g in 1:length(gene_to_run)){
    temp_gene <- gene_to_run[g]
    cat("\n", temp_gene)
    temp_df <- prepareGenomicDF(maf = maf, clinical = clinical, gene=temp_gene, panel_presence=panel_presence)
    
    surv_list <- list()
    for (i in 1:length(outcomes)){
      tmp_outcome <- outcomes[i]
      print(tmp_outcome)
      
      if (baseline == "oncogenic"){
        rundf <- temp_df[[1]][get(tmp_outcome) %in% c("OncoKB_Oncogenic", "Rescued_Oncogenic"), ]
      } else if (baseline == "benign") {
        rundf <- temp_df[[1]][get(tmp_outcome) %in% c("No_mutation", "Rescued_Benign"), ]
      } else {
        rundf <- temp_df[[1]]
      }
      
      print(rundf[, .N, by=tmp_outcome])
      outcome_level <- setdiff(as.character(rundf[, .N, by=tmp_outcome][N >= count_threshold, ][[1]]),
                               "Unknown")
      
      rundf <- rundf[rundf[[tmp_outcome]] %in% outcome_level, ]
      rundf[[tmp_outcome]] <- droplevels(rundf[[tmp_outcome]])
      
      tmp_weights <- tryCatch(calculateIPTWeights(lhs = tmp_outcome,
                                                  rhs = matching_vars,
                                                  df = rundf),
                              error=function(e) NULL)


      if (!is.null(tmp_weights)){
        #print(tmp_weights[[1]])
        if (baseline == "oncogenic"){
          rundf[[tmp_outcome]] <- relevel(rundf[[tmp_outcome]], ref="OncoKB_Oncogenic") 
        } else if (baseline == "benign") {
          rundf[[tmp_outcome]] <- relevel(rundf[[tmp_outcome]], ref="No_mutation") 
        }
        
        surv_list[[i]] <- weightedSurvival(df=rundf, 
                                           survform=paste0(survobjform, "~", tmp_outcome), 
                                           weighted=tmp_weights[[2]])
      } else { surv_list[[i]] <- NULL }
    }
    surv_list_snv[[g]] <- surv_list
  }
  
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
  coefs_df <- coefs_df[, estimate_sig := ifelse(pval_adj <= 0.2, hr, NA)]
  
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
                                                                  
                                                                   grepl("gene_of_interest", term), "Gene",
                                                                   default="Other")]
  
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


getpercent <- function(x) { return(round(x*100, 2)) }

# create a comparison df where all oncogenic + reclassified oncogenic mutations in all genes within a pathway count. 
# basically set up a fisher's test where it's known oncogenic vs. known oncogenic + reclassified oncogenic 
createReclassifiedGAMs <- function(tmp, method, gene_list){
  
  # pathway <- pathway_data$NRF2
  # method <- methods[1]
  
  # tmp <- Reduce(function(x,y) merge(x,y, by="SAMPLE_ID", all=T), pathway)
  tmpgam <- tmp[, c("SAMPLE_ID", grep(paste0(method, "_Rescue"), colnames(tmp), value=T)), with=F]
  
  tmpgam <- tmpgam[, grep(paste0(method, "_"), colnames(tmpgam), value=T, fixed=T) := lapply(.SD, function(x) ifelse(x %in% c("Rescued_Oncogenic"), 1, 0)), .SDcols=grep(paste0(method, "_"), colnames(tmpgam), value=T)]
  
  tmpgam <- setnames(tmpgam, colnames(tmpgam), gsub(paste0("_", method, "_Rescue"), "", colnames(tmpgam)))
  
  tmpgam <- setnames(tmpgam, "SAMPLE_ID", "Tumor_Sample_Barcode")
  
  tmpgam <- rbind(tmpgam, gam_oncogenic[Tumor_Sample_Barcode %in% tmpgam$Tumor_Sample_Barcode, ][, colnames(gam_oncogenic) %in% c("Tumor_Sample_Barcode", gene_list), with=F], use.names=T, fill=T)
  
  tmpgam[is.na(tmpgam)] <- 0
  
  tmpgam <- unique(tmpgam[, lapply(.SD, sum), by=.(Tumor_Sample_Barcode), .SDcols=gene_list])
  
  return(tmpgam)
  
}


testMutualExclusivity <- function(pathway, gene_list, side1, side2, pairwise=T, oncokb_included=T, methods){
  
  tmp <- Reduce(function(x,y) merge(x,y, by="SAMPLE_ID", all=T), pathway)
  print(colnames(tmp))
  
  gene_results <- list()
  for (g in 1:length(gene_list)){
    anchor_gene <- gene_list[g] # one gene at a time
    print(anchor_gene)
    other_genes <- gene_list # test all other genes
    
    tmp_anchor <- tmp[, unique(c("SAMPLE_ID", 
                                 grep("Hugo_Symbol", colnames(tmp), value=T),
                                 grep(anchor_gene, colnames(tmp), value=T))), 
                      with=F]
    
    anchor_column <- paste0(anchor_gene, "_oncogenic_numeric")
    
    method_results <- list()
    
    for (i in 1:length(methods)){
      rescue_method <- methods[i]
      print(rescue_method)
      
      if (side2 == "oncogenic") {
        if (rescue_method == "oncogenic"){
          tmp_anchor <- tmp[, anchor_gene_altered := ifelse(get(anchor_column) == 1, 1, 0)] 
          newgam <- gam_oncogenic[, colnames(gam_oncogenic) %in% c("Tumor_Sample_Barcode", other_genes), with=F]
          
        } else if (oncokb_included == T & rescue_method != "oncogenic") { # known + reclassified oncogenic vs known + reclassified oncogenic
          anchor_rescue <- paste0(anchor_gene,  "_", rescue_method, "_Rescue") 
          tmp_anchor <- tmp[, anchor_gene_altered := ifelse(get(anchor_rescue) ==  "Rescued_Oncogenic" | get(anchor_column) == 1, 1, 0)]
          newgam <- createReclassifiedGAMs(tmp, rescue_method, other_genes) 
          
        } else if (oncokb_included == F  & rescue_method != "oncogenic"){ # reclassified oncogenic vs known oncogenic
          anchor_rescue <- paste0(anchor_gene,  "_", rescue_method, "_Rescue") 
          tmp_anchor <- tmp[, anchor_gene_altered := ifelse(get(anchor_rescue)  == "Rescued_Oncogenic", 1, 0)]
          newgam <- gam_oncogenic[, colnames(gam_oncogenic) %in% c("Tumor_Sample_Barcode", other_genes), with=F]
        }
        
      } else if (side2 == "benign") {
        if (rescue_method == "oncogenic"){
          tmp_anchor <- tmp[, anchor_gene_altered := ifelse(get(anchor_column) == 0, 1, 0)] 
          newgam <- gam_oncogenic[, colnames(gam_oncogenic) %in% c("Tumor_Sample_Barcode", other_genes), with=F]
          
        } else if (oncokb_included == T & rescue_method != "oncogenic") {
          anchor_rescue <- paste0(anchor_gene,  "_", rescue_method, "_Rescue") 
          tmp_anchor <- tmp[, anchor_gene_altered := ifelse(get(anchor_rescue)  == "Rescued_Benign" | get(anchor_column) == 0, 1, 0)]
          newgam <- createReclassifiedGAMs(tmp, rescue_method, other_genes)
          
        } else if (oncokb_included == F  & rescue_method != "oncogenic"){
          anchor_rescue <- paste0(anchor_gene,  "_", rescue_method, "_Rescue") 
          tmp_anchor <- tmp[, anchor_gene_altered := ifelse(get(anchor_rescue) == "Rescued_Benign", 1, 0)]
          newgam <- gam_oncogenic[, colnames(gam_oncogenic) %in% c("Tumor_Sample_Barcode", other_genes), with=F]
        }
        
      }
      
      tmp_anchor <- merge(tmp_anchor,
                          newgam,
                          by.x="SAMPLE_ID", by.y="Tumor_Sample_Barcode")
      
      other_genes_included <- intersect(other_genes, colnames(tmp_anchor))
      
      tmp_anchor <- tmp_anchor[, other_gene_altered := as.integer(rowSums(.SD != 0) > 0),
                               .SDcols=other_genes_included]
      
      if (pairwise == T){
        tmp_anchor_subset <- lapply(other_genes_included, function(x) table(tmp_anchor[, c("anchor_gene_altered", x), with=F]))
        names(tmp_anchor_subset) <- other_genes_included
        tmp_anchor_subset <- tmp_anchor_subset[which(sapply(tmp_anchor_subset, function(x) identical(dim(x), dim(ref_matrix))) == T)]
        test_df <- lapply(tmp_anchor_subset, function(x) setDT(tidy(fisher.test(x, alternative = "two.sided"))))
        if (length(test_df) > 0){
          for (t in 1:length(test_df)){
            test_df[[t]] <- test_df[[t]][, `:=`(anchor_gene = anchor_gene, 
                                                test_gene = names(tmp_anchor_subset)[t],
                                                method = rescue_method)]
          }
          test_df <- rbindlist(test_df)
        } else { test_df <- NULL }
        
        
      } else {
        
        testxtab <- table(tmp_anchor[, .(anchor_gene_altered, other_gene_altered)])
        
        if (0 %in% testxtab) {
          testxtab <-  testxtab + 0.5
        } 
        
        print(testxtab)
        
        if (identical(dim(testxtab), dim(ref_matrix))){
          test_df <- setDT(tidy(fisher.test(testxtab, alternative = "two.sided")))
          test_df <- test_df[, `:=`(gene = anchor_gene, method = rescue_method)]
        } else { test_df <- NULL }
      }
      
      method_results[[i]] <- test_df
      
    }
    
    gene_results[[g]] <- rbindlist(method_results)
    
  }
  
  results_df <- rbindlist(gene_results)[, `:=`(logOR = log(estimate),
                                               logOR_lowCI = log(conf.low),
                                               logOR_highCI = log(conf.high),
                                               anchor_gene_mutations_included = side2)]
  
  return(results_df)
  
}


processGeneList <- function(pathway, maf, clinical){
  
  gene_list <- intersect(pathways_df[Pathway == pathway, ]$Gene, unique(all_missense_process_nsclc$Hugo_Symbol))
  
  tmp_gene_process <- list()
  
  for (g in 1:length(gene_list)){
    tmp_gene_process[[g]] <- prepareGenomicDF(maf = maf, clinical = clinical, gene=gene_list[g])[[1]]
    tmp_gene_process[[g]] <- tmp_gene_process[[g]][, grepl("numeric|Hugo_Symbol|HGVSp_Short|SAMPLE_ID|_Rescue|", colnames(tmp_gene_process[[g]])), with=F]
    tmp_gene_process[[g]] <- setnames(tmp_gene_process[[g]], colnames(tmp_gene_process[[g]]), paste0(gene_list[g], "_", colnames(tmp_gene_process[[g]])))
    tmp_gene_process[[g]] <- setnames(tmp_gene_process[[g]], paste0(gene_list[g], "_SAMPLE_ID"), "SAMPLE_ID")
  }
  
  return(tmp_gene_process)
  
}

resultDTProcess <- function(tmp_results, tmpp){
  tmp_results <- tmp_results[, pval.adj := p.adjust(p.value, method = "fdr")][, pathway := tmpp]
  tmp_count <- tmp_results[, length(which(pval.adj <= 0.1 & estimate < 1 | pval.adj <= 0.1 & is.infinite(estimate)))/length(tmp_gene_list), by=c("method", "pathway")]
  tmp_count <- tmp_count[, method := case.(grepl("oncogenic", method), "OncoKB",
                                           grepl("CADD", method), "CADD",
                                           grepl("Clinvar", method), "ClinVar",
                                           grepl("SIFT", method), "SIFT",
                                           grepl("PolyPhen", method), "PolyPhen",
                                           grepl("MA", method), "MutationAssessor",
                                           grepl("AM", method), "AlphaMissense",
                                           grepl("REVEL", method), "REVEL",
                                           default=method)]
  return(list(results=tmp_results, count=tmp_count))
}


resultListProcess <- function(results_df, count_df) {
  pathway_results_df <- rbindlist(results_df)
  
  pathway_count_df <- merge(pathway_results_df[, length(which(pval.adj <= 0.1 & estimate < 1 | pval.adj <= 0.1 & is.infinite(estimate))), by=c("method", "pathway")],
                            pathway_results_df[, max(length(unique(gene))), by=c("pathway")], by="pathway")[, V1 := V1.x/V1.y]
  
  
  pathway_count_df <- pathway_count_df[, method := ifelse(grepl("oncogenic", method), "OncoKB", method)]
  
  pathway_count_df <- pathway_count_df[, method := relevel(factor(method), ref="OncoKB")]
  
  pathway_count_df <- pathway_count_df[, pathway := sub("_", " ", pathway)]
  
  return(list(results=pathway_results_df, count=pathway_count_df))
}
