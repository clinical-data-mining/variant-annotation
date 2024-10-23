# Machine learning cancer driver mutation predictions are valid in real-world data

This is the code repository accompanying the manuscript *Machine learning cancer driver mutation predictions are valid in real-world data* (Tran et al., 2024). Description of the project is [here](https://www.abstractsonline.com/pp8/#!/20272/presentation/9170). Preprint is [here](https://www.biorxiv.org/content/10.1101/2024.03.31.587410v2).

## Introduction

This codebase consists of three R notebooks that perform benchmarking tasks for multiple variant effect predictor (VEP) tools, focusing on survival outcomes, pathway mutual exclusivity, and occurrence at binding residues using real-world cancer data. These notebooks are designed to replicate the findings described in the manuscript and can be run with new data with similar format as the sample data we provided in `data`. 

We do not provide scripts to annotate variants with VEPs here, as the annotations can be more easily incorporated into upstream bioinformatic workflows e.g. with Ensembl VEP, or using the platforms recommended by the authors of individual methods. For researchers interested in annotating their data with multiple VEPs using pre-calculated scores, we recommend [dbNSFP](http://database.liulab.science/dbNSFP) (Liu et al., 2023) and [OpenCRAVAT](https://opencravat.org/) (Pagel et al., 2020) as two open-source, frequently updated resources. 

## Requirements

All analyses for the manuscript were performed in `R 4.2.2` (2022-10-21) and the following package versions: data.table_1.14.6, survival_3.5-3, survminer_0.4.9, survey_4.1-1, WeightIt_0.13.1, ggplot2_3.4.4, forcats_1.0.0, ranger_0.15.1, tidytext_0.4.1, RColorBrewer_1.1-3, stringr_1.5.0, stringdist_0.9.10, Matrix_1.5-1, cowplot_1.1.1, cobalt_4.4.1, ggpubr_0.6.0, here_1.0.1, tidytable_0.9.2, corrplot_0.92, xtable_1.8-4       

Installation of all dependencies should not take more than 30 minutes.

## How to run 
Set up and load data: Load the clinical data, MAF file, and optional GAM using fread() from the data.table package. Analyses can be run from individual Rmd files. 

### Binding residue enrichment (`binding_reclassified_mutation.Rmd`)
- Fisher’s tests: Performing Fisher's exact tests to evaluate the association between reclassified mutations (pathogenic or benign) and whether they occur at binding residues.

- Visualization: Plotting the results for different genes and methods, highlighting the proportion of mutations predicted to be oncogenic at binding residues and the statistical significance of these predictions.

- Performance metrics: Calculating odds ratios (OR) for reclassified mutations to assess the ability of different VEPs to distinguish pathogenic  from benign mutations.

### Survival outcome analysis (`OS_reclassified_mutation.Rmd`)
- Run Cox PH models: Define the genes of interest and use the weightedSurvivalGeneList function to stratify outcomes based on mutational status. This will run a series of Cox PH models for each gene and compile the results.

- Plot results: Use ggplot2 to create a bubble plot that visualizes the hazard ratios (HRs) for reclassified pathogenic mutations versus no mutation, categorized by different VEPs.

- Run control models: Compare patients with reclassified benign mutations versus no mutation, using the same workflow as for oncogenic mutations.

- Calculate relative risk (RR): For a gene of interest, calculate the RR of survival between reclassified pathogenic and benign mutations using the results of the previous analyses in order to compare VEP performance in identifying prognostic reclassified pathogenic mutations.

### Pathway mutual exclusivity analysis (`pathway_mutual_exclusivity.Rmd`)
- Run Fisher's tests for mutual exclusivity between reclassified oncogenic mutations vs. all other known oncogenic mutations: Given a list of genes within each pathway (we provided the gene lists for oncogenic signaling pathways defined in Sanchez-Vega et al. in `data/Pathways.Rdata`), the notebook uses a custom testMutualExclusivity function to run a series of Fisher's tests to test for mutual exclusivity between reclassified pathogenic mutations in each gene with other known oncogenic mutations in other genes within the same pathway

- Run control tests: the same Fisher's tests are run to test for mutual exclusivity between reclassified benign mutations in each gene with other known oncogenic mutations in other genes within the same pathway

- Performance metrics: Calculating odds ratios (OR_binding) of reclassified oncogenic mutations, compared to reclassified benign mutations, being mutually exclusive with other known oncogenic mutations in the same pathway.

## Repository components

### data

--msk_impact_nsclc_clinical.csv: flat file containing clinical characteristics, including overall survival status and duration, from 7,965 patients with non-small cell lung cancer in the MSK-IMPACT cohort

--msk_impact_nsclc_missense_mutations.csv: flat file containing missense mutation data as well as annotations from OncoKB, ClinVar and 13 VEPs for genes altered >= 2% in 7,965 patients with NSCLC in the MSK-IMPACT cohort

--msk_impact_gam.csv: flat file containing the gene alteration matrix, which is a binary matrix with samples as rows and genes as columns. A cell takes value 1 if a gene is altered in a sample, and 0 otherwise.

--BioLiP_PPI_processed.RData: R datafile containing pre-processed [BioLiP](https://zhanggroup.org/BioLiP/index.cgi) and [PPI-HotspotDB](https://pubs.acs.org/doi/10.1021/acs.jcim.2c00025) databases

--Pathways.Rdata: R datafile containing gene lists of oncogenic signaling pathways in Sanchez-Vega et al. (2018)

### src

--OS_reclassified_mutation.Rmd: run weighted OS analysis stratified by reclassified mutational status

--pathway_mutual_exclusivity.Rmd: run mutual exclusivity analysis for oncogenic pathways in Sanchez-Vega et al. (2018)

--binding_reclassified_mutation.Rmd: run Fisher's tests to test for enrichment of reclassified oncogenic mutations at binding residues

--utils.R: custom functions underlying analyses

## Running instruction

The .Rmd files are R notebooks that can be open and executed using RStudio. Running all chunks of codes will return the result tables used to created plots in Figure 2B and 2C of the manuscript.

Expected run time is ~30 minutes on a Mac laptop with 16GB RAM.



