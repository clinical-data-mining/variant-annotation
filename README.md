# Machine learning predictions improve identification of real-world cancer driver mutations

This is the code repository accompanying the manuscript *Machine learning predictions improve identification of real-world cancer driver mutations*. Description of the project is [here](https://www.abstractsonline.com/pp8/#!/20272/presentation/9170). Pre-print forthcoming.

## Requirements

***

All analyses for the manuscript were performed in `R 4.2.2` (2022-10-21) and the following package versions:

forcats_1.0.0      ranger_0.15.1      tidytext_0.4.1
RColorBrewer_1.1-3 stringr_1.5.0      stringdist_0.9.10 
survey_4.1-1       Matrix_1.5-1       cowplot_1.1.1     
cobalt_4.4.1       WeightIt_0.13.1    survminer_0.4.9   
ggpubr_0.6.0       survival_3.5-3     here_1.0.1        
ggplot2_3.4.4      tidytable_0.9.2    corrplot_0.92     
xtable_1.8-4       data.table_1.14.6

## Components

***

### data

--msk_impact_nsclc_clinical.csv: flat file containing clinical characteristics, including overall survival status and time, from 7,965 patients with non-small cell lung cancer in the MSK-IMPACT cohort

--msk_impact_nsclc_maf_annotated.csv: MAF file containing mutation data and additional annotations from 11 VEPs for 7,965 patients with NSCLC in the MSK-IMPACT cohort

### src

--OS_reclassified_mutation.Rmd: run weighted OS analysis stratified by reclassified mutational status

--pathway_mutual_exclusivity.Rmd: run mutual exclusivity analysis for oncogenic pathways in Sanchez-Vega (2018)

--utils.R: functions underlying analyses



