---
  title: "Upload the data from all models and set FDR thresholds"
---
  
#rm(list=ls())

library(tidyverse)
library(dplyr)


##directories
input_file_direc='./003_association_testing_results/input_files/'
output_file_direc='./003_association_testing_results/results/'
output_figure_direc='./003_association_testing_results/figures/'

##upload models to be discarded, extracted via lambdas.R script
models_discard <- read.csv(file = paste(output_file_direc,'models_removed_based_on_lambda.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) 


####### estimation of empirical FDRs for TCGA-discovery
##upload SKAT-O results from randomization
#remove models which will be removed based on inflation factor
SKAT_res <- read.csv(file= paste(input_file_direc,'Discovery_TCGA_SKATO_random.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  dplyr::filter(!model_snps_cancer_pheno %in% models_discard$model_snps_cancer_pheno)


##make plot of the p value distributions of all combined
p_values_random <- SKAT_res %>%
  group_by(cancer_type) %>%
  mutate(LOG10pValue=-log10(pValue)) %>%
  mutate(rankP=rank(-LOG10pValue)) %>%
  ungroup()  

##save the FDR information
FDR_calculated <- c()

for(i in unique(p_values_random$cancer_type)){
  rm(p_values_random_cancer_type)
  rm(FDR_cancer_type) 
  
  p_values_random_cancer_type <- p_values_random %>%
    dplyr::filter(cancer_type==i)
  
  FDR_cancer_type <- data.frame(FDR= seq(0,1,0.01),
                                FDR_threshold= quantile(p_values_random_cancer_type$pValue, seq(0,1,0.01)),
                                stringsAsFactors = F)
  FDR_cancer_type$cancer_type <- i
  
  FDR_calculated <- FDR_calculated %>%
    rbind(FDR_cancer_type)
}

write.table(FDR_calculated,
            file= paste(output_file_direc,'TCGA_empirical_FDRs.txt',sep=''),
            quote=FALSE, sep='\t',row.names=FALSE,col.names = T)



####### estimation of empirical FDRs for PCAWG_Hartwig-validation
##upload SKAT-O results from randomization
#remove models which will be removed based on inflation factor
SKAT_res <- read.csv(file= paste(input_file_direc,'Validation_PCAWG_Hartwig_SKATO_random.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  dplyr::filter(!model_snps_cancer_pheno %in% models_discard$model_snps_cancer_pheno)


##make plot of the p value distributions of all combined
p_values_random <- SKAT_res %>%
  group_by(cancer_type) %>%
  mutate(LOG10pValue=-log10(pValue)) %>%
  mutate(rankP=rank(-LOG10pValue)) %>%
  ungroup()  

##save the FDR information
FDR_calculated <- c()

for(i in unique(p_values_random$cancer_type)){
  rm(p_values_random_cancer_type)
  rm(FDR_cancer_type) 
  
  p_values_random_cancer_type <- p_values_random %>%
    dplyr::filter(cancer_type==i)
  
  FDR_cancer_type <- data.frame(FDR= seq(0,1,0.01),
                                FDR_threshold= quantile(p_values_random_cancer_type$pValue, seq(0,1,0.01)),
                                stringsAsFactors = F)
  FDR_cancer_type$cancer_type <- i
  
  FDR_calculated <- FDR_calculated %>%
    rbind(FDR_cancer_type)
}

write.table(FDR_calculated,
            file= paste(output_file_direc,'PCAWG_Hartwig_empirical_FDRs.txt',sep=''),
            quote=FALSE, sep='\t',row.names=FALSE,col.names = T)


