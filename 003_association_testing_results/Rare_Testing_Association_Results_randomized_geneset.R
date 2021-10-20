---
  title: "Main results - number of hits - number of replicated hits - etc - empirical FDR of 1% and 2% via randomized geneset"
---
  
#rm(list=ls())
  
library(tidyverse)
library(dplyr)
library(wesanderson)
library(reshape2)


##directories
input_file_direc='./003_association_testing_results/input_files/'
output_file_direc='./003_association_testing_results/results/'
output_figure_direc='./003_association_testing_results/figures/'


##discard models based on lambda
models_discard <- read.csv(file = paste(output_file_direc,'models_removed_based_on_lambda_randomized_geneset.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) 



##upload results from burden
TCGA_burden_results <- read.csv(file = paste(input_file_direc,'Discovery_TCGA_burden_results_randomized_geneset.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE)  
colnames(TCGA_burden_results)[2] <- 'burden_test_pValue'
WGS_burden_results <- read.csv(file = paste(input_file_direc,'Validation_PCAWG_Hartwig_burden_results_randomized_geneset.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE)  
colnames(WGS_burden_results)[2] <- 'burden_test_pValue'

##upload results from SKAT
TCGA_SKAT_results <- read.csv(file = paste(input_file_direc,'Discovery_TCGA_SKATO_results_randomized_geneset.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) 
WGS_SKAT_results <- read.csv(file = paste(input_file_direc,'Validation_PCAWG_Hartwig_SKATO_results_randomized_geneset.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) 

########## empirical FDR of 1% ########## ########## ########## ########## ########## 
#empirical FDR of 1%; can be changed here if needed
TCGA_FDR <- read.csv(file = paste(output_file_direc,'TCGA_empirical_FDRs_random_geneset.txt',sep=''), head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  dplyr::filter(FDR==0.01)
WGS_FDR <- read.csv(file = paste(output_file_direc,'PCAWG_Hartwig_empirical_FDRs_random_geneset.txt',sep=''), head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  dplyr::filter(FDR==0.01)

##combine discovery SKAT results with burden and remove models based on inflation factors
TCGA_results <- TCGA_SKAT_results %>%
  left_join(TCGA_burden_results, by=c('gene','pheno','cancer_type','snps','model','pheno_name','model_snps','model_snps_cancer_pheno')) %>%
  dplyr::filter(!is.na(estimate)) %>%
  left_join(TCGA_FDR, by=c('cancer_type')) %>%
  mutate(cohort='TCGA')  %>%
  mutate(model_snps_cancer_pheno=paste(model,snps,cancer_type,pheno,sep='__')) %>%
  dplyr::filter(!model_snps_cancer_pheno %in% models_discard$model_snps_cancer_pheno) 

##combine validations SKAT results with burden and remove models based on inflation factors
WGS_results <- WGS_SKAT_results %>%
  left_join(WGS_burden_results, by=c('gene','pheno','cancer_type','snps','model','pheno_name','model_snps','model_snps_cancer_pheno')) %>%
  dplyr::filter(!is.na(estimate)) %>%
  left_join(WGS_FDR, by=c('cancer_type')) %>%
  mutate(cohort='WGS') %>%
  mutate(model_snps_cancer_pheno=paste(model,snps,cancer_type,pheno,sep='__')) %>%
  dplyr::filter(!model_snps_cancer_pheno %in% models_discard$model_snps_cancer_pheno)

###identify replicated hits
validated_FDR_cohort <- TCGA_results %>%
  rbind(WGS_results) %>%
  dplyr::filter((cohort == 'TCGA' & pValue < FDR_threshold) | (cohort == 'WGS' & pValue < FDR_threshold )) %>% 
  dplyr::select(gene,pheno,cancer_type,model,snps,estimate,cohort,pheno_name) %>%
  mutate(gene_pheno_cancer_type_snp_model=paste(gene,pheno,cancer_type,snps,model,sep='__')) %>%
  mutate(gene_cancer_type= paste(gene,cancer_type,sep='__')) %>%
  spread(cohort,estimate) %>%
  dplyr::filter(sign(TCGA) == sign(WGS))

#print replicated hits
print(paste(nrow(validated_FDR_cohort), ' hits replicating at an FDR of 1% using the randomized genelist',sep=''))
print(paste(length(unique(validated_FDR_cohort$gene)),' genes replicating at an FDR of 1% using the randomized genelist',sep=''))
print(paste(unique(validated_FDR_cohort$gene), collapse = ', '))



########## empirical FDR of 2% ########## ########## ########## ########## ########## 
#empirical FDR of 2%; can be changed here if needed
TCGA_FDR <- read.csv(file = paste(output_file_direc,'TCGA_empirical_FDRs_random_geneset.txt',sep=''), head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  dplyr::filter(FDR==0.02)
WGS_FDR <- read.csv(file = paste(output_file_direc,'PCAWG_Hartwig_empirical_FDRs_random_geneset.txt',sep=''), head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  dplyr::filter(FDR==0.02)

##combine discovery SKAT results with burden and remove models based on inflation factors
TCGA_results <- TCGA_SKAT_results %>%
  left_join(TCGA_burden_results, by=c('gene','pheno','cancer_type','snps','model','pheno_name','model_snps','model_snps_cancer_pheno')) %>%
  dplyr::filter(!is.na(estimate)) %>%
  left_join(TCGA_FDR, by=c('cancer_type')) %>%
  mutate(cohort='TCGA')  %>%
  mutate(model_snps_cancer_pheno=paste(model,snps,cancer_type,pheno,sep='__')) %>%
  dplyr::filter(!model_snps_cancer_pheno %in% models_discard$model_snps_cancer_pheno) 

##combine validations SKAT results with burden and remove models based on inflation factors
WGS_results <- WGS_SKAT_results %>%
  left_join(WGS_burden_results, by=c('gene','pheno','cancer_type','snps','model','pheno_name','model_snps','model_snps_cancer_pheno')) %>%
  dplyr::filter(!is.na(estimate)) %>%
  left_join(WGS_FDR, by=c('cancer_type')) %>%
  mutate(cohort='WGS') %>%
  mutate(model_snps_cancer_pheno=paste(model,snps,cancer_type,pheno,sep='__')) %>%
  dplyr::filter(!model_snps_cancer_pheno %in% models_discard$model_snps_cancer_pheno)

###identify replicated hits
validated_FDR_cohort <- TCGA_results %>%
  rbind(WGS_results) %>%
  dplyr::filter((cohort == 'TCGA' & pValue < FDR_threshold) | (cohort == 'WGS' & pValue < FDR_threshold )) %>% 
  dplyr::select(gene,pheno,cancer_type,model,snps,estimate,cohort,pheno_name) %>%
  mutate(gene_pheno_cancer_type_snp_model=paste(gene,pheno,cancer_type,snps,model,sep='__')) %>%
  mutate(gene_cancer_type= paste(gene,cancer_type,sep='__')) %>%
  spread(cohort,estimate) %>%
  dplyr::filter(sign(TCGA) == sign(WGS))

#print replicated hits
print(paste(nrow(validated_FDR_cohort), ' hits replicating at an FDR of 1% using the randomized genelist',sep=''))
print(paste(length(unique(validated_FDR_cohort$gene)),' genes replicating at an FDR of 1% using the randomized genelist',sep=''))
print(paste(unique(validated_FDR_cohort$gene), collapse = ', '))



