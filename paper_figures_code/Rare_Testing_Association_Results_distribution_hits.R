---
  title: "Plot distribution of hits components, models of inheritance, etc."
---
  
#rm(list=ls())
  
library(tidyverse)
library(dplyr)
library(wesanderson)
library(reshape2)
library(RColorBrewer)
library(ggrastr)
library(ggpubr)


##directories
input_file_direc='./input_files/'
output_file_direc='./results/'
output_figure_direc='./figures/'


##discard models based on lambda
models_discard <- read.csv(file = paste(output_file_direc,'models_removed_based_on_lambda.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) 


##upload results from burden
TCGA_burden_results <- read.csv(file = paste(input_file_direc,'Discovery_TCGA_burden_results.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE)  
colnames(TCGA_burden_results)[2] <- 'burden_test_pValue'

##upload results from SKAT
TCGA_SKAT_results <- read.csv(file = paste(input_file_direc,'Discovery_TCGA_SKATO_results.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) 

##upload empirical FDR estimates
#empirical FDR of 1%; can be changed here if needed
TCGA_FDR <- read.csv(file = paste(output_file_direc,'TCGA_empirical_FDRs.txt',sep=''), head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  dplyr::filter(FDR==0.01)

##combine discovery SKAT results with burden and remove models based on inflation factors
TCGA_results <- TCGA_SKAT_results %>%
  left_join(TCGA_burden_results, by=c('gene','pheno','cancer_type','snps','model','pheno_name','model_snps','model_snps_cancer_pheno')) %>%
  dplyr::filter(!is.na(estimate)) %>%
  left_join(TCGA_FDR, by=c('cancer_type')) %>%
  mutate(cohort='TCGA')  %>%
  mutate(model_snps_cancer_pheno=paste(model,snps,cancer_type,pheno,sep='__')) %>%
  dplyr::filter(!model_snps_cancer_pheno %in% models_discard$model_snps_cancer_pheno) 


##upload replicated hits
validated_FDR_cohort_1FDR <- read.csv(file = paste(output_file_direc,'validated_genes.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) 
validated_FDR_cohort_2FDR <- read.csv(file = paste(output_file_direc,'validated_genes_FDR2.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) 


#####plot distribution of hits and validated hits across several factors

##cancer type
TCGA_hits_cancer_type <- TCGA_results %>%
  dplyr::filter(pValue < FDR_threshold) %>%
  group_by(cancer_type) %>%
  summarise(n=n()) %>%
  ungroup()  %>%
  mutate(cohort='Discovered hits')
Validated_hits_cancer_type <- validated_FDR_cohort_1FDR %>%
  group_by(cancer_type) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(cohort='Validated 1% FDR')
Validated_hits_cancer_type_2FDR <- validated_FDR_cohort_2FDR %>%
  group_by(cancer_type) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(cohort='Validated 2% FDR')
hits_cancer_type <- TCGA_hits_cancer_type %>%
  rbind(Validated_hits_cancer_type) %>%
  rbind(Validated_hits_cancer_type_2FDR) %>%
  spread(cohort,n,is.na(0)) %>%
  melt(id.vars=c('cancer_type'))


##plots number of validated hits only
hits_cancer_type_validated <- hits_cancer_type %>%
  dplyr::filter(variable != 'Discovered hits') %>%
  mutate(cancer_type=factor(cancer_type,ordered = F,levels = hits_cancer_type$cancer_type[order(hits_cancer_type$value[hits_cancer_type$variable=='Validated 1% FDR'], decreasing = T)])) 

pdf(file= paste(output_figure_direc,'distribution_validated_1FDR_hits_cancer_type','.pdf',sep=''),
    width= 2,
    height = 2)  
ggplot(hits_cancer_type_validated,
       aes(x=cancer_type,y=value,fill=variable)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  theme(plot.subtitle = element_text(hjust = 0.5,size=8,color='black'),
        axis.text.x = element_text(size=8,color='black',angle=70,hjust=1),
        axis.text.y = element_text(size=8,color='black'),
        axis.title = element_text(size=8,color='black'),
        legend.position = "none") + 
  scale_fill_manual(values=c("Discovered hits" = "#F21A00","Validated 1% FDR"="#3B9AB2","Validated 2% FDR"="#71B3C2")) +
  labs(fill='') +
  ylab('# replicated hits') +
  xlab('Cancer type')
dev.off()

##pheno
TCGA_hits_pheno <- TCGA_results %>%
  dplyr::filter(pValue < FDR_threshold) %>%
  group_by(pheno_name) %>%
  summarise(n=n()) %>%
  ungroup()  %>%
  mutate(cohort='Discovered hits')
Validated_hits_pheno <- validated_FDR_cohort_1FDR %>%
  group_by(pheno_name) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(cohort='Validated 1% FDR')
Validated_hits_pheno_type_2FDR <- validated_FDR_cohort_2FDR %>%
  group_by(pheno_name) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(cohort='Validated 2% FDR')
hits_pheno <- TCGA_hits_pheno %>%
  rbind(Validated_hits_pheno) %>%
  rbind(Validated_hits_pheno_type_2FDR) %>%
  spread(cohort,n,is.na(0)) %>%
  melt(id.vars=c('pheno_name'))

##plots number of validated hits only
hits_pheno_validated <- hits_pheno %>%
  dplyr::filter(variable != 'Discovered hits') %>%
  mutate(pheno_name=factor(pheno_name,ordered = F,levels = hits_pheno$pheno_name[order(hits_pheno$value[hits_pheno$variable=='Validated 1% FDR'], decreasing = T)])) 

pdf(file= paste(output_figure_direc,'distribution_validated_1FDR_hits_pheno','.pdf',sep=''),
    width= 4,
    height = 2)  
ggplot(hits_pheno_validated,
       aes(x=pheno_name,y=value,fill=variable)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  theme(plot.subtitle = element_text(hjust = 0.5,size=8,color='black'),
        axis.text.x = element_text(size=8,color='black',angle=70,hjust=1),
        axis.text.y = element_text(size=8,color='black'),
        axis.title = element_text(size=8,color='black'),
        legend.position = "none") + 
  scale_fill_manual(values=c("Discovered hits" = "#F21A00","Validated 1% FDR"="#3B9AB2","Validated 2% FDR"="#71B3C2")) +
  labs(fill='')  +
  ylab('# replicated hits') +
  xlab('Somatic Mutation Components')
dev.off()

##snps
TCGA_hits_snps <- TCGA_results %>%
  dplyr::filter(pValue < FDR_threshold) %>%
  group_by(snps) %>%
  summarise(n=n()) %>%
  ungroup()  %>%
  mutate(cohort='Discovered hits')
Validated_hits_snps <- validated_FDR_cohort_1FDR %>%
  group_by(snps) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(cohort='Validated 1% FDR')
Validated_hits_snps_2FDR <- validated_FDR_cohort_2FDR %>%
  group_by(snps) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(cohort='Validated 2% FDR')
hits_snps <- TCGA_hits_snps %>%
  rbind(Validated_hits_snps) %>%
  rbind(Validated_hits_snps_2FDR) %>%
  spread(cohort,n,is.na(0)) %>%
  melt(id.vars=c('snps')) %>%
  mutate(snps=sub('_0.1perc','',snps)) %>%
  mutate(snps=sub('MTR_25thCentile','Misssense_MTR',snps)) %>%
  mutate(snps=sub('CCR_90orhigher','Misssense_CCR',snps)) 


##plots number of validated hits only
hits_snps_validated <- hits_snps %>%
  dplyr::filter(variable != 'Discovered hits') %>%
  mutate(snps=factor(snps,ordered = F,levels = hits_snps$snps[order(hits_snps$value[hits_snps$variable=='Validated 1% FDR'], decreasing = T)])) 

pdf(file= paste(output_figure_direc,'distribution_validated_1FDR_hits_snps','.pdf',sep=''),
    width= 1.5,
    height = 2)  
ggplot(hits_snps_validated,
       aes(x=snps,y=value,fill=variable)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  theme(plot.subtitle = element_text(hjust = 0.5,size=8,color='black'),
        axis.text.x = element_text(size=8,color='black',angle=70,hjust=1),
        axis.text.y = element_text(size=8,color='black'),
        axis.title = element_text(size=8,color='black'),
        legend.position = "none") + 
  scale_fill_manual(values=c("Discovered hits" = "#F21A00","Validated 1% FDR"="#3B9AB2","Validated 2% FDR"="#71B3C2")) +
  labs(fill='')  +
  ylab('# replicated hits') +
  xlab('RDGV SNP set')
dev.off()

##model
TCGA_hits_model <- TCGA_results %>%
  dplyr::filter(pValue < FDR_threshold) %>%
  group_by(model) %>%
  summarise(n=n()) %>%
  ungroup()  %>%
  mutate(cohort='Discovered hits')
Validated_hits_model <- validated_FDR_cohort_1FDR %>%
  group_by(model) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(cohort='Validated 1% FDR')
Validated_hits_model_2FDR <- validated_FDR_cohort_2FDR %>%
  group_by(model) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(cohort='Validated 2% FDR')
hits_model <- TCGA_hits_model %>%
  rbind(Validated_hits_model) %>%
  rbind(Validated_hits_model_2FDR) %>%
  spread(cohort,n,is.na(0)) %>%
  melt(id.vars=c('model'))

##plots number of validated hits only
hits_model_validated <- hits_model %>%
  dplyr::filter(variable != 'Discovered hits') %>%
  mutate(model=factor(model,ordered = F,levels = hits_model$model[order(hits_model$value[hits_model$variable=='Validated 1% FDR'], decreasing = T)])) 

pdf(file= paste(output_figure_direc,'distribution_validated_1FDR_hits_model','.pdf',sep=''),
    width= 1.2,
    height = 2)  
ggplot(hits_model_validated,
       aes(x=model,y=value,fill=variable)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  theme(plot.subtitle = element_text(hjust = 0.5,size=8,color='black'),
        axis.text.x = element_text(size=8,color='black',angle=70,hjust=1),
        axis.text.y = element_text(size=8,color='black'),
        axis.title = element_text(size=8,color='black'),
        legend.position = "none") + 
  scale_fill_manual(values=c("Discovered hits" = "#F21A00","Validated 1% FDR"="#3B9AB2","Validated 2% FDR"="#71B3C2")) +
  labs(fill='')  +
  ylab('# replicated hits') +
  xlab('Model of inheritance')
dev.off()












