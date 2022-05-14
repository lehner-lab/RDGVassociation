---
  title: "Count number of RDGVs across SNP sets and the two cohorts (TCGA-WES and PCAWG_Hartwig-WGS)"
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


##upload files with gnomAD info
counted_rare_variants_TCGA_gnomAD <- read.csv(file = './TCGA_variants_gnomAD.txt',head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  mutate(cohort='Discovery-WES\nTCGA')

counted_rare_variants_PCAGW_Hartwig_gnomAD <- read.csv(file = './PCAWG_Hartwig_variants_gnomAD.txt',head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  mutate(cohort='Validation-WGS\nPCAWG_Hartwig')

##upload files with gnomAD info
counted_rare_variants_TCGA_random <- read.csv(file = './TCGA_counted_individuals_rare.txt',head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  mutate(cohort='Discovery-WES\nTCGA')

counted_rare_variants_PCAGW_Hartwig_random <- read.csv(file = './PCAWG_Hartwig_counted_individuals_rare.txt',head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  mutate(cohort='Validation-WGS\nPCAWG_Hartwig')


##upload FDR1/2 genes
genes_1FDR <- read.csv(file= "./validated_genes.txt",head=T,sep ="\t",stringsAsFactors = F);
genes_2FDR <- read.csv(file= "./validated_genes_FDR2.txt",head=T,sep ="\t",stringsAsFactors = F) %>%
  dplyr::filter(!variant %in% genes_1FDR$variant);


##distribution across gene sets
counted_rare_variants_frequencies <- counted_rare_variants_TCGA_gnomAD %>%
  rbind(counted_rare_variants_PCAGW_Hartwig_gnomAD) %>%
  dplyr::filter(!is.na(controls_AF)) %>%
  mutate(RDGV_set= sub('MRT','MTR',RDGV_set)) %>%
  mutate(RDGV_set=sub('MTR_25thCentile','Misssense_MTR',RDGV_set)) %>%
  mutate(RDGV_set=sub('CCR_90orhigher','Misssense_CCR',RDGV_set)) %>%
  dplyr::filter(RDGV_set != 'Misssense_CCR' & RDGV_set != 'Misssense_MTR') %>%
  mutate(gene_set= if_else(Gene.refGene %in% c('MSH2','MLH1','PMS2','MSH6'),"dMMR genes",
                           if_else(Gene.refGene %in% c('BRCA2','BRCA1','PALB2','RAD51C'),"dHR genes",
                                   if_else(Gene.refGene %in% genes_1FDR$variant,"Replicated <= 1% FDR",
                                           if_else(Gene.refGene %in% genes_2FDR$variant,"Replicated 1-2% FDR","none"))))) %>%
  dplyr::filter(gene_set != 'none')  %>%
  group_by(gene_set,cohort,RDGV_set) %>%
  summarise(n_carriers=n(),
            gnomAD_controls_freq_nfe=sum(controls_AF_nfe, na.rm = T)) %>%
  ungroup() %>%
  mutate(carriers_freq= if_else(cohort=='Discovery-WES\nTCGA', n_carriers/6799, n_carriers/4683)) %>%
  dplyr::select(-n_carriers) %>%
  melt(id.vars=c('gene_set','cohort','RDGV_set')) 

##distribution across gene sets/random set
counted_rare_variants_frequencies_random <- counted_rare_variants_TCGA_random %>%
  rbind(counted_rare_variants_PCAGW_Hartwig_random) %>%
  dplyr::filter(!is.na(AF)) %>%
  mutate(RDGV_set= sub('MRT','MTR',RDGV_set)) %>%
  mutate(RDGV_set=sub('MTR_25thCentile','Misssense_MTR',RDGV_set)) %>%
  mutate(RDGV_set=sub('CCR_90orhigher','Misssense_CCR',RDGV_set)) %>%
  mutate(gene_set=sub('Replicated 2% FDR','Replicated 1-2% FDR',gene_set)) %>%
  mutate(gene_set=sub('Replicated 1% FDR','Replicated <= 1% FDR',gene_set)) %>%
  dplyr::filter(RDGV_set != 'Misssense_CCR' & RDGV_set != 'Misssense_MTR') %>%
  group_by(gene_set,cohort,RDGV_set,set) %>%
  summarise(n_carriers=n()) %>%
  ungroup() %>%
  dplyr::filter(set != 'real') %>%
  mutate(gene_set= sub(' matched','',gene_set)) %>%
  mutate(set='Cancer Dataset random length-matched genes') %>%
  mutate(variable=set) %>%
  mutate(value= if_else(cohort=='Discovery-WES\nTCGA', n_carriers/6799, n_carriers/4683)) %>%
  dplyr::select(-n_carriers, -set) 

##combine all
frequencies_variants <-counted_rare_variants_frequencies %>%
  rbind(counted_rare_variants_frequencies_random) %>%
  mutate(gene_set= factor(gene_set, levels =c('dHR genes','dMMR genes','Replicated <= 1% FDR','Replicated 1-2% FDR'))) %>%
  mutate(RDGV_set= factor(RDGV_set, levels =c('PTV','PTV_Misssense_CADD25','PTV_Misssense_CADD15','Misssense_CCR','Misssense_MTR'))) %>%
  mutate(variable=sub('gnomAD_controls_freq_nfe','gnomAD controls Non-Finnish European',variable)) %>%
  mutate(variable=sub('carriers_freq','Cancer Dataset',variable)) %>%
  mutate(variable= factor(variable, levels =c('gnomAD controls Non-Finnish European','Cancer Dataset','Cancer Dataset random length-matched genes'))) %>%
  dplyr::mutate(value=value*100) %>%
  mutate(cohort=sub('Validation-WGS\nPCAWG_Hartwig','Validation-WGS PCAWG_Hartwig',cohort)) %>%
  mutate(cohort=sub('Discovery-WES\nTCGA','Discovery-WES TCGA',cohort)) 


##plot frequencies
pdf(file= paste(output_file,'rare_variants_freq_gnomad','.pdf',sep=''),useDingbats=FALSE,
    width= 8,
    height = 6)  
ggplot(frequencies_variants,
       aes(x=gene_set, y=value, color=variable)) +
  geom_boxplot(outlier.shape = NA,position="dodge") + 
  geom_jitter(size=0.3,shape=16,
              position=position_jitterdodge()) +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, color='black', hjust=1, angle=70),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=8, color='black'),
        strip.text.x = element_text(size=8, color='black'),
        strip.text.y = element_text(size=8, color='black'),
        legend.title = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(size=8, color='black')) + 
  scale_color_manual(values = c("#3B9AB2","#F21A00","#EBCC2A")) +
  #scale_y_continuous(limits = c(0,45)) +
  ylab('% of individuals\ncarrying rare pLoF variant') +
  xlab('') +
  facet_wrap(cohort ~ RDGV_set, scales='free') +
  guides(color = guide_legend(nrow = 2)) 
dev.off()