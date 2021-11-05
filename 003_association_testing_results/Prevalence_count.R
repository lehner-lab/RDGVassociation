---
  title: "Count number of RDGVs across SNP sets and the two cohorts (TCGA-WES and PCAWG_Hartwig-WGS)"
---
  
#rm(list=ls())
  
library(tidyverse)
library(dplyr)
library(wesanderson)
library(reshape2)


##directories
input_file_direc='./input_files/'
output_file_direc='./results/'
output_figure_direc='./figures/'


##upload files in which number of RDGVs was counted - also in length matched randomly seected genes
counted_rare_variants_TCGA <- read.csv(file =paste(input_file_direc,'TCGA_counted_individuals_rare.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  mutate(cohort='Discovery-WES\nTCGA')

counted_rare_variants_PCAGW_Hartwig <- read.csv(file = paste(input_file_direc,'PCAWG_Hartwig_counted_individuals_rare.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  mutate(cohort='Validation-WGS\nPCAWG_Hartwig')


##distribution across snp sets
counted_rare_variants_plot <- counted_rare_variants_TCGA %>%
  rbind(counted_rare_variants_PCAGW_Hartwig) %>%
  mutate(variable= sub('MRT','MTR',variable)) %>%
  mutate(variable=sub('MTR_25thCentile','Misssense_MTR',variable)) %>%
  mutate(variable=sub('CCR_90orhigher','Misssense_CCR',variable)) %>%
  mutate(gene_set=sub('Replicated 2% FDR','Replicated 2% FDR only',gene_set)) %>%
  mutate(set_group= if_else(grepl('random',set), 'random length-matched genes','real')) %>%
  mutate(gene_set_simple= sub(' matched','',gene_set)) %>%
  mutate(gene_set_simple= factor(gene_set_simple, levels =c('dHR genes','dMMR genes',
                                                            'Replicated 1% FDR',
                                                            'Replicated 2% FDR only'))) %>%
  mutate(variable= factor(variable, levels =c('PTV','PTV_Misssense_CADD25','PTV_Misssense_CADD15','Misssense_CCR','Misssense_MTR'))) %>%
  mutate(set_group= factor(set_group, levels =c('real','random length-matched genes'))) %>%
  dplyr::filter(variable != 'Misssense_CCR' & variable != 'Misssense_MTR') #excluded since almost no hits found here


##plot distribution across snp sets - percentage
pdf(file= paste(output_figure_direc,'rare_variants_precentage_individuals','.pdf',sep=''),
    width= 6,
    height = 4.5)  
ggplot(counted_rare_variants_plot %>%
         dplyr::filter(variable != 'Misssense_CCR' & variable != 'Misssense_MTR'),
       aes(x=gene_set_simple, y=n_perc, color=set_group)) +
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
  scale_color_manual(values = wes_palette("Zissou1",n=2,type = 'continuous')) +
  ylab('% of individuals\ncarrying RDGV') +
  xlab('') +
  facet_grid(cohort ~ variable)
dev.off()