---
  title: "Plot distribution of hits across sample sizes"
---
  
#rm(list=ls())
  
library(tidyverse)
library(dplyr)
library(wesanderson)
library(reshape2)


###run code inside /003_association_testing_results/ folder

##directories
output_file_direc='./results/'
output_figure_direc='./figures/'


##upload replicated hits
replicated_hits_FDR1 <- read.csv(file = paste(output_file_direc,'validated_genes.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) 
replicated_hits_FDR2 <- read.csv(file = paste(output_file_direc,'validated_genes_FDR2.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) 

##plot number of replicated hits by cancer type sample size
cancer_type_replicated_FDR1 <- replicated_hits_FDR1 %>%
  group_by(cancer_type) %>%
  summarise(n_replicated_FDR1=n()) %>%
  ungroup() 

cancer_type_replicated_FDR2 <- replicated_hits_FDR2 %>%
  group_by(cancer_type) %>%
  summarise(n_replicated_FDR2=n()) %>%
  ungroup() 

sample_sizes <- data.frame(discovery_TCGA=c(323,405,253,684,410,445,434,373,199,6799,386,403,363),
                           validation_PCAWG_Hartwig=c(87,283,283,656,417,168,299,299,180,4683,443,370,431),
                           cancer_type=c('Bladder','Brain_glioma_low','Brain_glioma_multi','Breast','Colon_Rectum',
                                         'Kidney','Lung_ad','Lung_sq','Ovary','Pancan','Prostate','Skin','Stomach_Eso'),
                           stringsAsFactors = F) 

sample_size_replication <- sample_sizes %>%
  left_join(cancer_type_replicated_FDR1, by=c('cancer_type')) %>%
  left_join(cancer_type_replicated_FDR2, by=c('cancer_type'))
sample_size_replication$n_replicated_FDR1[is.na(sample_size_replication$n_replicated_FDR1)] <- 0
sample_size_replication$n_replicated_FDR2[is.na(sample_size_replication$n_replicated_FDR2)] <- 0

sample_size_replication_plot <- sample_size_replication %>%
  melt(id.vars=c('cancer_type','discovery_TCGA','validation_PCAWG_Hartwig'))  %>%
  mutate(replicated=variable,n_replicated=value) %>%
  dplyr::select(-variable,-value) %>%
  melt(id.vars=c('cancer_type','replicated','n_replicated')) %>%
  dplyr::filter(cancer_type != 'Pancan' ) %>%
  mutate(replicated= sub('n_replicated_FDR1','replicated FDR 1%',replicated)) %>%
  mutate(replicated= sub('n_replicated_FDR2','replicated FDR 2%',replicated)) 
#%>%
 # dplyr::filter(replicated=='replicated FDR 2%' & variable=='validation_PCAWG_Hartwig')
#cor(sample_size_replication_plot$n_replicated,sample_size_replication_plot$value, method=c('pearson')) ##to calculate pearsons

##pearson correlations
sample_size_replication_pearson <- data.frame(replicated=c("replicated FDR 1%","replicated FDR 2%","replicated FDR 1%","replicated FDR 2%"),
                                              variable=c("discovery_TCGA","discovery_TCGA","validation_PCAWG_Hartwig","validation_PCAWG_Hartwig"),
                                              pearson=c(0.7283584,0.7589889,0.7097435,0.7846405),
                                              stringsAsFactors = F) %>%
  mutate(pearson=round(pearson,2))

pdf(file= paste(output_figure_direc,'Replication_vs_sample_size','.pdf',sep=''),
    width= 8,
    height = 3.2)  
ggplot(sample_size_replication_plot,
       aes(x=value,y=n_replicated)) +
  geom_point_rast(aes(fill=cancer_type), shape=21) +
  geom_smooth(method = "lm", se = T) +
  theme_bw() +
  theme(plot.subtitle = element_text(hjust = 0.5,size=8,color='black'),
        axis.text.x = element_text(size=8,color='black'),
        axis.text.y = element_text(size=8,color='black'),
        axis.title = element_text(size=8,color='black'),
        legend.position = "right") + 
  ylab('number of replicated hits') +
  xlab('sample sizes of cancer types') +
  scale_fill_brewer(palette = "Paired") +
  facet_grid(replicated ~ variable, scales = 'free') +
  geom_text(data = sample_size_replication_pearson, 
            aes(x = -Inf, y = +Inf, label = pearson), 
            hjust   = -0.2,
            vjust   = +2,
            color='red',
            size=3) +
  scale_x_continuous(limits = c(60,700)) +
  guides(colour = guide_legend(ncol = 2)) +
  labs(fill='')
dev.off()


