---
  title: "Main results - hits at a FDR of 1% - replicating Figure 3"
---

library(tidyverse)
library(dplyr)
library(ggsignif)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(grid)
library(tclust)
library(corrplot)
library(ggpmisc)
library(ggpubr)
library(cowplot)
library(ggrastr)
library(wesanderson)


##info for table with all info
output_file='.'
setwd('.')

##upload models to be discarded
models_discard <- read.csv(file = models_discard_input,head=TRUE,sep ='\t',stringsAsFactors=FALSE) 

##mapping features
mapping_features <- data.frame(pheno=c(paste('IC',seq(1,15,1),sep=''),paste('VAE_',seq(1,14,1),sep='')),
                               pheno_name=c('Sig.17',
                                            'Sig.MMR2+ampli.',
                                            'dMMR_ICA',
                                            'dHR_ICA',
                                            'Deletions_ICA',
                                            'APOBEC_ICA',
                                            'Sig.18',
                                            'Ploidy',
                                            'Sig.11+19',
                                            'DBS2',
                                            'Sig.5_ICA',
                                            'Smoking_ICA',
                                            'Small indels 2bp',
                                            'UV_ICA',
                                            'Sig.8',
                                            'APOBEC_VAE2',
                                            'Deletions_VAE',
                                            'Sig.5_VAE',
                                            'Sig.1',
                                            'UV_VAE',
                                            'dHR_VAE1',
                                            'Mitochondria',
                                            'dHR_VAE2',
                                            'Smoking_VAE',
                                            'dMMR_VAE2',
                                            'APOBEC_VAE1',
                                            'X-hypermutation',
                                            'dMMR_VAE1',
                                            'Amplifications'),
                               stringsAsFactors = F)

##TCGA burden testing results to get estimates
TCGA_burden_results <- c()
for(i in list.files(pattern = '^TCGA.*burden.*results_withgenelist.txt')){
  print(i)
  
  TCGA_burden <- read.csv(file =i,head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
    dplyr::select(variant,pValue,estimate,pheno,cancer_type,snps,model)
  
  TCGA_burden_results <- TCGA_burden_results %>%
    rbind(TCGA_burden)
}
colnames(TCGA_burden_results)[2] <- 'burden_test_pValue'


##WGS burden testing results to get estimates
WGS_burden_results <- c()
for(i in list.files(pattern = '^WGS.*burden.*results_withgenelist.txt')){
  print(i)
  
  WGS_burden <- read.csv(file =i,head=TRUE,sep ='\t',stringsAsFactors=FALSE)  %>%
    dplyr::select(variant,pValue,estimate,pheno,cancer_type,snps,model)
  
  WGS_burden_results <- WGS_burden_results %>%
    rbind(WGS_burden)
}
colnames(WGS_burden_results)[2] <- 'burden_test_pValue'


##TCGA SKAT testing results to get estimates
TCGA_SKAT_results <- c()
for(i in list.files(pattern = '^TCGA.*SKAT.*results_withgenelist.txt')){
  print(i)
  
  TCGA_SKAT <- read.csv(file =i,head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
    dplyr::select(variant,pValue,pheno,n_carriers,n_variants,n_samples,cancer_type,snps,model,test,rho_SKAT_O,n_genes)
  
  TCGA_SKAT_results <- TCGA_SKAT_results %>%
    rbind(TCGA_SKAT)
}

##WGS SKAT testing results to get estimates
WGS_SKAT_results <- c()
for(i in list.files(pattern = '^WGS.*SKAT.*results_withgenelist.txt')){
  print(i)
  
  WGS_SKAT <- read.csv(file =i,head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
    dplyr::select(variant,pValue,pheno,n_carriers,n_variants,n_samples,cancer_type,snps,model,test,rho_SKAT_O,n_genes)
  
  WGS_SKAT_results <- WGS_SKAT_results %>%
    rbind(WGS_SKAT)
}

##upload FDR estimates
TCGA_FDR <- read.csv(file = './TCGA_empirical_FDRs.txt',head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  dplyr::filter(FDR==0.01)
WGS_FDR <- read.csv(file = './WGS_empirical_FDRs.txt',head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  dplyr::filter(FDR==0.01)


##TCGA combine everything
TCGA_results <- TCGA_SKAT_results %>%
  left_join(TCGA_burden_results, by=c('variant','pheno','cancer_type','snps','model')) %>%
  dplyr::filter(!is.na(estimate)) %>%
  mutate(snps= sub('MRT','MTR',snps)) %>%
  mutate(snps=sub('_0.1perc','',snps)) %>%
  mutate(snps=sub('MTR_25thCentile','Misssense_MTR',snps)) %>%
  mutate(snps=sub('CCR_90orhigher','Misssense_CCR',snps)) %>%
  mutate(cancer_type= sub('BRCA','Breast',cancer_type)) %>%
  mutate(cancer_type= sub('COAD_READ','Colon_Rectum',cancer_type)) %>%
  mutate(cancer_type= sub('GBM','Brain_glioma_multi',cancer_type)) %>%
  mutate(cancer_type= sub('LGG','Brain_glioma_low',cancer_type)) %>%
  mutate(cancer_type= sub('LUAD','Lung_ad',cancer_type)) %>%
  mutate(cancer_type= sub('LUSC','Lung_sq',cancer_type)) %>%
  mutate(cancer_type= sub('PRAD','Prostate',cancer_type)) %>%
  mutate(cancer_type= sub('SKCM','Skin',cancer_type)) %>%
  mutate(cancer_type= sub('STAD_ESCA_EAC','Stomach_Eso',cancer_type)) %>%
  mutate(cancer_type= sub('OV','Ovary',cancer_type)) %>%
  mutate(cancer_type= sub('BLCA','Bladder',cancer_type)) %>%
  mutate(cancer_type= sub('KIRC_KIRP','Kidney',cancer_type)) %>%
  dplyr::filter(cancer_type %in% c('Breast','Colon_Rectum','Brain_glioma_multi','Brain_glioma_low',
                                   'Lung_ad','Lung_sq','Pancan','Prostate','Skin','Stomach_Eso',
                                   'Ovary','Bladder','Kidney')) %>%
  left_join(TCGA_FDR, by=c('cancer_type')) %>%
  mutate(cohort='TCGA')  %>%
  mutate(model_snps_cancer_pheno=paste(model,snps,cancer_type,pheno,sep='__')) %>%
  dplyr::filter(!model_snps_cancer_pheno %in% models_discard$model_snps_cancer_pheno) 


##WGS combine everything
WGS_results <- WGS_SKAT_results %>%
  left_join(WGS_burden_results, by=c('variant','pheno','cancer_type','snps','model')) %>%
  dplyr::filter(!is.na(estimate)) %>%
  mutate(snps= sub('MRT','MTR',snps)) %>%
  mutate(snps=sub('_0.1perc','',snps)) %>%
  mutate(snps=sub('MTR_25thCentile','Misssense_MTR',snps)) %>%
  mutate(snps=sub('CCR_90orhigher','Misssense_CCR',snps)) %>%
  mutate(cancer_type= sub('BRCA','Breast',cancer_type)) %>%
  mutate(cancer_type= sub('COAD_READ','Colon_Rectum',cancer_type)) %>%
  mutate(cancer_type= sub('GBM','Brain_glioma_multi',cancer_type)) %>%
  mutate(cancer_type= sub('LGG','Brain_glioma_low',cancer_type)) %>%
  mutate(cancer_type= sub('LUAD','Lung_ad',cancer_type)) %>%
  mutate(cancer_type= sub('LUSC','Lung_sq',cancer_type)) %>%
  mutate(cancer_type= sub('PRAD','Prostate',cancer_type)) %>%
  mutate(cancer_type= sub('SKCM','Skin',cancer_type)) %>%
  mutate(cancer_type= sub('STAD_ESCA_EAC','Stomach_Eso',cancer_type)) %>%
  mutate(cancer_type= sub('OV','Ovary',cancer_type)) %>%
  mutate(cancer_type= sub('BLCA','Bladder',cancer_type)) %>%
  mutate(cancer_type= sub('KIRC_KIRP','Kidney',cancer_type)) %>%
  dplyr::filter(cancer_type %in% c('Breast','Colon_Rectum','Brain_glioma_multi','Brain_glioma_low',
                                   'Lung_ad','Lung_sq','Pancan','Prostate','Skin','Stomach_Eso',
                                   'Ovary','Bladder','Kidney')) %>%
  left_join(WGS_FDR, by=c('cancer_type')) %>%
  mutate(cohort='WGS') %>%
  mutate(model_snps_cancer_pheno=paste(model,snps,cancer_type,pheno,sep='__')) %>%
  dplyr::filter(!model_snps_cancer_pheno %in% models_discard$model_snps_cancer_pheno)


###identify hits
TCGA_WGS_results <- TCGA_results %>%
  rbind(WGS_results) 

validated_FDR_cohort <-
  TCGA_WGS_results %>%
  dplyr::filter((cohort == 'TCGA' & pValue < FDR_threshold) | (cohort == 'WGS' & pValue < FDR_threshold )) %>% 
  dplyr::select(variant,pheno,cancer_type,model,snps,estimate,cohort) %>%
  mutate(variant_pheno_cancer_type_snp_model=paste(variant,pheno,cancer_type,snps,model,sep='__')) %>%
  mutate(variant_pheno_cancer_type=paste(variant,pheno,cancer_type,sep='__')) %>%
  mutate(variant_cancer_type= paste(variant,cancer_type,sep='__')) %>%
  spread(cohort,estimate) %>%
  dplyr::filter(sign(TCGA) == sign(WGS))

#print number of individual genes
print(nrow(validated_FDR_cohort))
print(length(unique(validated_FDR_cohort$variant)))
print(paste(unique(validated_FDR_cohort$variant), collapse = ', '))


##plot validated hits FDR1%
replicated_hits_snps <-
  TCGA_WGS_results %>%
  mutate(variant_pheno_cancer_type_snp_model=paste(variant,pheno,cancer_type,snps,model,sep='__')) %>%
  dplyr::filter(variant_pheno_cancer_type_snp_model %in% validated_FDR_cohort$variant_pheno_cancer_type_snp_model) %>%
  mutate(gene_cancer_type=paste(variant,' in ',cancer_type,sep='')) %>%
  mutate(pheno=sub('depth1_','',pheno)) %>%
  left_join(mapping_features, by=c('pheno')) %>%
  group_by(gene_cancer_type,pheno_name) %>%
  dplyr::select(gene_cancer_type,pheno_name,snps) %>% 
  distinct() %>%
  mutate(dummy=1) %>%
  spread(snps,dummy,is.na(0)) %>%
  mutate(snp_set=paste(PTV,PTV_Misssense_CADD25,PTV_Misssense_CADD15,Misssense_MTR,sep='|')) %>%
  mutate(snp_set= sub(0,'o',snp_set)) %>%
  mutate(snp_set= sub(0,'o',snp_set)) %>%
  mutate(snp_set= sub(0,'o',snp_set)) %>%
  mutate(snp_set= sub(0,'o',snp_set)) %>%
  mutate(snp_set= sub(1,'x',snp_set)) %>%
  mutate(snp_set= sub(1,'x',snp_set)) %>%
  mutate(snp_set= sub(1,'x',snp_set)) %>%
  mutate(snp_set= sub(1,'x',snp_set)) %>%
  dplyr::select(gene_cancer_type,pheno_name,snp_set) 

replicated_hits_model <-
  TCGA_WGS_results %>%
  mutate(variant_pheno_cancer_type_snp_model=paste(variant,pheno,cancer_type,snps,model,sep='__')) %>%
  dplyr::filter(variant_pheno_cancer_type_snp_model %in% validated_FDR_cohort$variant_pheno_cancer_type_snp_model) %>%
  mutate(gene_cancer_type=paste(variant,' in ',cancer_type,sep='')) %>%
  mutate(pheno=sub('depth1_','',pheno)) %>%
  left_join(mapping_features, by=c('pheno')) %>%
  dplyr::select(gene_cancer_type,pheno_name,model) %>% 
  distinct() %>%
  mutate(dummy=1) %>%
  spread(model,dummy,is.na(0)) %>%
  mutate(model_set=paste(dominant,additive,recessive,sep='|')) %>%
  mutate(model_set= sub(0,'o',model_set)) %>%
  mutate(model_set= sub(0,'o',model_set)) %>%
  mutate(model_set= sub(0,'o',model_set)) %>%
  mutate(model_set= sub(1,'x',model_set)) %>%
  mutate(model_set= sub(1,'x',model_set)) %>%
  mutate(model_set= sub(1,'x',model_set)) %>%
  dplyr::select(gene_cancer_type,pheno_name,model_set) %>%
  mutate(model_translation=(if_else(model_set == 'x|o|o', 'dominant',
                                    if_else(model_set == 'o|x|o', 'additive',
                                            if_else(model_set == 'o|o|x', 'recessive',
                                                    if_else(model_set == 'x|x|o', 'dominant & additive',
                                                            if_else(model_set == 'o|x|x', 'additive & recessive',
                                                                    if_else(model_set == 'x|x|x', 'dominant & additive & recessive','none'))))))))

replicated_hits <-
  TCGA_WGS_results %>%
  mutate(variant_pheno_cancer_type_snp_model=paste(variant,pheno,cancer_type,snps,model,sep='__')) %>%
  dplyr::filter(variant_pheno_cancer_type_snp_model %in% validated_FDR_cohort$variant_pheno_cancer_type_snp_model) %>%
  mutate(gene_cancer_type=paste(variant,' in ',cancer_type,sep='')) %>%
  mutate(pheno=sub('depth1_','',pheno)) %>%
  left_join(mapping_features, by=c('pheno')) %>%
  mutate(estimate_direction = if_else(estimate >0, 'positive', 'negative')) %>%
  dplyr::select(gene_cancer_type,pheno_name,estimate_direction) %>% 
  distinct() %>%
  left_join(replicated_hits_snps, by=c('gene_cancer_type','pheno_name')) %>%
  left_join(replicated_hits_model, by=c('gene_cancer_type','pheno_name')) %>%
  mutate(model_translation=factor(model_translation, levels=c('dominant','dominant & additive',
                                                              'additive','additive & recessive',
                                                              'recessive','dominant & additive & recessive')))  %>%
  mutate(components='Somatic mutational components \n') %>%
  mutate(estimate_direction=as.factor(estimate_direction))


ggplot(replicated_hits,
       aes(x=pheno_name,y=gene_cancer_type)) +
  geom_tile(aes(fill=model_translation, color=estimate_direction), size=1) +
  geom_text(aes(label = snp_set),color = "black", size=1.8) +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, color='black', hjust=0, angle=70),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=10, color='black'),
        legend.text  = element_text(size=8, color='black'),
        legend.title  = element_text(size=8, color='black'),
        legend.position = 'bottom') +
  xlab('') +
  ylab('') +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(factor(replicated_hits$gene_cancer_type)))) +
  labs(fill='')  +
  guides(fill = guide_legend(ncol = 2)) +
  labs(color='estimate size direction')  +
  guides(color = guide_legend(ncol = 1)) +
  scale_fill_manual(values=wes_palette("Zissou1", 6,type = "continuous")) +
  scale_color_manual(values=c('#C93312','#899DA4')) +
  facet_grid(. ~ components )

##add number of carriers info
replicated_hits_carriers <-
  TCGA_WGS_results %>%
  mutate(variant_pheno_cancer_type_snp_model=paste(variant,pheno,cancer_type,snps,model,sep='__')) %>%
  dplyr::filter(variant_pheno_cancer_type_snp_model %in% validated_FDR_cohort$variant_pheno_cancer_type_snp_model) %>%
  mutate(gene_cancer_type=paste(variant,' in ',cancer_type,sep='')) %>%
  mutate(pheno=sub('depth1_','',pheno)) %>%
  left_join(mapping_features, by=c('pheno')) %>%
  group_by(gene_cancer_type,snps,model,n_carriers,cohort) %>%
  summarise() %>%
  ungroup() %>%
  mutate(carrier_type= if_else(model== 'additive' | model== 'dominant', 'RDGVs',
                               if_else(model=='recessive', 'RDGVs + LOH','none' ))) %>%
  mutate(cohort_snps= paste(cohort,snps,sep='__')) %>%
  distinct() %>%
  mutate(cohort=sub('TCGA','Discovery',cohort)) %>%
  mutate(cohort=sub('WGS','Validation',cohort))


ggplot(replicated_hits_carriers,
       aes(x=snps,y=gene_cancer_type)) +
  geom_tile(aes(fill=n_carriers)) +
  geom_text(aes(label = n_carriers),color = "black", size=2) +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, color='black', hjust=0, angle=70),
        strip.text = element_text(size=8, color='black'),
        axis.title = element_text(size=8, color='black'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.text  = element_text(size=8, color='black'),
        legend.title  = element_text(size=8, color='black'),
        legend.position = 'bottom') +
  xlab('') +
  ylab('') +
  scale_fill_gradientn(colours = c('grey','#F21A00')) + 
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(factor(replicated_hits$gene_cancer_type)))) +
  facet_grid(. ~ carrier_type + cohort) +
  labs(fill='# of individuals')



##rho distribution in hits between different variant sets
##TCGA hits
tcga_hits_view <- TCGA_results %>%
  dplyr::filter(pValue < FDR_threshold) %>%
  dplyr::select(variant,pheno,cancer_type,model,snps,pValue,estimate,n_carriers,rho_SKAT_O,burden_test_pValue)
colnames(tcga_hits_view)[6:10] <- paste('TCGA',colnames(tcga_hits_view)[6:10],sep='_')

##wgs all 
WGS_hits_view <- WGS_results %>%
  dplyr::select(variant,pheno,cancer_type,model,snps,pValue,estimate,n_carriers,rho_SKAT_O,burden_test_pValue,FDR_threshold)
colnames(WGS_hits_view)[6:10] <- paste('WGS',colnames(WGS_hits_view)[6:10],sep='_')

##combined TCGA+WGS
combined_TCGA_WGS <- tcga_hits_view %>%
  left_join(WGS_hits_view, by=c('variant','pheno','cancer_type','model','snps')) %>%
  mutate(TCGA_log10= -log10(TCGA_pValue)) %>%
  mutate(WGS_log10= -log10(WGS_pValue)) 

SKAT_rho_distribution_hits <- combined_TCGA_WGS %>%
  dplyr::filter(sign(TCGA_estimate) == sign(WGS_estimate))   %>%
  dplyr::filter(WGS_pValue < FDR_threshold) %>%
  dplyr::select(variant,pheno,cancer_type,model,snps,TCGA_rho_SKAT_O,WGS_rho_SKAT_O) %>%
  melt(id.vars=c('variant','pheno','cancer_type','model','snps')) %>%
  mutate(variable=sub('TCGA_rho_SKAT_O','Discovery',variable))  %>%
  mutate(variable=sub('WGS_rho_SKAT_O','Validation',variable)) %>%
  mutate(variable= factor(variable, levels =c('Discovery','Validation')))  %>%
  mutate(snps= factor(snps, levels =c('PTV','PTV_Misssense_CADD25',
                                      'PTV_Misssense_CADD15','Misssense_MTR'))) %>%
  mutate(variable=sub('Discovery','Replicated Hits\nin discovery',variable)) %>%
  mutate(variable=sub('Validation','Replicated Hits\nin validation',variable)) 

ggplot(SKAT_rho_distribution_hits,
       aes(x=value, fill=variable)) +
  geom_density(alpha=0.4)+  
  theme_bw() +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=8, color='black'),
        strip.text.x = element_text(size=8, color='black'),
        strip.text.y = element_text(size=8, color='black'),
        legend.title = element_blank(),
        legend.position = 'bottom') + #legend.position = 'none'
  scale_fill_manual(values=c("Replicated Hits\nin discovery" = "#899DA4","Replicated Hits\nin validation"="#C93312")) +
  ylab('density') +
  xlab(expression(rho))


##rho distribution different variants sets
ggplot(SKAT_rho_distribution_hits,
       aes(x=snps, y=value, fill=snps)) +
  geom_violin()+  
  geom_jitter(size=0.3,shape=16, height = 0, width=0.1) +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, color='black', hjust=1, angle=70),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=8, color='black'),
        strip.text.x = element_text(size=8, color='black'),
        strip.text.y = element_text(size=8, color='black'),
        legend.title = element_blank(),
        legend.position = 'none') + #legend.position = 'none'
  scale_fill_manual(values=wes_palette("Zissou1", type = "continuous")) +
  ylab(expression(rho)) +
  xlab('') +
  facet_grid(variable ~ model) +
  scale_y_continuous(limits = c(0,1))



