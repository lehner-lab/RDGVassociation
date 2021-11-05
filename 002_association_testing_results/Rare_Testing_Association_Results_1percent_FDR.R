---
  title: "Main results - number of hits - number of replicated hits - etc - empirical FDR of 1%"
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


##discard models based on lambda
models_discard <- read.csv(file = paste(output_file_direc,'models_removed_based_on_lambda.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) 

##upload results from burden
TCGA_burden_results <- read.csv(file = paste(input_file_direc,'Discovery_TCGA_burden_results.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE)  
colnames(TCGA_burden_results)[2] <- 'burden_test_pValue'
WGS_burden_results <- read.csv(file = paste(input_file_direc,'Validation_PCAWG_Hartwig_burden_results.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE)  
colnames(WGS_burden_results)[2] <- 'burden_test_pValue'

##upload results from SKAT
TCGA_SKAT_results <- read.csv(file = paste(input_file_direc,'Discovery_TCGA_SKATO_results.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) 
WGS_SKAT_results <- read.csv(file = paste(input_file_direc,'Validation_PCAWG_Hartwig_SKATO_results.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE) 

##upload empirical FDR estimates
#empirical FDR of 1%; can be changed here if needed
TCGA_FDR <- read.csv(file = paste(output_file_direc,'TCGA_empirical_FDRs.txt',sep=''), head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  dplyr::filter(FDR==0.01)
WGS_FDR <- read.csv(file = paste(output_file_direc,'PCAWG_Hartwig_empirical_FDRs.txt',sep=''), head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  dplyr::filter(FDR==0.01)

##combine discovery SKAT results with burden and remove models based on inflation factors
TCGA_results <- TCGA_SKAT_results %>%
  left_join(TCGA_burden_results, by=c('gene','pheno','cancer_type','snps','model','pheno_name','model_snps','model_snps_cancer_pheno')) %>%
  dplyr::filter(!is.na(estimate)) %>%
  left_join(TCGA_FDR, by=c('cancer_type')) %>%
  mutate(cohort='TCGA')  %>%
  mutate(model_snps_cancer_pheno=paste(model,snps,cancer_type,pheno,sep='__')) %>%
  dplyr::filter(!model_snps_cancer_pheno %in% models_discard$model_snps_cancer_pheno) 
#real_result <- nrow(TCGA_results) ##count number of tests
#real_result_comb <- TCGA_results %>% group_by(cancer_type,model) %>% summarise(n=n()) %>% ungroup() ##count number of tests across models and cancer types

##combine validations SKAT results with burden and remove models based on inflation factors
WGS_results <- WGS_SKAT_results %>%
  left_join(WGS_burden_results, by=c('gene','pheno','cancer_type','snps','model','pheno_name','model_snps','model_snps_cancer_pheno')) %>%
  dplyr::filter(!is.na(estimate)) %>%
  left_join(WGS_FDR, by=c('cancer_type')) %>%
  mutate(cohort='WGS') %>%
  mutate(model_snps_cancer_pheno=paste(model,snps,cancer_type,pheno,sep='__')) %>%
  dplyr::filter(!model_snps_cancer_pheno %in% models_discard$model_snps_cancer_pheno)
#real_result_WGS <- nrow(WGS_results) ##count number of tests
#real_result_comb_WGS <- WGS_results %>% group_by(cancer_type,model) %>% summarise(n=n()) %>% ungroup() ##count number of tests across models and cancer types

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
print(nrow(validated_FDR_cohort))
print(length(unique(validated_FDR_cohort$gene)))
print(paste(unique(validated_FDR_cohort$gene), collapse = ', '))
 
#write output
write.table(validated_FDR_cohort,
            file= paste(output_file_direc,'/validated_genes.txt',sep=''),
            quote=FALSE, sep='\t',row.names=FALSE,col.names = T)

##count hits in discovery cohort TCGA-WES
##how many hits
tcga_hits <- TCGA_results %>%
  dplyr::filter(pValue < FDR_threshold) 

print(paste(nrow(tcga_hits), ' genes-phenotype pairs reached 1% FDR in TCGA',sep=''))
print(paste(length(unique(tcga_hits$gene)),' genes reached 1% FDR in TCGA',sep=''))
print(unique(tcga_hits$gene))


##validation stats 
#how many hits got tested? how many did not get tested?
tcga_hits_validation <- TCGA_results %>%
  rbind(WGS_results) %>%  
  dplyr::filter((cohort == 'TCGA' & pValue < FDR_threshold) | cohort != 'TCGA') %>%
  dplyr::select(gene,pheno,cancer_type,model,snps,estimate,cohort) %>%
  spread(cohort,estimate) %>%
  dplyr::filter(!is.na(TCGA))
  
tcga_hits_validation_not_tested <- tcga_hits_validation %>%
  dplyr::filter(is.na(WGS))
print(paste(nrow(tcga_hits_validation_not_tested), ' hits in the discovery TCGA and were not tested in the validation PCAWG_Hartwig',sep=''))

tcga_hits_validation_tested <- tcga_hits_validation %>%
  dplyr::filter(!is.na(WGS)) 
#%>%
 # dplyr::filter(sign(TCGA) == sign(WGS)) #checking for effect size directions from burden test
print(paste(nrow(tcga_hits_validation_tested), ' hits in the discovery TCGA and were tested in the validation PCAWG_Hartwig',sep=''))



##count number replicated hits across gene-pheno-cancer_types tuples
validation_1FDR <- validated_FDR_cohort %>%
  group_by(gene,pheno,cancer_type) %>%
  summarise(n=n()) %>%
  ungroup() 

###depmap ordered plot
validated_hits_depmap <- WGS_results %>%
  dplyr::select(gene,pheno,pheno_name,cancer_type,model,snps,estimate,cohort) %>%
  mutate(gene_pheno_cancer_type_snp_model=paste(gene,pheno,cancer_type,snps,model,sep='__')) %>%
  dplyr::filter(gene_pheno_cancer_type_snp_model %in% validated_FDR_cohort$gene_pheno_cancer_type_snp_model) %>% #only validated FDR 1%
  group_by(gene,pheno,pheno_name,cancer_type) %>%
  summarise(mean_estimate=mean(estimate)) %>%
  ungroup() %>%
  mutate(gene_cancer_type=paste(gene,' in ',cancer_type,sep='')) %>%
  left_join(validation_1FDR, by=c('gene','pheno','cancer_type')) %>%
  mutate(pheno= factor(pheno, levels =c('VAE_8','VAE_6','IC4','IC3',
                                        'VAE_2','VAE_1','IC2','IC10',
                                        'IC12','IC9','IC1','IC14',
                                        'IC8', 'VAE_9','VAE_12','IC13',
                                        'VAE_10','VAE_13'))) %>%
  mutate(pheno_name=factor(pheno_name, levels=c('dHR_VAE1','dHR_VAE2','dHR_ICA','Sig.MMR2+ampli.',
                                                'dMMR_ICA','dMMR_VAE1','dMMR_VAE2','Small indels 2bp',
                                                'Deletions_VAE','APOBEC_VAE2','DBS2',
                                                'Smoking_ICA','Smoking_VAE',
                                                'Sig.11+19','Sig.17','UV_ICA',
                                                'Ploidy','X-hypermutation'))) %>%
  mutate(gene_cancer_type=factor(gene_cancer_type, levels=c('TRRAP in Ovary','MDN1 in Brain_glioma_multi','SMC2 in Skin','ACTL6A in Pancan',
                                                            'RAD51 in Pancan','TOP2A in Pancan','AQR in Pancan','NCAPD2 in Colon_Rectum',
                                                            'RBBP5 in Pancan','MTOR in Prostate','KANSL1 in Lung_ad','SETD1A in Breast',
                                                            'TTI2 in Prostate','DNMT1 in Pancan','MAD2L2 in Pancan','NCAPG2 in Pancan',
                                                            'BRCA2 in Pancan','BRCA2 in Breast','PALB2 in Pancan','BRCA1 in Pancan','BRCA1 in Breast',
                                                            'APC in Pancan','APC in Breast','REV3L in Pancan','HERC2 in Prostate','HERC2 in Ovary',
                                                            'NFRKB in Breast','POT1 in Pancan','RIF1 in Skin','MSH2 in Pancan','FANCC in Pancan',
                                                            'EXO1 in Skin','PAXIP1 in Skin','APEX1 in Pancan','MSH3 in Pancan','MLH1 in Pancan',
                                                            'CHD3 in Bladder','RECQL in Pancan','SMC1B in Pancan','ASCC2 in Pancan','HMG20B in Pancan',
                                                            'EYA2 in Prostate','PIK3C2B in Lung_sq','PER1 in Breast','TP53BP1 in Skin','EP300 in Skin'))) 


pdf(file= paste(output_figure_direc,'Tile_validated_hit_FDR1_depmap','.pdf',sep=''),
    width= 8,
    height = 5)  
ggplot(validated_hits_depmap,
       aes(x=gene_cancer_type,y=pheno_name)) +
  geom_tile(aes(fill=mean_estimate)) +
  geom_text(aes(label = n),color = "black") +
  #scale_fill_gradientn(colours = wes_palette("Zissou1", 10, type = "continuous")) + 
  scale_fill_gradient2(low ='#3B9AB2',
                       mid='white',
                       high='#F21A00',
                       midpoint = 0) + 
  theme_bw() +
  theme(axis.text.x = element_text(size=8, color='black', hjust=1, angle=70),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=10, color='black'),
        strip.text.x = element_text(size=10, color='black'),
        strip.text.y = element_text(size=10, color='black'),
        #legend.title = element_blank(),
        legend.position = 'bottom') +
  xlab('replicated gene-tumor pairs') +
  ylab('') +
  labs(fill='mean estimate burden test in validation cohort')
dev.off()


##rho distribution in hits between different SNP sets
##Discovery TCGA hits
tcga_hits_view <- TCGA_results %>%
  dplyr::filter(pValue < FDR_threshold) %>%
  dplyr::select(gene,pheno,cancer_type,model,snps,pValue,estimate,n_carriers,rho_SKAT_O,burden_test_pValue)
colnames(tcga_hits_view)[6:10] <- paste('TCGA',colnames(tcga_hits_view)[6:10],sep='_')

##Validation PCAWG_Hartwig all 
WGS_hits_view <- WGS_results %>%
  dplyr::select(gene,pheno,cancer_type,model,snps,pValue,estimate,n_carriers,rho_SKAT_O,burden_test_pValue,FDR_threshold)
colnames(WGS_hits_view)[6:10] <- paste('WGS',colnames(WGS_hits_view)[6:10],sep='_')

##combined TCGA+WGS
combined_TCGA_WGS <- tcga_hits_view %>%
  left_join(WGS_hits_view, by=c('gene','pheno','cancer_type','model','snps')) %>%
  mutate(TCGA_log10= -log10(TCGA_pValue)) %>%
  mutate(WGS_log10= -log10(WGS_pValue)) 

SKAT_rho_distribution_hits <- combined_TCGA_WGS %>%
  dplyr::filter(sign(TCGA_estimate) == sign(WGS_estimate))   %>%
  dplyr::filter(WGS_pValue < FDR_threshold) %>%
  dplyr::select(gene,pheno,cancer_type,model,snps,TCGA_rho_SKAT_O,WGS_rho_SKAT_O) %>%
  melt(id.vars=c('gene','pheno','cancer_type','model','snps')) %>%
  mutate(variable=sub('TCGA_rho_SKAT_O','Discovery',variable))  %>%
  mutate(variable=sub('WGS_rho_SKAT_O','Validation',variable)) %>%
  mutate(variable= factor(variable, levels =c('Discovery','Validation')))  %>%
  mutate(snps= factor(snps, levels =c('PTV','PTV_Misssense_CADD25',
                                      'PTV_Misssense_CADD15','Misssense_MTR'))) %>%
  mutate(variable=sub('Discovery','Replicated Hits\nin discovery',variable)) %>%
  mutate(variable=sub('Validation','Replicated Hits\nin validation',variable)) 


pdf(file= paste(output_figure_direc,'rho_distribution_hits','.pdf',sep=''),
    width= 3,
    height = 3)  
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
dev.off()


##rho distribution different SNP sets
pdf(file= paste(output_figure_direc,'rho_distribution_hits_violin','.pdf',sep=''),
    width= 4,
    height = 3.5)  
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
dev.off()


##plot number of replicated hits/number of retested hits across models 
validated_FDR_cohort_by_model <- validated_FDR_cohort %>% 
  group_by(model) %>%
  summarise(n_replicated=n()) %>%
  ungroup()

hits_retested <- tcga_hits_validation_tested %>% 
  group_by(model) %>%
  summarise(n_hits_retested=n()) %>%
  ungroup()

prop_model <- hits_retested %>%
  left_join(validated_FDR_cohort_by_model, by=c('model')) %>%
  mutate(prop=(n_replicated/n_hits_retested)*100)

pdf(file= paste(output_figure_direc,'distribution_hits_Proportion_model_replicated','.pdf',sep=''),
    width= 1.5,
    height = 2.5)  
ggplot(prop_model,
       aes(x=model,y=prop)) +
  geom_bar(position="dodge", stat="identity", fill='grey55') + #wes_palette("Zissou1", 10, type = "continuous")
  theme_bw() +
  theme(plot.subtitle = element_text(hjust = 0.5,size=8,color='black'),
        axis.text.x = element_text(size=8,color='black',angle=70,hjust=1),
        axis.text.y = element_text(size=8,color='black'),
        axis.title = element_text(size=8,color='black'),
        legend.position = "none") + 
  ylab('% of replicated associations \n out of re-tested hits') +
  xlab('Model of inheritance')
dev.off()

