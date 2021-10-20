---
  title: "Main results - number of hits - number of replicated hits - etc - empirical FDR of 2%"
---
  
#rm(list=ls())
  
library(tidyverse)
library(dplyr)
library(wesanderson)
library(reshape2)


###run code inside /003_association_testing_results/ folder

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
#empirical FDR of 2%; can be changed here if needed
TCGA_FDR <- read.csv(file = paste(output_file_direc,'TCGA_empirical_FDRs.txt',sep=''), head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  dplyr::filter(FDR==0.02)
WGS_FDR <- read.csv(file = paste(output_file_direc,'PCAWG_Hartwig_empirical_FDRs.txt',sep=''), head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
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
print(nrow(validated_FDR_cohort))
print(length(unique(validated_FDR_cohort$gene)))
print(paste(unique(validated_FDR_cohort$gene), collapse = ', '))
 
#write output
write.table(validated_FDR_cohort,
            file= paste(output_file_direc,'/validated_genes_FDR2.txt',sep=''),
            quote=FALSE, sep='\t',row.names=FALSE,col.names = T)

##count hits in discovery cohort TCGA-WES
##how many hits
tcga_hits <- TCGA_results %>%
  dplyr::filter(pValue < FDR_threshold) 

print(paste(nrow(tcga_hits), ' genes-phenotype pairs reached 2% FDR in TCGA',sep=''))
print(paste(length(unique(tcga_hits$gene)),' genes reached 2% FDR in TCGA',sep=''))
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
validation_2FDR <- validated_FDR_cohort %>%
  group_by(gene,pheno,cancer_type) %>%
  summarise(n=n()) %>%
  ungroup() 

###depmap ordered plot
validated_hits_2FDR <- WGS_results %>%
  dplyr::select(gene,pheno,pheno_name,cancer_type,model,snps,estimate,cohort) %>%
  mutate(gene_pheno_cancer_type_snp_model=paste(gene,pheno,cancer_type,snps,model,sep='__')) %>%
  dplyr::filter(gene_pheno_cancer_type_snp_model %in% validated_FDR_cohort$gene_pheno_cancer_type_snp_model) %>%
  group_by(gene,pheno,pheno_name,cancer_type) %>%
  summarise(mean_estimate=mean(estimate)) %>%
  ungroup() %>%
  mutate(gene_cancer_type=paste(gene,' in ',cancer_type,sep='')) %>%
  left_join(validation_2FDR, by=c('gene','pheno','cancer_type')) %>%
  mutate(pheno_name=factor(pheno_name, levels=c('dHR_VAE1','dHR_VAE2','dHR_ICA','Sig.MMR2+ampli.',
                                                'dMMR_ICA','dMMR_VAE1','dMMR_VAE2','Small indels 2bp',
                                                'Deletions_VAE','Deletions_ICA','APOBEC_VAE1','APOBEC_ICA','APOBEC_VAE2',
                                                'DBS2', 'Smoking_ICA','Smoking_VAE',
                                                'Sig.11+19','Sig.17','UV_ICA',
                                                'Ploidy','X-hypermutation','Mitochondria',
                                                'Sig.1','Sig.5_VAE'))) %>%
  mutate(gene_cancer_type=factor(gene_cancer_type, levels=c('NCAPD2 in Colon_Rectum','RBBP6 in Brain_glioma_multi','RBBP5 in Pancan','KANSL3 in Skin',
                                                            'KANSL1 in Lung_ad','SETD1A in Breast','MTOR in Prostate','MTOR in Stomach_Eso',
                                                            'ATR in Pancan','ATR in Skin','KMT5A','KMT5A in Pancan','RBBP8 in Pancan','TTI2 in Prostate',
                                                            'TELO2 in Kidney','MAD2L2 in Pancan','NCAPG2 in Pancan','DNMT1 in Pancan','VPS72 in Prostate',
                                                            'ELP2 in Breast','POT1 in Pancan','BRCA2 in Breast','BRCA2 in Pancan','PALB2 in Pancan',
                                                            'PALB2 in Breast','MUS81 in Pancan','FANCM in Brain_glioma_multi','APC in Pancan','APC in Breast',
                                                            'KMT2B in Lung_ad','RIF1 in Skin','PIAS1 in Pancan','SOS1 in Pancan','KMT2D in Stomach_Eso',
                                                            'HERC2 in Ovary','HERC2 in Prostate','SETD2 in Colon_Rectum','REV3L in Pancan','REV3L in Kidney',
                                                            'NFRKB in Breast','ASCC3 in Stomach_Eso','BRCA1 in Lung_ad','BRCA1 in Colon_Rectum','BRCA1 in Breast',
                                                            'BRCA1 in Pancan','PAXIP1 in Skin','SUPT20H in Pancan','EP300 in Breast','EP300 in Skin',
                                                            'COL7A1 in Pancan','PRMT7 in Brain_glioma_low','KMT2E in Skin','PHF8 in Pancan','EXO1 in Skin',
                                                            'MSH2 in Pancan','SETX in Pancan','FANCC in Pancan','KDM6B in Kidney','CHD3 in Bladder','CHD3 in Pancan',
                                                            'ANKRD28 in Kidney','RECQL in Pancan','NUDT7 in Pancan','AXIN2 in Pancan','MLH1 in Colon_Rectum',
                                                            'MLH1 in Pancan','PIF1 in Pancan','SMC1B in Pancan','APEX1 in Stomach_Eso','APEX1 in Pancan',
                                                            'MLH3 in Pancan','MSH3 in Pancan','MSL2 in Lung_sq','KDM1A in Pancan','TOP3B in Pancan',
                                                            'TP53BP1 in Skin','WRN in Lung_ad','WRN in Prostate','KDM6A in Pancan','JADE2 in Colon_Rectum',
                                                            'SMARCAL1 in Brain_glioma_multi','ASCC2 in Pancan','DOCK8 in Breast','HMG20B in Pancan',
                                                            'DIS3L2 in Ovary','PARP3 in Pancan','CCNA1 in Pancan','EYA2 in Prostate','ZRANB3 in Pancan',
                                                            'PIK3C2B in Lung_sq','PADI4 in Skin','PER1 in Breast','PER1 in Pancan','SMC2 in Skin',
                                                            'SMC1A in Pancan','ACTL6A in Pancan','TRRAP in Ovary','TRRAP in Breast','MDN1 in Colon_Rectum',
                                                            'MDN1 in Brain_glioma_multi','MDN1 in Pancan','RAD51 in Pancan','TIMELESS in Brain_glioma_multi',
                                                            'AQR in Pancan','TOP2A in Pancan','PRPF19 in Pancan')))

pdf(file= paste(output_figure_direc,'Tile_validated_hit_FDR2_depmap','.pdf',sep=''),
    width= 11.5,
    height = 7)  
ggplot(validated_hits_2FDR,
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


##plot hits identified across several cancer types
replicated_hits_cancer_types_list <- validated_FDR_cohort %>%
  dplyr::filter(cancer_type != 'Pancan') %>%
  group_by(gene,cancer_type) %>%
  summarise() %>%
  ungroup() %>%
  group_by(gene) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  dplyr::filter(n>1)

replicated_hits_cancer_types <- validated_FDR_cohort %>%
  dplyr::filter(cancer_type != 'Pancan') %>%
  dplyr::filter(gene %in% replicated_hits_cancer_types_list$gene) %>%
  group_by(gene,cancer_type,pheno_name) %>%
  summarise() %>%
  ungroup() %>%
  mutate(gene=factor(gene, levels=c('TRRAP','BRCA1','MDN1','EP300','MTOR','HERC2','WRN'))) %>%
  mutate(cancer_type= factor(cancer_type, levels =c('Stomach_Eso','Prostate','Lung_ad','Brain_glioma_multi','Colon_Rectum',
                                                    'Ovary','Skin','APEX1 in Stomach_Eso','Breast'))) %>%
  mutate(gene_cancer_type=paste(gene,' in ',cancer_type,sep='')) 


pdf(file= paste(output_figure_direc,'Hits_FDR2_cancer_type_distribution','.pdf',sep=''),
    width= 8,
    height = 3.5)  
ggplot(replicated_hits_cancer_types,
       aes(x=gene_cancer_type,y=pheno_name, fill=gene)) +
  geom_tile() +
  theme_bw() +
  scale_fill_manual(values=wes_palette("Zissou1", 7,type = "continuous")) +
  #scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(size=8, color='black', hjust=1, angle=70),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=10, color='black'),
        #legend.title = element_blank(),
        legend.position = 'right') +
  xlab('replicated gene-tumor pairs') +
  ylab('somatic component') +
  labs(fill='')
dev.off()







##plot number of replicated hits by cancer type sample size
cancer_type_replicated_FDR1 <- validated_FDR_cohort %>%
  group_by(cancer_type) %>%
  summarise(n_replicated_FDR1=n()) %>%
  ungroup() 

cancer_type_replicated_FDR2 <- FDR2_hits %>%
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
#cor(sample_size_replication_plot$n_replicated,sample_size_replication_plot$value, method=c('pearson'))

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


