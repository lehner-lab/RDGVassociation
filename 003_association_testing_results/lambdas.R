---
  title: "Estimation of inflation factors"
---
  
#rm(list=ls())

library(tidyverse)
library(dplyr)
library(ggrastr)


##directories
input_file='./003_association_testing_results/input_files/Discovery_TCGA_SKATO_results.txt'
output_file_direc='./003_association_testing_results/results/'
output_figure_direc='./003_association_testing_results/figures/'

table_name='SKAT_models'


##upload SKAT-O results
SKAT_res <- read.csv(file =input_file,head=TRUE,sep ='\t',stringsAsFactors=FALSE) 

##calculate lambdas
lambdas_results <- c()

for(i in unique(SKAT_res$model_snps_cancer_pheno)){
  print(i)
  rm(Inflation_pvalues)
  rm(Inflation_chisq)
  rm(lambda)
  rm(assc_test)
  
  assc_test <- SKAT_res %>%
    dplyr::filter(model_snps_cancer_pheno == rlang::sym(i))
  
  ##estimate genomic inflation factor lambda 
  Inflation_pvalues <- assc_test$pValue
  Inflation_chisq <- qchisq(1-Inflation_pvalues,1)
  lambda = round(median(Inflation_chisq)/qchisq(0.5,1),2)
  
  ##save lambda values in table
  rm(lambda_values)
  lambda_values <- data.frame(model_snps_cancer_pheno=i,
                              model_snps=unique(assc_test$model_snps),
                              model=unique(assc_test$model),
                              snps=unique(assc_test$snps),
                              cancer_type=unique(assc_test$cancer_type),
                              pheno=unique(assc_test$pheno),
                              lambda_qq=lambda,
                              n_genes=nrow(assc_test),
                              n_samples=unique(assc_test$n_samples),
                              stringsAsFactors = F)
  
  lambdas_results <- lambdas_results %>%
    rbind(lambda_values)
}#end looping through cancer pheno pairs

##plotting lambdas
lambdas_results_filtered <- lambdas_results %>%
  mutate(snps= factor(snps, levels =c('PTV','PTV_Misssense_CADD25','PTV_Misssense_CADD15','Misssense_MTR','Misssense_CCR'))) %>%
  mutate(cancer_type=sub('Brain_glioma_multi','Brain\nglioma\nmulti',cancer_type)) %>%
  mutate(cancer_type=sub('Brain_glioma_low','Brain\nglioma\nlow',cancer_type))  %>%
  mutate(cancer_type=sub('Stomach_Eso','Stomach\nEso',cancer_type))  %>%
  mutate(cancer_type=sub('Colon_Rectum','Colon\nRectum',cancer_type)) %>%
  dplyr::filter(n_genes >= 100) %>% #only lambdas where at least 100 genes were tested
  mutate(removed=if_else(lambda_qq >= 1.5, 'discarded','kept')) #inflation factors bigger or equal to 1.5 should be marked
  
pdf(file= paste(output_figure_direc,table_name,'_lambdas_distribution_VAE_ICA','.pdf',sep=''),
    width= 8,
    height =10)  
ggplot(lambdas_results_filtered, 
       aes(x=snps,y=lambda_qq, label=pheno, fill=n_genes)) +
  geom_jitter_rast(alpha=0.4, 
                   size=1,
                   height=0,
                   width=0.2,
                   aes(x=snps,
                       y=lambda_qq,
                       color=removed),
                   lambdas_results_filtered) + 
  geom_boxplot(outlier.shape = NA,
               color='black') + 
  scale_fill_gradientn(colours = rainbow(5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=8, angle = 70,hjust=1, color='black'),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=10, color='black'),
        strip.text.x = element_text(size=10, color='black'),
        strip.text.y = element_text(size=8, color='black'),
        legend.position = "right",
        legend.text = element_text(size=8)) +
  facet_grid(cancer_type ~ model) +
  xlab('') +
  ylab('inflation factor') +
  labs(fill='n genes tested') +
  labs(color='discarded?') +
  scale_colour_manual(values=c("discarded" = "#F21A00","kept" = "#3B9AB2"))
dev.off()


####write out the ones which will be removed
models_removed <- lambdas_results_filtered %>%
  dplyr::filter(removed == 'discarded')

write.table(models_removed,
            file= paste(output_file_direc,'/models_removed_based_on_lambda.txt',sep=''),
            quote=FALSE, sep='\t',row.names=FALSE,col.names = T)




