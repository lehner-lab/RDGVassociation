---
  title: "Plot for the 15 extracted ICs correlation/contribution with input features"
---
  
#rm(list=ls())
  
library(tidyverse)
library(dplyr)
library(ggplot2);
library(wesanderson);
library(corrplot)


##directories
input_file_direc='./input_files/'
output_file_direc='./results/'
output_figure_direc='./figures/'

##mapping features
mapping_features <- data.frame(variable=c(paste('IC',seq(1,15,1),sep=''),paste('VAE_',seq(1,14,1),sep='')),
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

##ica with contributions and correlations already calculated
contri_cor_selected <-  read.csv(file= paste(input_file_direc,'TCGA_Hartwig_PCAWG_least45of56_allforGWAS_14832samples_ICs_ContrCorr.txt',sep=''),
                                 head=T,sep ="\t", stringsAsFactors = F) %>%
  mutate(pheno= sub('transcription_bias','trx',pheno)) %>%
  mutate(pheno= sub('small_indel_','indel_',pheno)) %>%
  mutate(pheno= sub('deletions_10bp','deletions>=10bp',pheno)) %>%
  mutate(pheno= sub('MSI','MS',pheno)) %>%
  mutate(pheno= sub('mh_2bp_more','mh>=2bp',pheno)) %>%
  left_join(mapping_features, by=c('variable')) %>%
  mutate(pheno_compo_name=paste(variable,': ',pheno_name,sep=''))


##selection of components to plot
plotOut <- NULL
n=1
for(PC in c('IC1','IC2','IC3','IC4','IC5','IC6','IC7','IC8',
            'IC9','IC10','IC11','IC12','IC13','IC14','IC15')){
  print(PC)
  
  rm(contri_cor_selected_order)
  contri_cor_selected_order <-
    contri_cor_selected %>%
    dplyr::filter(variable == PC) %>%
    mutate(pheno=factor(pheno,ordered = F,levels = contri_cor_selected$pheno[order(contri_cor_selected$value[contri_cor_selected$variable==PC & contri_cor_selected$info== 'contribution'], decreasing = T)])) %>%
    dplyr::filter(pheno %in% contri_cor_selected$pheno[order(contri_cor_selected$value[contri_cor_selected$variable==PC & contri_cor_selected$info== 'contribution'], decreasing = T)][1:10]) #only the 10 with the highest rankings
  
  bla <-
    ggplot(contri_cor_selected_order,
           aes(x=pheno,y=value,fill=info)) +
    geom_bar(stat="identity",color='black') +
    geom_hline(yintercept = 0,linetype="dashed", color = "black") +
    theme_bw() +
    theme(axis.text.x = element_text(size=8,angle = 70,hjust=1, color = 'black'),
          axis.text.y = element_text(size=8, color = 'black'),
          axis.title = element_blank(),
          strip.text.y = element_text(size=7, color = 'black'),
          legend.position = 'none') +
    xlab('') +
    ylab('') +
    facet_grid(info ~ pheno_compo_name, scales = 'free') +
    scale_fill_manual(values=c("contribution" = "#bdbdbd", "correlation" = "white"))
  
  plotOut[[paste('bla',n,sep='_')]] <- bla
  rm(bla)
  n= n+1
}


##plot other variables
pdf(file= paste(output_figure_direc,'ICA_overview','.pdf',sep=''),
    width= 8.2,
    height = 10.8)   
annotate_figure(
  ggarrange(plotOut$bla_1,plotOut$bla_2,plotOut$bla_3,
          plotOut$bla_4,plotOut$bla_5,plotOut$bla_6,
          plotOut$bla_7,plotOut$bla_8,plotOut$bla_9,
          plotOut$bla_10,plotOut$bla_11,plotOut$bla_12,
          plotOut$bla_13,plotOut$bla_14,plotOut$bla_15,
          nrow=4,ncol=4,common.legend = F,align='hv'),
  bottom = text_grob("somatic components", color = "black")
  )
dev.off()

