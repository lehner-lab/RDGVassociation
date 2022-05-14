---
  title: "Plot overview of all extracted components"
---
  
#rm(list=ls())
  
library(tidyverse)
library(dplyr)
library(ggplot2);
library(wesanderson);
library(gplots);
library(RColorBrewer)
library(ggrastr)
library(ggpubr)

##directories
input_file_direc='./input_files/'
output_file_direc='./results/'
output_figure_direc='./figures/'

###inputs
IC_features_input=''
VAE_features_input=''

##upload the IC features
IC_features <-  read.csv(file= IC_features_input, head=T,sep ="\t", stringsAsFactors = F) %>%
  dplyr::filter(info=='correlation')

##upload the VAE features
VAE_features <-  read.csv(file= VAE_features_input,head=T,sep ="\t", stringsAsFactors = F) 

##combine VAE+IC stuff
IC_VAE_features <- IC_features %>%
  rbind(VAE_features) %>%
  mutate(pheno= sub('transcription_bias','trx',pheno)) %>%
  mutate(pheno= sub('small_indel_','indel_',pheno)) %>%
  mutate(pheno= sub('deletions_10bp','deletions>=10bp',pheno)) %>%
  mutate(pheno= sub('MSI','MS',pheno)) %>%
  mutate(pheno= sub('mh_2bp_more','mh>=2bp',pheno))  %>%
  mutate(variable= sub('depth1_','',variable)) %>%
  mutate(variable= factor(variable, levels =c('VAE_12','VAE_5','IC14','VAE_13','IC3',
                                              'VAE_2','IC5','IC2','IC15','VAE_6','IC8',
                                              'VAE_14', 'IC12','VAE_9','IC1', 'VAE_8',
                                              'IC4','IC9','IC13','VAE_3','IC11','VAE_4',
                                              'VAE_7','IC7','IC10','VAE_10','VAE_11',
                                              'IC6','VAE_1'))) %>%
  mutate(pheno= factor(pheno, levels =c('mh>=2bp','mh_1bp','deletions>=10bp','indel_6to10bp','indel_2to5bp_nonMS',
                                        'indel_1bp_nonMS','amp_100to1000kb','amp_bigger1000kb','amp_10to100kb',
                                        'amp_1to10kb','ploidy','WGD','del_1to10kb','del10to100kb','del_bigger100kb',
                                        'DBS2','Ref.Sig.4','ID3', 'Ref.Sig.22','mt_total','ID4','Ref.Sig.11','Ref.Sig.19',
                                        'Ref.Sig.MMR2','Fork','trx_TtoA','trx_CtoA','trx_TtoC','Xhyper','indel_1bp_MS',
                                        'ID2','indel_2to5bp_MS','Ref.Sig.MMR1','Ref.Sig.1','DBS9','DNase','CTCF',
                                        'trx_CtoG','trx_TtoG','trx_CtoT','Ref.Sig.33','Ref.Sig.18','DBS4','Ref.Sig.17',
                                        'Ref.Sig.5','Ref.Sig.30','Ref.Sig.3','Ref.Sig.8','Ref.Sig.2','Ref.Sig.13',
                                        'H3K36me3','Expression','RT','ID8','DBS1','Ref.Sig.7'))) 


ggplot(IC_VAE_features,
       aes(x=variable,y=pheno)) +
  geom_tile(aes(fill=value)) +
  scale_fill_gradientn(colours = c('#3B9AB2','white','#F21A00')) + 
  theme_bw() +
  theme(axis.text.x = element_text(size=8, color='black', angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=10, color='black'),
        strip.text.x = element_text(size=10, color='black'),
        strip.text.y = element_text(size=10, color='black'),
        legend.position = 'bottom') +
  xlab('') +
  ylab('') +
  labs(fill='pearson correlation')


##heatmap in order to order accordingly
IC_VAE_features_heatmap <- IC_VAE_features %>%
  dplyr::select(variable,pheno,value) %>%
  spread(variable,value,is.na(0))

IC_VAE_features_heatmap_matrix <- as.matrix(IC_VAE_features_heatmap[,2:ncol(IC_VAE_features_heatmap)])
row.names(IC_VAE_features_heatmap_matrix) <- IC_VAE_features_heatmap$pheno

heatmap.2(IC_VAE_features_heatmap_matrix,
          trace='none',
          mar=c(11,5),
          col = wes_palette("Zissou1", 10, type = "continuous"),
          cexRow = 0.5)




