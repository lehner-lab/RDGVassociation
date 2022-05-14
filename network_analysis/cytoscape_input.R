---
  title: "Input to visualize stuff with cytoscape"
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


####cytoscape input nodes
##upload validated genes FDR1
validated_genes_FDR1 <-  read.csv(file = '../002_association_testing_results/results/validated_genes.txt',head=TRUE,sep ='\t',stringsAsFactors=FALSE)  %>%
  dplyr::select(gene,pheno) %>%
  distinct() %>%
  mutate(gene_pheno= paste(gene,pheno,sep='_')) 

##upload validated genes FDR2
validated_genes_FDR2 <- read.csv(file = '../002_association_testing_results/results/validated_genes_FDR2.txt',head=TRUE,sep ='\t',stringsAsFactors=FALSE)  %>%
  dplyr::select(gene,pheno) %>%
  distinct() %>%
  mutate(gene_pheno= paste(gene,pheno,sep='_')) %>%
  dplyr::filter(!gene_pheno %in% validated_genes_FDR1$gene_pheno)

##FDR1 hits
validated_genes_plot_FDR1 <- validated_genes_FDR1 %>%
  left_join(mapping_features, by=c('pheno')) %>%
  group_by(gene) %>%
  mutate(fraction=1/length(gene)) %>%
  ungroup() %>%
  dplyr::select(gene,pheno_name,fraction) %>%
  spread(pheno_name,fraction,is.na(0)) %>%
  mutate(validation=if_else(gene %in% validated_genes_FDR1$gene, 'FDR1','FDR2'))

#write output
write.table(validated_genes_plot_FDR1,
            file= paste(output_file_direc,'network_validated_hits_FDR1.txt',sep=''),
            quote=FALSE, sep='\t',row.names=FALSE,col.names = T)


##combine hits from FDR1 and FDR2
validated_genes_plot_FDR1_FDR2 <- validated_genes_FDR1 %>%
  rbind(validated_genes_FDR2) %>%
  left_join(mapping_features, by=c('pheno')) %>%
  group_by(gene) %>%
  mutate(fraction=1/length(gene)) %>%
  ungroup() %>%
  dplyr::select(gene,pheno_name,fraction) %>%
  spread(pheno_name,fraction,is.na(0)) %>%
  mutate(validation=if_else(gene %in% validated_genes_FDR1$gene, 'FDR1','FDR2'))

#write output
write.table(validated_genes_plot_FDR1_FDR2,
            file= paste(output_file_direc,'network_validated_hits_FDR1_FDR2.txt',sep=''),
            quote=FALSE, sep='\t',row.names=FALSE,col.names = T)


####cytoscape input network
##upload STRING
i='STRING_combinedScore_atleast0.8_converted_genelist_interactions.txt'
network_name=sub('_converted_genelist_interactions.txt','',i)
  
#interactions file FDR1
gene_interactions_FDR1 <- read.csv(file= paste(input_file_direc,i,sep=''), head=TRUE,sep ='\t',stringsAsFactors=FALSE)%>%
  group_by(ID = paste0(pmax(symbol1, symbol2), pmin(symbol1, symbol2))) %>%
  slice(1) %>% ##remove duplicated pair scenarios
  ungroup() %>%
  mutate(gene=symbol1) %>%
  dplyr::filter(symbol1 %in% validated_genes_FDR1$gene & symbol2 %in% validated_genes_FDR1$gene) %>% #only interactions from the validated geneset
  mutate(pp='pp')

#write output
write.table(gene_interactions_FDR1,
            file= paste(output_file_direc,network_name,'_network_validated_hits_FDR1.txt',sep=''),
            quote=FALSE, sep='\t',row.names=FALSE,col.names = T)


#interactions file FDR2 and FDR1
gene_interactions_FDR2 <-  read.csv(file= paste(input_file_direc,i,sep=''), head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  group_by(ID = paste0(pmax(symbol1, symbol2), pmin(symbol1, symbol2))) %>%
  slice(1) %>% ###remove duplicated pair scenarios
  ungroup() %>%
  mutate(gene=symbol1) %>%
  dplyr::filter(symbol1 %in% validated_genes_plot_FDR1_FDR2$gene & symbol2 %in% validated_genes_plot_FDR1_FDR2$gene) %>% #only interactions from the validated geneset
  mutate(pp='pp')

#write output
write.table(gene_interactions_FDR2,
            file= paste(output_file_direc,network_name,'_network_validated_hits_FDR1_FDR2.txt',sep=''),
            quote=FALSE, sep='\t',row.names=FALSE,col.names = T)
  

  