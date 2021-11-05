---
  title: "Plot overlpa of hits at 1% FDR via venn diagramm"
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


##upload replicated hits at an empirical FDR of 1%
validated_genes <- read.csv(file = paste(output_file_direc,'validated_genes.txt',sep=''),head=TRUE,sep ='\t',stringsAsFactors=FALSE)  %>%
  dplyr::select(gene,pheno,cancer_type,model,snps)

# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}


###overlap hits between vae vs ica features
validated_genes_pheno_overlap <- validated_genes %>%
  mutate(pheno_type= if_else(grepl('VAE',pheno),'VAE','ICA')) %>%
  group_by(pheno_type,gene) %>%
  summarise() %>%
  ungroup() %>%
  mutate(dummy=1) %>%
  spread(pheno_type,dummy,is.na(0)) 

x_pheno <- list(
  ICA = validated_genes_pheno_overlap$gene[validated_genes_pheno_overlap$ICA==1],
  VAE = validated_genes_pheno_overlap$gene[validated_genes_pheno_overlap$VAE==1] 
)

pdf(file= paste(output_figure_direc,'venn_validated_hits_pheno_overlap','.pdf',sep=''), #_FDR2
    width= 1.8,
    height = 1.5)  
display_venn(x_pheno,
             fill = c("#3B9AB2","#F21A00"),
             cat.cex = .8,
             cex = .8,
             lwd=c(1,1))
dev.off()


###overlap hits between snp sets
validated_genes_snps_overlap <- validated_genes %>%
  mutate(variable_group=paste(gene,cancer_type,pheno,model,sep='_')) %>%
  group_by(snps,variable_group) %>%
  summarise() %>%
  ungroup() %>%
  mutate(snps=sub('_0.1perc','',snps)) %>%
  mutate(snps=sub('MTR_25thCentile','Misssense_MTR',snps)) %>%
  mutate(snps=sub('CCR_90orhigher','Misssense_CCR',snps)) %>%
  mutate(dummy=1) %>%
  spread(snps,dummy,is.na(0)) 

x_snps <- list(
  PTV_Misssense_CADD15 = validated_genes_snps_overlap$variable_group[validated_genes_snps_overlap$PTV_Misssense_CADD15==1] ,
  PTV_Misssense_CADD25 = validated_genes_snps_overlap$variable_group[validated_genes_snps_overlap$PTV_Misssense_CADD25==1] ,
  Misssense_MTR = validated_genes_snps_overlap$variable_group[validated_genes_snps_overlap$Misssense_MTR==1],
  PTV = validated_genes_snps_overlap$variable_group[validated_genes_snps_overlap$PTV==1]
)

pdf(file= paste(output_figure_direc,'venn_validated_hits_snps_overlap','.pdf',sep=''), #_FDR2
    width= 2.5,
    height = 2)   
display_venn(x_snps,
             fill = c("#3B9AB2","#9EBE91","#E4B80E","#F21A00"),
             cat.cex = .8,
             cex = .8,
             lwd=c(1,1,1,1))
dev.off()


###overlap hits between models 
validated_genes_models_overlap <- validated_genes %>%
  mutate(variable_group=paste(gene,cancer_type,pheno,snps,sep='_')) %>%
  group_by(model,variable_group) %>%
  summarise() %>%
  ungroup() %>%
  mutate(dummy=1) %>%
  spread(model,dummy,is.na(0)) 

x_models <- list(
  additive = validated_genes_models_overlap$variable_group[validated_genes_models_overlap$additive==1],
  dominant = validated_genes_models_overlap$variable_group[validated_genes_models_overlap$dominant==1],
  recessive = validated_genes_models_overlap$variable_group[validated_genes_models_overlap$recessive==1] 
)

pdf(file= paste(output_figure_direc,'venn_validated_hits_models_overlap','.pdf',sep=''), #_FDR2
    width= 2.5,
    height = 2)   
display_venn(x_models,
             fill = c("#3B9AB2","#EBCC2A","#F21A00"),
             cat.cex = .8,
             cex = .8,
             lwd=c(1,1,1))
dev.off()

