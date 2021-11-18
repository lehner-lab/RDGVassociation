---
  title: "Network analysis - either with String or HumanNet - controlling overall for number of interactions each gene has"
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
output_figure_direc='./figures/'


##upload the validated genes
validated_genes <- read.csv(file = '../002_association_testing_results/results/validated_genes.txt',head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  dplyr::select(gene,pheno,cancer_type,model,snps) %>%
   dplyr::filter(!gene %in% c('MSH2','MLH1','BRCA2','BRCA1','PALB2'))
length(unique(validated_genes$gene))

validated_genes_FDR2 <- read.csv(file = '../002_association_testing_results/results/validated_genes_FDR2.txt', head=TRUE,sep ='\t',stringsAsFactors=FALSE) %>%
  dplyr::select(gene,pheno,cancer_type,model,snps)  %>%
   dplyr::filter(!gene %in% c('MSH2','MLH1','BRCA2','BRCA1','PALB2'))
length(unique(validated_genes_FDR2$gene))


###loop through different nets
for(i in list.files(path =input_file_direc,pattern = 'interactions.txt')){
  print(i)
  
  network_name=paste(sub('_converted_genelist_interactions.txt','',i),'_removed_known_genes',sep='')

  #upload file
  rm(gene_interactions)
  gene_interactions <- read.csv(file = paste(input_file_direc,i,sep=''), head=TRUE,sep ='\t',stringsAsFactors=FALSE)  %>%
    group_by(ID = paste0(pmax(symbol1, symbol2), pmin(symbol1, symbol2))) %>%
    slice(1) %>%
    ungroup() ###remove duplicated pair scenarios
  
  ############### part 1 total number of interactions between FDR1 genes and randomization check  ############### 
  
  ##counting the total number of interactions of each gene within gene list
  rm(total_interactions)
  total_interactions <- data.frame(interaction= c(gene_interactions$symbol1,gene_interactions$symbol2),stringsAsFactors=FALSE) %>%
    group_by(interaction) %>%
    summarise(total_n_interactions=n()) %>%
    ungroup()

  ##counting the total number of interactions from validated genes with each other
  rm(interactions_validated)
  interactions_validated <- gene_interactions %>%
    dplyr::filter(symbol1 %in% validated_genes$gene & symbol2 %in% validated_genes$gene)
  print(paste(nrow(interactions_validated), ' interaction between validated genes in total',sep=''))
  
  ##counting the total number of interactions of each gene in validated set
  rm(total_interactions_validated)
  total_interactions_validated <- total_interactions %>%
    dplyr::filter(interaction %in% validated_genes$gene)
  
  
  ##histogram of number of gene interactions per gene - all
  pdf(file= paste(output_figure_direc,network_name,'_hist_total_interactions_genes','.pdf',sep=''),
      width=3,
      height = 1.5)  
  print(ggplot(total_interactions,
         aes(x=total_n_interactions)) +
    geom_histogram(fill='#3B9AB2',alpha=1,binwidth = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(size=8, color='black'),
          axis.text.y = element_text(size=8, color='black'),
          axis.title = element_text(size=8, color='black'),
          legend.position = 'none') +
    xlab('Numer of interactions per gene in geneset') +
    ylab('Count'))
  dev.off()
  
  ##histogram of number of gene interactions per gene - validated genes
  pdf(file= paste(output_figure_direc,network_name,'_hist_total_interactions_genes_validated','.pdf',sep=''),
      width=3,
      height = 1.5)  
  print(ggplot(total_interactions_validated,
               aes(x=total_n_interactions)) +
          geom_histogram(fill='#3B9AB2',alpha=1,binwidth = 1) +
          geom_vline(xintercept=quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[2:10], col = 'red',size=0.4, linetype='dashed') +
          theme_bw() +
          theme(axis.text.x = element_text(size=8, color='black'),
                axis.text.y = element_text(size=8, color='black'),
                axis.title = element_text(size=8, color='black'),
                legend.position = 'none') +
          xlab('Numer of interactions per gene in geneset \nFDR 1% validated genes') +
          ylab('Count'))
  dev.off()
  
  ##split number of interactions in validated geneset into 10 equal sized bins
  rm(total_interactions_validated_10bins)
  total_interactions_validated_10bins <- total_interactions_validated %>%
    mutate(bin= if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[2],'bin1',
                        if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[3],'bin2',
                                if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[4],'bin3',
                                        if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[5],'bin4',
                                                if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[6],'bin5',
                                                        if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[7],'bin6',
                                                                if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[8],'bin7',
                                                                        if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[9],'bin8',
                                                                                if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[10],'bin9','bin10'))))))))))

  
  
  ##now also split the total geneset into 10 based on the thresholds of above
  rm(total_interactions_10bins)
  total_interactions_10bins <- total_interactions %>%
    dplyr::filter(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[11]) %>%
    mutate(bin= if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[2],'bin1',
                        if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[3],'bin2',
                                if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[4],'bin3',
                                        if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[5],'bin4',
                                                if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[6],'bin5',
                                                        if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[7],'bin6',
                                                                if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[8],'bin7',
                                                                        if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[9],'bin8',
                                                                                if_else(total_n_interactions <= quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[10],'bin9','bin10'))))))))))

  
  ##randomization
  number_edges_random <- c()
  
  for(n_random in 1:1000){
    rm(random_gene_list_controlled_interactions)
    rm(random_interactions)
    rm(random_egdes)
    
    ##create a random list of genes, sample from each bin same number of genes
    set.seed(n_random)
    random_gene_list_controlled_interactions <- data.frame(
      genes=c(sample(total_interactions_10bins$interaction[total_interactions_10bins$bin=='bin1'],nrow(total_interactions_validated_10bins[total_interactions_validated_10bins$bin=='bin1',])),
              sample(total_interactions_10bins$interaction[total_interactions_10bins$bin=='bin2'],nrow(total_interactions_validated_10bins[total_interactions_validated_10bins$bin=='bin2',])),
              sample(total_interactions_10bins$interaction[total_interactions_10bins$bin=='bin3'],nrow(total_interactions_validated_10bins[total_interactions_validated_10bins$bin=='bin3',])),
              sample(total_interactions_10bins$interaction[total_interactions_10bins$bin=='bin4'],nrow(total_interactions_validated_10bins[total_interactions_validated_10bins$bin=='bin4',])),
              sample(total_interactions_10bins$interaction[total_interactions_10bins$bin=='bin5'],nrow(total_interactions_validated_10bins[total_interactions_validated_10bins$bin=='bin5',])),
              sample(total_interactions_10bins$interaction[total_interactions_10bins$bin=='bin6'],nrow(total_interactions_validated_10bins[total_interactions_validated_10bins$bin=='bin6',])),
              sample(total_interactions_10bins$interaction[total_interactions_10bins$bin=='bin7'],nrow(total_interactions_validated_10bins[total_interactions_validated_10bins$bin=='bin7',])),
              sample(total_interactions_10bins$interaction[total_interactions_10bins$bin=='bin8'],nrow(total_interactions_validated_10bins[total_interactions_validated_10bins$bin=='bin8',])),
              sample(total_interactions_10bins$interaction[total_interactions_10bins$bin=='bin9'],nrow(total_interactions_validated_10bins[total_interactions_validated_10bins$bin=='bin9',])),
              sample(total_interactions_10bins$interaction[total_interactions_10bins$bin=='bin10'],nrow(total_interactions_validated_10bins[total_interactions_validated_10bins$bin=='bin10',]))),
      stringsAsFactors = F
    )
    
    ## extract these genes and count number of interactions between them
    random_interactions <- gene_interactions %>%
      dplyr::filter(symbol1 %in% random_gene_list_controlled_interactions$genes & symbol2 %in% random_gene_list_controlled_interactions$genes)
    print(paste(nrow(random_interactions), ' interaction between random genes in total',sep=''))  
    
    ##histogram of number of gene interactions per gene - random set to check that everything looks fine
    pdf(file= paste(output_figure_direc,network_name,'_hist_total_interactions_genes_random','.pdf',sep=''),
        width=3,
        height = 1.5)  
    print(ggplot(total_interactions %>%
                   dplyr::filter(interaction %in% random_gene_list_controlled_interactions$genes),
                 aes(x=total_n_interactions)) +
            geom_histogram(fill='grey20',alpha=1,binwidth = 1) +
            geom_vline(xintercept=quantile(total_interactions_validated$total_n_interactions,seq(0,1,0.1))[2:10], col = 'red',size=0.4, linetype='dashed') +
            theme_bw() +
            theme(axis.text.x = element_text(size=8, color='black'),
                  axis.text.y = element_text(size=8, color='black'),
                  axis.title = element_text(size=8, color='black'),
                  legend.position = 'none') +
            xlab('Numer of interactions per random gene in geneset') +
            ylab('Count'))
    dev.off()
    
    ##count
    random_egdes <- data.frame(random=n_random, n_interactions=nrow(random_interactions), stringsAsFactors = F)
    
    ##combine all into one df
    number_edges_random <- number_edges_random %>%
      rbind(random_egdes)
  }#end of randomization for edges
  
  ##histo number of edges real vs random
  pdf(file= paste(output_figure_direc,network_name,'_hist_interactions_random_vs_real','.pdf',sep=''),
      width=3,
      height = 1.5)  
  print(ggplot(number_edges_random,
         aes(x=n_interactions)) +
    geom_histogram(color='black',fill='#3B9AB2',alpha=1,binwidth = 1) +
    geom_vline(xintercept=nrow(interactions_validated), col = 'red',size=1.5)+
    theme_bw() +
    theme(axis.text.x = element_text(size=8, color='black'),
          axis.text.y = element_text(size=8, color='black'),
          axis.title = element_text(size=8, color='black'),
          legend.position = 'none') +
    xlab('Number of interactions in random subset of the tested gene set') +
    ylab('Count'))
  dev.off()
  
  #write output
  write.table(number_edges_random %>%
                mutate(number_observed=nrow(interactions_validated)),
              file= paste(output_file_direc,network_name,'_interactions_random_vs_real.txt',sep=''),
              quote=FALSE, sep='\t',row.names=FALSE,col.names = T)
  
  
  ############### part 2 count number of FDR2 genes interacting with at least one FDR1 gene and randomization check ############### 
  validates_genesFDR2_noFDR1 <- validated_genes_FDR2[!validated_genes_FDR2$gene %in% validated_genes$gene,]
  genesFDR2_noFDR1 <- unique(validates_genesFDR2_noFDR1$gene)
  length(genesFDR2_noFDR1)
  
  ##get distribution of interactions of the FDR2 only genes (excluding the ones found at FDR1)
  rm(total_interactions_validated_FDR2only)
  total_interactions_validated_FDR2only <- total_interactions %>%
    dplyr::filter(interaction %in% genesFDR2_noFDR1)
  
  ##get distribution of interactions excluding the FDR1% genes from the geneist
  rm(total_interactions_noFDR1)
  total_interactions_noFDR1 <- total_interactions %>%
    dplyr::filter(!interaction %in% validated_genes$gene)
  
  ##count how many of them have a direct interaction with at least one FDR1 gene
  rm(validates_genesFDR1_FDR2_interaction)
  validates_genesFDR1_FDR2_interaction <- gene_interactions %>% 
    dplyr::filter((symbol1 %in% validates_genesFDR2_noFDR1$gene & symbol2 %in% validated_genes$gene)|
                    (symbol1 %in% validated_genes$gene & symbol2 %in% validates_genesFDR2_noFDR1$gene))  
  
  interactions_FDR2_with_FDR1 <- length(genesFDR2_noFDR1[genesFDR2_noFDR1 %in% unique(c(validates_genesFDR1_FDR2_interaction$symbol1,validates_genesFDR1_FDR2_interaction$symbol2))])
  print(paste(interactions_FDR2_with_FDR1,' FDR2 genes (out of ',nrow(total_interactions_validated_FDR2only), ') directly interacting with a FDR1 gene',sep=''))
  
  ##histogram of number of gene interactions per gene - all without FDR1% genes
  pdf(file= paste(output_figure_direc,network_name,'_hist_total_interactions_genes_noFDR1','.pdf',sep=''),
      width=3,
      height = 1.5)  
  print(ggplot(total_interactions_noFDR1,
               aes(x=total_n_interactions)) +
          geom_histogram(fill='#3B9AB2',alpha=1,binwidth = 1) +
          theme_bw() +
          theme(axis.text.x = element_text(size=8, color='black'),
                axis.text.y = element_text(size=8, color='black'),
                axis.title = element_text(size=8, color='black'),
                legend.position = 'none') +
          xlab('Numer of interactions per gene in geneset') +
          ylab('Count'))
  dev.off()
  
  ##histogram of number of gene interactions per gene - validated genes
  pdf(file= paste(output_figure_direc,network_name,'_hist_total_interactions_genes_validated_FDR2only','.pdf',sep=''),
      width=3,
      height = 1.5)   
  print(ggplot(total_interactions_validated_FDR2only,
               aes(x=total_n_interactions)) +
          geom_histogram(fill='#3B9AB2',alpha=1,binwidth = 1) +
          geom_vline(xintercept=quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[2:10], col = 'red',size=0.4, linetype='dashed') +
          theme_bw() +
          theme(axis.text.x = element_text(size=8, color='black'),
                axis.text.y = element_text(size=8, color='black'),
                axis.title = element_text(size=8, color='black'),
                legend.position = 'none') +
          xlab('Numer of interactions per gene in geneset \nFDR 2% validated genes without FDR 1% genes') +
          ylab('Count'))
  dev.off()
  
  ##split number of interactions in validated geneset into 10 equal sized bins, based on FDR2% only genes
  rm(total_interactions_validated_10bins_FDR2only)
  total_interactions_validated_10bins_FDR2only <- total_interactions_validated_FDR2only %>%
    mutate(bin= if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[2],'bin1',
                        if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[3],'bin2',
                                if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[4],'bin3',
                                        if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[5],'bin4',
                                                if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[6],'bin5',
                                                        if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[7],'bin6',
                                                                if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[8],'bin7',
                                                                        if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[9],'bin8',
                                                                                if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[10],'bin9','bin10'))))))))))
 
  
  
  ##now also split the total geneset into 10 based on the thresholds of above - without FDR1 genes
  rm(total_interactions_noFDR1_10bins)
  total_interactions_noFDR1_10bins <- total_interactions_noFDR1 %>%
    dplyr::filter(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[11]) %>%
    mutate(bin= if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[2],'bin1',
                        if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[3],'bin2',
                                if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[4],'bin3',
                                        if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[5],'bin4',
                                                if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[6],'bin5',
                                                        if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[7],'bin6',
                                                                if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[8],'bin7',
                                                                        if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[9],'bin8',
                                                                                if_else(total_n_interactions <= quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[10],'bin9','bin10'))))))))))

  
  ##randomization
  number_interactions_with_FDR1_gene <- c()
  
  for(n_random in 1:1000){
    rm(random_gene_list_controlled_interactions)
    rm(FDR1_interaction_random)
    rm(FDR1_interaction_random_count)
    rm(FDR1_interaction_random_count_df)
    
    ##create a random list of genes
    set.seed(n_random)
    random_gene_list_controlled_interactions <- data.frame(
      genes=c(sample(total_interactions_noFDR1_10bins$interaction[total_interactions_noFDR1_10bins$bin=='bin1'],nrow(total_interactions_validated_10bins_FDR2only[total_interactions_validated_10bins_FDR2only$bin=='bin1',])),
              sample(total_interactions_noFDR1_10bins$interaction[total_interactions_noFDR1_10bins$bin=='bin2'],nrow(total_interactions_validated_10bins_FDR2only[total_interactions_validated_10bins_FDR2only$bin=='bin2',])),
              sample(total_interactions_noFDR1_10bins$interaction[total_interactions_noFDR1_10bins$bin=='bin3'],nrow(total_interactions_validated_10bins_FDR2only[total_interactions_validated_10bins_FDR2only$bin=='bin3',])),
              sample(total_interactions_noFDR1_10bins$interaction[total_interactions_noFDR1_10bins$bin=='bin4'],nrow(total_interactions_validated_10bins_FDR2only[total_interactions_validated_10bins_FDR2only$bin=='bin4',])),
              sample(total_interactions_noFDR1_10bins$interaction[total_interactions_noFDR1_10bins$bin=='bin5'],nrow(total_interactions_validated_10bins_FDR2only[total_interactions_validated_10bins_FDR2only$bin=='bin5',])),
              sample(total_interactions_noFDR1_10bins$interaction[total_interactions_noFDR1_10bins$bin=='bin6'],nrow(total_interactions_validated_10bins_FDR2only[total_interactions_validated_10bins_FDR2only$bin=='bin6',])),
              sample(total_interactions_noFDR1_10bins$interaction[total_interactions_noFDR1_10bins$bin=='bin7'],nrow(total_interactions_validated_10bins_FDR2only[total_interactions_validated_10bins_FDR2only$bin=='bin7',])),
              sample(total_interactions_noFDR1_10bins$interaction[total_interactions_noFDR1_10bins$bin=='bin8'],nrow(total_interactions_validated_10bins_FDR2only[total_interactions_validated_10bins_FDR2only$bin=='bin8',])),
              sample(total_interactions_noFDR1_10bins$interaction[total_interactions_noFDR1_10bins$bin=='bin9'],nrow(total_interactions_validated_10bins_FDR2only[total_interactions_validated_10bins_FDR2only$bin=='bin9',])),
              sample(total_interactions_noFDR1_10bins$interaction[total_interactions_noFDR1_10bins$bin=='bin10'],nrow(total_interactions_validated_10bins_FDR2only[total_interactions_validated_10bins_FDR2only$bin=='bin10',]))),
      stringsAsFactors = F
    )

    
    ##histogram of number of gene interactions per gene  - random set to check that everything looks fine
    pdf(file= paste(output_figure_direc,network_name,'_hist_total_interactions_genes_random_FDR2only','.pdf',sep=''),
        width=3,
        height = 1.5)  
    print(ggplot(total_interactions_noFDR1 %>%
                   dplyr::filter(interaction %in% random_gene_list_controlled_interactions$genes),
                 aes(x=total_n_interactions)) +
            geom_histogram(fill='grey20',alpha=1,binwidth = 1) +
            geom_vline(xintercept=quantile(total_interactions_validated_FDR2only$total_n_interactions,seq(0,1,0.1))[2:10], col = 'red',size=0.4, linetype='dashed') +
            theme_bw() +
            theme(axis.text.x = element_text(size=8, color='black'),
                  axis.text.y = element_text(size=8, color='black'),
                  axis.title = element_text(size=8, color='black'),
                  legend.position = 'none') +
            xlab('Numer of interactions per random gene in geneset') +
            ylab('Count'))
    dev.off()
    
    ##count for these genes how many have at least one interaction with a FDR1 gene
    FDR1_interaction_random <- gene_interactions %>% 
      dplyr::filter((symbol1 %in% random_gene_list_controlled_interactions$genes & symbol2 %in% validated_genes$gene)|
                      (symbol1 %in% validated_genes$gene & symbol2 %in% random_gene_list_controlled_interactions$genes))  
    
    FDR1_interaction_random_count <- length(random_gene_list_controlled_interactions$genes[random_gene_list_controlled_interactions$genes %in% unique(c(FDR1_interaction_random$symbol1,FDR1_interaction_random$symbol2))])
    FDR1_interaction_random_count_df <- data.frame(random=n_random, n_interactions=FDR1_interaction_random_count, stringsAsFactors = F)
    
    print(paste(FDR1_interaction_random_count, ' genes (out of ',nrow(total_interactions_validated_FDR2only),') interacting with at least one FDR1 gene',sep=''))  
      
    ##combine all
    number_interactions_with_FDR1_gene <- number_interactions_with_FDR1_gene %>%
      rbind(FDR1_interaction_random_count_df)
  }#end of randomization for edges
  
  
  ##histo number of random set of genes interacting with at least one FDR1 gene 
  pdf(file= paste(output_figure_direc,network_name,'_hist_interaction_FDR2_with_FDR1_random_vs_real','.pdf',sep=''),
      width=3,
      height = 1.5)   
  print(ggplot(number_interactions_with_FDR1_gene,
               aes(x=n_interactions)) +
          geom_histogram(color='black',fill='#F21A00',alpha=1,binwidth = 1) +
          geom_vline(xintercept=interactions_FDR2_with_FDR1, col = 'red',size=1.5)+
          theme_bw() +
          theme(axis.text.x = element_text(size=8, color='black'),
                axis.text.y = element_text(size=8, color='black'),
                axis.title = element_text(size=8, color='black'),
                legend.position = 'none') +
          xlab(paste('Number of random genes (out of ',nrow(total_interactions_validated_FDR2only),') interacting directly \n with at least one gene validated at FDR 1%',sep='')) +
          ylab('Count'))
  dev.off()
  
  
  #write output
  write.table(number_interactions_with_FDR1_gene %>%
                mutate(number_observed=interactions_FDR2_with_FDR1),
              file= paste(output_file_direc,network_name,'_interaction_FDR2_with_FDR1_random_vs_real.txt',sep=''),
              quote=FALSE, sep='\t',row.names=FALSE,col.names = T)

}#end looping through generated network files
