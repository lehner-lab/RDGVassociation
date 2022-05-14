---
  title: "Extract independent components from somatic mutational data"
---
  
#rm(list=ls())
  
library(tidyverse)
library(dplyr)
library(fastICA)
library(cluster)
library(ggplot2);
library(wesanderson);
library(RColorBrewer)
library(ggrastr)
library(ggpubr)

##directories
input_file_direc='./input_files/'
output_file_direc='./results/'
output_figure_direc='./figures/'

####upload combined phenotypes output
#sample x phenotype matrix
#scaled, samples with many NAs removed, NAs with median replaced

Combined_phenotypes_cohorts <-  read.csv(file= paste(input_file_direc,"/TCGA_Hartwig_PCAWG_samples_vs_phenotypes_withSigs_least45of56_allforGWAS_14832samples_table.txt",sep=''),
                                         head=T,sep ="\t")
input_name='TCGA_Hartwig_PCAWG_least45of56_allforGWAS_14832samples'


################### perform fastICA on matrix ####################  
#ICA changes with every run
#run it n times
#run it with different component extractions
#run different component extractions with different n cluster extractions

ICA_n=200; #number of repeating runs
component_max=30; #maximum extracting components
clustering_max=50; #maximum extracting clusters with k-medoids
total_fields=(clustering_max-1)*(component_max-1); #total number of runs

##save output
rm(ICA_silhouette_results)
ICA_silhouette_results = data.frame(component=rep(NA,total_fields),
                                    k=rep(NA,total_fields),
                                    min_silh_cluster=rep(NA,total_fields),
                                    second_low_silh_cluster=rep(NA,total_fields),
                                    average_silh_cluster=rep(NA,total_fields), 
                                    average_silh=rep(NA,total_fields), 
                                    stringsAsFactors = F);

n=1;                                
for(n_components in 2:component_max){
  print(paste('starting with extracting ',n_components,' components',sep=''))
  
  ica_results <- c() #results for each k run
  ica_score_results <- c() #results scores for each component extraction run
  
  for(i in 1:ICA_n){
    #set random number generator
    set.seed(sample(1:1000000,1))
    
    rm(ica_on_phenotypes)
    ica_on_phenotypes <- fastICA(Combined_phenotypes_cohorts[,2:dim(Combined_phenotypes_cohorts)[2]], n.comp=n_components)
    
    rm(ica_loadings_matrix)
    ica_loadings_matrix <- ica_on_phenotypes$A
    rownames(ica_loadings_matrix) <- paste('ICA',as.character(seq(1,nrow(ica_loadings_matrix),1)),as.character(i),sep='_')
    
    ica_results <- ica_results %>%
      rbind(ica_loadings_matrix)
  
    ##save score matrix as well
    rm(ica_score_matrix)
    ica_score_matrix <- ica_on_phenotypes$S
    colnames(ica_score_matrix) <- paste('ICA',as.character(seq(1,nrow(ica_loadings_matrix),1)),as.character(i),sep='_')
    
    ica_score_results <- ica_score_results %>%
      cbind(ica_score_matrix)
    
  }#finished running all randomizations for one component extraction
  
  write.table(ica_score_results,
              file=paste(output_file_direc,"output_ICA/fastICA",'_',input_name,'_',as.character(n_components),'components_','score_matrix_',as.character(ICA_n),'nonresampled_replicates.txt',sep=''),
              col.names = T, row.names = F, sep='\t', quote = F);
  
  ###do the k-means clustering based on medoids with pam, k as number of components
  for(i_k in 2:clustering_max){
    
    #calculate optimal number of clusters
    rm(clustering_output)
    clustering_output <- pam(ica_results, i_k, metric = "euclidean", stand = FALSE)
    
    ICA_silhouette_results$component[n] <- n_components
    ICA_silhouette_results$k[n] <- i_k
    ICA_silhouette_results$min_silh_cluster[n] <- min(clustering_output$silinfo$clus.avg.widths)
    ICA_silhouette_results$second_low_silh_cluster[n] <- sort(clustering_output$silinfo$clus.avg.widths)[2]
    ICA_silhouette_results$average_silh_cluster[n] <- mean(clustering_output$silinfo$clus.avg.widths)
    ICA_silhouette_results$average_silh[n] <- clustering_output$silinfo$avg.width
    
    
    ###saving the output
    save(clustering_output, 
         file=paste(output_file_direc,"output_ICA/fastICA",'_',input_name,'_',as.character(n_components),'components_',as.character(i_k),'kmedoid_clusters_',as.character(ICA_n),'nonresampled_replicates.RData',sep='')
    )
    
    ###saving silhouette coefficients output
    cluster_medoid_numbers <- data.frame(cluster_medoid=row.names(clustering_output$medoids),stringsAsFactors = F) %>%
      mutate(cluster=as.character(row_number()))
    
    sh_cluster <- data.frame(clustering_output$silinfo$widths) %>%
      mutate(cluster=as.character(cluster)) %>%
      left_join(cluster_medoid_numbers,by=c('cluster'))

    write.table(sh_cluster,
                file=paste(output_file_direc,"output_ICA/fastICA",'_',input_name,'_',as.character(n_components),'components_',as.character(i_k),'kmedoid_clusters_',as.character(ICA_n),'nonresampled_replicates_silhouette.txt',sep=''),
                col.names = T, row.names = F, sep='\t', quote = F);
    
    n= n+1
  } #finish run for each clustering
}; #finish run for each component


####matrix min_silh_cluster
pdf(file=paste(output_figure_direc,"fastICA",'_',input_name,'_',as.character(ICA_n),'nonresampled_replicates_min_silhouette_vs_k.pdf',sep=''),
    width= 8,
    height = 7)    
ggplot(ICA_silhouette_results, aes(component, k)) +
  geom_tile(aes(fill = min_silh_cluster)) +
  geom_text(aes(label=round(min_silh_cluster,2)),
            size=2.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=10,  color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black'),
        strip.text.x = element_text(size=10, color='black'),
        strip.text.y = element_text(size=10, color='black'),
        legend.title = element_text(size=10, color='black'),
        legend.position = 'bottom') +
  scale_x_discrete(breaks=seq(2,30,1), limits=factor(seq(1,30+1,1))) +
  scale_y_discrete(breaks=seq(2,50,1), limits=factor(seq(1,50+1,1))) +
  scale_fill_gradientn(colours = wes_palette("Zissou1")) +
  xlab('independent component') +
  ylab('extracted number of clusters') +
  labs(fill='minimum silhouette index cluster')
dev.off();

write.table(ICA_silhouette_results,
            file=paste(output_file_direc,"output_ICA/fastICA",'_',input_name,'_ICA_run_silhouette_results.txt',sep=''),
            col.names = T, row.names = F, sep='\t', quote = F)


##plot only the k=2*n_ICA
sh_ICA_plot <- ICA_silhouette_results %>%
  dplyr::filter(k == 2*component) %>% 
  dplyr::select(-average_silh) %>%
  melt(id.vars=c('k','component'))  %>%
  mutate(variable= sub('min_silh_cluster','minimum\nsilh_cluster',variable))  %>%
  mutate(variable= sub('second_low_silh_cluster','second minimum\nsilh_cluster',variable))  %>%
  mutate(variable= sub('average_silh_cluster','average\nsilh_cluster',variable))  %>%
  mutate(variable= factor(variable, levels =c('average\nsilh_cluster','minimum\nsilh_cluster','second minimum\nsilh_cluster')))


pdf(file= paste(output_figure_direc,'fastICA_silhouette_average_k_2n','.pdf',sep=''),
    width= 8,
    height = 4)  
ggplot(sh_ICA_plot, aes(component, value)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(size=10,  color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black'),
        strip.text.x = element_text(size=10, color='black'),
        strip.text.y = element_text(size=10, color='black'),
        legend.title = element_blank(),
        legend.position = 'none') +
  xlab('independent component')  +
  ylab('silhouette index')  +
  facet_grid(variable ~ .)
dev.off()

##plot the correlation between the ICs with 15 component extractions and 30-medoid clustering 

#upload the best run
load(paste(input_file_direc,'/fastICA_TCGA_Hartwig_PCAWG_least45of56_allforGWAS_14832samples_15components_30kmedoid_clusters_200nonresampled_replicates.RData',sep=''))

###taking clustered mediods
ica_loadings_matrix <- t(clustering_output$medoids)

###look at pearson between medoids to find mirrored ones
ica_loadings_pearson <- ica_loadings_matrix %>%
  cor(method = c("pearson"), use = "pairwise.complete.obs") #pairwise pearson wird berechnent, du kannst dir jeden command anscheun mit ?cor
colnames(ica_loadings_pearson) <- paste('IC_medoid_',seq(1,length(colnames(ica_loadings_pearson))),sep='')
rownames(ica_loadings_pearson) <- paste('IC_medoid_',seq(1,length(colnames(ica_loadings_pearson))),sep='')

colors_plot <- wes_palette("Zissou1", 10, type = "continuous")
res1 <- cor.mtest(ica_loadings_pearson, conf.level=.99) #significance test

pdf(file= paste(output_figure_direc,'fastICA_pearsonComponents','.pdf',sep=''),
    width= 9,
    height = 5) 
corrplot(ica_loadings_pearson,
         is.corr=T,
         tl.cex= .8,
         tl.col = "black",
         col=colors_plot,
         title='',
         #tl.pos = "lt",
         mar=c(0,0,2,0), #margins fürs pdf plotten
         method= "square",
         type = "upper",
         cl.cex = 1,
         p.mat = res1$p, sig.level = 0.05, pch.cex = .3,
         insig = "blank",
         na.label = "o"
)
dev.off() 


