---
  title: "Plot VAE results"
---
  
#rm(list=ls())
  
library(tidyverse)
library(dplyr)
library(ggplot2);
library(wesanderson);
library(corrplot)
library(RColorBrewer)
library(ggrastr)
library(ggpubr)

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

##plot parameter sweep results
vae_parameter_sweep <-  read.csv(file= paste(input_file_direc,'/VAE_TCGA_Hartwig_PCAWG_least45of56_allforGWAS_14832samples_parameters_test_summary.txt',sep=''),
                                 head=T,sep=' ', stringsAsFactors = F) 


##process to plot parameter sweep
vae_parameter_sweep_plot <- vae_parameter_sweep %>%
  mutate(ln_mean_sumActivity_mu_layer= log(mean_sum_activity_mean_encoded_layer + exp(-10))) %>% #0 will become -10
  mutate(ln_validation_loss= log(validation_loss)) %>%
  dplyr::select(ln_mean_sumActivity_mu_layer,ln_validation_loss,r_mean_reconstruction,mean_max_cor_golden_ICs,num_components,learning_rate,batch_size,epochs,kappa,depth) %>%
  reshape2::melt(id.vars=c('num_components','learning_rate','batch_size','epochs','kappa','depth')) %>%
  mutate(epochs= factor(epochs)) %>%
  mutate(batch_size= factor(batch_size)) %>%
  mutate(num_components= factor(num_components)) %>%
  mutate(kappa= factor(kappa)) %>%
  mutate(depth= factor(depth)) %>%
  mutate(variable = sub('mean_max_cor_golden_ICs','correlation\n with ICs',variable)) %>%
  mutate(variable = sub('r_mean_reconstruction','correlation\n with input',variable)) %>%
  mutate(variable=factor(variable,ordered = F,levels = c('correlation\n with ICs','ln_mean_sumActivity_mu_layer','ln_validation_loss','correlation\n with input'))) %>%
  dplyr::filter(variable == 'correlation\n with input' | variable == 'correlation\n with ICs')

##plot VAE parameter sweep
mycolors <- colorRampPalette(brewer.pal(5, "Set1"))(length(unique(vae_parameter_sweep_plot$batch_size)))

pdf(file= paste(output_figure_direc,'VAE_parameter_sweep','.pdf',sep=''),
    width= 8,
    height = 9)  
ggplot(data=vae_parameter_sweep_plot , 
       aes(x=epochs, y=value)) + 
  geom_point_rast(aes(shape=kappa,fill=batch_size), size=2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=8, angle = 70,hjust=1, color='black'),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=10, color='black'),
        strip.text.x = element_text(size=8, color='black'),
        strip.text.y = element_text(size=8, color='black')) +
  facet_grid(cohort ~ tissue) +
  scale_shape_manual(values=c(21, 22, 23, 24)) +
  guides(fill = guide_legend(override.aes=list(shape=c(21)))) +
  scale_fill_manual(values = mycolors) +
  facet_grid(variable + learning_rate + depth ~ num_components, scales='free') +
  ylab('Pearson correlation') +
  xlab('number of epochs')
dev.off()


##plot results from 5 random intis with optimized conditions
vae_5_random_inits <-  read.csv(file= paste(input_file_direc,'VAE_TCGA_Hartwig_PCAWG_least45of56_allforGWAS_14832samples_5randomInit_summary.txt',sep=''),
                                 head=T,sep=' ', stringsAsFactors = F) 

##process for plotting
vae_5_random_inits_plot <- vae_5_random_inits %>%
  dplyr::select(mean_max_cor_golden_ICs,num_components,learning_rate,batch_size,epochs,kappa,depth,seed) %>%
  reshape2::melt(id.vars=c('num_components','learning_rate','batch_size','epochs','kappa','depth','seed')) %>%
  mutate(variable = sub('mean_max_cor_golden_ICs','correlation\n with ICs',variable)) %>%
  mutate(num_components = factor(num_components))


##plot 5 inits
pdf(file= paste(output_figure_direc,'VAE_5_random_inits','.pdf',sep=''),
    width= 8,
    height = 2.5)  
ggplot(data=vae_5_random_inits_plot , 
       aes(x=num_components, y=value)) + 
  geom_boxplot(outlier.shape = NA,
               color='grey') + 
  geom_jitter_rast(size=0.5,shape=16, height = 0, width=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=8, angle = 70,hjust=1, color='black'),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=10, color='black'),
        strip.text.x = element_text(size=8, color='black'),
        strip.text.y = element_text(size=8, color='black')) +
  facet_grid(. ~ depth, scales='free') +
  ylab('average correlation with \nselected independent components') +
  xlab('number of components')
dev.off()


##upload VAEs for both depth with highest IC correlations
vae_depth1 <- read.csv(file= paste(input_file_direc,'VAE_mse_plus_KL_14_components_0.0005_lr_200_batchsize_200_epochs_0.5_kappa_1_beta_28_hiddenDim_1_depth_7869_seed_latent_mean_encoded.txt',sep=''),
                       head=T,sep='\t', stringsAsFactors = F) 
colnames(vae_depth1)[1:length(colnames(vae_depth1))-1] <- paste('VAE_depth1_',seq(1,length(colnames(vae_depth1))-1),sep='')


vae_depth2 <- read.csv(file= paste(input_file_direc,'VAE_mse_plus_KL_14_components_0.0005_lr_200_batchsize_200_epochs_0.5_kappa_1_beta_28_hiddenDim_2_depth_7304_seed_latent_mean_encoded.txt',sep=''),
                       head=T,sep='\t', stringsAsFactors = F) 
colnames(vae_depth2)[1:length(colnames(vae_depth2))-1] <- paste('VAE_depth2_',seq(1,length(colnames(vae_depth2))-1),sep='')

##also upload ICA resuts for comparison
ica_results <- read.csv(file= paste(input_file_direc,'TCGA_Hartwig_PCAWG_least45of56_allforGWAS_14832samples_IC_values_table.txt',sep=''),
                       head=T,sep='\t', stringsAsFactors = F) 

##correlation between VAEs
vae_depth <- vae_depth1 %>%
  left_join(vae_depth2, by = c('sample_short')) %>%
  dplyr::select(-sample_short) %>% 
  cor(method = c("pearson"), use = "pairwise.complete.obs")

colors_plot <- wes_palette("Zissou1", 10, type = "continuous")
res1 <- cor.mtest(vae_depth, conf.level=.99) 

pdf(file= paste(output_figure_direc,'VAE_component_correlation','.pdf',sep=''),
    width= 9,
    height = 6) 
corrplot(vae_depth,
         is.corr=T,
         tl.cex= .8,
         tl.col = "black",
         col=colors_plot,
         title='',
         mar=c(0,0,2,0), #margins fürs pdf plotten
         method= "square",
         type = "upper",
         cl.cex = 1,
         p.mat = res1$p, sig.level = 0.05, pch.cex = .3,
         insig = "blank",
         na.label = "o"
)
dev.off() 


##correlation depth 1 with ICA
vae_depth1_cor <- vae_depth1 %>%
  left_join(ica_results, by = c('sample_short')) %>%
  dplyr::select(-sample_short) %>% 
  cor(method = c("pearson"), use = "pairwise.complete.obs")
res_depth1 <- cor.mtest(vae_depth1_cor, conf.level=.99) # Signifikanztest

pdf(file= paste(output_figure_direc,'VAE_component_correlation_ICs','.pdf',sep=''),
    width= 9,
    height = 6) 
corrplot(vae_depth1_cor,
         is.corr=T,
         tl.cex= .8,
         tl.col = "black",
         col=colors_plot,
         title='',
         mar=c(0,0,2,0), #margins fürs pdf plotten
         method= "square",
         type = "upper",
         cl.cex = 1,
         p.mat = res_depth1$p, sig.level = 0.05, pch.cex = .3,
         insig = "blank",
         na.label = "o")
dev.off()


##show the input somatic features correlating with VAE components
input_features <- read.csv(file= paste(input_file_direc,'TCGA_Hartwig_PCAWG_samples_vs_phenotypes_withSigs_least45of56_allforGWAS_14832samples_table.txt',sep=''),
                           head=T,sep='\t', stringsAsFactors = F) 

#check correlation with input features
vae_correlation_depth1 <- vae_depth1 %>%
  left_join(input_features, by=c('sample_short')) %>%
  dplyr::select(-sample_short)
vae_features_correlation1 <- t(cor(vae_correlation_depth1[1:14],vae_correlation_depth1[(14+1):(dim(vae_correlation_depth1)[2])]))

cor <-
  data.frame(vae_features_correlation1) %>%
  mutate(pheno= row.names(vae_features_correlation1), info='correlation') %>%
  reshape2::melt(id.vars=c('pheno','info'))  %>%
  mutate(pheno= sub('transcription_bias','trx',pheno)) %>%
  mutate(pheno= sub('small_indel_','indel_',pheno)) %>%
  mutate(pheno= sub('deletions_10bp','deletions>=10bp',pheno)) %>%
  mutate(pheno= sub('MSI','MS',pheno)) %>%
  mutate(pheno= sub('mh_2bp_more','mh>=2bp',pheno)) %>%
  mutate(variable=sub('depth1_','',variable)) %>%
  mutate(info=factor(info,ordered = F,levels = c('correlation'))) %>%
  left_join(mapping_features, by=c('variable')) %>%
  mutate(pheno_compo_name=paste(variable,':\n',pheno_name,sep=''))

plotOut <- NULL
n=1
for(VAE in c('VAE_1','VAE_2','VAE_3','VAE_4',
            'VAE_5','VAE_6','VAE_7','VAE_8',
            'VAE_9','VAE_10','VAE_11','VAE_12',
            'VAE_13','VAE_14')){
  
  print(VAE)
  
  rm(contri_cor_selected_order)
  contri_cor_selected_order <-
    cor %>%
    dplyr::filter(variable == VAE) %>%
    #dplyr::filter(info == 'correlation') %>%
    mutate(pheno=factor(pheno,ordered = F,levels = cor$pheno[order(abs(cor$value[cor$variable==VAE & cor$info== 'correlation']), decreasing = T)])) %>%
    dplyr::filter(pheno %in% cor$pheno[order(abs(cor$value[cor$variable==VAE & cor$info== 'correlation']), decreasing = T)][1:10]) #only the 20 with the highest rankings
  
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
    scale_fill_manual(values=c("correlation" = "white"))
  
  plotOut[[paste('bla',n,sep='_')]] <- bla
  rm(bla)
  n= n+1
}

##### #plot other variables
pdf(file= paste(output_figure_direc,'VAE_overview','.pdf',sep=''),
    width= 8.2,
    height = 9)   
annotate_figure(
  ggarrange(plotOut$bla_1,plotOut$bla_2,plotOut$bla_3,
          plotOut$bla_4,plotOut$bla_5,plotOut$bla_6,
          plotOut$bla_7,plotOut$bla_8,plotOut$bla_9,
          plotOut$bla_10,plotOut$bla_11,plotOut$bla_12,
          plotOut$bla_13,plotOut$bla_14,
          nrow=4,ncol=4,common.legend = F,align='hv'),
  left = text_grob("Pearson correlation", color = "black", rot = 90),
  bottom = text_grob("somatic components", color = "black")
)
dev.off()
