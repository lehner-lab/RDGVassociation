# Association of Rare Damaging Germline Variants with Somatic Mutational Components 
Scripts to extract components via ICA & VAE, to perform network analysis and to replicate manuscript figures: [The impact of rare germline variants on human somatic mutation processes](https://www.biorxiv.org/content/10.1101/2021.11.14.468508v1).

## Getting Started
* Scripts can be all run in R (at least version 3.5.0)
* Required R packages: tidyverse, dplyr, reshape2, ggplot2, wesanderson, gplots, cluster, fastICA, corrpot, ggrepel, ggrastr, ggpubr, RColorBrewer

## Content
* Description of the respective scripts can be found in each sub-directory

### somatic_component_extraction
* Scripts to extract somatic components via Independent Component Analysis (ICA) or Variational Autoencoder (VAE) neural network 
* Variational Autoencoder in Python, recommended running in Singularity Environment (more info in /somatic_component_extraction/)

### paper_figures_code
* All results from association study between rare damaging germline variant and somatic components using [SKAT-O](https://www.cell.com/ajhg/fulltext/S0002-9297(12)00316-3) are provided
* Scripts to replicate manuscript figures

### network_analysis
* Scripts to perform network analysis as described in the manuscript

## Interactive visualization of association results via Shiny app
https://mischanvp.shinyapps.io/rare_association_shiny/
