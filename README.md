# RDGVassociation
Scripts and data to replicate paper figures.
Scripts to extract components (either via ICA or VAE).

## Getting Started
* Run all scripts in their respective directories (otherwise need to change directories in scripts)
* Variational Autoencoder in Python, recommended to run in Singularity Environment (more info in /001_somatic_component_extraction/)
* Other scripts can be all run in R (at least version 3.5.0)
* Required R packages: tidyverse, dplyr, reshape2, ggplot2, wesanderson, gplots, cluster, fastICA, corrpot, ggrepel

## Content
### 001_somatic_component_extraction
* Scripts to extract somatic components via Independent Component Analysis (ICA) or Variational Autoencoder (VAE) neural network (input files provided)
* Scripts to replicate manuscript figures (output files from somatic component extraction provided)

### 002_association_testing_results
* All resuts from association study between rare damaging germline variants (RDGVs) and somatic components using [SKAT-O](https://www.cell.com/ajhg/fulltext/S0002-9297(12)00316-3) provided
* Scripts to replicate manuscript figures

### 003_network_analysis_results
* Scripts to perform network analysis as described in manuscript
* Scripts to replicate manuscript figures

## Visualization of association results via Shiny APP
https://mischanvp.shinyapps.io/rare_association_shiny/

