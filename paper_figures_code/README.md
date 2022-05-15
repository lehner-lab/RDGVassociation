# Code for main manuscript figures

## Description of Scripts
* **lambdas.R:** Estimation of inflation factors for all tests for the selected gene set (enriched for DNA repair genes).
* **lambdas_randomized_geneset.R:** Estimation of inflation factors for all tests for the randomized gene set.
* **Empirical_FDR_thresholds.R:** Estimation of False Discovery Rates across cancer_types and pan-cancer for selected gene set.
* **Empirical_FDR_thresholds_randomized_geneset.R:** Estimation of False Discovery Rates across cancer_types and pan-cancer for the randomized gene set.
* **Rare_Testing_Association_Results_replicated_genes.R:** Main results and stats at a FDR of 1%.
* **Rare_Testing_Association_Results_distribution_hits.R:** Overview of distribution of hits across components, rare pLoF (putative loss of function) variant sets, models of inheritance, and cancer types.
* **Rare_Testing_Association_Results_randomized_geneset.R:** Number of replicated hits when using the randomized set of genes. 2nd approach of calculating an upper limit of the False Discovery Rate.
* **venn_overlap_hits.R:**  Plotting Venn diagrams to visualize the number of hits in different groups. 
* **Counted_variants_gnomAD.R:**  Calculating frequencies of rare pLoF (putative loss of function) variants across different set of genes and variants sets within the cancer cohort and gnomAD.
