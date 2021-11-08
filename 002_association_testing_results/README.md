# Results of RDGV Association Study

## Description of Scripts
* **lambdas.R:** Estimation of inflation factors for all tests for selected gene set (enriched for DNA repair genes).
* **lambdas_randomized_geneset.R:** Estimation of inflation factors for all tests for the randomized gene set.
* **Empirical_FDR_thresholds.R:** Estimation of False Discovery Rates across cancer_types and pan-can for selected gene set.
* **Empirical_FDR_thresholds_randomized_geneset.R:** Estimation of False Discovery Rates across cancer_types and pan-can for randomized gene set.
* **Rare_Testing_Association_Results_1percent_FDR.R:** Main results and stats at a FDR of 1%.
* **Rare_Testing_Association_Results_2percent_FDR.R:** Main results and stats at a FDR of 2%.
* **Rare_Testing_Association_Results_distribution_hits.R:** Overview of distribution of hits across components, RDGV sets, models of inheritance and cancer types.
* **Rare_Testing_Association_Results_randomized_geneset.R:** Number of replicated hits when using the randomized set of genes. 2nd approach of calculating an upper limit of the False Discovery Rate.
* **Distribution_hits_sample_size.R:** Plotting the number of replicated hits across cancer types vs. sample size.
* **venn_overlap_hits.R:**  Venn diagramms for overview of replicated hits. 
* **Prevalence_count.R:**  Counting number of RDGVs across different set of genes and RDGV sets.
