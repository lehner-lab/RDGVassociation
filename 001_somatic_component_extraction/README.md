# Extraction of Somatic Mutational Components

## Description of Scripts

### Extraction of Independent Components
* **ICA_extraction.R:** Extraction of Independent Components from 2 to 30. Runs over 24h and not recommended to be run locally. Final plots with ICA extraction using 15 components.
* **ICA_overview.R:** Plotting correlations/contributions of input somatic features with ICs.


### Extractions of Variational Autoencoder derived Components
* **VAE_tensorflow1.py:** Running Variational Autoencoder with dataset spit (90 % training and 10% test) to find optimal parameters. Recommended to be run in Singularity environment.
* **VAE_tensorflow1_nosplit_alldata.py:** Running on complete dataset after finding optimal parameters. 
* **VAE_overview.R:** Plotting results from VAE: parameter sweep, optimal parameters, overview of selected components, correlation with selected ICs.


## Running Variational Autoencoder Neural Network in Singularity Environment
* Download Docker image as Singularity image:
```
singularity pull docker://mvpandapaw/tensorflow1.15.5_gpu_jupyter_moredependencies:v1
```
* Underying Dockerfile can be found in the directory /Dockerfile/
* Scripts can then be run as follows:
```
singularity exec tensorflow1.15.5_gpu_jupyter_moredependencies-v1.simg python script.py
```
