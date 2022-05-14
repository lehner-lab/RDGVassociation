# Extraction of Somatic Mutational Components

## Description of Scripts

### Extraction of Independent Components
* **ICA_extraction.R:** Extraction of Independent Components from 2 to 30. Runs over 24h and it is not recommended to be run locally. Final plots with ICA extraction using 15 components.

### Extractions of Variational Autoencoder derived Components
* **VAE_tensorflow1.py:** Running Variational Autoencoder with dataset split (90 % training and 10% test) to find optimal parameters. Recommended to be run in Singularity environment. Example:
```
singularity exec tensorflow1.15.5_gpu_jupyter_moredependencies-v1.simg python VAE_tensorflow1.py --learning_rate 0.0005 --batch_size 200 --epochs 200 --num_components 14 --dataset_training 'TCGA_Hartwig_PCAWG_least45of56_14832samples_90split_balanced.txt' --dataset_test 'TCGA_Hartwig_PCAWG_least45of56_14832samples_10split_balanced.txt' --kappa 0.5 --depth 1
```
* **VAE_tensorflow1_nosplit_alldata.py:** Running on complete dataset after finding optimal parameters. Example:
```
singularity exec tensorflow1.15.5_gpu_jupyter_moredependencies-v1.simg python VAE_tensorflow1_nosplit_alldata.py --learning_rate 0.0005 --batch_size 200 --epochs 200 --num_components 14 --dataset_training 'TCGA_Hartwig_PCAWG_least45of56_14832samples_nosplit_balanced.txt' --kappa 0.5 --depth 1
```

### Final Set of Components
* **Components_overview.R:** Plotting heatmap with ICs and VAE-derived components.


## Running Variational Autoencoder Neural Network in Singularity Environment
* Download Docker image as Singularity image:
```
singularity pull docker://mvpandapaw/tensorflow1.15.5_gpu_jupyter_moredependencies:v1
```
* Underlying Dockerfile can be found in the directory /Dockerfile/
* Scripts can then be run as follows:
```
singularity exec tensorflow1.15.5_gpu_jupyter_moredependencies-v1.simg python script.py
```
