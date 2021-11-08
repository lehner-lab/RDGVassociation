# Extraction of Somatic Mutational Components

## Description of Scripts

### Extraction of Independent Components
### Extractions of Variational Autoencoder derived Components


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
