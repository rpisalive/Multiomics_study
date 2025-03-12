#!/bin/bash

#SBATCH --partition=shared
#SBATCH --ntasks-per-node=10
#SBATCH --mem=20G

for file in *h5ad; do cellbender remove-background --input $file --output ../clean_adata/$(basename $file .h5ad)denoised.h5ad --total-droplets-included 50000
# --cuda; done
