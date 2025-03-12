#!/bin/bash

#SBATCH --partition=high_mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=120:00:00                 #Maximum time for job to run
#SBATCH --mem=200G                        #Amount of memory per node
#SBATCH --job-name=DEIMoS_add

# Run Python script for the selected group
python /parallel_scratch/mp01950/multiomics/deimos_trial/MS2_deconvolution.py
