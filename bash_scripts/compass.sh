#!/bin/bash

#SBATCH --partition=shared 
#SBATCH --nodes=2 
#SBATCH --ntasks-per-node=10 
#SBATCH --mem=20G

compass --data /users/mp01950/GSE242894/raw_data/compass_input/add.tsv --output-dir /users/mp01950/GSE242894/raw_data/compass_result --temp-dir /parallel_scratch/mp01950/compass_tem --species homo_sapiens
