#!/bin/bash

#SBATCH --partition=high_mem             #Selecting “shared” Queue
#SBATCH --output=/parallel_scratch/mp01950/multiomics/deimos_trial/logs/alignment_%A_%a.out
#SBATCH --error=/parallel_scratch/mp01950/multiomics/deimos_trial/logs/alignment_%A_%a.err
#SBATCH --ntasks=1                       # 1 task per job array element
#SBATCH --cpus-per-task=4                # 4 CPU cores per task
#SBATCH --time=120:00:00                 # Max job runtime
#SBATCH --mem=350G                       # Requesting 250GB per task
#SBATCH --array=0-6                       # Job array for 7 groups (0-based index)
#SBATCH --job-name=DEIMoS_align

# Set environment variables for threading if needed
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

# List of experiment groups
groups=("QC" "DO50" "PO50" "PIV10" "DIV10" "PO0" "DO0")

# Select the group based on SLURM_ARRAY_TASK_ID
GROUP=${groups[$SLURM_ARRAY_TASK_ID]}

# Run Python script for the assigned group
python /parallel_scratch/mp01950/multiomics/deimos_trial/Alignment_agglomerative_clustering_parallel.py "$GROUP"
