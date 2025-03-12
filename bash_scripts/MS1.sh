#!/bin/bash

#SBATCH --partition=shared             #Selecting “shared” Queue
#SBATCH --nodes=2                       #No of nodes to run job
#SBATCH --ntasks-per-node=5            #No of cores to use per node
#SBATCH --time=48:00:00                 #Maximum time for job to run
#SBATCH --mem=15G                        #Amount of memory per node
#SBATCH --job-name=MS1

# Bind the 'metar' virtual environment into the Apptainer container
apptainer exec --bind /users/mp01950/.conda/envs/metar:/mnt/metar /parallel_scratch/mp01950/MetaboAnalystR_container/metaboanalystr.sif Rscript /parallel_scratch/mp01950/multiomics/MetaboAnalystR_trial/MS1.R

