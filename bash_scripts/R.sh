#!/bin/bash

#SBATCH --partition=high_mem             #Selecting “shared” Queue
#SBATCH --nodes=2                       #No of nodes to run job
#SBATCH --ntasks-per-node=5            #No of cores to use per node
#SBATCH --time=36:00:00                 #Maximum time for job to run
#SBATCH --mem=15G                        #Amount of memory per node

Rscript /users/mp01950/bash_temp/Norm_and_Inte_without_154826.R
