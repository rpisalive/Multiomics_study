#!/bin/bash

#SBATCH --partition=shared             #Selecting “shared” Queue
#SBATCH --nodes=2                       #No of nodes to run job
#SBATCH --ntasks-per-node=5            #No of cores to use per node
#SBATCH --time=72:00:00                 #Maximum time for job to run
#SBATCH --mem=15G                        #Amount of memory per node
#SBATCH --job-name=jupyter

# Define paths
CONTAINER=/parallel_scratch/mp01950/MetaboAnalystR_container/metaboanalystr.sif
METAR_ENV=/users/mp01950/.conda/envs/metar

# Bind and launch Jupyter Notebook
apptainer exec --bind ${METAR_ENV}:/metar \
               --bind /parallel_scratch/mp01950:/workspace \
               --env LANG=C.UTF-8,LC_ALL=C.UTF-8 \
               ${CONTAINER} \
               bash -c "export PATH=/metar/bin:\$PATH && \
                        export PYTHONPATH=/metar/lib/python3.9/site-packages:\$PYTHONPATH && \
                        /metar/bin/jupyter notebook --no-browser --ip=0.0.0.0 --port=8888 --notebook-dir=/workspace"

