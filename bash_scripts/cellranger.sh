#!/bin/bash

cellranger count --id=GSM7774439 \
   --fastqs=/parallel_scratch/mp01950/rawfastq \
   --sample=RHP7294,RHP7295,RHP7296,RHP7297 \
   --transcriptome=/parallel_scratch/mp01950/refdata-gex-GRCh38-2020-A \
   --output-dir=/parallel_scratch/mp01950/GSM7774439
