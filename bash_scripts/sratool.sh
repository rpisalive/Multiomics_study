#!/bin/bash

for line in $(cat SRR_Acc_List.txt); do fastq-dump --gzip --split-files ; done
