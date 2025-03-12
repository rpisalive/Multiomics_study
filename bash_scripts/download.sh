#!/bin/bash

# Base URL parts
base_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/GSE154826%5Famp%5Fbatch%5FID%5F"

# Loop through each directory
for dir in GSE154826_batch_*; do
	    if [[ -d "$dir" ]]; then
		            # Extract the number from the directory name
			            number=$(echo "$dir" | grep -oP '(?<=_batch_)\d+')
				            
				            # Construct the full URL
					            url="${base_url}${number}.tar.gz"
						            
						            # Download the file
							            wget "$url"
								        fi
								done
