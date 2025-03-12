#!/bin/bash

# Loop through each .tar.gz file
for file in GSE154826_amp_batch_ID_*.tar.gz; do
    # Extract the number from the file name
    number=$(echo "$file" | grep -oP '(?<=_ID_)\d+(?=\.tar\.gz)')
    
    # Construct the corresponding directory name
    dir="GSE154826_batch_$number"
    
    # Check if the directory exists
    if [[ -d "$dir" ]]; then
        # Extract the .tar.gz file into the corresponding directory
        tar -xzvf "$file" -C "$dir"
    else
        echo "Directory $dir does not exist. Skipping extraction for $file."
    fi
done
