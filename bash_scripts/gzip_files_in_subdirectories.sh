#!/bin/bash

# Path to the parent directory containing the 77 directories
parent_directory="/parallel_scratch/mp01950/raw_data/Lung"

# Iterate over each subdirectory in the parent directory
for dir in "$parent_directory"/*; do
  if [ -d "$dir" ]; then
    echo "Processing directory: $dir"
    
    # Iterate over each file in the current directory
    for file in "$dir"/*; do
      if [ -f "$file" ]; then
        echo "Gzipping file: $file"
        gzip "$file"
      fi
    done
  fi
done

echo "All files have been gzipped."
