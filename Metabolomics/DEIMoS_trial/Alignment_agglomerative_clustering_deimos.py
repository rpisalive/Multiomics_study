#!/usr/bin/env python
# coding: utf-8

# In[28]:


import argparse
import os
import dask.dataframe as dd
import pandas as pd
from dask.distributed import Client
import deimos
import numpy as np
import sys
import gc


# In[29]:


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run DEIMoS Alignment on SLURM Array")
    parser.add_argument("--group", type=str, required=True, help="Experimental group to process")
    
    # Parse arguments
    args = parser.parse_args()

    # Debug output: print the received group argument
    print(f"Received group: {args.group}")

    # Check if the group is valid
    group_to_process = args.group
    print(f"Processing group: {group_to_process}")

    # Define valid groups
    GROUPS = ["QC", "DO50", "PO50", "PIV10", "DIV10", "PO0", "DO0"]

    if group_to_process not in GROUPS:
        print(f"Group {group_to_process} not found. Exiting.")
        sys.exit(1)

    # Initialize Dask client for distributed processing
    client = Client(processes=False)  # Ensure Dask does not spawn processes

    # Define input/output paths
    os.chdir('/parallel_scratch/mp01950/multiomics/')
    out_path = os.getcwd() + "/deimos_trial/alignment_result/"

    # Get the full paths of all .h5 files
    h5_files = sorted([os.path.join(os.getcwd() + "/deimos_trial/peak_detected", f) for f in os.listdir(os.getcwd() + "/deimos_trial/peak_detected") if f.endswith(".h5")])

    # Load the data into a dictionary using DEIMoS
    data = {}
    for h5 in h5_files:
        filename = h5.split("/")[-1].replace("_peaks.h5", "")
        data[filename] = deimos.load(h5, key='ms1', columns=['scan', 'mz', 'retention_time', 'intensity', 'persistence', 'mz_weighted', 'retention_time_weighted'])

    # Define group mapping based on the identifier
    group_mapping = {
        "QC": "QC",
        "S12": "DO50",
        "S9": "PO50",
        "S19": "PIV10",
        "S22": "DIV10",
        "S1": "PO0",
        "S7": "DO0"
    }

    # Initialize dictionary to store grouped DataFrames
    grouped_data = {group: [] for group in set(group_mapping.values())}

    # Iterate over the data dictionary to assign samples to groups
    for filename, df in data.items():
        parts = filename.split("_")  # Split by underscores
        print(f"Filename: {filename}")
        print(f"Parts: {parts}")
        identifier = parts[1]  # The identifier is the second part of the filename

        # Determine the correct group based on the identifier
        group = group_mapping.get(identifier)
        print(f"Group: {group}")
        if group:
            df = df.copy()  # Avoid modifying original data
            df["sample"] = filename  # Assign the filename to a new column
            grouped_data[group].append(df)

    # Filter for the selected experimental group (if specified)
    if group_to_process in grouped_data:
        grouped_data = {group_to_process: grouped_data[group_to_process]}
    else:
        print(f"Group {group_to_process} not found.")
        return

    # Concatenate DataFrames within the selected group
    for group in grouped_data:
        grouped_data[group] = dd.from_pandas(pd.concat(grouped_data[group], ignore_index=True), npartitions=10)

    # Assign sample_idx for clustering
    all_samples = set()
    for group, df in grouped_data.items():
        all_samples.update(df["sample"].unique())

    sample_mapping = {sample: i for i, sample in enumerate(sorted(all_samples))}

    for group, df in grouped_data.items():
        df["sample_idx"] = df["sample"].map(sample_mapping)

    # Perform agglomerative clustering for each group using DEIMoS
    def align_group_data(df):
        if 'drift_time' in df.columns:
            dims = ['mz', 'retention_time', 'drift_time']
            tol = [2e-05, 0.03, 0.3]
            relative = [True, True, False]
        else:
            dims = ['mz', 'retention_time']
            tol = [2e-05, 0.3]
            relative = [True, False]

        return deimos.alignment.agglomerative_clustering(df, dims=dims, tol=tol, relative=relative)

    # Apply alignment for the selected experimental group
    aligned_grouped_data = {group: df.map_partitions(align_group_data) for group, df in grouped_data.items()}

    # Compute and save results
    for group, df in aligned_grouped_data.items():
        print(f"Processing {group}...")
        result_df = df.compute()
        output_file = os.path.join(out_path, f"{group}_clustered.csv")
        result_df.to_csv(output_file, index=False)
        del result_df
        gc.collect()
        print(f"Saved {output_file}")


# In[30]:


# Ensure the script runs only when executed directly
if __name__ == "__main__":
    main()


# In[ ]:




