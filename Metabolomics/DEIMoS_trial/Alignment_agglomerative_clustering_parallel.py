#!/usr/bin/env python
# coding: utf-8

# In[16]:


import numpy as np
import deimos
import os
import pandas as pd
import sys
from sklearn.cluster import AgglomerativeClustering
from scipy.sparse import coo_matrix
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import cKDTree


# In[11]:


# Set working directory
os.chdir('/parallel_scratch/mp01950/multiomics/')
out_path = os.path.join(os.getcwd(), "deimos_trial/alignment_result/")


# In[12]:


# Get input argument (group name)
if len(sys.argv) != 2:
    print("Usage: python Alignment_agglomerative_clustering.py <group_name>")
    sys.exit(1)

group_name = sys.argv[1]
print(f"Processing group: {group_name}")


# In[14]:


# Function to assign groups based on the filename pattern
def assign_group_from_filename(filename):
    if 'QC' in filename:
        return 'QC'
    elif 'S12' in filename:
        return 'DO50'
    elif 'S19' in filename:
        return 'PO50'
    elif 'S22' in filename:
        return 'PIV10'
    elif 'S7' in filename:
        return 'DIV10'
    elif 'S1' in filename:
        return 'PO0'
    elif 'S9' in filename:
        return 'DO0'
    else:
        return None  # Return None if no group found

# Get all relevant H5 files
all_h5_files = sorted([os.path.join(os.getcwd(), "deimos_trial/peak_detected", f)
                        for f in os.listdir(os.getcwd() + "/deimos_trial/peak_detected") if f.endswith(".h5")])

# Filter files for the given group based on filename
group_files = [f for f in all_h5_files if assign_group_from_filename(os.path.basename(f)) == group_name]

if not group_files:
    print(f"No files found for group {group_name}, exiting...")
    sys.exit(1)


# In[15]:


# Load data only for the specified group
data = {}
for h5 in group_files:
    filename = os.path.basename(h5).replace("_peaks.h5", "")
    data[filename] = deimos.load(h5, key='ms1',
                                 columns=['scan', 'mz', 'retention_time', 'intensity', 'persistence',
                                          'mz_weighted', 'retention_time_weighted'])

# Create DataFrame for the specified group
df_list = []
for filename, df in data.items():
    identifier = filename.split("_")[1]
    group = assign_group_from_filename(filename)  # Assign group using the function
    if group == group_name:
        df = df.copy()
        df["sample"] = filename
        df_list.append(df)

if not df_list:
    print(f"Error: No data found for group '{group_name}'")
    sys.exit(1)

df = pd.concat(df_list, ignore_index=True)
del data, df_list  # Free memory

# Optimize data types
df["retention_time"] = df["retention_time"].astype("float32")
df["mz"] = df["mz"].astype("float32")
df["intensity"] = df["intensity"].astype("float32")

# Assign sample_idx
sample_mapping = {sample: i for i, sample in enumerate(sorted(df["sample"].unique()))}
df["sample_idx"] = df["sample"].map(sample_mapping)


# In[ ]:


# Calculate persistence ratio and filter based on the threshold
def filter_by_persistence_ratio(df, threshold=0.8):
    # Add persistence ratio column
    df['persistence_ratio'] = df['persistence'] / df['intensity']
    
    # Apply threshold to filter rows based on persistence_ratio
    filtered_df = deimos.threshold(df, by='persistence_ratio', threshold=threshold)
    
    return filtered_df


# In[ ]:


# Filter the DataFrame before alignment
print(f"Filtering peaks based on persistence ratio for {group_name}...")
df_filtered = filter_by_persistence_ratio(df, threshold=0.5)


# In[ ]:


# Function to create connectivity matrix (prevents within-sample clustering)
def create_connectivity_matrix(df, tol=0.1):
    sample_indices = df["sample_idx"].values
    rt_values = df["retention_time"].values.reshape(-1, 1)

    # Construct KDTree for efficient neighborhood search
    tree = cKDTree(rt_values)
    pairs = tree.query_pairs(r=tol)

    if not pairs:
        return coo_matrix(([], ([], [])), shape=(len(df), len(df))).tocsr()

    row, col = zip(*pairs)

    # Ensure inter-sample connectivity only
    mask = [sample_indices[r] != sample_indices[c] for r, c in zip(row, col)]
    row = [r for r, m in zip(row, mask) if m]
    col = [c for c, m in zip(col, mask) if m]

    # Mirror the indices to maintain symmetry
    row.extend(col)
    col.extend(row)

    data = np.ones(len(row), dtype=np.uint8)
    return coo_matrix((data, (row, col)), shape=(len(df), len(df))).tocsr()


# In[ ]:


# Agglomerative clustering function with DEIMoS specifications
def agglomerative_alignment_subset(df, tol=0.1):
    if df.empty:
        return df

    rt_values = df["retention_time"].values.reshape(-1, 1)

    # Use Chebyshev (L-infinity norm) distance
    distance_matrix = squareform(pdist(rt_values, metric="chebyshev").astype(np.float32))

    # Create connectivity matrix to prevent within-sample clustering
    connectivity = create_connectivity_matrix(df, tol=tol)

    # Perform clustering with complete linkage
    clustering = AgglomerativeClustering(
        n_clusters=None,
        metric="precomputed",
        linkage="complete",
        distance_threshold=tol,
        connectivity=connectivity
    )

    labels = clustering.fit_predict(distance_matrix)
    df['cluster'] = labels

    # Assign median retention time to each cluster
    cluster_medians = df.groupby('cluster')['retention_time'].median().to_dict()
    df['rt_aligned'] = df['cluster'].map(cluster_medians)

    return df


# In[ ]:


# Bin-based alignment function
def bin_based_alignment(df, tol=0.1, bin_width=0.5):
    if df.empty:
        return df

    df_sorted = df.sort_values("retention_time").reset_index(drop=True)
    min_rt, max_rt = df_sorted["retention_time"].min(), df_sorted["retention_time"].max()
    bins = np.arange(min_rt, max_rt + bin_width, bin_width)

    clustered_dfs = []
    for i in range(len(bins) - 1):
        bin_df = df_sorted[(df_sorted["retention_time"] >= bins[i]) & (df_sorted["retention_time"] < bins[i + 1])]
        if bin_df.empty:
            continue
        
        clustered_bin = agglomerative_alignment_subset(bin_df, tol=tol)
        clustered_dfs.append(clustered_bin)

    return pd.concat(clustered_dfs, ignore_index=True) if clustered_dfs else df_sorted.copy()


# In[ ]:


# Perform alignment
print(f"Running bin-based alignment for {group_name}...")
aligned_df = bin_based_alignment(df_filtered, tol=0.1, bin_width=0.5)


# In[ ]:


# Save results
output_file = os.path.join(out_path, f"{group_name}.csv")
aligned_df.to_csv(output_file, index=False, float_format="%.5f")
print(f"Finished processing {group_name}. Results saved to {output_file}")


# In[ ]:





# In[ ]:





# In[ ]:




