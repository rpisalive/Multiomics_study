#!/usr/bin/env python
# coding: utf-8

# In[1]:


import deimos
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import ast


# In[2]:


os.chdir('/parallel_scratch/mp01950/multiomics/')
out_path = os.getcwd()+"/deimos_trial/iso_detected/"


# In[3]:


# Get the full paths of all .h5 files
h5_files = sorted([os.path.join(os.getcwd()+"/deimos_trial/peak_detected/", f) for f in os.listdir(os.getcwd()+"/deimos_trial/peak_detected/") if f.endswith(".h5")])
#h5_files = h5_files[0:1]
h5_files


# In[4]:


#This function is created to adjust the length of mz_iso & intensity_iso
#The deimos.isotopes.detect function would collapse mz_iso and intensity_iso of the same value, leading to inconsistent length and subsequent problem in adduct detection.
def update_isotopes(isotopes, ms1_peaks):
    for i, row in isotopes.iterrows():
        idx_iso_list = row['idx_iso']  # List of index numbers referring to ms1_peaks
        mz_iso_list = row['mz_iso']  # List of m/z values
        intensity_iso_list = row['intensity_iso']  # List of intensity values
        
        updated_mz_iso = []
        updated_intensity_iso = []
        
        for j, idx in enumerate(idx_iso_list):
            if idx in ms1_peaks.index:
                mz_value = ms1_peaks.loc[idx, 'mz']
                intensity_value = ms1_peaks.loc[idx, 'intensity']
                
                if j < len(mz_iso_list) and j < len(intensity_iso_list):
                    # Check if values match
                    if mz_iso_list[j] == mz_value and intensity_iso_list[j] == intensity_value:
                        updated_mz_iso.append(mz_iso_list[j])
                        updated_intensity_iso.append(intensity_iso_list[j])
                    else:
                        # Replace with correct values from ms1_peaks
                        updated_mz_iso.append(mz_value)
                        updated_intensity_iso.append(intensity_value)
                else:
                    # Append missing values
                    updated_mz_iso.append(mz_value)
                    updated_intensity_iso.append(intensity_value)
        
        # Update the row in isotopes DataFrame
        isotopes.at[i, 'mz_iso'] = updated_mz_iso
        isotopes.at[i, 'intensity_iso'] = updated_intensity_iso
    
    return isotopes


# In[5]:


for h5 in h5_files:
    ms1_peaks = deimos.load(h5, key='ms1', columns = ['scan', 'mz', 'retention_time', 'intensity', 'persistence', 'mz_weighted', 'retention_time_weighted'])
    ms2_peaks = deimos.load(h5, key='ms2', columns = ['scan', 'mz', 'retention_time', 'intensity', 'persistence', 'mz_weighted', 'retention_time_weighted'])
    # Partition the data
    partitions = deimos.partition(ms1_peaks, size=2000, overlap=5.1)
    # Map isotope detection over partitions
    isotopes = partitions.map(deimos.isotopes.detect,
                              dims=['mz', 'retention_time'],
                              tol=[0.1, 0.15],
                              delta=1.003355,
                              max_isotopes=5,
                              max_charge=1,
                              max_error=50E-6)
    isotopes = update_isotopes(isotopes, ms1_peaks)
    # Ensure idx in isotopes is an integer index (convert to match ms1_peaks index type)
    isotopes['idx'] = isotopes['idx'].astype(int)
    # Set index of ms1_peaks to ensure proper mapping
    ms1_peaks_merged = ms1_peaks.copy()
    # Merge isotope columns into ms1_peaks by aligning on index
    ms1_peaks_merged = ms1_peaks_merged.join(isotopes.set_index('idx')[['charge', 'multiple', 'dx', 'mz_iso', 'intensity_iso', 'idx_iso', 'error', 'decay', 'n']], how='left')
    float_columns = ['charge','n']
    list_columns = ['multiple', 'dx', 'mz_iso', 'intensity_iso', 'idx_iso', 'error', 'decay']
    for col in float_columns:
        ms1_peaks_merged[col] = pd.to_numeric(ms1_peaks_merged[col], errors="coerce")
    for col in list_columns:
        ms1_peaks_merged[col] = ms1_peaks_merged[col].astype(str)
    name = h5.split('/')[-1].replace('_peaks.h5', '')    
    deimos.save(out_path + name + '_isopeaks.h5', ms1_peaks_merged, key='ms1', mode='w')
    deimos.save(out_path + name + '_isopeaks.h5', ms2_peaks, key='ms2', mode='a')  


# In[ ]:




