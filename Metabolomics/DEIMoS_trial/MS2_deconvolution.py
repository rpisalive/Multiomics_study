#!/usr/bin/env python
# coding: utf-8

# In[1]:


import deimos
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns


# In[2]:


os.chdir('/parallel_scratch/mp01950/multiomics/')
output_path = os.getcwd()+"/deimos_trial/MS2_decon/"


# In[3]:


MS2_ce = 4


# In[4]:


# Get the full paths of all raw .h5 files (i.e. peaks not detected)
raw_h5_files = sorted([os.path.join(os.getcwd()+"/deimos_trial/mzML2h5/", f) for f in os.listdir(os.getcwd()+"/deimos_trial/mzML2h5/") if f.endswith(".h5")])
raw_h5_files = raw_h5_files[0:1]
raw_h5_files
# Get the full paths of all .peaks detected .h5 files
h5_files = sorted([os.path.join(os.getcwd()+"/deimos_trial/peak_detected/", f) for f in os.listdir(os.getcwd()+"/deimos_trial/peak_detected/") if f.endswith(".h5")])
h5_files = h5_files[0:1]
h5_files


# In[5]:


for raw_file, file in zip(raw_h5_files, h5_files):  # Use zip() to iterate both lists in parallel
    ms1 = deimos.load(raw_file, key='ms1')
    ms2 = deimos.load(raw_file, key='ms2')

    if 'drift_time' in ms1.columns:
        dt = True
        columns_to_select = ['scan', 'mz', 'retention_time', 'drift_time', 'intensity']
    else:
        dt = False
        columns_to_select = ['scan', 'mz', 'retention_time', 'intensity']

    # Keep only existing columns
    existing_columns = [col for col in columns_to_select if col in ms1.columns]
    ms1 = ms1[existing_columns]
    ms2 = ms2[existing_columns]

    # Load peak files with the correct column structure
    peak_columns = existing_columns + ['persistence']
    ms1_peaks = deimos.load(file, key='ms1', columns=peak_columns)
    ms2_peaks = deimos.load(file, key='ms2', columns=peak_columns)

    # Apply thresholding
    ms1 = deimos.threshold(ms1, threshold=500)
    ms2 = deimos.threshold(ms2, threshold=500)
    ms1_peaks = deimos.threshold(ms1_peaks, threshold=5E3)
    ms2_peaks = deimos.threshold(ms2_peaks, threshold=5E2)

    #Indexing resetting only required for subsetting
    # Reset index
    #ms1_peaks = ms1_peaks.reset_index(drop=False)
    #ms2_peaks = ms2_peaks.reset_index(drop=False)
    # Rename index column
    #ms1_peaks = ms1_peaks.rename(columns={'index': 'original_index'})
    #ms2_peaks = ms2_peaks.rename(columns={'index': 'original_index'})

    # Deconvolution
    decon = deimos.deconvolution.MS2Deconvolution(ms1_peaks, ms1, ms2_peaks, ms2)

    if not dt:
        decon.construct_putative_pairs(dims=['retention_time'], low=[-0.1], high=[0.1], ce=MS2_ce, require_ms1_greater_than_ms2=True)
        decon.configure_profile_extraction(dims=['mz', 'retention_time'], low=[-200E-6, -0.1], high=[600E-6, 0.1], relative=[True, False])
        res = decon.apply(dims='retention_time', resolution=0.01)
    else:
        decon.construct_putative_pairs(dims=['drift_time', 'retention_time'], low=[-0.12, -0.1], high=[1.4, 0.1], ce=MS2_ce, require_ms1_greater_than_ms2=True, error_tolerance=0.12)
        decon.configure_profile_extraction(dims=['mz', 'drift_time', 'retention_time'], low=[-200E-6, -0.05, -0.1], high=[600E-6, 0.05, 0.1], relative=[True, True, False])
        res = decon.apply(dims=['drift_time', 'retention_time'], resolution=[0.01, 0.01])

        # Filtering results
    filtered_res = res.loc[res['retention_time_score'] > 0.9].groupby(
        by=[x for x in res.columns if x.endswith('_ms1')],
        as_index=False
    ).agg(list).sort_values(by='persistence_ms1', ascending=False)

    # Saving output
    name = file.split('/')[-1].replace('_peaks.h5', '')
    res.to_csv(output_path + name + '_5E4full.csv', index=True)
    filtered_res.to_csv(output_path + name + '_5E4filtered.csv', index=True)


# In[8]:


#Deconvolution would not work with just 1 scan if there is only one unique retention time in MS1 and MS2
#ms1_subset = ms1_subset[ms1_subset['scan'] == 5]
#ms2_subset = ms2_subset[ms2_subset['scan'] == 5]
#ms1_peaks_subset = ms1_peaks[ms1_peaks['scan'] == 5]
#ms2_peaks_subset = ms2_peaks[ms2_peaks['scan'] == 5]


# In[ ]:




