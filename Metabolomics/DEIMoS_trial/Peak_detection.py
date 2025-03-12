#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
import deimos
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
os.chdir('/parallel_scratch/mp01950/multiomics/')
mzML_dir = os.getcwd()+"/Brain_trauma/MS2/"
mzML_files = os.listdir(mzML_dir)
out_path = os.getcwd()+"/deimos_trial/peak_detected/"


# In[7]:


for h5_file in h5_files:
    #Data frame loading
    ms1 = deimos.load(h5_file, key='ms1')
    ms2 = deimos.load(h5_file, key='ms2')
    if 'drift_time' in ms1.columns:
        dimension = ['mz', 'drift_time', 'retention_time']
        smooth_radius = [0, 1, 0]
        pp_radius = [2, 10, 4]
    else:
        dimension = ['mz', 'retention_time']
        smooth_radius = [0, 1]
        pp_radius = [1, 5]
    ##Peak detection
    # Build factors from raw data
    factors1 = deimos.build_factors(ms1, dims = dimension)
    factors2 = deimos.build_factors(ms2, dims = dimension)
    # Nominal threshold
    ms1 = deimos.threshold(ms1, threshold=500)
    ms2 = deimos.threshold(ms2, threshold=500)
    # Build index
    index1 = deimos.build_index(ms1, factors1)
    index2 = deimos.build_index(ms2, factors2)
    # Smooth data
    ms1 = deimos.filters.smooth(ms1, index=index1, dims = dimension, radius = smooth_radius, iterations=7)
    ms2 = deimos.filters.smooth(ms2, index=index2, dims = dimension, radius = smooth_radius, iterations=7)
    # Perform peak detection
    ms1_peaks = deimos.peakpick.persistent_homology(ms1, index=index1, dims=dimension, radius = pp_radius)
    ms2_peaks = deimos.peakpick.persistent_homology(ms2, index=index2, dims=dimension, radius = pp_radius)
    # Sort by persistence
    ms1_peaks = ms1_peaks.sort_values(by='persistence', ascending=False).reset_index(drop=True)
    ms2_peaks = ms2_peaks.sort_values(by='persistence', ascending=False).reset_index(drop=True)
    name = h5_file.rsplit("/", 1)[-1].rsplit(".h5", 1)[0]
    # Save ms1 to new file
    deimos.save(out_path+name+'_peaks.h5', ms1_peaks, key='ms1', mode='w')
    # Save ms2 to same file
    deimos.save(out_path+name+'_peaks.h5', ms2_peaks, key='ms2', mode='a')


# In[ ]:





# In[ ]:




