#!/usr/bin/env python
# coding: utf-8

# In[4]:


import deimos
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.stats import pearsonr
import json


# In[5]:


os.chdir('/parallel_scratch/mp01950/multiomics/')
out_path = os.getcwd()+"/deimos_trial/adduct_detected/"


# In[6]:


# Get the full paths of all .h5 files
h5_files = sorted([os.path.join(os.getcwd()+"/deimos_trial/iso_detected/", f) for f in os.listdir(os.getcwd()+"/deimos_trial/iso_detected/") if f.endswith(".h5")])
#h5_files = h5_files[0:1]
h5_files


# In[7]:


adduct_shifts = {
    "[M+H]+": 1.0073,      # Protonation (H+)
    "[M+Na]+": 22.9892,    # Sodium adduct (Na+)
    "[M+K]+": 38.9632,     # Potassium adduct (K+)
    "[M+NH4]+": 18.0338,   # Ammonium adduct (NH4+)
    "[M+H-H2O]+": -18.0106 # Loss of H2O
}


# In[8]:


# Set tolerances
mz_tolerance = 0.01  # m/z tolerance for matching
rt_tolerance = 0.1  # Retention time tolerance in seconds
iso_intensity_threshold = 0.7  # Minimum correlation threshold for isotopic similarity


# In[9]:


def process_iso_column(iso_str):
    if isinstance(iso_str, str) and iso_str.strip():
        try:
            return [float(x) for x in iso_str.strip("[]").split(',') if x.strip()]
        except ValueError:
            return []
    return []


# In[10]:


def is_constant(array):
    """Check if an array is empty or has constant values."""
    return len(array) == 0 or np.all(array == array[0])


# In[11]:


def normalize_intensities(intensities):
    """Normalize intensities by dividing by the sum of intensities."""
    total = np.sum(intensities)
    return intensities / total if total > 0 else intensities  # Avoid division by zero


# In[12]:


def align_isotopic_peaks(parent_iso_mz, parent_iso_intensity, adduct_iso_mz, adduct_iso_intensity, shift, ppm_tolerance=40):
    aligned_parent_intensities = []
    aligned_adduct_intensities = []
    
    for parent_mz, parent_intensity in zip(parent_iso_mz, parent_iso_intensity):
        expected_adduct_mz = parent_mz + shift
        ppm_tol = expected_adduct_mz * ppm_tolerance / 1e6  # Calculate ppm tolerance in Da
        
        closest_idx = np.argmin(np.abs(adduct_iso_mz - expected_adduct_mz))
        adduct_mz_shifted = adduct_iso_mz[closest_idx]

        if np.abs(adduct_mz_shifted - expected_adduct_mz) < ppm_tol:
            aligned_adduct_intensities.append(adduct_iso_intensity[closest_idx])
            aligned_parent_intensities.append(parent_intensity)
    
    return aligned_parent_intensities, aligned_adduct_intensities


# In[13]:


for h5 in h5_files:
    print(f"Processing file: {h5}")
    ms1 = deimos.load(h5, key='ms1')
    ms2 = deimos.load(h5, key='ms2')  # Load MS2 data
    
    print("Processing isotopic columns...")
    ms1['mz_iso'] = ms1['mz_iso'].apply(process_iso_column)
    ms1['intensity_iso'] = ms1['intensity_iso'].apply(process_iso_column)

    ms1['adducts'] = ['' for _ in range(len(ms1))]
    ms1['add_idx'] = ['' for _ in range(len(ms1))]

    for index, row in ms1.iloc[0:len(ms1)].iterrows():
        mz_value = row['mz']
        rt_value = row['retention_time']

        adducts_list = []
        add_idx_list = []

        print(f"Processing row {index} with m/z: {mz_value} and RT: {rt_value}...")

        for adduct, shift in adduct_shifts.items():
            expected_adduct_mz = mz_value + shift
            matching_rows = ms1[
                (ms1['mz'] >= expected_adduct_mz - mz_tolerance) & 
                (ms1['mz'] <= expected_adduct_mz + mz_tolerance) & 
                (abs(ms1['retention_time'] - rt_value) <= rt_tolerance)
            ]
            
            if matching_rows.empty:
                continue
            
            for match_idx in matching_rows.index:
                matched_row = ms1.loc[match_idx]

                if not row['mz_iso'] or not matched_row['mz_iso']:
                    continue

                parent_iso_mz = np.array(row['mz_iso'])
                adduct_iso_mz = np.array(matched_row['mz_iso'])
                parent_iso_intensity = np.array(row['intensity_iso'])
                adduct_iso_intensity = np.array(matched_row['intensity_iso'])

                print(f"Aligning isotopic peaks for {adduct}...")
                parent_aligned_intensities, adduct_aligned_intensities = align_isotopic_peaks(
                    parent_iso_mz, parent_iso_intensity, adduct_iso_mz, adduct_iso_intensity, shift
                )

                if len(parent_aligned_intensities) == 0 or len(adduct_aligned_intensities) == 0:
                    continue

                # Normalize intensities before correlation
                parent_aligned_intensities = normalize_intensities(np.array(parent_aligned_intensities))
                adduct_aligned_intensities = normalize_intensities(np.array(adduct_aligned_intensities))

                print(f"Checking intensity patterns for {adduct} at index {match_idx}...")

                if is_constant(parent_aligned_intensities) or is_constant(adduct_aligned_intensities):
                    print(f"Skipping correlation due to constant intensity values for {adduct}.")
                    continue

                corr, _ = pearsonr(parent_aligned_intensities, adduct_aligned_intensities)
                print(f"Normalized intensity correlation for {adduct} = {corr:.3f}")

                if corr >= iso_intensity_threshold:
                    adducts_list.append(adduct)
                    add_idx_list.append(str(match_idx))
                    print(f"Assigned adduct {adduct} at index {match_idx}.")

        ms1.at[index, 'adducts'] = ', '.join(adducts_list)
        ms1.at[index, 'add_idx'] = ', '.join(add_idx_list)

    ms1['mz_iso'] = ms1['mz_iso'].apply(json.dumps)
    ms1['intensity_iso'] = ms1['intensity_iso'].apply(json.dumps)

    name = h5.split('/')[-1].replace('_isopeaks.h5', '')
    print(f"Saving updated file: {name}_isoaddpeaks.h5")
    deimos.save(out_path + name + '_isoaddpeaks.h5', ms1, key='ms1', mode='w')
    deimos.save(out_path + name + '_isoaddpeaks.h5', ms2, key='ms2', mode='a')


# In[ ]:




