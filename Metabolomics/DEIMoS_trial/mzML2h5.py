#!/usr/bin/env python
# coding: utf-8

# In[9]:


import os
import deimos
import numpy as np
import matplotlib.pyplot as plt
os.chdir('/parallel_scratch/mp01950/multiomics/')
mzML_dir = os.getcwd()+"/Brain_trauma/MS2/"
mzML_files = os.listdir(mzML_dir)
out_path = os.getcwd()+"/deimos_trial/mzML2h5/"


# In[10]:


mzML_files


# In[12]:


accessions = deimos.get_accessions(mzML_dir+mzML_files[0])
accessions


# In[13]:


if any(key in accessions for key in ['drift time', 'ion mobility drift time']):
    dimension = ['mz', 'drift_time', 'retention_time', 'intensity']
else:
    dimension = ['mz', 'retention_time', 'intensity']

print(dimension)


# In[14]:


#%%time
for files in mzML_files:
    name = files.split(".mzML")[0]  # Splitting at ".mzML" and taking the first part
    if 'drift_time' in dimension:
        data = deimos.load(mzML_dir+files, accession={'retention_time': 'MS:1000016','drift_time': 'MS:1002476'})
    else:
        data = deimos.load(mzML_dir+files, accession={'retention_time': 'MS:1000016'})
    
    # Swapping retention_time if the retention times of ms1 and ms2 contrdict
    # If the scan number(s) is missing in either data frame, do nothing
    # Find common scan numbers in both MS1 and MS2
    common_scans = set(data["ms1"]["scan"]).intersection(set(data["ms2"]["scan"]))

    # Extract retention_time mappings only for common scans
    ms2_rt_dict = {scan: rt for scan, rt in zip(data["ms2"]["scan"], data["ms2"]["retention_time"]) if scan in common_scans}
    ms1_rt_dict = {scan: rt for scan, rt in zip(data["ms1"]["scan"], data["ms1"]["retention_time"]) if scan in common_scans}

    # Swap retention times only for common scans
    data["ms1"]["retention_time"] = data["ms1"].apply(lambda row: ms2_rt_dict.get(row["scan"], row["retention_time"]), axis=1)
    data["ms2"]["retention_time"] = data["ms2"].apply(lambda row: ms1_rt_dict.get(row["scan"], row["retention_time"]), axis=1)
    
    # Save ms1 to new file
    deimos.save(out_path+name+'.h5', data['ms1'], key='ms1', mode='w')
    # Save ms2 to same file
    deimos.save(out_path+name+'.h5', data['ms2'], key='ms2', mode='a')


# In[ ]:




