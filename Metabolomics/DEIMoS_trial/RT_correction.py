#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pyopenms as oms
import os


# In[2]:


os.chdir('/parallel_scratch/mp01950/multiomics/')
out_path = os.getcwd()+"/Brain_trauma/RT_corrected/"
mzML_dir = os.getcwd()+"/Brain_trauma/Original/"
mzML_files = os.listdir(mzML_dir)


# In[3]:


mzML_files


# In[6]:


for file in mzML_files:
    # Load the mzML file
    exp = oms.MSExperiment()
    oms.MzMLFile().load("Brain_trauma/Original/"+file, exp)
    # Sort scans into MS1 and MS2
    ms1_scans = [s for s in exp if s.getMSLevel() == 1]
    ms2_scans = [s for s in exp if s.getMSLevel() == 2]
    # Find scans that have both MS1 and MS2 pairs
    num_swaps = 0
    min_length = min(len(ms1_scans), len(ms2_scans))  # Only iterate over pairs
    for i in range(min_length):
        # Ensure they have the same retention time before swapping
        rt1, rt2 = ms1_scans[i].getRT(), ms2_scans[i].getRT()
        ms1_scans[i].setRT(rt2)
        ms2_scans[i].setRT(rt1)
        num_swaps += 1
    # Create a new experiment and add spectra back in the original order
    corrected_exp = oms.MSExperiment()
    for spec in exp:  # Iterate over the original experiment spectra
        if spec.getMSLevel() == 1 and ms1_scans:
            corrected_exp.addSpectrum(ms1_scans.pop(0))  # Add modified MS1 spectrum
        elif spec.getMSLevel() == 2 and ms2_scans:
            corrected_exp.addSpectrum(ms2_scans.pop(0))  # Add modified MS2 spectrum
        else:
            corrected_exp.addSpectrum(spec)  # Add unchanged spectrum (if any)
    # Save the corrected mzML file
    name = file.split(".mzML")[0]  # Splitting at ".mzML" and taking the first part
    oms.MzMLFile().store(out_path + name + "rc.mzML", corrected_exp)


# In[47]:





# In[ ]:




