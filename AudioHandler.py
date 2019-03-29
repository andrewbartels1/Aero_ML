#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'inline')

import csv
import os
import numpy as np
import glob
import matplotlib
import numpy as np
import matplotlib.pyplot as plt  


class AudioHandler():
    
def __init__(self):
    mic_text = 'test.txt'
    
# microphone sensitivity
mic_sensitivity = np.array([50.2, 49.3, 53.5, 52.9, 52.8, 48.7, 47.0])
Pref = 20e-6

def FileSort(self, dir1, dir2):
    '''Sorts the files into the different types of tests
    '''
    self.TestSectionCalibartion_data = glob.glob(dir1, '?.dat')
    self.MomentumWakeSheilding_data =  glob.glob(dir2, '?.dat')

def MicCalibartion(self):
    '''Microphone calibration here
    '''
    
def main():
    TestSectionCalibartion_dir = 'C:/Users/andre/Dropbox (CSU Fullerton)/EGME597_AB/RAWDATA/TSC'
    MomentumWakeSheilding_dir = 'C:/Users/andre/Dropbox (CSU Fullerton)/EGME597_AB/RAWDATA/MWS'
    

if __name__ == "__main__":
    main()
