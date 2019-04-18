# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 14:30:24 2019

@author: andre
"""
# =============================================================================
# This script is meant to preprocess all the data and label it to be fed into 
# the neural network.
# =============================================================================

import sys
#Adds dev_util git files found useful into python path (temporarily)
sys.path.append("~\\dev\\dev_utils")
from tic_toc import tic,toc
tic()

import os
import h5py
import glob
import numpy as np
import pandas as pd
import scipy.interpolate as spinterp
from scipy import signal
from sklearn.decomposition import FastICA, PCA

# Set home dir
from os.path import expanduser
home = expanduser("~")


filename = 'TSC202.hdf5'
folderPath = 'C:/Users/andre/Dropbox (CSU Fullerton)/EGME597_AB/ML_DATA/'
fileType = '*.hdf5'


#    fileList = glob.glob(folderPath + fileType)
file = folderPath + filename
f = h5py.File(file)

key_list = list(f)

RawPressure = f[key_list[0]][:]
CleanedAcousticPressure = f[key_list[1]][:]
SmSpectrum = f[key_list[2]][:]

#take file names, decode and put into list
field = f[key_list[3]][:]
field = [i.decode('utf-8') for i in field]

data = [RawPressure, CleanedAcousticPressure, SmSpectrum]


# =============================================================================
# list[dicts_of_all_hdf5_files]
#
#dicts_of_all_hdf5_files = {
#key(filename): list_of_nupy_arrays[ Raw Pressure, cleaned acoustic noise...]
#}
# =============================================================================

#data_dict = dict(zip(field, data))








#print('data labeled and processed')