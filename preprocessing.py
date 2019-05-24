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
import matplotlib
import glob
import numpy as np
import pandas as pd
import scipy as sp
import scipy.interpolate as spinterp
from sklearn import preprocessing 
from scipy import signal
from sklearn.decomposition import FastICA, PCA
import matplotlib.pyplot as plt

# Set home dir
from os.path import expanduser
home = expanduser("~")

#   Specify constants
sampleRate  = 48000                        #  DAQ Sample rate (S/sec)
NS          = 48000                        #  Number of samples
fn          = sampleRate/2                 #  maximum resolvoble frequency
NFFT        = 2**12                        #  4096 point FFT
NF          = NFFT/2                       #  No. point for powerspecturm
NR          = 4                            #  No. of runs
NP          = 4                            #  No. of positions
NM          = 7                            #  No. of microphones
Pref        = 20e-6                        #  Reference pressure
sensitivity = [50.2, 49.3, 53.5, 52.9, 52.8, 48.7, 47.0] # microphone sensitivity
c = 343
MicArrayElements = 7
arraySpacing = 0.00858 # Array spacing (cm)

filename = 'MWS203.hdf5'
folderPath = home + '/Dropbox (CSU Fullerton)/EGME597_AB/ML_DATA/'
fileType = '*.hdf5'


#    fileList = glob.glob(folderPath + fileType)
file = folderPath + filename
f = h5py.File(file, 'r')

key_list = list(f)

CleanedAcousticPressure = f[key_list[0]][:]
RawPressure = f[key_list[1]][:]
SmSpectrum = f[key_list[2]][:]
field = f[key_list[3]][:]
field = [i.decode('utf-8') for i in field]

# Now, the next thing on the list is to creat a csv file with all the data and 
# label it with wall noise, or noise of object.

             key_list[0]: list(CleanedAcousticPressure),
             'noiseFlag': np.ones(CleanedAcousticPressure.shape[0])
        }

df = pd.DataFrame(data_dict, columns=key_list)

print(df[key_list[0]].shape)

# Next thing on the list... BSS on each of the different Cleaned pressures and 
# Labeling all the points. First thing that needs to be done is the 

def normalize_complex_arr(a):
    a_oo = a - a.real.min() - 1j*a.imag.min() # origin offsetted
    return a_oo/np.abs(a_oo).max()


#Normalize audio data
Std_CleanedAcousticPressure = preprocessing.normalize(CleanedAcousticPressure)

# Separate out all signals from the raw microphone pressures


# =============================================================================
# list[dicts_of_all_hdf5_files]
#
#dicts_of_all_hdf5_files = {
#key(filename): list_of_nupy_arrays[ Raw Pressure, cleaned acoustic noise...]
#}
# =============================================================================

#data_dict = dict(zip(field, data))

# Plots spectogram with lots of overlap 
fig, ax = plt.subplots(nrows=1)
#plt.specgram(CleanedAcousticPressure[0],Fs=sampleRate, noverlap=200)
Pxx, freqs, bins, im = ax.specgram(CleanedAcousticPressure[0], Fs=sampleRate, noverlap=200)
plt.xlabel('Samples')
plt.ylabel('Normalized Frequency')
plt.show()

#print('data labeled and processed')
toc()