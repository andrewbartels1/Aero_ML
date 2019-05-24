# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 11:25:21 2019

@author: andrewbartels1
"""
import sys
#Adds dev_util git files found useful into python path (temporarily)
sys.path.append("~\\dev\\dev_utils")
from tic_toc import tic,toc
tic()

import os
import glob
import numpy as np
import pandas as pd
import acoular as ac
import scipy.interpolate as spinterp
from scipy import signal
from sklearn.decomposition import FastICA, PCA
import matplotlib.pyplot as plt

# Set home dir
from os.path import expanduser
home = expanduser("~")

# Get lists of all the different subfolders into lists

path = home + '/Dropbox (CSU Fullerton)/EGME597_AB' 
fileType = '*.dat'

os.chdir(path)

folderPath = 'RAWDATA/MWS/'

# Gets a list of all the files in all directories (600 of them on 4-9-19)
dat_fileList = glob.glob(folderPath + fileType, recursive=True)

 

#   Specify constants
sampleRate  = 48000                        #  DAQ Sample rate (S/sec)
NS          = 48000                        #  Number of samples
fn          = sampleRate/2                         #  maximum resolvoble frequency
NFFT        = 2**12                         #  4096 point FFT
NF          = NFFT/2                       #  No. point for powerspecturm
NR          = 4                            #  No. of runs
NP          = 4                            #  No. of positions
NM          = 7                            #  No. of microphones
Pref        = 20e-6                        #  Reference pressure
sensitivity = [50.2, 49.3, 53.5, 52.9, 52.8, 48.7, 47.0] # microphone sensitivity
c = 343
MicArrayElements = 7
arraySpacing = 0.00858 # Array spacing (cm)

#Physical locations of each X,Y of each mic
MicArray = np.linspace(0, MicArrayElements*arraySpacing, MicArrayElements)

#Defining frequency range
frequency = np.linspace(0,NF,2048)
frequency = sampleRate*frequency/NFFT


# Setup and read in all the microphone sensitivity text files
calDirectory = home +  '/Dropbox (CSU Fullerton)/EGME597_AB/CALIBRATION/' 
micCalibrationList = glob.glob(calDirectory + '*.txt')   

for i in range(0, len(dat_fileList)):
    file_names = [dat_fileList[i].split('/')]

calibrationArray = np.array([])


# Take out all the mic calibration data
for i in range(0,NM):
    temp = np.loadtxt(micCalibrationList[i])
    fm = temp[:,0]
    Sm = temp[:,1]
    micResponses =  spinterp.pchip_interpolate(fm,Sm,frequency) #Same as pchip interp1 MATLAB
    calibrationArray = np.append(calibrationArray, micResponses, axis=0)
calibrationArray = calibrationArray.reshape(len(sensitivity),frequency.shape[0]).T

# Get all the files without the extension
file_name_list = []

for i in range(0, len(dat_fileList)):
    temp = os.path.splitext(dat_fileList[i])[0]
    temp1 = temp.split('/')
    temp2 = temp1[1].split('\\')
    temp2 = temp2[1]
    file_name_list.append(temp2)
   
# Import all the data
Pressure = []

for i in range(0, len(dat_fileList)):
    data = pd.read_csv(path +'\\'+ dat_fileList[i], delimiter='\s+',
                       header=None)
    temp = (1000*data.values)/(np.asarray(sensitivity)*float(Pref))  #  convert voltage to pressure and normalize by Pref
    Pressure.append(temp)

# Compute ICA
ica = FastICA(n_components=1)
S_ = ica.fit_transform(Pressure[0])  # Reconstruct signals
A_ = ica.mixing_  # Get estimated mixing matrix


plt.figure()
models = [Pressure[0], S_]
names = ['Observations (mixed signal)',
         'True Sources',
         'ICA recovered signals']
colors = ['red', 'steelblue', 'orange', 'green', 'blue', 'yellow', 'orange']

for ii, (model, name) in enumerate(zip(models, names), 1):
    plt.subplot(8, 1, ii)
    plt.title(name)
    for sig, color in zip(model.T, colors):
        plt.plot(sig, color=color)

plt.subplots_adjust(0.09, 0.04, 1.00, 0.94, 0.26, 0.46)
plt.show()


toc()
