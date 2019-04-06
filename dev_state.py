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
import csv
import os
import glob
import numpy as np
import matplotlib.pyplot as plt  
import tensorflow as tf
import keras as ks



# Get lists of all the different subfolders into lists
path ='C:\\Users\\andrewbartels1\\Dropbox (CSU Fullerton)\\EGME597_AB' 
fileType = '*.dat'

os.chdir(path)

folderPath = 'RAWDATA\\**\\'

dat_fileList = glob.glob(folderPath + fileType, recursive=True)

#   Specify constants
Fs          = 48000                        #  DAQ Sample rate (S/sec)
NS          = 48000                        #  Number of samples
fn          = Fs/2                         #  maximum resolvoble frequency
NFFT        = 2^12                         #  4096 point FFT
NF          = NFFT/2                       #  No. point for powerspecturm
NR          = 4                            #  No. of runs
NP          = 4                            #  No. of positions
NM          = 7                            #  No. of microphones
Pref        = 20e-6                        #  Reference pressure
sensitivity = [50.2, 49.3, 53.5, 52.9, 52.8, 48.7, 47.0] # microphone sensitivity
c = 343
MicArrayElements = 7



calDirectory = 'C:\\Users\\andrewbartels1\\Dropbox (CSU Fullerton)\\EGME597_AB\\CALIBRATION\\' 
micCalibrationList = glob.glob(calDirectory + '*.txt')   

print(len(micCalibrationList))

#f = load(, delimiter=';')












































toc()
