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
import acoular as ac

# Set home dir
from os.path import expanduser
home = expanduser("~")

# Get lists of all the different subfolders into lists

path = home + '/Dropbox (CSU Fullerton)/EGME597_AB' 
fileType = '*.dat'

os.chdir(path)

folderPath = 'RAWDATA/**/'

# Gets a list of all the files in all directories (600 of them on 4-9-19)
dat_fileList = glob.glob(folderPath + fileType, recursive=True)

 

#   Specify constants
Fs          = 48000                        #  DAQ Sample rate (S/sec)
NS          = 48000                        #  Number of samples
fn          = Fs/2                         #  maximum resolvoble frequency
NFFT        = 2**12                         #  4096 point FFT
NF          = NFFT/2                       #  No. point for powerspecturm
NR          = 4                            #  No. of runs
NP          = 4                            #  No. of positions
NM          = 7                            #  No. of microphones
Pref        = 20e-6                        #  Reference pressure
sensitivity = [50.2, 49.3, 53.5, 52.9, 52.8, 48.7, 47.0] # microphone sensitivity
c = 343
MicArrayElements = 7
arraySpacing = 0.00858 # Array spacing (in meters?, maybe inches?)

MicArray = np.linspace(0, MicArrayElements*arraySpacing, MicArrayElements)

#Defining frequency range
frequency = np.linspace(0,NF,2048)
frequency = Fs*frequency/NFFT


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
calibrationArray = calibrationArray.reshape(len(micCalibrationList),frequency.shape[0]).T

# Get all the files without the extension
file_name_list = []

for i in range(0, len(dat_fileList)):
    temp = os.path.splitext(dat_fileList[i])[0]
    temp1 = temp.split('\\')
    temp2 = temp1[2]
    file_name_list.append(temp2)
   
# Import all the data
Pressure = []

for i in range(0, len(dat_fileList)-590):
    data = pd.read_csv(path +'\\'+ dat_fileList[i], delimiter='\s+',
                       header=None)
# ASK THE PROF WHY HE ONLY IS USING 4 MICS HERE
    temp = (1000*data.values[:,1:5])/(np.asarray(sensitivity[1:5])*Pref)  #  convert voltage to pressure and normalize by Pref
    Pressure.append(temp)

time = np.linspace(0,NS,NS)


micGeoFile = 'C:/Users/andre/dev/Aero_ML/mic_array.xml'
mg = ac.MicGeom( from_file=micGeoFile )

# =============================================================================
# This is the MATLAB way to do it.
# =============================================================================
#%   Create array object
#    array      = phased.ULA('NumElements',NE,'ElementSpacing',d);   %  Define array
#    microphone = phased.OmnidirectionalMicrophoneElement('FrequencyRange',...
#                 [20 20e3]);                                        %  Mic frequency response
#    collector  = phased.WidebandCollector('Sensor',array,'SampleRate',48e3,...
#                 'PropagationSpeed',c,'ModulatedInput',false);      %  Array properties
#    sigang     = [-1 -.5 0 .5 1; 0 0 0 0 0];
#    beamformer = phased.TimeDelayBeamformer('SensorArray',array,...
#                 'SampleRate',48e3,'PropagationSpeed',c,'Direction',[0; 0]);
#%
#%   ---------------------------------------------------------------------------
#%
#%   COMPUTE POWER SPECTRUM
#%   --------------------------------------------------------------------------- 
#%   Loop through cases   
#    for i = 1:NP
#        for j = 1:NR+1
#            rsig = collector(V(:,:,j,i),sigang);
#            y = beamformer(rsig)/NE;
#            Y = fft(y,NFFT)/NFFT;
#            P(:,j,i) = abs(Y(1:NF)).^2;
#            S(:,j,i) = 10*log10(P(:,j,i))-mean(R,2);
#            P(:,j,i) = 10.^(S(:,j,i)/10); 
#        end
#        for j = 2:NR+1
#            P(:,j,i) = abs(P(:,j,i) - P(:,1,i));
#        end
#        p(:,i) = 10*log10( mean(P(:,2:end,i),2) );
#        s(:,i) = smooth(f,p(:,i),0.005,'rloess' );
#        B      = 10*log10( P(:,1,i));
#        b(:,i) = smooth(f,B,0.005,'rloess' );
#    end
#%











# 










































toc()
