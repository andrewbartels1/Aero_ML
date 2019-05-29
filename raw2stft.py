# -*- coding: utf-8 -*-
"""
Created on Fri May 24 10:40:34 2019

@author: andre
"""

import os
import sys
import h5py
import numpy as np
import pandas as pd

from numpy.lib import stride_tricks
from matplotlib import pyplot as plt

# Do some path and import setup
from os.path import expanduser
home = expanduser("~")
sys.path.append("~\\dev\\dev_utils")
from tic_toc import tic,toc
# =============================================================================
# The main purpose of this script is to generate a csv file and data directory 
# with the 2 raw signals and then the 'answer' or bramformed and processed signal
# from the hdf5 files processed in MATLAB
# =============================================================================
tic()
""" short time fourier transform of audio signal """
def stft(sig, frameSize, overlapFac=0.5, window=np.hanning):
    win = window(frameSize)
    hopSize = int(frameSize - np.floor(overlapFac * frameSize))
    
    # zeros at beginning (thus center of 1st window should be for sample nr. 0)
    samples = np.append(np.zeros(int(np.floor(frameSize/2.0))), sig)    
    
    # cols for windowing
    cols = np.ceil((len(samples) - frameSize) / float(hopSize)) + 1
    
    # zeros at end (thus samples can be fully covered by frames)
    samples = np.append(samples, np.zeros(frameSize))
    frames = stride_tricks.as_strided(samples, shape=(int(cols), frameSize), strides=(samples.strides[0]*hopSize, samples.strides[0])).copy()
    frames *= win
    
    return np.fft.rfft(frames)  

""" scale frequency axis logarithmically """    
def logscale_spec(spec, sr=48000, factor=20.):
    
    timebins, freqbins = np.shape(spec)
    scale = np.linspace(0, 1, freqbins) ** factor
    scale *= (freqbins-1)/max(scale)
    scale = np.unique(np.round(scale))
    
    # create spectrogram with new freq bins
    newspec = np.complex128(np.zeros([timebins, len(scale)]))
    for i in range(0, len(scale)):
        if i == len(scale)-1:
            newspec[:,i] = np.sum(spec[:,int(scale[i]):], axis=1)
        else:
            newspec[:,i] = np.sum(spec[:,int(scale[i]):int(scale[i+1])], axis=1)
    
    # list center freq of bins
    allfreqs = np.abs(np.fft.fftfreq(freqbins*2, 1./sr)[:freqbins+1])
    freqs = []
    for i in range(0, len(scale)):
        if i == len(scale)-1:
            freqs += [np.mean(allfreqs[int(scale[i]):])]
        else:
            freqs += [np.mean(allfreqs[int(scale[i]):int(scale[i+1])])]
    
    return newspec, freqs

"""Read in all the text files and place them in a csv"""
def read_raw_microphone_data(folderPath, fileType='.dat'):
    fileList = [f for f in os.listdir(folderPath) if f.endswith(fileType)] 
    fileList.sort()
    
    DatList = [np.loadtxt(os.path.join(folderPath, dat_file)) for dat_file in fileList]
    DatList = np.array(DatList)    
    DatListExist = True
    print('Raw .dat files read in')
    return fileList, DatList, DatListExist

#    
#def datfiles_to_csv():

""" plot spectrogram"""
def plotstft(rawdata, samplerate=48000, binsize=2**10, plotpath=None, colormap="jet"): 
    
    s = stft(rawdata, binsize)
    
    sshow, freq = logscale_spec(s, factor=1.0, sr=samplerate)
    ims = 20.*np.log10(np.abs(sshow)/10e-6) # amplitude to decibel
    
    timebins, freqbins = np.shape(ims)
    
    plt.figure(figsize=(15, 7.5))
    plt.imshow(np.transpose(ims), origin="lower", aspect="auto", cmap=colormap, interpolation="none")
    plt.colorbar()

    plt.xlabel("time (s)")
    plt.ylabel("frequency (hz)")
    plt.xlim([0, timebins-1])
    plt.ylim([0, freqbins])

    xlocs = np.float32(np.linspace(0, timebins-1, 5))
    plt.xticks(xlocs, ["%.02f" % l for l in ((xlocs*len(rawdata)/timebins)+(0.5*binsize))/samplerate])
    ylocs = np.int16(np.round(np.linspace(0, freqbins-1, 10)))
    plt.yticks(ylocs, ["%.02f" % freq[i] for i in ylocs])
    
    if plotpath:
        plt.savefig(plotpath, bbox_inches="tight", dpi=1200)
    else:
        plt.show()
        
    plt.clf()

def saveSpecPlots(DatListExist, MicNum,PlotDirectory='RawSTFTPlots'):
    if DatListExist == True:
        if not os.path.exists(PlotDirectory):
            os.mkdir(PlotDirectory)
            print("Directory " , PlotDirectory ,  " Created ")
    
        else:    
            print("Directory " , PlotDirectory ,  " already exists")
        
    for i in range(0, len(DatList)):
        name = fileList[i][:9]+ '_{}_{}'.format(MicNum, i)
        plotName =  os.getcwd() + '\\'+ PlotDirectory +'\\' + name
        print(plotName)
        plt.specgram(DatList[i, :, MicNum], cmap='jet', Fs=48000)
        plt.title(name)
        plt.savefig(plotName)
        plt.close()
# =============================================================================
# The next thing to do would be to import all the processed .hdf5 files from MATLAB 
# and sync up the raw first and last microphone (raw) signal with the processed
# (finalized) signal where is is then batch fed into the model with (research this)
# csv/pandas datafrome? or numpy array?
# =============================================================================

def hdf5_to_spec(hdf5fullpath, fileType='.hdf5'):
    '''This takes the hdf5 file that is produced in MATLAB and imports the field (file name)
    Cleaned Acoustic Pressure, Raw Pressure, and the SmSpectrum (smoothed spectrum)
    then plots onto a spectogram and saves as a png'''
    hdfFileList = [f for f in os.listdir(hdf5fullpath) if f.endswith(fileType)]
    df_list = []
    for i in range(0,len(hdfFileList)):
        print(i)
        print(hdfFileList[i])
        f = h5py.File(hdf5fullpath+ '\\' + hdfFileList[i], 'r')
        key_list = list(f)
        print(key_list)
        
        CleanedAcousticPressure = f[key_list[0]][:]
        RawPressure = f[key_list[1]][:]
        SmSpectrum = f[key_list[2]][:]
        field = f[key_list[3]][:]
        field = [i.decode('utf-8') for i in field]
        
        # Now, the next thing on the list is to creat a csv file with all the data and 
        # label it with wall noise, or noise of object.
        data_dict = {key_list[3]: field,
                     key_list[0]: list(CleanedAcousticPressure),
                     key_list[1]: list(RawPressure),
                     key_list[2]: list(SmSpectrum),
                     }
        keyList1 = [key_list[3],key_list[0],key_list[1],key_list[2]]
        print(keyList1)
        df = pd.DataFrame(data_dict, columns=keyList1)
        print(df[key_list[0]].shape)
        df_list.append(df)
    return hdfFileList, df_list



## INPUT COMMANDS BELOW
folderPath = home + '\\Dropbox (CSU Fullerton)\\EGME597_AB\\ML_DATA\\RAWDATA'

hdf5fullpath = 'C:\\Users\\andre\\Dropbox (CSU Fullerton)\\EGME597_AB\\ML_DATA'

fileList, DatList, DatListExist = read_raw_microphone_data(folderPath)
saveSpecPlots(DatListExist, MicNum=1)
hdfFileList, df_list = hdf5_to_spec(hdf5fullpath)

toc()















