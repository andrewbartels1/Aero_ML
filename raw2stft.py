# -*- coding: utf-8 -*-
"""
Created on Fri May 24 10:40:34 2019

@author: andre
"""

import sys
import os
#Adds dev_util git files found useful into python path (temporarily)
sys.path.append("~\\dev\\dev_utils")
from tic_toc import tic,toc
import numpy as np
from matplotlib import pyplot as plt
import h5py
import glob
from numpy.lib import stride_tricks
import pandas as pd
from os.path import expanduser
home = expanduser("~")

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
#    path = home + '/Dropbox (CSU Fullerton)/EGME597_AB' 
#    os.chdir(path)
#    
    return fileList, DatList

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
    plt.xticks(xlocs, ["%.02f" % l for l in ((xlocs*len(samples)/timebins)+(0.5*binsize))/samplerate])
    ylocs = np.int16(np.round(np.linspace(0, freqbins-1, 10)))
    plt.yticks(ylocs, ["%.02f" % freq[i] for i in ylocs])
    
    if plotpath:
        plt.savefig(plotpath, bbox_inches="tight")
    else:
        plt.show()
        
    plt.clf()
    
folderPath = home + '\\Dropbox (CSU Fullerton)\\EGME597_AB\\RAWDATA\\MWS'
fileList, DatList = read_raw_microphone_data(folderPath)

plotstft(DatList[0][:,1])
plt.specgram(DatList[0][:,1], cmap='jet', Fs=48000)
toc()

# =============================================================================
# The next thing to do would be to import all the processed .hdf5 files from MATLAB 
# and sync up the raw first and last microphone (raw) signal with the processed
# (finalized) signal where is is then batch fed into the model with (research this)
# csv/pandas datafrome? or numpy array?
# =============================================================================



