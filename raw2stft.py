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
import copy
from scipy.io import wavfile
from scipy.signal import butter, lfilter
import scipy.ndimage
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

"""Read in all the text files and hdf5 files"""
def read_raw_microphone_data(folderPath, fileType='.dat'):
    
    fileList = [f for f in os.listdir(folderPath) if f.endswith(fileType)] 
    fileList.sort()
    
    DatList = [np.loadtxt(os.path.join(folderPath, dat_file)) for dat_file in fileList]
    DatList = np.array(DatList)    
    DatListExist = True
    print('Raw .dat files read in')
    return fileList, DatList, DatListExist

def read_hdf5(hdf5fullpath, fileType='.hdf5'):
    '''This takes the hdf5 file that is produced in MATLAB and imports the field (file name)
    Cleaned Acoustic Pressure, Raw Pressure, and the SmSpectrum (smoothed spectrum)
    then plots onto a spectogram and saves as a png'''
    hdfFileList = [f for f in os.listdir(hdf5fullpath) if f.endswith(fileType)]
    df_list = []
    for i in range(0,len(hdfFileList)):
        f = h5py.File(hdf5fullpath+ '\\' + hdfFileList[i], 'r')
        key_list = list(f)
#        print(key_list)
        
        CleanedAcousticPressure = f[key_list[0]][:]
        RawPressure = f[key_list[1]][:]
        SmSpectrum = f[key_list[2]][:]
        field = f[key_list[3]][:]
        field = [i.decode('utf-8') for i in field]
        
        # Now, the next thing on the list is to creat a csv file with all the data and 
        # label it with wall noise, or noise of object. I could make the allocation from
        # an hdf5 to dict general but there should really only be the below fields for 
        # this specific research project.
        data_dict = {key_list[3]: field,
                     key_list[0]: list(RawPressure),
                     key_list[1]: list(CleanedAcousticPressure),
                     key_list[2]: list(SmSpectrum),
                     }
        keyList1 = [key_list[3],key_list[0],key_list[1],key_list[2]]
#        print(keyList1)
        df = pd.DataFrame(data_dict, columns=keyList1)
#        print(df[key_list[0]].shape)
        df_list.append(df)
    return hdfFileList, df_list, key_list
#    
#def datfiles_to_csv():


def folder_helper(dirName, parentPath):
    if not os.path.exists(os.path.join(parentPath, dirName)):
        os.makedirs(os.path.join(parentPath, dirName))
#        print("Directory " , dirName ,  " Created ")
        return dirName
    else:    
#        print("Directory " , dirName ,  " already exists") 
        return dirName
    
    
def saveSpecPlots(Data, MicNum, parentPath, PlotDirectory=None, key_list=None,
                  hdf5plotFolder=None, fileList=None):
    ''' The main reason this function is so complicated looking is it is 
    is equipped to handle the read in hdf5 file or the raw dat files'''
    if isinstance(Data, np.ndarray):
        
       temp =  folder_helper(PlotDirectory, parentPath)
       for i in range(1, len(Data)):
           name = fileList[i][:9]+ '_{}_{}'.format(MicNum, i)
           plotName =  os.getcwd() + '\\'+ PlotDirectory +'\\' + name
           plt.specgram(Data[i, :, MicNum], cmap='jet', Fs=48000)
           plt.title(name)
           plt.xlabel('Time [S]')
           plt.ylabel('Frequency [Hz]')
           plt.savefig(plotName)
           plt.close()
            
    elif  isinstance(Data, list):
        
        for hdfFiles in range(0,len(df_list)):
#            Take the pd.Datafame out of the list
            data = df_list[hdfFiles]
#           Make the new plot directory the key list
        
#           Make the maine and sub directories
            subDIR = folder_helper(hdf5plotFolder, parentPath)
            print(subDIR)
            for folder in key_list:
                if folder == 'field':
                    continue
                else:
                   subsubDir =  folder_helper(folder, subDIR)
                   
                   for ind, row in data.iterrows():
                        name = data['field'][ind]
                        plotName = name[:9] + '_{}'.format(ind)
                        fileName = os.path.join(parentPath, subDIR, subsubDir)+ '\\' + plotName
                        plt.specgram(data[folder][ind], cmap='jet', Fs=48000)
                        plt.title(plotName)
                        plt.xlabel('Time [S]')
                        plt.ylabel('Frequency [Hz]')
                        plt.savefig(fileName)
                        plt.close()
    else:
        raise Exception('saveSpecPlots accepts 3D np.ndarray or a list of pd.core.frame.DataFromes only')
            
        

# =============================================================================
# The next thing to do would be to import all the processed .hdf5 files from MATLAB 
# and sync up the raw first and last microphone (raw) signal with the processed
# (finalized) signal where is is then batch fed into the model with (research this)
# csv/pandas datafrome? or numpy array?
# =============================================================================

## INPUT COMMANDS BELOW
# Raw data into spec

'''Next thing to do is saveSpecPlots for the hdf5 file... need to figure out how to
save cleaned spectrum. maybe put each type into a different folder and then create a 
generate csv with the _1 and _6 with the target in another column.
'''



# Most of the Spectrograms and Inversion are taken from: https://gist.github.com/kastnerkyle/179d6e9a88202ab0a2fe

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def overlap(X, window_size, window_step):
    """
    Create an overlapped version of X
    Parameters
    ----------
    X : ndarray, shape=(n_samples,)
        Input signal to window and overlap
    window_size : int
        Size of windows to take
    window_step : int
        Step size between windows
    Returns
    -------
    X_strided : shape=(n_windows, window_size)
        2D array of overlapped X
    """
    if window_size % 2 != 0:
        raise ValueError("Window size must be even!")
    # Make sure there are an even number of windows before stridetricks
    append = np.zeros((window_size - len(X) % window_size))
    X = np.hstack((X, append))

    ws = window_size
    ss = window_step
    a = X

    valid = len(a) - ws
    nw = (valid) // ss
    out = np.ndarray((int(nw),ws),dtype = a.dtype)

    for i in range(int(nw)):
        # "slide" the window along the samples
        start = i * ss
        stop = start + ws
        out[i] = a[start : stop]

    return out


def stft(X, fftsize=128, step=65, mean_normalize=True, real=False,
         compute_onesided=True):
    """
    Compute STFT for 1D real valued input X
    """
    if real:
        local_fft = np.fft.rfft
        cut = -1
    else:
        local_fft = np.fft.fft
        cut = None
    if compute_onesided:
        cut = fftsize // 2
    if mean_normalize:
        X -= X.mean()
    
    X = overlap(X, fftsize, step)
    
    size = fftsize
    win = 0.54 - .46 * np.cos(2 * np.pi * np.arange(size) / (size - 1))
    X = X * win[None]
    X = local_fft(X)[:, :cut]
    return X

def pretty_spectrogram(d,log = True, thresh= 5, fft_size = 512, step_size = 64):
    """
    creates a spectrogram
    log: take the log of the spectrgram
    thresh: threshold minimum power for log spectrogram
    """
    specgram = np.abs(stft(d, fftsize=fft_size, step=step_size, real=False, compute_onesided=True))
  
    if log == True:
        specgram /= specgram.max() # volume normalize to max 1
        specgram = np.log10(specgram) # take log
        specgram[specgram < -thresh] = -thresh # set anything less than the threshold as the threshold
    else:
        specgram[specgram < thresh] = thresh # set anything less than the threshold as the threshold
    
    return specgram

# Also mostly modified or taken from https://gist.github.com/kastnerkyle/179d6e9a88202ab0a2fe
def invert_pretty_spectrogram(X_s, log = True, fft_size = 512, step_size = 512/4, n_iter = 10):
    
    if log == True:
        X_s = np.power(10, X_s)

    X_s = np.concatenate([X_s, X_s[:, ::-1]], axis=1)
    X_t = iterate_invert_spectrogram(X_s, fft_size, step_size, n_iter=n_iter)
    return X_t

def iterate_invert_spectrogram(X_s, fftsize, step, n_iter=10, verbose=False):
    """
    Under MSR-LA License
    Based on MATLAB implementation from Spectrogram Inversion Toolbox
    References
    ----------
    D. Griffin and J. Lim. Signal estimation from modified
    short-time Fourier transform. IEEE Trans. Acoust. Speech
    Signal Process., 32(2):236-243, 1984.
    Malcolm Slaney, Daniel Naar and Richard F. Lyon. Auditory
    Model Inversion for Sound Separation. Proc. IEEE-ICASSP,
    Adelaide, 1994, II.77-80.
    Xinglei Zhu, G. Beauregard, L. Wyse. Real-Time Signal
    Estimation from Modified Short-Time Fourier Transform
    Magnitude Spectra. IEEE Transactions on Audio Speech and
    Language Processing, 08/2007.
    """
    reg = np.max(X_s) / 1E8
    X_best = copy.deepcopy(X_s)
    for i in range(n_iter):
        if verbose:
            print("Runnning iter %i" % i)
        if i == 0:
            X_t = invert_spectrogram(X_best, step, calculate_offset=True,
                                     set_zero_phase=True)
        else:
            # Calculate offset was False in the MATLAB version
            # but in mine it massively improves the result
            # Possible bug in my impl?
            X_t = invert_spectrogram(X_best, step, calculate_offset=True,
                                     set_zero_phase=False)
        est = stft(X_t, fftsize=fftsize, step=step, compute_onesided=False)
        phase = est / np.maximum(reg, np.abs(est))
        X_best = X_s * phase[:len(X_s)]
    X_t = invert_spectrogram(X_best, step, calculate_offset=True,
                             set_zero_phase=False)
    return np.real(X_t)

def invert_spectrogram(X_s, step, calculate_offset=True, set_zero_phase=True):
    """
    Under MSR-LA License
    Based on MATLAB implementation from Spectrogram Inversion Toolbox
    References
    ----------
    D. Griffin and J. Lim. Signal estimation from modified
    short-time Fourier transform. IEEE Trans. Acoust. Speech
    Signal Process., 32(2):236-243, 1984.
    Malcolm Slaney, Daniel Naar and Richard F. Lyon. Auditory
    Model Inversion for Sound Separation. Proc. IEEE-ICASSP,
    Adelaide, 1994, II.77-80.
    Xinglei Zhu, G. Beauregard, L. Wyse. Real-Time Signal
    Estimation from Modified Short-Time Fourier Transform
    Magnitude Spectra. IEEE Transactions on Audio Speech and
    Language Processing, 08/2007.
    """
    size = int(X_s.shape[1] // 2)
    wave = np.zeros((X_s.shape[0] * step + size))
    # Getting overflow warnings with 32 bit...
    wave = wave.astype('float64')
    total_windowing_sum = np.zeros((X_s.shape[0] * step + size))
    win = 0.54 - .46 * np.cos(2 * np.pi * np.arange(size) / (size - 1))

    est_start = int(size // 2) - 1
    est_end = est_start + size
    for i in range(X_s.shape[0]):
        wave_start = int(step * i)
        wave_end = wave_start + size
        if set_zero_phase:
            spectral_slice = X_s[i].real + 0j
        else:
            # already complex
            spectral_slice = X_s[i]

        # Don't need fftshift due to different impl.
        wave_est = np.real(np.fft.ifft(spectral_slice))[::-1]
        if calculate_offset and i > 0:
            offset_size = size - step
            if offset_size <= 0:
                print("WARNING: Large step size >50\% detected! "
                      "This code works best with high overlap - try "
                      "with 75% or greater")
                offset_size = step
            offset = xcorr_offset(wave[wave_start:wave_start + offset_size],
                                  wave_est[est_start:est_start + offset_size])
        else:
            offset = 0
        wave[wave_start:wave_end] += win * wave_est[
            est_start - offset:est_end - offset]
        total_windowing_sum[wave_start:wave_end] += win
    wave = np.real(wave) / (total_windowing_sum + 1E-6)
    return wave

def xcorr_offset(x1, x2):
    """
    Under MSR-LA License
    Based on MATLAB implementation from Spectrogram Inversion Toolbox
    References
    ----------
    D. Griffin and J. Lim. Signal estimation from modified
    short-time Fourier transform. IEEE Trans. Acoust. Speech
    Signal Process., 32(2):236-243, 1984.
    Malcolm Slaney, Daniel Naar and Richard F. Lyon. Auditory
    Model Inversion for Sound Separation. Proc. IEEE-ICASSP,
    Adelaide, 1994, II.77-80.
    Xinglei Zhu, G. Beauregard, L. Wyse. Real-Time Signal
    Estimation from Modified Short-Time Fourier Transform
    Magnitude Spectra. IEEE Transactions on Audio Speech and
    Language Processing, 08/2007.
    """
    x1 = x1 - x1.mean()
    x2 = x2 - x2.mean()
    frame_size = len(x2)
    half = frame_size // 2
    corrs = np.convolve(x1.astype('float32'), x2[::-1].astype('float32'))
    corrs[:half] = -1E30
    corrs[-half:] = -1E30
    offset = corrs.argmax() - len(x1)
    return offset


### Parameters ###
fft_size = int(2**12) # window size for the FFT
step_size = int(fft_size/16) # distance to slide along the window (in time)
spec_thresh = 3 # threshold for spectrograms (lower filters out more noise)
lowcut = 1000 # Hz # Low cut for our butter bandpass filter
highcut = 48000/2 # Hz # High cut for our butter bandpass filter


sampleRate  = 48000                        #  DAQ Sample rate (S/sec)
NS          = 48000                        #  Number of samples
fn          = sampleRate/2                 #  maximum resolvoble frequency
NFFT        = 2**12                        #  4096 point FFT
NF          = NFFT/2                       #  No. point for powerspecturm
Pref        = 20e-6                        #  Reference pressure
c = 343                                    # Speed of sound 
MicArrayElements = 7
arraySpacing = 0.00858 # Array spacing (cm)

folderPath = home + '\\Dropbox (CSU Fullerton)\\EGME597_AB\\ML_DATA\\RAWDATA'
parentPath = os.getcwd()
fileList, DatList, DatListExist = read_raw_microphone_data(folderPath)
saveSpecPlots(DatList, MicNum=1, PlotDirectory='RawSTFTPlots', parentPath=parentPath,
              fileList=fileList)
#
saveSpecPlots(DatList, MicNum=6, PlotDirectory='RawSTFTPlots', parentPath=parentPath,
              fileList=fileList)

# hdf5 file processing to stft in separate folders
hdf5fullpath = home + '\\Dropbox (CSU Fullerton)\\EGME597_AB\\ML_DATA'

hdfFileList, df_list, key_list = read_hdf5(hdf5fullpath)

saveSpecPlots(df_list, MicNum=1, parentPath=parentPath, key_list=key_list,
              hdf5plotFolder='ProcessedSTFTPlots')
saveSpecPlots(df_list, MicNum=6, parentPath=parentPath, key_list=key_list,
              hdf5plotFolder='ProcessedSTFTPlots')

print('STFT plots generated!')

wav_spectrogram = pretty_spectrogram(DatList[20,:,1], fft_size = fft_size, 
                                     step_size = step_size, log = True, thresh = spec_thresh)

fig, ax = plt.subplots(nrows=1,ncols=1, figsize=(20,4))
cax = ax.matshow(np.transpose(wav_spectrogram), interpolation='nearest', aspect='auto', cmap=plt.cm.jet, origin='lower')
fig.colorbar(cax)
plt.title('Raw Spectrogram')

recovered_audio_orig = invert_pretty_spectrogram(wav_spectrogram, fft_size = fft_size,
                                            step_size = step_size, log = True, n_iter = 10)



toc()


