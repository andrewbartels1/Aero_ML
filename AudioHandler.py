# -*- coding: utf-8 -*-



get_ipython().run_line_magic('matplotlib', 'inline')

import csv
import os
import glob
import numpy as np
#import matplotlib.pyplot as plt  
#import tensorflow as tf
#import keras as ks


class AudioHandler():
    
    def __init__(self, dir1, dir2):
        self.dir1 = dir1
        self.dir2 = dir1
        self.mic_text = 'test.txt'
        self.Fs          = 48000;                        #  DAQ Sample rate (S/sec)
        self.NS          = 48000;                        #  Number of samples
        self.fn          = self.Fs/2;                         #  maximum resolvoble frequency
        self.NFFT        = 2**12;                        #  4096 point FFT
        self.NF          = self.NFFT/2;                       #  No. point for powerspecturm
        self.NR          = 4;                            #  No. of runs
        self.NP          = 4;                            #  No. of positions
        self.NM          = 7;                            #  No. of microphones
        self.Pref        = 20e-6;                        #  Reference pressure
        self.sensitivity = [50.2, 49.3, 53.5, 52.9, 52.8, 48.7, 47.0] # microphone sensitivity

    # microphone sensitivity
    mic_sensitivity = np.array([50.2, 49.3, 53.5, 52.9, 52.8, 48.7, 47.0])
    Pref = 20e-6
    
    def FileSort(self):
        '''Sorts the files into the different types of tests
        '''
        self.TestSectionCalibartion_data = glob.glob(self.dir1, '?.dat')
        self.MomentumWakeSheilding_data =  glob.glob(self.dir2, '?.dat')
    
    def MicCalibartion(self):
        '''Microphone calibration here
        '''
    
def main():
    path = 'C:\\Users\\andrewbartels1\\Dropbox (CSU Fullerton)\\EGME597_AB'
    os.chdir(path)
    TestSectionCalibartion_dir = 'EGME597_AB/RAWDATA/TSC'
    MomentumWakeSheilding_dir = 'EGME597_AB/RAWDATA/MWS'
#    calibration_dir = 'EGME597_AB/CALIBRATION/sn442'
    AH = AudioHandler(TestSectionCalibartion_dir, MomentumWakeSheilding_dir)
    AH.FileSort()

if __name__ == "__main__":
    main()



