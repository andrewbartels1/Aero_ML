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
                     key_list[0]: list(CleanedAcousticPressure),
                     key_list[1]: list(RawPressure),
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
        print("Directory " , dirName ,  " Created ")
    else:    
        print("Directory " , dirName ,  " already exists") 
    
def saveSpecPlots(Data, MicNum, parentPath, PlotDirectory=None, key_list=None,
                  hdf5plots=None, fileList=None):
    ''' The main reason this function is so complicated looking is it is 
    is equipped to handle the read in hdf5 file or the '''
    if isinstance(Data, np.ndarray):
        folder_helper(PlotDirectory, parentPath)
        for i in range(1, len(Data)):
            name = fileList[i][:9]+ '_{}_{}'.format(MicNum, i)
            plotName =  os.getcwd() + '\\'+ PlotDirectory +'\\' + name
            plt.specgram(Data[i, :, MicNum], cmap='jet', Fs=48000)
            plt.title(name)
            plt.savefig(plotName)
            plt.close()
    elif  isinstance(Data, pd.core.frame.DataFrame):
        PlotDirectory = key_list;
        folder_helper(hdf5plots, parentPath)
        subdir = parentPath +'\\' + hdf5plots
        for folder in range(0, len(PlotDirectory)):
            folder_helper(PlotDirectory[folder], subdir)
        for j in range(0, len(folders)):
            for i in range(0, len(Data)):
                name = fileList[i][:9]+ '_{}_{}'.format(MicNum, i)
                plotName =  os.getcwd() + '\\'+ j +'\\' + name
                plt.specgram(Data[i, :, MicNum], cmap='jet', Fs=48000)
                plt.title(name)
                plt.savefig(plotName)
                plt.close()
    else:
            raise Exception('saveSpecPlots accepts 3D np.ndarray or pd.core.frame.DataFromes only')
            
        

# =============================================================================
# The next thing to do would be to import all the processed .hdf5 files from MATLAB 
# and sync up the raw first and last microphone (raw) signal with the processed
# (finalized) signal where is is then batch fed into the model with (research this)
# csv/pandas datafrome? or numpy array?
# =============================================================================

## INPUT COMMANDS BELOW
# Raw data into spec
folderPath = home + '\\Dropbox (CSU Fullerton)\\EGME597_AB\\ML_DATA\\RAWDATA'

fileList, DatList, DatListExist = read_raw_microphone_data(folderPath)
saveSpecPlots(DatList, MicNum=6, PlotDirectory='RawSTFTPlots', parentPath=folderPath,
              fileList=fileList)
saveSpecPlots(DatList, MicNum=1, PlotDirectory='RawSTFTPlots', parentPath=folderPath,
              fileList=fileList)

# hdf5 file processing to stft in separate folders
hdf5fullpath = home + '\\Dropbox (CSU Fullerton)\\EGME597_AB\\ML_DATA'
parentPath = os.getcwd()
hdfFileList, df_list, key_list = read_hdf5(hdf5fullpath)
data = df_list[0]
saveSpecPlots(data, MicNum=1, parentPath=parentPath, key_list=key_list,
              hdf5plots='ProcessedSTFTPlots')

print('STFT plots generated!')
toc()
'''Next thing to do is saveSpecPlots for the hdf5 file... need to figure out how to
save cleaned spectrum. maybe put each type into a different folder and then create a 
generate csv with the _1 and _6 with the target in another column.
'''














