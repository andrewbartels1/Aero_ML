%   ----------------------------------------------------------------------------
%   This script a training harness can be made to produce everything needed start to
%   finish.
%
%  Inputs:
%  1) File directory of microphone .dat files (must have 7 mics and be in the order of:
%     ...R0.dat
%     ...R1.dat
%     ...R2.dat
%     ...R3.dat
%     ...R4.dat
%  such that R0 is the background noise of the test section and R1-4 are
%  different speeds etc. with the noise or point source on.
%  Outputs:
%
%
%  DEPENDENCIES:
%  MATLAB 2019a (mainly for MelSpectrogram in the LSTM_feature_preprocessing.m)
%  glob.m
%  MATLAB dependancies:
%  audio_system_toolbox
%  distrib_computing_toolbox
%  image_toolbox
%  matlab 2019a
%  neural_network_toolbox
%  phased_array_system_toolbox
%  signal_blocks
%  signal_toolbox
%  statistics_toolbox
% %    Written By: Andrew Bartels
% %    Written on: 8-1-2019
% %    Modified By: Andrew Bartels
% %    Modified on: 8-27-2019
%     Original base code written by: Andrew Bartels
%     contact info: andrewbartels1@gmail.com or ~@csu.fullerton.edu
%
%
%   ----------------------------------------------------------------------------
%
%   INITIALIZE CODE
%   ----------------------------------------------------------------------------
clc
clear all
close all
tic


%-------------------------------------------------------------------------------
% Setup all the save files and input directories
%-------------------------------------------------------------------------------

% Put the machine's home dir here!
home = pwd;

% if this is pulled right from the git without data this should be true or errors will happen
fromGit = true; 

% Setup file directories
MWS_dir = '/RAW/MWS/';
ACC_dir = '/RAW/ACC/';
BOTH_dir = '/RAW/BOTH/';

calDirectory = '/CALIBRATION/sn442';                  %  calibration directory
 
% If you have a new folder put the directory of data to feed the LSTM, put
% it here.

MWS_dir = strcat(home,MWS_dir);
ACC_dir = strcat(home,ACC_dir);
BOTH_dir = strcat(home,BOTH_dir);
Calibration_dir = strcat(home,calDirectory);

% Setup preprocessing save file names
PreprocessingSaveFile = 'CleanedDataForNN_';

% Setup model save files
SaveFile = 'simple_LSTM_justmelFreq_coeffs2';



%----------------------------------------------------------------------------
% Preprocess and train on only the old MWS data
%----------------------------------------------------------------------------

% Set Case to use here:
Case = 'BOTH';

switch(Case)
    case 'MWS'
        SaveFile = strcat(SaveFile,'_',Case);
        PreprocessingSaveFile = strcat(PreprocessingSaveFile,'_',Case);
        disp('MWS being used')
    case 'ACC'
        SaveFile = strcat(SaveFile,'_',Case);
        PreprocessingSaveFile = strcat(PreprocessingSaveFile,'_',Case);
        disp('ACC being used')
    case 'BOTH'
        SaveFile = strcat(SaveFile,'_',Case);
        PreprocessingSaveFile = strcat(PreprocessingSaveFile,'_',Case);
        disp('Both being used')
    otherwise
        error('MWS, ACC, and BOTH are the only cases available')
end

% Preprocess all the old data (of the same test)
% ONCE ALL THE SAVEFILES ARE DONE COMMENT THIS OUT!!!!
% THIS SHOULD ONLY BE RUN A COUPLE OF TIMES. IT TAKES A LONG TIME

run('LSTM_preprocessing.m')




% Rename save variable used for preprocessing
ModelSaveFile = SaveFile;


% Run and train the network. If you only want to pre-process the data,
% comment out the line below.
run('LSTM_model.m')


clear Case

% MATLAB dependancies:
% audio_system_toolbox
% distrib_computing_toolbox
% image_toolbox
% matlab 2019a
% neural_network_toolbox
% phased_array_system_toolbox
% signal_blocks
% signal_toolbox
% statistics_toolbox
