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
% %    Written on: 7-16-2019
% %    Modified By: Andrew Bartels
% %    Modified on: 7-15-2019
%     Original base code written by: Andrew Bartels
%     contact info: andrewbartels1@gmail.com or ~@csu.fullerton.edu
%
%  Run time of this entire script is roughly:
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
home = 'C:\Users\andrewbartels1';

% Setup file directories
MWS_dir = '\Dropbox (CSU Fullerton)\EGME597_AB\ML_DATA\MWS\';
ACC_dir = '\Dropbox (CSU Fullerton)\EGME597_AB\ML_DATA\ACC\';
BOTH_dir = '\Dropbox (CSU Fullerton)\EGME597_AB\ML_DATA\BOTH\';

calDirectory = '\Dropbox (CSU Fullerton)\EGME597_AB\CALIBRATION\sn442';                  %  calibration directory
 
% If you have a new folder put the directory (MINUS HOME) here.

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
Case = 'ACC';

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

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% DEFINE THE MODEL STRUCTURE BELOW.

% two ecamples are given here, a C-BiLSTM and a simple LSTM.
% -------------------------------------------------------------------------
% DEEP biLSTM
% -------------------------------------------------------------------------
tempLayers = [
    sequenceInputLayer(InputSize,"Name","sequence","Normalization","zerocenter")
    sequenceFoldingLayer("Name","seqfold")];
layers = addLayers(layers,tempLayers);

tempLayers = [
    convolution2dLayer([3 3],3,"Name","conv_1","Padding","same")
    batchNormalizationLayer("Name","batchnorm_1")
    dropoutLayer(0.5,"Name","dropout_1")
    convolution2dLayer([5 5],3,"Name","conv_3","Padding","same")
    batchNormalizationLayer("Name","batchnorm_2")
    eluLayer(1,"Name","elu")];
layers = addLayers(layers,tempLayers);

tempLayers = [
    sequenceUnfoldingLayer("Name","sequnfold")
    flattenLayer("Name","flatten")
    lstmLayer(500,"Name","bilstm1")
    fullyConnectedLayer(numResponses,"Name","fc_1")
    %     fullyConnectedLayer(numResponses,"Name","fc_2")
    dropoutLayer(0.5,"Name","dropout_2")
    fullyConnectedLayer(numResponses,"Name","fc_3")
    regressionLayer("Name","Output")];
layers = addLayers(layers,tempLayers);

layers = connectLayers(layers,"seqfold/out","conv_1");
layers = connectLayers(layers,"seqfold/miniBatchSize","sequnfold/miniBatchSize");
layers = connectLayers(layers,"elu","sequnfold/in");

% -------------------------------------------------------------------------
% Simple LSTM
% -------------------------------------------------------------------------
% tempLayers = [
%     sequenceInputLayer(InputSize,"Name","sequence")
%     flattenLayer("Name","flatten")
%     lstmLayer(2000,"Name","lstm1")
%     dropoutLayer(0.5,"Name","dropout_1")
%     fullyConnectedLayer(250,"Name","fc_1")
%     fullyConnectedLayer(numResponses,"Name","fc_2")
%     regressionLayer("Name","Output")];
%
% layers = addLayers(layers,tempLayers);

% Rename save variable used for preprocessing
ModelSaveFile = SaveFile;


% Run and train the network. If you only want to pre-process the data,
% comment out the line below.
run('LSTM_feature_preprocessing.m')


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
