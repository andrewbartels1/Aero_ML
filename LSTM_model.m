%   ----------------------------------------------------------------------------
%   This script computes all the features needed to be fed into the LSTM network
%   and saves it in another .mat file (the reason for multiple .mat files
%   are so the development process is faster and takes less time to
%   compute. Once the pipeline is all setup, all the scripts can be combine
%   or a training harness can be made to produce everything needed start to
%   finish.
%   
%  Input file:
%       1) .mat file with 2 mic and 7 mic processed signals.
%
% %    Written By: Andrew Bartels
% %    Written on: 1-1-2019
% %    Modified By: Andrew Bartels
% %    Modified on: 7-15-2019
%     Original base code written by: Salvador Mayoral
%   Run time of this script is roughly 380 seconds or 6ish minutes
%   ----------------------------------------------------------------------------
%
%   INITIALIZE CODE
%   ----------------------------------------------------------------------------
clc
close all

tic

% Load in the precalculated 2 and 7 mic signals
load(PreprocessingSaveFile)

disp('Preprocessing being used:')
PreprocessingSaveFile

% Produce time stamp here and insert into savefile
datestamp = char(datetime('now'));
date = datestamp(1:11);
time = strcat(datestamp(13:14),'_',datestamp(16:17),'_',datestamp(19:end));

ModelSaveFile = strcat(strcat(date,'_',time), ModelSaveFile);

% take the real out of smooth so MFCCs can deal with it.
SPL_smooth1 = real(SPL_smooth);

% Plot some mel Specs for pretty pictures...
randNum = randi(size(SPL_smooth1,2),1);

for i=1:size(Data,3)    
[f0(:,i), idx(:,i)] = pitch(Data(:,1,10), SampleRate);
end

disp('pitch calculated');

FigHandle_01 = figure('Position', [100, 150, 350, 290]);
subplot(3,1,1)
title('raw signal input')
plot(Data(:,1,randNum))
title('Raw Microphone Data')
ylabel('Amplitude')

subplot(3,1,2)
title('pitch plot')
plot(idx(:,randNum),f0(:,randNum))
ylabel('Pitch (Hz)')
xlabel('Sample Number')

subplot(3,1,3)
semilogx(f,SPL_smooth1(:,randNum))
xlabel('freq')
ylabel('SPL (dB)')
title('SPL vs. Freq')
saveas(FigHandle_01,'2_vs_7 MelSpec','jpeg')


%   ----------------------------------------------------------------------------
%
%   Mel Frequency Cepstral Coefficients (MFCC) 
%   ---------------------------------------------------------------------------
% Extract features going into the network

% TODO: Split up and train models the over lap. i.e. have one model predict
% 1:150, 125:175 etc. and use the compilation of models to see what is
% best. OR look up better model architectures like a super-resolution model
% with repeating and skipping layers.

% Also revise paper and see about publishing at a converence etc after much
% better result.

tic
cepFeatures  = cepstralFeatureExtractor('InputDomain', 'Frequency','SampleRate', SampleRate,...
                                          'NumCoeffs', 40);
for i=1:size(Data,3)
[coeffs_7(:, i),delta7(:,i),deltaDelta7(:,i)] = cepFeatures(SPL_smooth1(:,i));
[coeffs_2(:, i),delta2(:,i),deltaDelta2(:,i)] = cepFeatures(SPL1(:,i));                                                                 

[MelSpecs_7(:, i), centFreq7(:,i)] = melSpectrogram(SPL_smooth1(:,i),SampleRate, ...
                                     'WindowLength',2048,...
                                     'OverlapLength',1024, ...
                                     'FFTLength',4096, ...
                                     'NumBands',41, ...
                                     'FrequencyRange',[62.5,fn]);
[MelSpecs_2(:, i), centFreq2(:,i)] = melSpectrogram(SPL1(:,i),SampleRate, ...
                                     'WindowLength',2048,...
                                     'OverlapLength',1024, ...
                                     'FFTLength',4096, ...
                                     'NumBands',41, ...
                                     'FrequencyRange',[62.5,fn]);
% Spectrograms(:,:,i) = 
FigHandle_01 = figure('Position', get(0, 'Screensize'),'Visible','off');
SpecCellArray{1,i} = real(spectrogram(SPL1(:,i))); %images, just in arrays
end
toc


% 

FigHandle_02 = figure('Position', [100, 150, 350, 290]);
subplot(2,1,1)
melSpectrogram(SPL1(:,randNum), SampleRate)
title('2 mic MFCC')
subplot(2,1,2)
melSpectrogram(SPL_smooth1(:,randNum), SampleRate)
title('7 mic MFCC')
saveas(FigHandle_02,'2_vs_7 MelSpec','jpeg')

disp('MCFFS calculated')



%% Get all the needed data into GPU and setup the network

% Check the status of the onboard GPU
GPU = gpuDevice(1);
reset(GPU);

% Get all of the extracted features and put it in an NxM matrix inside a
% cell. Such that N are columns of features and M are the total number of
% training files.
for i=1:size(Data,3)
sequences{i,1}(:,1) = MelSpecs_2(:,i);
sequences{i,1}(:,2) = coeffs_2(:,i);

% Target is the cell array of 7 beamformed microphones
target{i,1} = SPL_smooth1(:,i);
end

% Do any data augmentation to audio stuff here.
% Possibly try this if spectrogram images are being used.

% Setup the train/val/test split!
% cv = cvpartition(size(sequences,1),'HoldOut',0.3);
% idx = cv.test;
TrainSplit = round(size(sequences,1)*.7)
TestSplit = round(size(sequences,1)*.3)

% Separate to training and test data
seqTrain = sequences(1:TrainSplit,:);
seqTest  = sequences(1:TestSplit,:);
tarTrain = target(1:TrainSplit,:);
tarTest  = target(1:TestSplit,:);

% Setup inputs into the network
InputSize = [size(sequences{1,1}),1];

numResponses = size(target{1},1);
clc
diary on
layers = layerGraph();

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



% Plot the network architecture
FigHandle_04 = figure('Position', [100, 150, 350, 290]);
plot(layers) 
saveas(FigHandle_04,ModelSaveFile,'jpeg')
clear FigHandle_04

%%Setup training options
miniBatchSize = 28;
numObservations = size(seqTrain,1);
numIterationsPerEpoch = floor(numObservations / miniBatchSize);

options = trainingOptions('adam', ...
    'MiniBatchSize',miniBatchSize, ...
    'InitialLearnRate',1e-4, ...
    'ExecutionEnvironment','auto', ...
    'MaxEpochs', 10000,...
    'GradientThreshold',1, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{seqTest,tarTest}, ...
    'ValidationFrequency',numIterationsPerEpoch, ...
    'LearnRateSchedule','none', ...
    'LearnRateDropFactor',0.2, ...
    'LearnRateDropPeriod',50, ...
    'Plots','training-progress', ...
    'Verbose',1);

% Train the RNN
close all
trainedNet = trainNetwork(seqTrain,tarTrain,layers,options);



save(ModelSaveFile)
diary(strcat(ModelSaveFile,'.txt'))
disp('NETWORK TRAINED')
toc