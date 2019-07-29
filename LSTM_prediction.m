%   ----------------------------------------------------------------------------
%   This script takes all of the different models trained and imports all
%   their data and plots them. Lastly it makes a full signal prediction
%   from 2 trained microphones.
%   Cases:
%       - MWS202, clean airfoil 0 AOA, 10 - 40 in 10 m/s
%       - MWS203, turbulent airfoil 0 AOA, 10 - 40 in 10 m/s
%       - MWS102, clean cylinder 0 AOA, 10 - 40 in 10 m/s
%       - MWS103, turbulent cylinder 0 AOA, 10 - 40 in 10 m/s
%       - TSC101,  Kevlar only (one screen)
%       - TSC202, Two Kevlar screens and speed varies from 0 - 40 in 5 m/s
%       - TSC203, Two Kevlar screens, wind tunnel off
%       - MLD,    Two Kevlar screens anechoic chambers all sides 5-25 m/s
%
%   All the inputs are managed by 'LSTM_harness.m'
%   Required inputs:
%   
%    1) micResponse ('MWS' or 'ACC')
%    
%
%
% %    Written By: Andrew Bartels
% %    Written on: 1-1-2019
% %    Modified By: Andrew Bartels
% %    Modified on: 7-22-2019
%     Original base code written by: Salvador Mayoral
%   ----------------------------------------------------------------------------
%
%   INITIALIZE CODE
%   ----------------------------------------------------------------------------
clear all 
close all
clc

tic
matDirectory = 'C:\Users\andrewbartels1\Dropbox (CSU Fullerton)\EGME597_AB\MATLAB\'
ModelsFileList = glob(strcat(matDirectory,'*.mat'));

%Import all the data here

% Run just this section and select the correct model to load. Double-click
% and open the ModelsFileList variable.
FileNumber = 10;

disp('Model being loaded:')
ModelsFileList{FileNumber}
%%
load(ModelsFileList{FileNumber},'trainedNet');
load('CleanedDataForNN_notsmoothed__ACC.mat');

dirLength = numel(ModelsFileList{FileNumber});
% Make a dir to save all the prediction images in
DirName = ModelsFileList{FileNumber}(1:dirLength-7)
mkdir(DirName);
% take the real out of smooth so MFCCs can deal with it.
SPL_smooth1 = real(SPL_smooth);

% Plot some mel Specs for pretty pictures...
randNum = randi(size(SPL_smooth1,2),1);

for i=1:size(Data,3)    
[f0(:,i), idx(:,i)] = pitch(Data(:,1,10), SampleRate);
end


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
                                 
FigHandle_01 = figure('Position', get(0, 'Screensize'),'Visible','off');
SpecCellArray{1,i} = real(spectrogram(SPL1(:,i))); %images, just in arrays
end

% Get all of the extracted features and put it in an NxM matrix inside a
% cell. Such that N are columns of features and M are the total number of
% training files.
for i=1:size(Data,3)
% sequences{i,1}(:,1) = MelSpecs_2(:,i);
% sequences = SpecCellArray';
sequences{i,1}(:,1) = coeffs_2(:,i);
% sequences{i,1}(:,3) = f0(1:size(MelSpecs_2,1),i);
% sequences{i,1}(:,4) = delta2(:,i);
% sequences{i,1}(:,5) = deltaDelta2(:,i);
% sequences{i,1}(:,6) = centFreq2(:,i);
target{i,1} = SPL_smooth1(:,i);
end
tic
for i=1:size(sequences,1)
YPred(:,i) = predict(trainedNet,sequences{i},'MiniBatchSize',1);
% [trainedNetUpdated, YPred(:,i)] = predictAndUpdateState(trainedNet,sequences{i});

% The model obviously understands the trend of things, but doesn't
% understand magnitudes are super important, so the prediction is corrected
% to move it up or down to start at the same point at the input sequence.

% calculates the distance between the first prediction point and the mean
% of the first 3 points of the target signal.
% Distance(i) = YPred(1,i) - mean(SPL1(1:3,i));
% 
% % and then moves the line up or down by Distance
% if Distance(i) >= 0 
%     YPred(:,i) = YPred(:,i) - abs(Distance(i));
% else
%     YPred(:,i) = YPred(:,i) + abs(Distance(i));
% end
% 
% test(i) = YPred(1,i) - mean(SPL1(1:3,i));

end

disp('the prediction time is:')
toc



disp('/n /n Prediction made, now plotting...')
% % Save in a loop to look through later in a .jpg
for i=1:size(sequences,1)
i
FigHandle_01 = figure('Position', [100, 150, 350, 290],'Visible','off');
semilogx(f,[YPred(:,i), target{i}, SPL1(:,i)]);
legend('Prediction','Target','Input')
title(sprintf('prediction %d',i))
xlabel('frequency (Hz)')
ylabel('Sound Pressure Level (SPL)')

axis([0 4000 0 130])
saveas(FigHandle_01,fullfile(DirName,sprintf('prediction %d',i)),'jpeg')
clear FigHandle_04
end


% Get all the text files and see which model is the best so far.

% then predict with that model.

toc









