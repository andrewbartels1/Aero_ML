%   ----------------------------------------------------------------------------
%   This script takes all of the raw microphone data from a folder, and
%   beamforms all microphones into one beamformed, smoothed signal of all 7 microphones,
%   and one sig
%   Cases:
%       - MWS202, clean airfoil 0 AOA, 10 - 40 in 10 m/s
%       - MWS203, turbulent airfoil 0 AOA, 10 - 40 in 10 m/s
%       - MWS102, clean cylinder 0 AOA, 10 - 40 in 10 m/s
%       - MWS103, turbulent cylinder 0 AOA, 10 - 40 in 10 m/s
%       - TSC101, Kevlar only (one screen)
%       - TSC202, Two Kevlar screens and speed varies from 0 - 40 in 5 m/s
%       - TSC203, Two Kevlar screens, wind tunnel off 
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
% clc
close all

tic
%   Specify constants
SampleRate  = 48000;                        %  DAQ Sample rate (S/sec)
NS          = 48000;                        %  Number of samples
fn          = SampleRate/2;                 %  maximum resolvoble frequency
NFFT        = 2^12;                         %  4096 point FFT
NF          = NFFT/2;                       %  No. point for powerspecturm
NR          = 4;                            %  No. of runs
NP          = 4;                            %  No. of positions
NM          = 7;                            %  No. of microphones
Pref        = 20e-6;                        %  Reference pressure
sensitivity = [50.2 49.3 53.5 52.9 52.8 48.7 47.0]; % microphone sensitivity
c = 343;                                    % Speed of sound (m/s)
d  = 0.00858;                               % spacing of microphone array (cm)


%   ----------------------------------------------------------------------------
%
%   IMPORT AND SORT RAWDATA
%   ---------------------------------------------------------------------------
%   Define microphone frequency responces files
micResponse_MWS  = {'50' '56' '51' '54' '55' '49' '79'};    %  microphone freq response (in the proper order written)
micResponse_ACC  = {'50' '51' '54' '49' '79' '55' '56'};    %  microphone freq response (in the proper order written)

% this just makes using the LSTM_harness.m easier to use.
switch(Case)
    case 'MWS'
        micResponse = micResponse_MWS;
        disp('MWS being used')
    case 'ACC'   
        micResponse = micResponse_ACC;
        disp('ACC being used')
    case 'BOTH'
        micResponse = micResponse_ACC;
        disp('Both being used')
    otherwise
        error('MWS or ACC are the only cases available')
end
        

R = zeros(NF,NM);                           %  allocate mic response
f = SampleRate*(0:NF-1)/NFFT;                       %  define frequency range

%   Import microphone frequency responces
for i = 1:NM
    filename = strcat(Calibration_dir,micResponse{i},'.txt');
    A = dlmread(filename);   %  import presssure rawdata
    fm = A(:,1);
    Sm = A(:,2);
    clear A
    R1(:,i) = interp1(fm,Sm,f','pchip'); %ends up being 2048x7
    R = mean(R1,2);
end

%%   Define data files
switch(Case)
    case 'MWS'
        datDirectory = MWS_dir;
    case 'ACC'
        datDirectory = ACC_dir;
    case 'BOTH'
        datDirectory = BOTH_dir;
    otherwise
        error('MWS or ACC are the only cases available')
end

% Get all the files and make a file list. 
datDirectory
fileList = glob(strcat(datDirectory,filesep,'*.dat'))

%   ----------------------------------------------------------------------------
%
%   PREALLOCATE MEMORY
%   ---------------------------------------------------------------------------
% Elapsed time is 286.441785 seconds. (without preallocation)
% Preallocate a bunch of arrays for faster performance
% this is done now since fileList gives many of the sizing.

% Data = zeros(SampleRate, 6, size(fileList,1));
acousticPressureTwoMics = zeros(NF,size(fileList,1));
SPL1 = zeros(NF,size(fileList,1));


% Pull out the data from the glob filelist (rev 2 way)
for i = 1:size(fileList,1)
    filename_test = fileList{i};
    A2 = dlmread(filename_test);            %  import mic rawdata
    
    %  convert voltage to pressure and normalize by Pref and 
    Data(:,:,i) = 1000*A2(:,2:7)./(sensitivity(size(A2(:,2:7),2))*Pref);  %slam all the data into a 3-D array (a little easier to work with) 
    
    %%%Depreciated version of storing data%%%
    % temp_string = char(fileList{i});
    % field{i} = temp_string(76:end-4);
    % V1{1, i} = struct(field{i},data); %feed into struct for processing pipeline (a little outdated and painful)
    % If structs are used later this is the best way to index throught them
    % V1{i}.(field{i}) with i being the final length
end

disp('All data read in.')

NumMics = (size(Data,2));

% This raises an error if 7 microphones are fed in since this will cause a
% lot of issues
if NumMics ~= 6
    message = 'Total microphones read in is not 7!'
    error(message)
end


%   ----------------------------------------------------------------------------
%
%   Delay-and-Sum Beamforming
%   ---------------------------------------------------------------------------
% Setup the beamforming arrays and sizing

%% Setup Beamforming for just 2 mics
NumMics1 = 2;
array1 = phased.ULA('NumElements',NumMics1,'ElementSpacing',d); % define arrat
microphone1 = phased.OmnidirectionalMicrophoneElement(...
    'FrequencyRange',[20 20e3]); % frequency response
collector1 = phased.WidebandCollector('Sensor',array1,'SampleRate',24e3,...
    'PropagationSpeed',c,'ModulatedInput',false);
sigang1 = zeros(2, NumMics1);
beamformer1 = phased.TimeDelayBeamformer('SensorArray',array1,...
    'SampleRate',SampleRate,'PropagationSpeed',c,'Direction',[0; 0]);

% acousticPressureTwoMics beamforming
for i = 1:size(Data,3)
%     REMINDER: Data is a 3D vector of NORMALIZED PRESSURES!
    rsig1 = collector1(Data(:,1:2,i),sigang1);
    amplituedBeamformedSig1 = beamformer1(rsig1); % 4096 point time domain signal

    % Converting to the frequency domain for some subtractions
    FreqBeamSig1 = fft(amplituedBeamformedSig1/NFFT,NFFT); 
    acousticPressureTwoMics(:,i) = 20*log10(abs(FreqBeamSig1(1:NF))/NumMics1); % SPL in dB form 2048 spectral conversion
    
    % taking out the microphone sensitivity profile corrected for freq response  
    SPL1(:,i) = acousticPressureTwoMics(:,i)-R;
    
    % Put it back to a normalized pressure    
%     CleanedAcousticPressureTwoMics(:,i) = real(10.^(SPL1(2:end,i)/10));
end

%% Setup Beamforming for all 7 mics

array = phased.ULA('NumElements',NumMics,'ElementSpacing',d); % define arrat
microphone = phased.OmnidirectionalMicrophoneElement(...
    'FrequencyRange',[20 20e3]); % frequency response
collector = phased.WidebandCollector('Sensor',array,'SampleRate',24e3,...
    'PropagationSpeed',c,'ModulatedInput',false);
sigang = zeros(2, NumMics);
beamformer = phased.TimeDelayBeamformer('SensorArray',array,...
    'SampleRate',SampleRate,'PropagationSpeed',c,'Direction',[0; 0]);

for i = 1:size(Data,3)
    rsig = collector(Data(:,:,i),sigang);
    amplituedBeamformedSig = beamformer(rsig);
%     Converting to the frequency domain for some subtractions
    FreqBeamSig = fft(amplituedBeamformedSig/NFFT,NFFT);
    acousticPressure(:,i) = 20*log10(abs(FreqBeamSig(1:NF))/NumMics);
    
    % taking out the microphone sensitivity profile
    SPL_uncleaned(:,i) = acousticPressure(:,i)-R;
end

clear SPL
% Subtract out noise from wind tunnel (no point source in all R0 files)
% Dr. Fowler said to subtract the mean of the noise signal from the
% beamformed signal (check with Mayoral.
% Setup counters for subtracting the R0 (noise) from the other signals
% (OUTSIDE THE LOOP)
    j  = 1;  
    k0 = 2;
    k1 = 3;
    k2 = 4;
    k3 = 5;
SPL_cleaned = zeros(size(SPL_uncleaned,1),size(SPL_uncleaned,2));

for j = 1:5:size(Data,3)
     
    SPL_cleaned(:,k0) = 10*log10(10.^(SPL_uncleaned(:,k0)/10) - 10.^(SPL_uncleaned(:,j)/10));
    SPL_cleaned(:,k1) = 10*log10(10.^(SPL_uncleaned(:,k1)/10) - 10.^(SPL_uncleaned(:,j)/10));
    SPL_cleaned(:,k2) = 10*log10(10.^(SPL_uncleaned(:,k2)/10) - 10.^(SPL_uncleaned(:,j)/10));
    SPL_cleaned(:,k3) = 10*log10(10.^(SPL_uncleaned(:,k3)/10) - 10.^(SPL_uncleaned(:,j)/10));
    
%     These are just counters that help do the above signal - noise maths
%     (don't worry unless there are a different number  of runs (aka max of
%     .dat files is not in sets of R0-R5)
    k0 = k0 + 5;
    k1 = k1 + 5;
    k2 = k2 + 5;
    k3 = k3 + 5;
end

% Take out all the NaNs (result of R0 being subracted from itself)
SPL_cleaned(:,1:5:end) = [];
SPL1(:,1:5:end) = [];
Data(:,:,1:5:end) = [];

% smoothing is apparently the slowest part (takes upwards of 68/315 seconds
% of script runtime.
for j = 1:size(SPL_cleaned,2)
    SPL_smooth(:,j) = smooth(f,SPL_cleaned(:,j),0.0055,'rloess' );
    SPL1(:,j) = smooth(f,SPL1(:,j),0.005,'rloess' );
end

% Plot example to ensure things make sense
randNum = randi(size(SPL1,2),1)

FigHandle_01 = figure('Position', [100, 150, 350, 290]);
semilogx(f,[SPL_cleaned(:,randNum),SPL1(:,randNum), SPL_smooth(:,randNum)],'LineWidth', 1)
legend('SPL_{denoised}','2 mic SPL_{smoothed}','7 mic SPL_{smoothed}', 'Location', 'Southwest')
title('2 mics Beamformed vs 7 mics smoothed & Beamformed', 'FontSize', 8.5)
xlabel('frequency (Hz)')
ylabel('Sound Pressure Level (dB)')
saveas(FigHandle_01,'example_signal_output','jpeg')

% Save a bunch of stuff as an intermediary file... It's probably a fairly
% large file size, so it's best to delete after reading in but won't since
% Dropbox space isn't an issue.
save(PreprocessingSaveFile, 'Data','SPL1','SPL_smooth','SPL_cleaned','SampleRate','NS','fn','NFFT',...
                'NF','NR','NP','NM','Pref','sensitivity','c','d', 'f', 'fileList')                            
disp('Data ready for LSTM network')


toc