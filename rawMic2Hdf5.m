%   ----------------------------------------------------------------------------
%   This script takes all of the raw microphone data from a folder, and
%   beamforms all microphones into one summed signal, where it is then
%   saved to an hdf5 file to be processed and fed to a neural network in
%   python
%   Cases:
%       - MWS202, clean airfoil 0 AOA, 10 - 40 in 10 m/s
%       - MWS203, turbulent airfoil 0 AOA, 10 - 40 in 10 m/s
%       - MWS102, clean cylinder 0 AOA, 10 - 40 in 10 m/s
%       - MWS103, turbulent cylinder 0 AOA, 10 - 40 in 10 m/s
%       - TSC101, Kevlar only (one screen)
%       - TSC202, Two Kevlar screens and speed varies from 0 - 40 in 5 m/s
%       - TSC203, Two Kevlar screens, wind tunnel off
%
%
% Input files:
%  1) Dir of raw data (from cases above) with varying speed
%
%
%
%
% %   Written By: Salvador Mayoral
% %   Written on: ???
% %    Modified By: Andrew Bartels
% %    Modified on: 4-2-2019
%
%   Run time of this script is roughly 41 seconds
%   ----------------------------------------------------------------------------
%
%   INITIALIZE CODE
%   ----------------------------------------------------------------------------
clc
clear all
close all


tic
%   Specify constants
SampleRate  = 48000;                        %  DAQ Sample rate (S/sec)
NS          = 48000;                        %  Number of samples
fn          = SampleRate/2;                         %  maximum resolvoble frequency
NFFT        = 2^12;                         %  4096 point FFT
NF          = NFFT/2;                       %  No. point for powerspecturm
NR          = 4;                            %  No. of runs
NP          = 4;                            %  No. of positions
NM          = 7;                            %  No. of microphones
Pref        = 20e-6;                        %  Reference pressure
sensitivity = [50.2 49.3 53.5 52.9 52.8 48.7 47.0]; % microphone sensitivity
c = 343;                                    % Speed of sound (m/s)
d  = 0.00858;                               % spacing of microphone array (cm)
%
%   ----------------------------------------------------------------------------
%
%   IMPORT AND SORT RAWDATA
%   ---------------------------------------------------------------------------
%   Define microphone frequency responces files
micResponse7  = {'50' '56' '51' '54' '55' '49' '79'};    %  microphone freq response
micResponse8  = {'50' '56' '51' '54' '55' '49' '79' '56'};    %  microphone freq response
calDirectory = '../CALIBRATION/sn442';                  %  calibration directory

R = zeros(NF,NM);                           %  allocate mic response
f = SampleRate*(0:NF-1)/NFFT;                       %  define frequency range

%   Import microphone frequency responces
for i = 1:NM
    filename = strcat(calDirectory,micResponse7{i},'.txt');
    A = dlmread(filename);   %  import presssure rawdata
    fm = A(:,1);
    Sm = A(:,2);
    clear A
    R(:,i) = interp1(fm,Sm,f','pchip');
end
%
%   ---------------------------------------------------------------------------
%%   Define data files
casename     = 'ALL';                    %  casename for data-processing
datDirectory = 'C:\Users\andrewbartels1\Dropbox (CSU Fullerton)\EGME597_AB\ML_DATA\RAWDATA\';
fileList = glob(strcat(datDirectory,'*.dat'));
%    This is all the things to parse the file list into
splices = {'V','P','R','.dat'};

fileList1 = fileList;

% Pull out the data from the glob filelist (rev 2 way)
for i = 1:size(fileList1,1)
    filename_test = fileList1{i};
    A2 = dlmread(filename_test);            %  import mic rawdata
    data_sizing(:,i) = size(A2);
    temp_string = char(fileList1{i});
    field{i} = temp_string(76:end-4);
    data = 1000*A2./(sensitivity*Pref); %  convert voltage to pressure and normalize by Pref
    V1{1, i} = struct(field{i},data);
end
disp('All data read in.')
NumMics = (size(data,2));
%
%   ----------------------------------------------------------------------------
%
%   Delay-and-Sum Beamforming
%   ---------------------------------------------------------------------------
% Setup the beamforming arrays and sizing
array = phased.ULA('NumElements',NumMics,'ElementSpacing',d); % define arrat
microphone = phased.OmnidirectionalMicrophoneElement(...
    'FrequencyRange',[20 20e3]); % frequency response
collector = phased.WidebandCollector('Sensor',array,'SampleRate',24e3,...
    'PropagationSpeed',c,'ModulatedInput',false);
sigang = zeros(2, NumMics);
beamformer = phased.TimeDelayBeamformer('SensorArray',array,...
    'SampleRate',SampleRate,'PropagationSpeed',c,'Direction',[0; 0]);

%%  Beamform each series of arrays of microphones %%
% initialize OUTSIDE THE ARRAY. (DON"T BE DUMB LIKE ME)
j=1;                                                        %counter for subracting pressures
    k0 = 2;
    k1 = 3;
    k2 = 4;
    k3 = 5;
    A = 1:5:size(V1,2);
    
for i = 1:size(V1,2)
  
    rsig = collector(V1{i}.(field{i}),sigang);
    amplituedBeamformedSig = beamformer(rsig)/NumMics;
    FreqBeamSig = fft(amplituedBeamformedSig,NFFT)/NFFT;
    acousticPressure(:,i) = abs(FreqBeamSig).^2;
    
    % taking out the microphone sensitivity profile
    R1 = interp1(R,1:size(acousticPressure),'pchip');
    SPL(:,i) = 10*log10(acousticPressure(:,i))-mean(R1,2);
    CAcousticPressure(:,i) = 10.^(SPL(:,i)/10);
end
for j = 1:5:size(V1,2)
    CleanedAcousticPressure(:,k0) = abs(CAcousticPressure(:,k0) - CAcousticPressure(:,j));
    CleanedAcousticPressure(:,k1) = abs(CAcousticPressure(:,k1) - CAcousticPressure(:,j));
    CleanedAcousticPressure(:,k2) = abs(CAcousticPressure(:,k2) - CAcousticPressure(:,j));
    CleanedAcousticPressure(:,k3) = abs(CAcousticPressure(:,k3) - CAcousticPressure(:,j));
    k0 = k0 + 5;
    k1 = k1 + 5;
    k2 = k2 + 5;
    k3 = k3 + 5;
end
f1 = interp1(f,1:size(R1,1),'pchip');
for j = 1:size(SPL,2)
    SPL(:,j) = 10*log10(CleanedAcousticPressure(:,j));
    SPL_smooth(:,j) = smooth(f1,SPL(:,j),0.005,'rloess' );
end

% original smoothing function
% B      = 10*log10(Pressure1(i,:));
% b(:,i) = smooth(f,B,0.005,'rloess' );

disp('numbers crunched')

%%   ----------------------------------------------------------------------------
% Write each structure to an hdf5 file %%
% where each field is a file with a beamformed pressure, and smoothed
% Items to export:
%  field (name of each file)
%  amplitudeBeamformedSig (raw beamformed pressure)
%  raw_spectrum (fft(abs()^2) of the raw pressure
%  s (smoothed
%   ----------------------------------------------------------------------------
hdfFilename = strcat(casename,'.hdf5')

hdf5write(hdfFilename, '/field', field);
hdf5write(hdfFilename, '/RawPressure', CleanedAcousticPressure,'writemode','append');
hdf5write(hdfFilename, '/RawSpectrum', SPL,'writemode','append');
hdf5write(hdfFilename, '/SmSpecturm', SPL_smooth,'writemode','append');

% h5create(hdfFilename,'/',[10 20])
% mydata = rand(10,20);
% h5write('myfile.h5', '/DS1', mydata)

disp('beamforming and file writing complete')
%     This ensures that all the data arrays are the same size and will
%     raise an error if they aren't


% Code from: An Introduction to Field Analysis Techniques: The Power
% Spectrum and Coherence (2013) and this essentially sees at what
% frequencies the signals are coherent. (important assumption to prove)
% K = size(x,1); %Define the number of trials.
% N = size(x,2); %Define the number of indices
% per trial.
% dt = t(2)-t(1); %Define the sampling interval.
% T = t(end); %Define the duration of data.
% Sxx = zeros(K,N); %Create variables to save the
% spectra.
% Syy = zeros(K,N);
% Sxy = zeros(K,N);
% for k=1:K %Compute the spectra for each
% trial.
%  Sxx(k,:) = 2*dt^2/T * fft(x(k,:)) .* conj(fft(x(k,:)));
%  Syy(k,:) = 2*dt^2/T * fft(y(k,:)) .* conj(fft(y(k,:)));
%  Sxy(k,:) = 2*dt^2/T * fft(x(k,:)) .* conj(fft(y(k,:)));
% end
% Sxx = Sxx(:,1:N/2+1); %Ignore negative frequencies.
% Syy = Syy(:,1:N/2+1);
% Sxy = Sxy(:,1:N/2+1);
% Sxx = mean(Sxx,1); %Average the spectra across
% trials.
% Syy = mean(Syy,1);
% Sxy = mean(Sxy,1);
% cohr = abs(Sxy) ./ (sqrt(Sxx) .* sqrt(Syy));
%  %Compute the coherence.
% df = 1/max(T); %Determine the frequency
% resolution.
% fNQ = 1/ dt / 2; %Determine the Nyquist
% frequency.
% faxis = (0:df:fNQ); %Construct frequency axis.
% plot(faxis, real(cohr)); %Plot the results
% xlim([0 50]); ylim([0 1]) %Set the axes limits
% xlabel('Frequency [Hz]') %Label axes.
% ylabel('Coherence [ ]')






toc