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
    datDirectory = 'C:\Users\andre\Dropbox (CSU Fullerton)\EGME597_AB\ML_DATA\RAWDATA\';              
    fileList = glob(strcat(datDirectory,'*.dat'));
%    This is all the things to parse the file list into
    splices = {'V','P','R','.dat'};
    
%     Initialize Arrays for speed
%     filenames_TF = contains(fileList,casename);
%     fileList1 = fileList(filenames_TF);
    fileList1 = fileList; 
%     parsed_file_list = cell(size(fileList1));
%     Velocity_list = cell(size(fileList1));
%     TSC_list = cell(size(fileList1));
%     Position_list = cell(size(fileList1));
%     Run_list = cell(size(fileList1));
%     
%     %     Pull every .dat file from fileList
%     for i=1:size(fileList1)
%         parsed_file_list{i} = strsplit(fileList1{i},splices);
%         TSC_list(i) = parsed_file_list{i}(1,4);
%         Velocity_list(i) = parsed_file_list{i}(1,5);
%         Position_list(i) = parsed_file_list{i}(1,6);
%         Run_list(i) = parsed_file_list{i}(1,7);       
%     end
%     
% display('The Velocities are:')
% Velocities = unique(Velocity_list)
% TSCs = unique(TSC_list);
% Max_runs = unique(Run_list);
% Number_Positions = unique(Position_list);

% Pull out the data from the glob filelist (rev 2 way)
for i = 1:size(fileList1,1)
    filename_test = fileList1{i};
    A2 = dlmread(filename_test);            %  import mic rawdata
    data_sizing(:,i) = size(A2);
    temp_string = char(fileList1{i});
    field{i} = temp_string(67:end-4);
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
for i = 1:size(V1,2)
rsig = collector(V1{i}.(field{i}),sigang);
amplituedBeamformedSig(:,i) = beamformer(rsig)/NumMics;
FreqBeamSig = fft(amplituedBeamformedSig,NFFT)/NFFT;
acousticPressure(:,i) = abs(FreqBeamSig(1:NF)).^2;
SPL(:,i) = 10*log10(acousticPressure(:,i))-mean(R,2);
CleanedAcousticPressure(:,i) = 10.^(SPL(:,i)/10);
% This is subtracting out the background noise. Trying to have the machine
% learning do this part. Will do this to label nonnoise features of the
% dataset
% Pressure(:,i) = abs(Pressure(:,i) - Pressure(:,1,i));

SPL(:,i) = 10*log10(CleanedAcousticPressure(:,i));
SPL_smooth(:,i) = smooth(f,SPL(:,i),0.005,'rloess' );

% Not sure what this is.
% B      = 10*log10(Pressure1(i,:));
% b(:,i) = smooth(f,B,0.005,'rloess' );
end
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

hdf5write(hdfFilename, '/field', field)
hdf5write(hdfFilename, '/RawPressure', SPL,'writemode','append')
hdf5write(hdfFilename, '/RawSpectrum', CleanedAcousticPressure,'writemode','append')
hdf5write(hdfFilename, '/SmSpecturm', SPL_smooth,'writemode','append')

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