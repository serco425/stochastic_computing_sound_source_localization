close all

addpath(fullfile(pwd,'RIR Erlangen'));
addpath(fullfile(pwd,'misc'));
addpath(fullfile(pwd,'voicebox'));
addpath(fullfile(pwd,'stochastic_computing'))

% Init
%------------------------------------------------
% General settings
interference = 0; % Interfering speaker/noise
diffuse = 0;
param.fs = 15600;%12000;       % Sampling frequency fft
param.c_speed = 345.8;  % Sound velocity (24°C = 345.8m/s)
reverberation_ref = 0.0;
RevValues = reverberation_ref;%[0,0.3,0.6];     % Reverberation (param.T60 = RevValues(ind_rev))
                    
% Frame settings                -> 0.032 s window length ->
param.window_length = 2.^nextpow2(0.016*param.fs);  % Window length in s
param.frameshift = param.window_length/2;
param.nfft = 4*param.window_length;
param.window_type = 'hann';

% Acoustic settings
param.SPL_int = 50;         % Sound pressure level of the interference in 1m distance
param.SPL = 60;             % Sound pressure level of the speaker in 1m distance
param.sensitivity = -34;    % default -34 Reference sensitivity of the microphones at SPL = 94dB
SnrValues = 75;%linspace(60,90,6); % SNR of the microphones (param.sensor_snr = SnrValues(i))

% Geometry settings
array_pos = [1 2.5 1.6];        % Center position of the microphone array (m)
room_geo = [5 5 3];             % Room dimensions (m) [5 3 3]
AzimuthValues = 0:5:180;
source_distance = [0,90,1];    % Distance of the source to the array center [azimuth, polarangle, radius]
 

param.plot = 0;                 % Display setup
% Missmatch settings
param.gain_std = 0;             % Standard Deviation of the microphone gains
param.d_std = 0;                % Standard Deviation of the microphone displacement [m]
param.tf_mems = [0,0,0,0,0,0,0,0,0,0];
array_geo = load(fullfile(pwd,'linear_2.dat'));
target_type = 'speaker';
% ------------------------------------------------
% Load source signal
% ------------------------------------------------
%[s, fs_sig] = audioread('swab7a.wav');   % s2_swwp2s.wav
%s = resample(s,param.fs, fs_sig);

 % Read all .wav files in folder 'WAVs'
 filenames =  dir('WAVs/*.wav');
 s = zeros(0,0);
 for ind = 1:length(filenames)
     [signal, fs_sig] = audioread( ['WAVs/' filenames(ind).name] );
     signal = resample(signal,param.fs, fs_sig);
     s = cat(1,s,signal);
 end

len = length(s);
param.len = len;
s = s(1:len);
n = (0:len-1)';
t = n/param.fs;

%-------------------------------------
% generate standard settings for all algorithms
%-------------------------------------

degree_steps = 5; %5 degree seach steps for azimuth search based algos
stdvalues = generate_std(array_geo);
stdvalues.azimuth_ref = AzimuthValues;
std.values.reverberation_ref = reverberation_ref;
stdvalues.fs = param.fs;
stdvalues.c = param.c_speed;
stdvalues.windowlength = param.window_length;
stdvalues.nfft = param.nfft;
stdvalues.PhiAngles = 0:degree_steps*pi/180:pi-degree_steps*pi/180; %resolution of search based algorithms
stdvalues.ThetaAngles = linspace(0,0,1); %only azimuth estimation

% % Calculate noisy output signals
%------------------------------------------------
param.signal_type = target_type;

X = [];
for ind_azimuth = 1:length(AzimuthValues)
        source_distance(1) = AzimuthValues(ind_azimuth);
        [X(:,:,ind_azimuth), x_ref] = generateTestData(s, source_distance, array_geo, array_pos, room_geo ,param);
end

% VAD
pp.pr = 0.999;
x_ref_dummy = x_ref;
x_ref_dummy(x_ref_dummy == 0) = [];
x_ref_vadsohn_values = logical(vadsohn(x_ref_dummy,stdvalues.fs,'a',pp));
x_ref_va = x_ref_dummy(x_ref_vadsohn_values);
X = X(x_ref_vadsohn_values,:,:,:); % short version

save('fs_15k6_5degree_winkelabhängigkeit_extensive.mat');

