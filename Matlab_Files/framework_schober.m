%{
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
param.fs = 12000;       % Sampling frequency fft
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

% Geometry settingss
array_pos = [1 2.5 1.6];        % Center position of the microphone array (m)
room_geo = [5 5 3];             % Room dimensions (m) [5 3 3]
azimuth_ref = 35;
source_distance = [azimuth_ref,90,1];    % Distance of the source to the array center [azimuth, polarangle, radius]


param.plot = 0;                 % Display setup
% Missmatch settings
param.gain_std = 0;             % Standard Deviation of the microphone gains
param.d_std = 0;                % Standard Deviation of the microphone displacement [m]
param.tf_mems = [0,0,0,0,0,0,0,0,0,0];
array_geo = load(fullfile(pwd,'linear_2.dat'));
target_type = 'speaker';
%------------------------------------------------
% Load source signal
%------------------------------------------------
% [s, fs_sig] = audioread('s2_swwp2s_short.wav');   % s2_swwp2s.wav
% s = resample(s,param.fs, fs_sig);

% Read all .wav files in folder 'WAVs'
filenames =  dir('WAVs/*.wav');
s = zeros(0,0);
for ind = 1:length(filenames);
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
stdvalues.azimuth_ref = azimuth_ref;
std.values.reverberation_ref = reverberation_ref;
stdvalues.fs = param.fs;
stdvalues.c = param.c_speed;
stdvalues.windowlength = param.window_length;
stdvalues.nfft = param.nfft;
stdvalues.PhiAngles = 0:degree_steps*pi/180:pi-degree_steps*pi/180; %resolution of search based algorithms
stdvalues.ThetaAngles = linspace(0,0,1); %only azimuth estimation

% Calculate noisy output signals
%------------------------------------------------
param.signal_type = target_type;
X = [];
for ind_rev = 1:length(RevValues)
    param.T60 = RevValues(ind_rev);% Reverberation time
    for ind = 1:length(SnrValues);
        param.sensor_snr = SnrValues(ind);     % SNR of the microphones
        [X(:,:,ind,ind_rev), x_ref] = generateTestData(s, source_distance, array_geo, array_pos, room_geo ,param);
    end
end

% -------------------------------------
% VAD
% -------------------------------------

pp.pr = 0.998;
x_ref_dummy = x_ref;
x_ref_dummy(x_ref_dummy == 0) = [];
x_ref_vadsohn_values = logical(vadsohn(x_ref_dummy,stdvalues.fs,'a',pp));
x_ref_va = x_ref_dummy(x_ref_vadsohn_values);
X = X(x_ref_vadsohn_values,:,:,:); % short version

% save('fs_12k.mat')
%}

% 
load('fs_15k6_5degree_winkelabhängigkeit_extensive.mat');
stdvalues.windowlength = 256;
scale = 200;
X = X*scale;

% -------------------------------------
% Initializing
total_blocks = floor(length(X)/stdvalues.windowlength);
phi_CC_interpol = zeros(total_blocks,length(SnrValues),length(RevValues));

% -------------------------------------
optional.phat = false;
optional.interpolate = true;

for index_round = 1:1
for index_azimuth = 1:length(AzimuthValues)
    
for index_frame = 1:total_blocks
    X_Frame = squeeze(X(1+(index_frame-1)*stdvalues.windowlength:(index_frame)*stdvalues.windowlength,:,index_azimuth)); %1.
    %-------------------------------------
    % Cross Correlation with Interpolation
    %-------------------------------------
    [ phi_CC_interpol(index_frame,index_azimuth,index_round), ~ ] = GCC_PHAT( X_Frame, stdvalues, optional );

end
end
end

mean_phi_CC_interpol = real(mean(phi_CC_interpol));

figure();plot(mean_phi_CC_interpol);