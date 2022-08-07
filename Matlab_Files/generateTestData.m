function [X, x_ref, h] = generateTestData(s, source_distance ,array_geo, array_pos, room_geo ,param)
%generateTestData Multichannel delay generator
%  X = generateTestData(s,source_dir,array_geo,param) produces delayed and
%  distorted copies of the source signal s.
%
%  Output:  X ........ matrix containing the delayed and distorted channel
%                      signals. Channels are organized in columns.
%           x_ref .... reference signal at the array origin
%           h ....... matrix containing the impuls responses of all channels
%           
%  Input:   s ................. source signal (mono)
%           source_distance ... vector containing the direction and distance 
%                               of the source [azimuth, elevation, radius]
%                               relative to the microphone array (°,°,m)
%           array_geo ......... matrix containing the array geometry. Each
%                               microphone position is described by the row
%                               vector [x, y, z, az, el] (m,m,m,°,°)
%           array_pos ......... origin of the microphone array
%           room_geo .......... room dimensions [x,y,z] (m)
%           param ............. struct containing following parameters:
%            .fs - sampling frequency (Hz)
%            .c_speed - sound velocity (m/s)        {Default: 345.8}
%            .T60 - reverberation time T60 (s)      {Default: 0}
%            .plot - 3D room setup plot             {Default: 0}
%            .SPL - sound pressure level (dB) of the source at 1m distance
%                                                   {Default: 60}
%            .sensitivity - Reference sensitivity of the microphones at SPL = 94dB
%                                                   {Default: -34}
%            .sensor_snr - SNR of the microphones   {Default: 63}

%% Init

s = s(:);                               % Guarantee column vector

if ~isfield(param,'c_speed')
    param.c_speed = 345.8;
end
if ~isfield(param,'T60')
    param.T60 = 0;
end
if ~isfield(param,'plot')
    param.plot = 0;                     % Do not plot setup
end

if ~isfield(param,'SPL')
    param.SPL = 60;
end

if ~isfield(param,'SPL')
    param.sensitivity = -34;
end

if ~isfield(param,'SPL')
    param.sensor_snr = 63;
end

phi = source_distance(1) * (pi/180);    % Source directions in radians
theta = source_distance(2) * (pi/180);
radius = source_distance(3);

fs = param.fs;
c = param.c_speed;                      % Sound velocity (m/s)
beta = param.T60;                       % Reverberation time (s)

M = size(array_geo,1);                  % Number of microphones

%% Calculate output signals

rp = array_geo(:,1:3) + repmat(array_pos,M,1);                      % Receiver positions [x y z] (m)
sp = array_pos + ...                                                % Source position [x y z] (m)
    radius*[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
L = room_geo;                                                       % Room dimensions [x y z] (m)


if ~isfield(param,'diffuse')

    if beta == 0    % T60 = 0;
        ir_length = 128;                                                % Length of the impulse responses for non-reverberant environment
        h = rir_generator(c, fs, rp, sp, L, beta, ir_length, [], [], [], [], 0)';
    else
        h = rir_generator(c, fs, rp, sp, L, beta, [], [], [], [], [], 1)';
    end
    
    s = s(1:param.len);
    X = fftfilt(h,s,length(h));                                         % Fast convolution

    % Virtual noisy (reverberant) reference signal at the array origin

else
    X = gen_diffuse(s,M,param);
    X = X(1:param.len,:);
end

X = addMismatch(X, param);

% Virtual clean reference signal at the array origin
h_ref = rir_generator(c, fs, array_pos, sp, L, 0, 128, [], [], [], [], 0)';
x_ref = fftfilt(h_ref,s,length(h_ref));

%% Normalization

SPL_ref = 94;       % Reference sound pressure level at 1kHz (1m distance) in dBSPL
f_ref = 1000;       % Reference sinus at 1kHz and 94 dBSPL
a_ref = 10.^(param.sensitivity/20);
sinus_ref = a_ref*sin(2*pi*f_ref/fs*(0:size(X,1)-1)');

SPL_mic = param.SPL + 20*log10(1/radius); % SPL at the reference microphone
SPL_diff = SPL_mic-SPL_ref;
g_SPL_diff = 10.^(SPL_diff/20);

if strcmp(param.signal_type,'noise')
    pp.pr = 0;
else
    pp.pr = 0.999;
end

x_ref_dummy = x_ref;
x_ref_dummy(x_ref_dummy == 0) = [];
x_ref_va = x_ref_dummy(logical(vadsohn(x_ref_dummy,fs,'a',pp)));
% x_ref_va = x_ref(logical(vadsohn(x_ref(10001:end),fs,'a',pp)));
g_SPL94 = rms(sinus_ref)/rms(x_ref_va);
if isfield(param,'diffuse')
    g_SPL94 = rms(sinus_ref)/rms(X(:,1));
end

X = g_SPL_diff*g_SPL94*X;
x_ref = g_SPL_diff*g_SPL94*x_ref;

% Add sensor noise
if isfield(param,'sensor_snr') && ~isfield(param,'interference')
    N = randn(size(X));             % Sensor noise
    NA = filterA(N(:,1), fs, 0);    % A-weighted sensor noise
%     g_SNR = 10.^(param.sensor_snr/20)*sum(abs(NA.^2))/sum(abs(sinus_ref.^2));
    se = norm(sinus_ref)^2;
    nsc = se/(10^(param.sensor_snr/10));
    ne = norm(N(:,1))^2;
    N = sqrt(nsc/ne)*N;
%     g_SNR = 10.^(param.sensor_snr/10)*norm(N(:,1)).^2/norm(sinus_ref).^2;
%     N = N/sqrt(g_SNR);
%     NA = NA/sqrt(g_SNR);
    
    X = X + N;
end 

%% 3D plot of the above setup:

if param.plot
    figure
    plot3(sp(1),sp(2),sp(3),'ro-','markersize',4); hold on;                             % Source position
    plot3(rp(:,1),rp(:,2),rp(:,3),'ko','markersize',4,'markerfacecolor',ones(1,3)*.6);  % Receiver positions
    axis equal; axis([0 room_geo(1) 0 room_geo(2) 0 room_geo(3)]);                      % Room dimensions
    box on; xlabel('x-axis (m)'); ylabel('y-axis (m)'); zlabel('z-axis (m)');
end


function X_mismatch = addMismatch(X,param)

param.M = size(X,2);

fs = param.fs;
wl = param.window_length;
frameshift = param.frameshift;
nfft = param.nfft;

k = [0:nfft/2 -nfft/2+1:-1]';
omega = 2*pi*fs/nfft*k;
param.f = omega/(2*pi);

[Xf, ~, param] = create_input_buffer(X,param);
param = calcMismatch(param.gain_std, param.d_std, param.M, param);
H = (1 + param.gain_dev).*exp(1i*omega*param.place_dev(1,:)/param.c_speed);

for channel = 1:param.M
    mems_tf(:,channel) = mems_tf_var(param.tf_mems(channel),param.f);
    Xf_delayed(:,:,channel) = Xf(:,:,channel).*repmat(H(:,channel).*mems_tf(:,channel),1,param.num_frames);
%     Xf_delayed(:,:,channel) = Xf(:,:,channel);
end

X_buf = ifft(Xf_delayed,nfft);
X_buf = real(X_buf);

for channel = 1:param.M
    X_oa(:,channel) = overlap_add(X_buf(:,:,channel),param);
end

X_mismatch = X_oa(1:size(X,1),:);

function [Xf_buf, X_buf, param] = create_input_buffer(X,param)

wl = param.window_length;
frameshift = param.frameshift;
nfft = param.nfft;
M = param.M;

window = getWindow(param.window_type,wl);

for channel = 1:M
    X_temp = buffer(X(:,channel),wl,wl-frameshift,'nodelay');
    param.num_frames = size(X_temp,2);
    X_temp = X_temp.*repmat(window,1,param.num_frames);
    X_temp = [zeros(wl/2,param.num_frames); X_temp; zeros(wl/2,param.num_frames)];
    X_buf(:,:,channel) =  X_temp;
    Xf_buf(:,:,channel) = fft(X_buf(:,:,channel),nfft);
end

function y = overlap_add(x_buf,param)

wl = param.window_length;
frameshift = param.frameshift;
nfft = param.nfft;
num_frames = param.num_frames;

y = zeros(nfft + frameshift*(num_frames-1),1);

for frame_number = 1:num_frames
    y((frame_number-1)*(frameshift)+1:(frame_number-1)*frameshift+nfft) = ...
        y((frame_number-1)*(frameshift)+1:(frame_number-1)*frameshift+nfft) ...
        + x_buf(:,frame_number);
end

y = y((wl/2+1):(end-wl/2));

function window = getWindow(window_type,wl)
    
if strcmp(window_type,'hann')
    window = hann(wl,'periodic');
elseif strcmp(window_type,'rect')
    window = ones(wl,1);
end
