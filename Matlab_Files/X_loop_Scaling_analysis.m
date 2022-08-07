clear all;
% load('fs_12k.mat','X','stdvalues','RevValues','SnrValues');
load('fs_12k_winkelabhängigkeit_extensive.mat','X','AzimuthValues','stdvalues','RevValues','SnrValues');
addpath(fullfile(pwd,'stochastic_computing'))

% -------------------------------------
% Initializing
stdvalues.windowlength = 128; %overwrite the default windowlength
max_input = max(X,[],'all');
ScaleValues = [1/max_input/10,1/max_input/5, 1/max_input, 1/max_input*5,1/max_input*10]; 
% X(X<0) = 0;
stochastic_mul_binary_add.type = 'MultiOrApproxBinary'; %to add up multiplication values
stochastic_mul_binary_add.pwm = 'LTC6992'; % model the LTC
full_stochastic.type = 'RotateSum'; % do thefault rotation
full_stochastic.pwm = 'LTC6992'; % model the LTC
Prime_select = [16,15];
LSB = 2^-4;
maxLag = 3;


maxLag_real = (2*stdvalues.radiusmean)/stdvalues.c*stdvalues.fs;
total_blocks = floor(length(X)/stdvalues.windowlength);
curPrime = Prime_select;
stochastic_mul_binary_add.curPrime = curPrime;
full_stochastic.curPrime = curPrime;

%results are stored here
phi_binary_quant = zeros(total_blocks,length(ScaleValues), length(SnrValues),length(RevValues));
phi_stochastic_full_wave = zeros(total_blocks,length(ScaleValues), length(SnrValues),length(RevValues));
phi_stochastic_none = zeros(total_blocks,length(ScaleValues),length(SnrValues),length(RevValues));
phi_stochastic_multiOrApproxBinary = zeros(total_blocks,length(ScaleValues),length(SnrValues),length(RevValues));
phi_stochastic_half_wave = zeros(total_blocks,length(ScaleValues), length(SnrValues),length(RevValues));

max_sum_result = zeros(total_blocks,length(ScaleValues));

for revindex = 1:length(RevValues)  % Simulating for different reveberation
for azimuthindex = 1:length(AzimuthValues)  % Simulating for different Signal to noise Ratios
for scale_index = 1:length(ScaleValues) 

    cur_scale = ScaleValues(scale_index);
    
for index_frame = 1:total_blocks
X_Frame = squeeze(X(1+(index_frame-1)*stdvalues.windowlength:(index_frame)*stdvalues.windowlength,:,azimuthindex,revindex)); %1.
X_Frame = X_Frame*cur_scale;    
    %-------------------------------------
    % Binary Correlation
    %------------------------------------- Reference
    sum_result = My_Corr(myquant(X_Frame(:,1)',LSB),myquant(X_Frame(:,2)',LSB),maxLag);
    [yc, lag] = max(sum_result);
    max_sum_result(index_frame,scale_index) = yc;
    if lag~=1 && lag ~= length(sum_result) %avoid endpoints
        yl = sum_result(lag-1);
        yr = sum_result(lag+1);

        delta=(yl-yr)./(2.*(yl-2*yc+yr));
        TapsLag = maxLag+1-lag-delta;
    elseif lag == 1 % unpossible delay >180°
        TapsLag = maxLag_real;
    elseif lag ==(length(sum_result)) % unpossible delay <0°
        TapsLag = -maxLag_real;
    end 
    
    
    if max(sum_result) ~= 0 %no voice detected
        phi_binary_quant(index_frame,scale_index,azimuthindex,revindex) = 90-asin(TapsLag*1/stdvalues.fs*stdvalues.c/(2*stdvalues.radiusmean))*180/pi;
    else %store invalid when nothing detected
        phi_binary_quant(index_frame,scale_index,azimuthindex,revindex) = -180;
    end
    %-------------------------------------
    % Clearing Negative Half Wave
    %-------------------------------------    
    X_Frame_pos = X_Frame;
    X_Frame_pos(X_Frame_pos<0) = 0;

    %-------------------------------------
    % Stochastic Computation Section
    %-------------------------------Best achievable with Stochastic OR
    sum_result = unipolarPWMCorr(X_Frame_pos(:,1)',X_Frame_pos(:,2)',curPrime(1),curPrime(2),maxLag,stochastic_mul_binary_add);
    [yc, lag] = max(sum_result);
    if lag~=1 && lag ~= length(sum_result) %avoid endpoints
        yl = sum_result(lag-1);
        yr = sum_result(lag+1);

        delta=(yl-yr)./(2.*(yl-2*yc+yr));
        TapsLag = maxLag+1-lag-delta;
    elseif lag == 1 % unpossible delay >180°
        TapsLag = maxLag_real;
    elseif lag ==(length(sum_result)) % unpossible delay <0°
        TapsLag = -maxLag_real;
    end 

    
    if max(sum_result) ~= 0 %no voice detected
        phi_stochastic_multiOrApproxBinary(index_frame,scale_index,azimuthindex,revindex) = 90-asin(TapsLag*1/stdvalues.fs*stdvalues.c/(2*stdvalues.radiusmean))*180/pi;
    else %store invalid when nothing detected
        phi_stochastic_multiOrApproxBinary(index_frame,scale_index,azimuthindex,revindex) = -180;
    end
    %-------------------------------------
    % Stochastic Computation Section
    %------------------------------- Stochastic Implementation
    sum_result = unipolarPWMCorr(X_Frame_pos(:,1)',X_Frame_pos(:,2)',curPrime(1),curPrime(2),maxLag,full_stochastic);
    [yc, lag] = max(sum_result);
    if lag~=1 && lag ~= length(sum_result) %avoid endpoints
        yl = sum_result(lag-1);
        yr = sum_result(lag+1);

        delta=(yl-yr)./(2.*(yl-2*yc+yr));
        TapsLag = maxLag+1-lag-delta;
    elseif lag == 1 % unpossible delay >180°
        TapsLag = maxLag_real;
    elseif lag ==(length(sum_result)) % unpossible delay <0°
        TapsLag = -maxLag_real;
    end 

    
    if max(sum_result) ~= 0 %no voice detected
       temp = 90-asin(TapsLag*1/stdvalues.fs*stdvalues.c/(2*stdvalues.radiusmean))*180/pi;
        phi_stochastic_half_wave(index_frame,scale_index,azimuthindex,revindex) = 90-asin(TapsLag*1/stdvalues.fs*stdvalues.c/(2*stdvalues.radiusmean))*180/pi;
    else %store invalid when nothing detected
        phi_stochastic_half_wave(index_frame,scale_index,azimuthindex,revindex) = -180;
    end
end
end
end
end

% save('X_loop_winkel_scale.mat');
% save('X_loop_4bit.mat');
% % all valid results are between 0 and 180°
% % all invalid results = -180;