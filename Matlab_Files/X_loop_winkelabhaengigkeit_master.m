clear all;
load('fs_15k6_5degree_winkelabhängigkeit_extensive.mat','X','AzimuthValues','stdvalues','RevValues','SnrValues');
addpath(fullfile(pwd,'stochastic_computing'))
% do for filtered and unfiltered 
% decide if i want to have the n and n-1 quantization effect

% -------------------------------------
% Initializing
stdvalues.windowlength = 128; %overwrite the default windowlength
%LTC_Transfer.type = 'RotateSum'; % do thefault rotation
LTC_Transfer.pwm = 'LTC6992'; % model the LTC
LTC_Transfer.type = 'CyclicBuffer';
LTC_Transfer.rotate = 254;

NO_Transfer = LTC_Transfer;
NO_Transfer.pwm = 'None'; %anything except LTC6992 works

Prime_select = [16,15];
LSB_unipolar = 2^-4;
maxLag = 3;


LSBValues = [2^-2,2^-3,2^-4,2^-7]; % That is 3Bit 4Bit 5Bit 8 Bit
%Also count the Sign bit!
%wenn ich so niedrig gehe muss ich schon drauf achten, dasse ein input auf
%4 und der andere auf 3 quantisiert wird...
Round_Setting = {'truncate','round'}; 
%Round_Setting = {'round'};
truncate_index = find(contains(Round_Setting,'truncate'));
round_index = find(contains(Round_Setting,'round'));

% %A scaling of 300 limits the input signal to +-1 - However
Scaling = 300; %150
X = X*Scaling; %A scaling of 300 limits the input signal to +-1 - However
% the maximum of the cross crorelation results is 12 and the median is 0.3
%##########################################################################################################################
%##########################################################################################################################
%##########################################################################################################################
%X_filtered = X;
%aliasing_frequency = stdvalues.c/2/stdvalues.IndPairDists;
%for ind = 1:size(X,3)
%   X_filtered(:,:,ind) =  lowpass(squeeze(X(:,:,ind)),aliasing_frequency,stdvalues.fs);
%end
%
%X = X_filtered;
%%%%
%
%##########################################################################################################################
%##########################################################################################################################
%##########################################################################################################################
total_blocks = round(floor(length(X)/stdvalues.windowlength));


%results are stored here
LSBValues_Ref_yc_mat = zeros(total_blocks,length(AzimuthValues),length(Round_Setting),length(LSBValues));
LSBValues_Ref_lag_mat = zeros(total_blocks,length(AzimuthValues),length(Round_Setting),length(LSBValues));

%results are stored here
Double_yc_mat = zeros(total_blocks,length(AzimuthValues),length(Round_Setting));
Double_yl_mat = zeros(total_blocks,length(AzimuthValues),length(Round_Setting));
Double_yr_mat = zeros(total_blocks,length(AzimuthValues),length(Round_Setting));

Double_lag_mat = zeros(total_blocks,length(AzimuthValues),length(Round_Setting));

%results are stored here
Unipolar_LTC6992_yc_mat = zeros(total_blocks,length(AzimuthValues),length(Round_Setting));
Unipolar_LTC6992_lag_mat = zeros(total_blocks,length(AzimuthValues),length(Round_Setting));

%results are stored here
Unipolar_NO_Transfer_yc_mat = zeros(total_blocks,length(AzimuthValues),length(Round_Setting));
Unipolar_NO_Transfer_lag_mat = zeros(total_blocks,length(AzimuthValues),length(Round_Setting));

Unipolar_counter_lag_mat = zeros(total_blocks,length(AzimuthValues),length(Round_Setting));
Unipolar_counter_yc_mat = zeros(total_blocks,length(AzimuthValues),length(Round_Setting));

for index_round = 1:length(Round_Setting)
for index_azimuth = 1:length(AzimuthValues)
    
for index_frame = 1:total_blocks
X_Frame = squeeze(X(1+(index_frame-1)*stdvalues.windowlength:(index_frame)*stdvalues.windowlength,:,index_azimuth)); %1.
    %-------------------------------------
    % Binary Correlation - Bit Dependence
    %------------------------------------- 
    for lsbindex = 1:length(LSBValues)
        sum_result = xcorr(myquant(X_Frame(:,1)',LSBValues(lsbindex),Round_Setting{index_round}),myquant(X_Frame(:,2)',LSBValues(lsbindex),Round_Setting{index_round}),maxLag);    
        [yc, lag] = max(sum_result);
        LSBValues_Ref_yc_mat(index_frame,index_azimuth,index_round,lsbindex) = yc;
        LSBValues_Ref_lag_mat(index_frame,index_azimuth,index_round,lsbindex) = lag;
    end
    
   
    %-------------------------------------
    % Double Precision
    %-------------------------------------     
    sum_result = xcorr(X_Frame(:,1)',X_Frame(:,2)',maxLag);
    [yc, lag] = max(sum_result); %there will always be only 1 maximum
    Double_yc_mat(index_frame,index_azimuth,index_round) = yc;  %same for truncate or round
    if lag~=1 && lag ~= length(sum_result) %avoid endpoints
        Double_yl_mat(index_frame,index_azimuth,index_round) = sum_result(lag-1);  %same for truncate or round
        Double_yr_mat(index_frame,index_azimuth,index_round) = sum_result(lag+1);  %same for truncate or round
    else
        Double_yl_mat(index_frame,index_azimuth,index_round) = nan;  %same for truncate or round
        Double_yr_mat(index_frame,index_azimuth,index_round) = nan;  %same for truncate or round        
    end    
    Double_lag_mat(index_frame,index_azimuth,index_round) = lag;
     
    %-------------------------------------
    % Clearing Negative Half Wave
    %-------------------------------------    
    X_Frame_pos = X_Frame;
    X_Frame_pos(X_Frame_pos<0) = 0;

    %-------------------------------------
    % Exact Quantized Computation with positive half waves / Study effect
    % of Half wave rectifier
    %-------------------------------------      
    sum_result = xcorr(myquant(X_Frame_pos(:,1)',LSB_unipolar,Round_Setting{index_round}),myquant(X_Frame_pos(:,2)',LSB_unipolar,Round_Setting{index_round}),maxLag);
    [yc, lag] = max(sum_result);
    Unipolar_counter_yc_mat(index_frame,index_azimuth,index_round) = yc;
    Unipolar_counter_lag_mat(index_frame,index_azimuth,index_round) = lag;
 
    %-------------------------------------
    % Stochastic Computation Section
    %------------------------------- Stochastic Implementation
    if Scaling==150 && index_round==round_index %For more speed!
        % Inupts are quantisized internally
        sum_result = unipolarPWMCorr(X_Frame_pos(:,1)',X_Frame_pos(:,2)',Prime_select(1),Prime_select(2),maxLag,LTC_Transfer); %quantization happens internally
        [yc, lag] = max(sum_result);
        Unipolar_LTC6992_yc_mat(index_frame,index_azimuth,index_round) = yc;
        Unipolar_LTC6992_lag_mat(index_frame,index_azimuth,index_round) = lag;
        
        
        sum_result = unipolarPWMCorr(X_Frame_pos(:,1)',X_Frame_pos(:,2)',Prime_select(1),Prime_select(2),maxLag,NO_Transfer); %quantization happens internally
        [yc, lag] = max(sum_result);
        Unipolar_NO_Transfer_yc_mat(index_frame,index_azimuth,index_round) = yc;
        Unipolar_NO_Transfer_lag_mat(index_frame,index_azimuth,index_round) = lag;
    end

end
end
end

clear X X_filtered
save('fs_15k6_5degree_winkelabhängigkeit_extensive_results.mat')
%save('fs_15k6_5degree_winkelabhängigkeit_extensive_scaling_150_results.mat')

