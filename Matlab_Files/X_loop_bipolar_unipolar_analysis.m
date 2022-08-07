close all;

%% Load and Xloop
load('fs_15k6_5degree_winkelabhängigkeit_extensive.mat','X','AzimuthValues','stdvalues','RevValues','SnrValues');
addpath('./stochastic_computing')

% -------------------------------------
% Initializing
Round_Setting = {'round'};
index_round = 1;
stdvalues.windowlength = 128; %overwrite the default windowlength
stochastic_mul_binary_add.type = 'MultiOrApproxBinary'; %to add up multiplication values
stochastic_mul_binary_add.pwm = 'LTC6992'; % model the LTC
full_stochastic.type = 'RotateSum'; % do thefault rotation
full_stochastic.pwm = 'LTC6992'; % model the LTC
bipolar.pwm = 'none';%'LTC6992';
bipolar.type = 'none';
Prime_select = [16,15];
LSB = 2^-3; %+1 sign bit! so in total 4_bit
maxLag = 3;
LSB_unipolar = 2^-4;
% The bipolar implementation calculates 0.5(x+1) and therefore reduces the
% precision before generating 16-bit (B=4) bit-streams
% With relatively prime periode lengths we get a systematic error dut to
% 0.5 cannot be represented with odd bitstream length!
% use clock division method in stead of relative prime method (comparable
% hardware)
% 
% with the LTC6992 we would need to shift it to 0.5 volts
% maybe there are pwm modules that accept a variable voltage range [-1,1]
% or [-5,5], similar to the ADC
%
% Because I don't have an implementation of clock division - i can select
% sobol in the bipolarPWMCorr function.

X = X*150; %scaling

maxLag_real = (2*stdvalues.radiusmean)/stdvalues.c*stdvalues.fs;
total_blocks = floor(length(X)/stdvalues.windowlength);
total_blocks = round(total_blocks/40); %for faster simulation
curPrime = Prime_select;
stochastic_mul_binary_add.curPrime = curPrime;
full_stochastic.curPrime = curPrime;

%results are stored here
Bipolar_LTC6992_yc_mat = zeros(total_blocks,length(AzimuthValues));
Bipolar_LTC6992_lag_mat = zeros(total_blocks,length(AzimuthValues));

%results are stored here
Unipolar_LTC6992_yc_mat = zeros(total_blocks,length(AzimuthValues));
Unipolar_LTC6992_lag_mat = zeros(total_blocks,length(AzimuthValues));

%results are stored here
LSBValues4Bit_yc_mat = zeros(total_blocks,length(AzimuthValues));
LSBValues4Bit_lag_mat = zeros(total_blocks,length(AzimuthValues));




for index_azimuth = 1:length(AzimuthValues)  % Simulating for different Signal to noise Ratios
    
for index_frame = 1:total_blocks
X_Frame = squeeze(X(1+(index_frame-1+1000)*stdvalues.windowlength:(index_frame+1000)*stdvalues.windowlength,:,index_azimuth)); %1.
    %-------------------------------------
    % Binary Correlation
    %------------------------------------- Reference
    sum_result = xcorr(myquant(X_Frame(:,1)',LSB,Round_Setting{index_round}),myquant(X_Frame(:,2)',LSB,Round_Setting{index_round}),maxLag);    
    [yc, lag] = max(sum_result);
    LSBValues4Bit_yc_mat(index_frame,index_azimuth) = yc;
    LSBValues4Bit_lag_mat(index_frame,index_azimuth) = lag;

    % Stochastic    ! BEI DEN QUANTISIERTEN NICHT AUserhalb quantisieren -
    % die werden eh intern quantisiert.

    sum_result = bipolarPWMCorr(X_Frame(:,1)',X_Frame(:,2)',curPrime(1),curPrime(2),maxLag,bipolar);
    [yc, lag] = max(sum_result);
    Bipolar_LTC6992_yc_mat(index_frame,index_azimuth) = yc;
    Bipolar_LTC6992_lag_mat(index_frame,index_azimuth) = lag;  


    %-------------------------------------
    % Clearing Negative Half Wave
    %-------------------------------------    
    X_Frame_pos = X_Frame;
    X_Frame_pos(X_Frame_pos<0) = 0;

    %-------------------------------------
    % Stochastic Computation Section
    %------------------------------- Stochastic Implementation
    sum_result = unipolarPWMCorr(X_Frame_pos(:,1)',X_Frame_pos(:,2)',curPrime(1),curPrime(2),maxLag,full_stochastic);
    [yc, lag] = max(sum_result);

    Unipolar_LTC6992_yc_mat(index_frame,index_azimuth) = yc;
    Unipolar_LTC6992_lag_mat(index_frame,index_azimuth) = lag;
    

end

end
clear X
save('fs_15k6_winkelabhaengigkeit_bipolar_unipolar_results.mat')


%% winkelabhaengigkeit_master mit average
pos_frame = 1;
pos_azimuth = 2;
pos_lsbvalues= 4;
round_dimensions = 1;
% to be compatible with the other simulations repmat the rounding
% simulation to also account for not having a truncation sim
%Unipolar_LTC6992_lag_mat=repmat(Unipolar_LTC6992_lag_mat,1,1,2);
%Unipolar_LTC6992_yc_mat=repmat(Unipolar_LTC6992_yc_mat,1,1,2);
%Bipolar_LTC6992_lag_mat=repmat(Bipolar_LTC6992_lag_mat,1,1,2);
%Bipolar_LTC6992_yc_mat=repmat(Bipolar_LTC6992_yc_mat,1,1,2);

Unipolar_LTC6992_winkelabhaengigkeit = zeros(total_blocks,length(AzimuthValues),round_dimensions); % LTCTransfer Function
Bipolar_LTC6992_winkelabhaengigkeit = zeros(total_blocks,length(AzimuthValues),round_dimensions); % LTCTransfer Function
LSBValues4Bit_winkelabhaengigkeit = zeros(total_blocks,length(AzimuthValues),round_dimensions); % LTCTransfer Function


mapfunction = @MapAndInterpolateAndAverage;
%mapfunction = @MapAndAverage;
maximum_of_mapping = 15;
init_average = (sin((AzimuthValues-90)*pi/180)*0.5+0.5)*maximum_of_mapping;
% VAD Setting
%Calculate the VAD threshold. Depends on scale of X (currently 300)
%Setting it too high would result in only few estimation (slow adapt).
%Too high = many "weak estimations".
vad_threshold = 0.3;    %not sure if this is the best (~0.3)
vad_threshold_or_based = vad_threshold/4; %To account for lower scaling

for index_azimuth = 1:length(AzimuthValues) 

bipolar_LTC6992_prev_estimation_Ref = ones(1,round_dimensions)*init_average(index_azimuth);
unipolar_LTC6992_prev_estimation_Ref = ones(1,round_dimensions)*init_average(index_azimuth);
LSBValues4Bit_prev_estimation_Ref = ones(1,round_dimensions)*init_average(index_azimuth);


for index_frame = 1:total_blocks
 
    %bipolar inputs
    new_lag = squeeze(LSBValues4Bit_lag_mat(index_frame, index_azimuth,:))';
    new_yc = squeeze(LSBValues4Bit_yc_mat(index_frame, index_azimuth,:))';
    LSBValues4Bit_prev_estimation_Ref = mapfunction(LSBValues4Bit_prev_estimation_Ref, new_yc, new_lag, vad_threshold);   


    %bipolar inputs
    new_lag = squeeze(Bipolar_LTC6992_lag_mat(index_frame, index_azimuth,:))';
    new_yc = squeeze(Bipolar_LTC6992_yc_mat(index_frame, index_azimuth,:))';
    bipolar_LTC6992_prev_estimation_Ref = mapfunction(bipolar_LTC6992_prev_estimation_Ref, new_yc, new_lag, vad_threshold);     

    
    %unipolar inputs
    new_lag = squeeze(Unipolar_LTC6992_lag_mat(index_frame, index_azimuth,:))';
    new_yc = squeeze(Unipolar_LTC6992_yc_mat(index_frame, index_azimuth,:))';
    unipolar_LTC6992_prev_estimation_Ref = mapfunction(unipolar_LTC6992_prev_estimation_Ref, new_yc, new_lag, vad_threshold_or_based);      
        
    Bipolar_LTC6992_winkelabhaengigkeit(index_frame,index_azimuth,:) = bipolar_LTC6992_prev_estimation_Ref;
    Unipolar_LTC6992_winkelabhaengigkeit(index_frame,index_azimuth,:) = unipolar_LTC6992_prev_estimation_Ref;
    LSBValues4Bit_winkelabhaengigkeit(index_frame,index_azimuth,:) = LSBValues4Bit_prev_estimation_Ref;
    
end
end

Unipolar_LTC6992_winkelabhaengigkeit_mean   	= squeeze(mean( Unipolar_LTC6992_winkelabhaengigkeit,pos_frame));
Bipolar_LTC6992_winkelabhaengigkeit_mean   	= squeeze(mean( Bipolar_LTC6992_winkelabhaengigkeit,pos_frame));
LSBValues4Bit_winkelabhaengigkeit_mean   	= squeeze(mean( LSBValues4Bit_winkelabhaengigkeit,pos_frame));


%% Plotting to verify


ideal = AzimuthValues-90;
cur_data = [LSBValues4Bit_winkelabhaengigkeit_mean',Bipolar_LTC6992_winkelabhaengigkeit_mean',Unipolar_LTC6992_winkelabhaengigkeit_mean'];


data_function = real(asin((cur_data-7.5)./7.5).*180./pi);
data_function = [data_function]-repmat(ideal',1,size(data_function,2));


width = 3.5;    %allowed inches in ieee access
height = 3.0;
font_size = 10;  %footnotesize
fig = figure('Units','inches',...
'Position',[0 0 width height], ...
'PaperPositionMode','auto');
repmat_ideal = repmat(ideal',1,size(data_function,2));
h=plot(repmat_ideal,data_function);
hold on;
xlabel_str = 'Source Location ($^{\circ}$)';
ylabel_str = 'Mean Estimation Error ($^{\circ}$)';
set(gca,'YColor', 'k');
%ylim(gca,[-20,20]);
xlim(gca,[-90,90]);
%xticks([-90,-45,0,45,90])
plot(get(gca,'xlim'),[0 0],'k--');
%yticks([0,45,90,135,180])
legend(h,{'LSBValues4Bit','Bipolar','Unipolar'})
xlabel(xlabel_str,...
'Interpreter','latex')
ylabel(ylabel_str,...
'Interpreter','latex')