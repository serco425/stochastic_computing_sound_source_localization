load('35_degree_short.mat','X','stdvalues','RevValues','SnrValues');

addpath(fullfile(pwd,'stochastic_computing'))

% -------------------------------------
% Initializing
total_blocks = floor(length(X)/stdvalues.windowlength);
LSBValues = [2^-4,2^-5,2^-6, 2^-7,2^-8]; 
stochastic_mul_binary_add.type = 'MultiOrApproxBinary'; %to add up multiplication values
maxLag = 4;

phi_binary_quant = zeros(total_blocks,length(LSBValues), length(SnrValues),length(RevValues));
phi_stochastic_full_wave = zeros(total_blocks,length(LSBValues), length(SnrValues),length(RevValues));
phi_stochastic_none = zeros(total_blocks,length(LSBValues),length(SnrValues),length(RevValues));
phi_stochastic_multiOrApproxBinary = zeros(total_blocks,length(LSBValues),length(SnrValues),length(RevValues));
phi_stochastic_half_wave = zeros(total_blocks,length(LSBValues), length(SnrValues),length(RevValues));

for revindex = 1:length(RevValues)
for snrindex = 1:length(SnrValues)  
 
    LSB = LSBValues(lsbindex);
    Prime_select = [[251,241];[127,98];[61,53];[31,23];[16,13]];
    if LSB == 2^-8
        curPrime = Prime_select(1,:);
    elseif LSB == 2^-7
        curPrime = Prime_select(2,:);
    elseif LSB == 2^-6
        curPrime = Prime_select(3,:);
    elseif LSB == 2^-5
        curPrime = Prime_select(4,:);
    elseif LSB == 2^-4
        curPrime = Prime_select(5,:);
    else
        clear curPrime; % throw an error
    end
    stochastic_mul_binary_add.curPrime = curPrime;
    full_stochastic.curPrime = curPrime;
    full_stochastic.type = 'RotateSum'; % do thefault rotation
    
for index_frame = 1:total_blocks
X_Frame = squeeze(X(1+(index_frame-1)*stdvalues.windowlength:(index_frame)*stdvalues.windowlength,:,snrindex,revindex)); %1.
    %-------------------------------------
    % Binary Correlation
    %------------------------------------- Reference
    sum_result = My_Corr(myquant(X_Frame(:,1)',LSB),myquant(X_Frame(:,2)',LSB),maxLag);
    [yc, lag] = max(sum_result);
    if lag~=1 && lag ~= length(sum_result) %avoid endpoints
        yl = sum_result(lag-1);
        yr = sum_result(lag+1);

        delta=(yl-yr)./(2.*(yl-2*yc+yr));  
    else
        delta = 0;
    end 
    TapsLag = maxLag+1-lag-delta;
    phi_binary_quant(index_frame,lsbindex,snrindex,revindex) = 90-asin(TapsLag*1/stdvalues.fs*stdvalues.c/(2*stdvalues.radiusmean))*180/pi;
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
    else
        delta = 0;
    end 
    TapsLag = maxLag+1-lag-delta;
    phi_stochastic_multiOrApproxBinary(index_frame,lsbindex,snrindex,revindex) = 90-asin(TapsLag*1/stdvalues.fs*stdvalues.c/(2*stdvalues.radiusmean))*180/pi;
    %-------------------------------------
    % Stochastic Computation Section
    %------------------------------- Stochastic Implementation
    sum_result = unipolarPWMCorr(X_Frame_pos(:,1)',X_Frame_pos(:,2)',curPrime(1),curPrime(2),maxLag,full_stochastic);
    [yc, lag] = max(sum_result);
    if lag~=1 && lag ~= length(sum_result) %avoid endpoints
        yl = sum_result(lag-1);
        yr = sum_result(lag+1);

        delta=(yl-yr)./(2.*(yl-2*yc+yr));  
    else
        delta = 0;
    end 
    TapsLag = maxLag+1-lag-delta;
    phi_stochastic_half_wave(index_frame,lsbindex,snrindex,revindex) = 90-asin(TapsLag*1/stdvalues.fs*stdvalues.c/(2*stdvalues.radiusmean))*180/pi;

end
end
end
end