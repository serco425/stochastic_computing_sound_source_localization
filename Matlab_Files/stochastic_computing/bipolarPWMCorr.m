function [sum_result] = bipolarPWMCorr(x,y,x_pwm_len,y_pwm_len,maxLag,optional)
%bipolarPWMCorr Stochastic correlation using xnor and mux
%   Inputs: Binary vectors x and y of length vecLength [example 128, 258, 512]
%           pwm_len: should be relative prime Output Len =
%               x_pwm_len*y_pwm_len
%           maxLag: similar to xcorr function - the center +-maxLag will be
%               calculated
%   Output: Returns binary vector of 2*maxLag+1 with center bins

% optional generate sobel rand here
if ~exist('optional','var') %
    optional.type = 'none';
end

if strcmp(optional.type,'MultiMuxAdd')
    add_len = x_pwm_len*y_pwm_len*LCM_ad_multiple;
    total_len = add_len*(2*maxLag+1);
    Bitwidth = nextpow2(length(x));
    p = sobolset(Bitwidth,'Skip',1e3);
    sobolseq = (net(p,total_len))';
    scale = 2.^(0:Bitwidth-1);
    sobolseq_round = round(sobolseq);
    
    scaled_selvec = sobolseq_round.*transpose(scale);
    selvec = sum(scaled_selvec);
end


x_inp = [x zeros(1,maxLag)];
y_inp = [y zeros(1,maxLag)];

sum_result = zeros(1,2*maxLag+1);
vecLength = length(x);

% SOBOL
sobol = 0;
if sobol ==1
    p = sobolset(2);
    N_sobol_length = 256;
    sobol_rand = net(p,N_sobol_length);


    x_sobol = sobolVec(myquant(1/2*(1+x_inp),2^(-log2(x_pwm_len))),sobol_rand(:,1))';
    y_sobol = sobolVec(myquant(1/2*(1+y_inp),2^(-log2(x_pwm_len))),sobol_rand(:,2))';

    %mul_sobol = zeros(vecLength,N_sobol_length);    
else
    x_sc = pwmVec(x_inp,x_pwm_len,optional);
    y_sc = pwmVec(y_inp,y_pwm_len,optional);
    %mul_sc = zeros(vecLength,size(x_sc,2)*size(y_sc,2));
end

for curLag = -maxLag:+maxLag
   
    if curLag > 0
        shift_x = curLag;
        shift_y = 0;
    else   
        shift_x = 0;
        shift_y = -curLag;
    end


    if sobol ==1
        
        mul_sobol =  scXNOrMul(x_sobol((1+shift_x):(shift_x+vecLength),:),   y_sobol((1+shift_y):(shift_y+vecLength),:),N_sobol_length);

        %for ind = 1:vecLength
        %    x_sobol(ind,:) = sobolVec(myquant(1/2*(1+x_inp(ind+shift_x)),2^(-log2(x_pwm_len))),sobol_rand(:,1));
        %    y_sobol(ind,:) = sobolVec(myquant(1/2*(1+y_inp(ind+shift_y)),2^(-log2(x_pwm_len))),sobol_rand(:,2));
        %end
           mul = mul_sobol;
    else
        
        mul_sc = scXNOrMul(x_sc((1+shift_x):(shift_x+vecLength),:),   y_sc((1+shift_y):(shift_y+vecLength),:),x_pwm_len*y_pwm_len);
        %for ind = 1:vecLength
        %    mul_sc(ind,:) = scXNOrMul(x_sc((ind+shift_x),:),y_sc(ind+shift_y,:),x_pwm_len*y_pwm_len);
        %end
        mul = mul_sc;
    end

    if strcmp(optional.type,'MultiMuxAdd') 
        sum_sc = scSobolMultiAdd(mul,add_len,selvec(1+(curLag+maxLag)*add_len:(1+curLag+maxLag)*add_len))./128;
        sum_result(curLag+maxLag+1) = Bipolar2Binary(sum_sc);
    else %for example counter
        % add by using a counter
        sum_result(curLag+maxLag+1) = sum(mul,'all')/(0.5*size(mul,2))-128; % equals x1y1/240*2-1 + x2y2/240*2-1 ... basically individuall bipolar to binary conversions
    end


end

end

