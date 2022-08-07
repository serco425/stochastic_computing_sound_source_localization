function [cur_estimation] = MapAndInterpolateAndAverage(prev_estimation, new_yc, new_lag, vad_threshold , new_yl,new_yr)


if nargin <4
    %VAD_Threshold = 16; % values below 32 will not be processed [0,256];
    %vad_threshold = VAD_Threshold/(2^8);
    vad_threshold = 16/2^8;
    disp('vad_threshold_not_set');
end


if nargin <5
    new_yl = nan;
    new_yr = nan;
end


if ~(isempty(find(isnan(new_yl)))|| isempty(find(isnan(new_yr))))
       cur_mapped = 2.5*(new_lag-1);
else
        delta=(new_yl-new_yr)./(2.*(new_yl-2*new_yc+new_yr));
        %we don't have to use switch mapping we can multiply [0,6] x 2.5    
        cur_mapped = 2.5.*((new_lag-1)+delta);               
end

    %y[n]=αx[n]+(1−α)y[n−1]
alpha = 2^-4;
cur_estimation = prev_estimation; %default assignement;
cur_estimation(new_yc >= vad_threshold) = reshape(alpha.*cur_mapped(new_yc >= vad_threshold),1,[])  + reshape((1-alpha).*prev_estimation(new_yc >= vad_threshold),1,[]); % update those that are above threshold

%cur_estimation(new_yc < vad_threshold) = prev_estimation(new_yc < vad_threshold); % keep those taht are smaller
end

%{
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
%}