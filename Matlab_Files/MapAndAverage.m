function [cur_estimation] = MapAndAverage(prev_estimation, new_yc, new_lag)

VAD_Threshold = 16; % values below 32 will not be processed [0,256];
VAD_Threshold_QFormat = VAD_Threshold/(2^8);
alpha = 2^-4;
cur_mapped_vec = linspace(0,15,7);

%y[n]=αx[n]+(1−α)y[n−1]
cur_mapped = cur_mapped_vec(new_lag);

cur_estimation = prev_estimation; %default assignement;
cur_estimation(new_yc >= VAD_Threshold_QFormat) = alpha*cur_mapped(new_yc >= VAD_Threshold_QFormat)  + (1-alpha)*prev_estimation(new_yc >= VAD_Threshold_QFormat); % update those that are above threshold

%cur_estimation(new_yc < VAD_Threshold_QFormat) = prev_estimation(new_yc < VAD_Threshold_QFormat); % keep those taht are smaller
end