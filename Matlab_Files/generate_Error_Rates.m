function [ Error_Rates ] = generate_Error_Rates( Angle_Mat, degrees, optional )
%Angle_Error_Rate The error rate will be used to evaluate performance of estimation techniques. In general, the error rate is the percentage of estimation with error greater or equal to some threshold, plotted as function of that threshold.
%% Inputs: 
%%% 1) Angle_Mat is (Frame_Index x SNR x Reverberation)
%%% 2) degrees Error Rate will be calculated for all values in degrees
%% Outputs:
%%% 1) Error_Rates: (Error_Rate for degrees x SNR x Reverberation)
%%

MatSize = size(Angle_Mat);
MatSize(1) = length(degrees);
Error_Rates = zeros(MatSize);

if (nargin == 2) 
    optional.offset = 0;
end


for ind_dim2 = 1:size(Error_Rates,2)
    
for ind_dim3 = 1:size(Error_Rates,3)    

for ind_rate = 1:length(degrees);
    Error_Rates(ind_rate,ind_dim2,ind_dim3) = length(  Angle_Mat(abs(Angle_Mat(:,ind_dim2,ind_dim3))-optional.offset  >=degrees(ind_rate),ind_dim2,ind_dim3)  )/size(Angle_Mat,1);
end

end

end


end

