function [ minindphi, minindtheta ] = findMinimum( Error_Square)
%% Returns the index of minimum of Error_Square
%% Inputs: 
%%% 1) Error_Square Matrix to find minimum
%% Outputs:
%%% 2) minindphi: azimuth
%%% 3) theta: elevation
%% 
% Error_Sqare
% Phi   x   Theta
% Zeile x   Spalte
% ylabelx   xlabel

threshhold = 10^(-10);

[MinValTheta, MinIndTheta] = min(Error_Square, [], 1);
[MinValPhi, MinIndPhi] = min(Error_Square, [], 2);

[minvaltheta, minindtheta] = min(MinValTheta);
[minvalphi, minindphi] = min(MinValPhi);

end




