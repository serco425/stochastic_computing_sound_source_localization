function [ minindphi, minindtheta ] = findMaximum( Error_Square)
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


[MinValTheta, MinIndTheta] = max(Error_Square, [], 1);
[MinValPhi, MinIndPhi] = max(Error_Square, [], 2);

[minvaltheta, minindtheta] = max(MinValTheta);
[minvalphi, minindphi] = max(MinValPhi);

end




