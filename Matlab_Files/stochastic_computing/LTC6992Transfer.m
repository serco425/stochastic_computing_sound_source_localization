function [out] = LTC6992Transfer(mat)
%LTC6992Transfer Takes an input matrix or vector and does the same
% nonlinear transformation than a LTC6992 does before generating pwm
%   (val-0.1).*1.25; and with saturation

out = (mat-0.1).*1.25;
out(out>1) = 1;
out(mat<0.1) = 0;   %zero all that were below 0.1

end

