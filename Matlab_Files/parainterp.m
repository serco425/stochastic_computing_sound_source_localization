% parabolische Interpolation 
%https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html
%https://dspguru.com/dsp/howtos/how-to-interpolate-fft-peak/
%
% Wahre Position des globalen Maximums
%
% INPUT
% yl,y,yr ... 3 Punkt Eingangsfolge um das lokale Maximum
% 
% OUTPUT
% delta  relativer Versatz des wahren Max. zum detektierten Max.

function delta=parainterp(yl,y,yr)

delta=(yl-yr)./(2.*(yl-2*y+yr)); 
%pmax=p+delta; %