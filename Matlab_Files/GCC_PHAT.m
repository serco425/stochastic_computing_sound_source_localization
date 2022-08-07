function [ phi, theta ] = GCC_PHAT( X_Frame , stdvalues, optional)
%% DOA Estimation by Minimization of the RMS TDOA Errors
% Reference 2000-DiBiaseThesis.pdf
% Page 59/122
%% Inputs: 
%%% 1) X_frame is a frame of Data (Window_Length x Microphones)
%%% 2) stdvalues contains information about Array Geometry, search resolution,
%%% sampling rate
%%% 3) optional.phat = true (default) states if phat weighing should be
%%% used
%% Outputs:
%%% 2) phi: azimuth
%%% 3) theta: elevation
%% 
theta = 0;

Window = hann(length(X_Frame),'periodic'); % Window constant
X_Frame = X_Frame.*Window; % applying window


if (nargin == 2) 
    optional.phat = true;
    optional.interpolate = true;
end

%init
Xf_Frame = fft(X_Frame,length(X_Frame));

%Calculation
Cross_Spectrum = Xf_Frame(:,1).*conj(Xf_Frame(:,2));
%-------------------------------------
% Phat_Weighing in Frequency Domain
%-------------------------------------
       if optional.phat
           Psi = (abs(Cross_Spectrum)+eps);
       else
           Psi = 1;
       end
GCC_PHAT = ifftshift(ifft((Cross_Spectrum./(Psi))));
MaxLag = stdvalues.maxlagcircular7;
[yc, lag] = max(GCC_PHAT(end/2-MaxLag:end/2+MaxLag));
lag = lag+(length(GCC_PHAT)/2-MaxLag);
if optional.interpolate
    yl = GCC_PHAT(lag-1);
    yr = GCC_PHAT(lag+1);

    delta=(yl-yr)./(2.*(yl-2*yc+yr));  
else
    delta = 0;
end

TapsLag = (length(GCC_PHAT))/2-lag-delta+1;

%% Time Delay to Angle trough trigonometric function

constant_val = 1/stdvalues.fs*stdvalues.c/(2*stdvalues.radiusmean);
phi = 90-asin(TapsLag*constant_val)*180/pi;

  
end

