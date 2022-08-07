function [ phi, theta ] = ErrorCriterion( X_Frame, stdvalues , optional)
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

%% copy begin
% Cross_Spectrum = zeros(stdvalues.nfft,stdvalues.indpairs);
% GCC_PHAT = zeros(stdvalues.nfft,stdvalues.indpairs);
% Xf_Frame = fft(X_Frame,stdvalues.nfft);
% ind = 1;
% for i1 = 1:stdvalues.m
% 
%     for i2 = i1+1:stdvalues.m
%         
%        Cross_Spectrum(:,ind) = Xf_Frame(:,i1).*conj(Xf_Frame(:,i2));
%         %-------------------------------------
%         % Phat_Weighing in Frequency Domain
%         %-------------------------------------
%        Psi = (abs(Cross_Spectrum(:,ind))+eps);
%   
%        GCC_PHAT(:,ind) = ifftshift(ifft((Cross_Spectrum(:,ind)./(Psi))));
%        ind = ind +1;
%        
%     end
% end
% 
% [yc, lag] = max(GCC_PHAT);
% yl = GCC_PHAT(lag-1);
% yr = GCC_PHAT(lag+1);
% 
% delta=(yl-yr)./(2.*(yl-2*yc+yr)); 
% TapsLag2 = (length(GCC_PHAT))/2-lag-delta+1;
% 
% constant_val = 1/stdvalues.fs*stdvalues.c/(2*stdvalues.radiusmean);
% phi = asin(TapsLag2*constant_val)*180/pi
%copy end


if (nargin == 2) 
    optional.phat = true;
end


Cross_Spectrum = zeros(stdvalues.nfft,stdvalues.indpairs);
Generalized_Cross_Correlation = zeros(stdvalues.nfft,stdvalues.indpairs);
Xf_Frame = fft(X_Frame,stdvalues.nfft);

ind = 1;
for i1 = 1:stdvalues.m

    for i2 = i1+1:stdvalues.m
        
       Cross_Spectrum(:,ind) = Xf_Frame(:,i1).*conj(Xf_Frame(:,i2));
        %-------------------------------------
        % Phat_Weighing
        %-------------------------------------
       if (optional.phat == true)
           Psi = (abs(Cross_Spectrum(:,ind))+eps);
       else
           Psi = 1;
       end
       Generalized_Cross_Correlation(:,ind) = ifftshift(ifft((Cross_Spectrum(:,ind)./(Psi))));
       ind = ind +1;
       
    end
end
  

[MaxValue, MaxIndex] = max(Generalized_Cross_Correlation,[],1);

ylvec = Generalized_Cross_Correlation(sub2ind(size(Generalized_Cross_Correlation),MaxIndex-1,1:length(MaxIndex)));
yvec =  MaxValue;
yrvec = Generalized_Cross_Correlation(sub2ind(size(Generalized_Cross_Correlation),MaxIndex+1,1:length(MaxIndex)));
delta_time = parainterp(ylvec,yvec,yrvec);

TapsLag = MaxIndex+delta_time-2*stdvalues.windowlength-1;

Error_Square = zeros(length(stdvalues.PhiAngles)-1,length(stdvalues.ThetaAngles)-1);

for iPhi = 1:length(stdvalues.PhiAngles)-1
    
   
    for iTheta = 1:length(stdvalues.ThetaAngles)
        
       
        Error_Square(iPhi,iTheta) = ErrorFar(stdvalues.PhiAngles(iPhi),stdvalues.ThetaAngles(iTheta),TapsLag,stdvalues);
        
        
    end
    
    
end

[phiind, thetaind] = findMinimum (Error_Square);

phi = stdvalues.PhiAngles(phiind)*180/pi;
theta = stdvalues.ThetaAngles(thetaind)*180/pi;

end

