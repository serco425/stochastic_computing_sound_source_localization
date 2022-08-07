function [ phi, theta ] = ClosedFormEquationError( X_Frame, stdvalues, optional )
%% DOA Estimation by Minimization Equation Error of Pytagorean Theorem
% Moves Microphone m2 into origin and minimizes for each m2mi i = 1...6
% Pair the Pytagorean theorem angle
% Currently 1D Azimuth Estimation
% Short version of "Closed-Form Least-Squares Source Location Estimation
% from Range Difference Measurements"
%% Inputs: 
%%% 1) X_frame is a frame of Data (Window_Length x Microphones)
%%% 2) stdvalues contains information about Array Geometry, search resolution,
%%% sampling rate
%% Outputs:
%%% 1) phi: azimuth
%%% 2) theta: CURRENTLY ALWAYS 0
%% 

if (nargin == 2) 
    optional.phat = true;
end


Cross_Spectrum = zeros(stdvalues.nfft,stdvalues.indpairs);
Generalized_Cross_Correlation = zeros(stdvalues.nfft,stdvalues.m-1);
Xf_Frame = fft(X_Frame,stdvalues.nfft);

ind = 1;
for i1 = 1:1 
%m2-m3
%m2-m4
%m2-m5
%m2-m6
%m2-m7

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

% R = transpose(stdvalues.IndPairDists(1:stdvalues.m-1));
[MaxValue, MaxIndex] = max(Generalized_Cross_Correlation(:,1:stdvalues.m-1),[],1);
% prepare values for parainterp

ylvec = Generalized_Cross_Correlation(sub2ind(size(Generalized_Cross_Correlation),MaxIndex-1,1:5));
yvec = MaxValue;
yrvec = Generalized_Cross_Correlation(sub2ind(size(Generalized_Cross_Correlation),MaxIndex+1,1:5));
delta_time = parainterp(ylvec,yvec,yrvec);

d = -transpose((MaxIndex+delta_time-2*stdvalues.windowlength-1)*stdvalues.c/stdvalues.fs); % - 1
S = transpose(stdvalues.Array_Geo(:,2:end)) - repmat(transpose(stdvalues.Array_Geo(:,1)),stdvalues.m-1,1); % move point m2
%delta = R.^2-d.^2;
Sw_star = pinv(transpose(S)*S)*transpose(S);
% Ps = S*Sw_star;
% Ps_bot = eye(length(Ps))-Ps;
% Rs_tild = (transpose(d)*Ps_bot*Ps_bot*delta)/(2*transpose(d)*Ps_bot*Ps_bot*d);
Xs = -Sw_star*d;

[phi,theta,r] = cart2sph(Xs(1),Xs(2),Xs(3));
phi = phi*180/pi - 180/pi * atan(5*10^(-4)); % Winkelkorrektor (1 meter entfernung, da nicht UCA Array)
theta = theta*180/pi;
end

