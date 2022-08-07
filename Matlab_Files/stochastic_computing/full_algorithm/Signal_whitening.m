close all;

load('simulation_data.mat','X_Frame','stdvalues')

scale = 100; % magnify the signal. (Will be done Analog)
%first sequence
x=X_Frame(:,1)'.*scale;
%second sequence
y=X_Frame(:,2)'.*scale;

inLength = length(x);

maxLag = 4;
%% Whitening in Time Domain LMS Adaptive
rho = 128;
order = 10;

W = zeros(order,1);
L =  zeros(order,1); %spalten vector
R = zeros(order,1);
Lout = zeros(1,length(x));
Rout = zeros(1,length(y));

Lout(1) = x(1);
Rout(1) = y(1);

for ind = 2:length(x)                  
    L = [x(ind-1); L(1:end-1)];
    R = [y(ind-1); R(1:end-1)];
    
    Lout(ind) = x(ind)-L'*W;
    Rout(ind) = y(ind)-R'*W; %to keep phase relation
    
    W = W + rho*Lout(ind)*L;
end

x_post = Lout;
y_post = Rout;

%% Use Matlab corsscor function
xcorr_res = abs(xcorr(x_post,y_post,inLength/2));

%Interpolation
[yc, lag_time] = max(xcorr_res);
yl = xcorr_res(lag_time-1);
yr = xcorr_res(lag_time+1);

delta=(yl-yr)./(2.*(yl-2*yc+yr)); 
TapsLag_time = (length(xcorr_res)-1)/2-lag_time-delta+1;

constant_val = 1/stdvalues.fs*stdvalues.c/(2*stdvalues.radiusmean);
phi_time = asin(TapsLag_time*constant_val)*180/pi; % needs 90°



%% Reference of GCC-PHAT in Frequency Domain (optimal)
X_Frame = [x;y]';

%init
Xf_Frame = fft(X_Frame,inLength);

%Calculation
Cross_Spectrum = Xf_Frame(:,1).*conj(Xf_Frame(:,2));
%-------------------------------------
% Phat_Weighing in Frequency Domain
%-------------------------------------
Psi = abs(Cross_Spectrum+eps);
GCC_PHAT = ifftshift(ifft((Cross_Spectrum./(Psi))));

%Interpolation
[yc, lag] = max(GCC_PHAT);
yl = GCC_PHAT(lag-1);
yr = GCC_PHAT(lag+1);

delta=(yl-yr)./(2.*(yl-2*yc+yr)); 
TapsLag = (length(GCC_PHAT))/2-lag-delta+1;

constant_val = 1/stdvalues.fs*stdvalues.c/(2*stdvalues.radiusmean);
phi = asin(TapsLag*constant_val)*180/pi; % needs 90°

%% Plotting
figure();
sgtitle('GCC-PHAT time-frequency comparison');
subplot(511);
plot(GCC_PHAT(1+end/2-inLength/2:end/2+inLength/2));
legend('Frequency domain whitening')
subplot(512)
plot(xcorr_res);
legend('Time domain whitening filter')
subplot(513);
plot(inLength/2-maxLag:inLength/2+maxLag,GCC_PHAT(inLength/2-maxLag:inLength/2+maxLag));
legend('Reference freq-domain zoom')
subplot(514)
plot(inLength/2-maxLag:inLength/2+maxLag,xcorr_res(inLength/2-maxLag:inLength/2+maxLag));
legend('Time domain zoom')
subplot(515)
without_phat = abs(xcorr(x,y,inLength/2));
plot(inLength/2-maxLag:inLength/2+maxLag,without_phat(inLength/2-maxLag:inLength/2+maxLag));
legend('Without PHAT zoom')
