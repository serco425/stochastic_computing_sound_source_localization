close all;
clear all;

load('simulation_data.mat','X','stdvalues');


inLength = 512; % process 512 samples
scale = 100; % magnify the signal. (Will be done Analog)
%first sequence
x=X(1:inLength,1)'.*scale;
%second sequence
y=X(1:inLength,2)'.*scale;


%% whitening
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
% x_post = Lout; %input (with or without whitening)
% y_post = fliplr(Rout); % elements flipped left to right


%% Cross Correlation
x_post = x; %input (with or without whitening)
y_post = fliplr(y);% elements flipped left to right

% maxlag = inLength/2;
maxlag =4;

A=zeros(1,2*maxlag+1);
ind = 1;
for i=inLength-maxlag:inLength+maxlag   %the taps to be calculated
    A(ind)=0;
    for j=1:inLength   %length of the one signal
        if(((i-j+1)>0) && (i-j+1)<inLength)
            A(ind)=A(ind)+(x_post(j)*y_post(i-j+1));
        else
        end
    end
    ind = ind+1;
end

%% parabolic Interpolation
[yc, lag] = max(A);
yl = A(lag-1);
yr = A(lag+1);

delta=(yl-yr)./(2.*(yl-2*yc+yr)); 
TapsLag = maxlag+1-lag-delta;

constant_val = 1/stdvalues.fs*stdvalues.c/(2*stdvalues.radiusmean);
phi = asin(TapsLag*constant_val);


%% Reference of GCC-PHAT in Frequency Domain (optimal)
X_Frame = [x;y]';
stdvalues.indpairs = 1;
stdvalues.m = 2;
stdvalues.nfft = inLength;
optional.phat = true;

Cross_Spectrum = zeros(stdvalues.nfft,stdvalues.indpairs);
GCC_PHAT = zeros(stdvalues.nfft,stdvalues.indpairs);
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
       GCC_PHAT(:,ind) = ifftshift(ifft((Cross_Spectrum(:,ind)./(Psi))));
       ind = ind +1;
       
    end
end
%% Time Domain Reference Function
xcorr_res = xcorr(x,y,maxlag);

figure();
subplot(411);
plot(A);
subplot(412);
plot(xcorr_res);
subplot(413);
plot(xcorr_res-A);
subplot(414);
plot(GCC_PHAT(inLength/2-maxlag:inLength/2+maxlag));
