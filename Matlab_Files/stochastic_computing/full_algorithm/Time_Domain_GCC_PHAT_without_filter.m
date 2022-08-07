% Drüben einbauen und performance für verschiedene Bitwidths analysieren


clear all;
load('simulation_data.mat','X','stdvalues');
inLength = 256; % [128, 256, 512] precision vs complexity
scale =300; % magnify the signal. (Will be done Analog)
LSB = 2^-6;
%first microphone
x=myquant(X(1:inLength,1)'.*scale,LSB);
%second microphone
y=myquant(X(1:inLength,2)'.*scale,LSB);

% corr = myquant(xcorr(x,y,5),LSB);
% figure();plot(corr);

maxLag = 5;

%positive half-wave only
x(x<0) = 0;
y(y<0) = 0;

corr = xcorr(x,y,maxLag);
% figure();plot(corr);
max(corr)

% max(corr)
%% Stochastic performance

Prime_select = [[251,241];[127,98];[61,53];[31,23]];
if LSB == 2^-8
    curPrime = Prime_select(1,:)
elseif LSB == 2^-7
    curPrime = Prime_select(2,:)
elseif LSB == 2^-6
    curPrime = Prime_select(3,:)
else
    curPrime = Prime_select(4,:)
end
optional.curPrime = curPrime;

%without options
stochastic = unipolarPWMCorr(x,y,curPrime(1),curPrime(2),maxLag); 

%without options
optional.type = 'MultiOrApproxBinary';
stochastic_uncorrelated = unaryPWMCorr(x,y,curPrime(1),curPrime(2),maxLag,optional);

optional.type = 'RotateSum';
stochastic_rotate_sum = unaryPWMCorr(x,y,curPrime(1),curPrime(2),maxLag,optional);


figure();
subplot(141);
plot(corr)
legend('XCORR Matlab');
subplot(142);
plot(stochastic_uncorrelated);
legend('And with Binary-Add-Approx');
subplot(143);
plot(stochastic);
legend('Or and And');
subplot(144);
plot(stochastic_rotate_sum);
legend('IDB Rotate Sum');



%% Testing Math operation


Iinputs = 3;
num1 = 0.1;
num2 = 0.2;
num3 = 0.3;

length = 2^(inputs*-log2(LSB));
%select a random part of sobol sequence
p = sobolset(3);
sobol_rand = net(p,length);

vec1 = sobolVec(num1,sobol_rand(:,1));
vec2 = sobolVec(num2,sobol_rand(:,2));
vec3 = sobolVec(num3,sobol_rand(:,3));

res = Unary2Binary(vec1 | vec2 | vec3)

% binary = myquant(num1*num2*num3,LSB)


%% Or With three Inputs
%num1+num2+num3-num1*num2-num1*num3-num2*num3+num1*num2*num3

%% Multiply with predifined Window
Window = hann(inLength,'periodic'); % Window constant
x = x.*Window';
y = y.*Window';


x = 1:10;
y = 1:10;
inLength = length(x);

%% Cross Correlation
x_post = x;
y_post = fliplr(y);
maxlag =100;

A=zeros(1,2*maxlag+1);
ind = 1;
for i=inLength-maxlag:inLength+maxlag   %the taps to be calculated
    A(ind)=0;
    for j=1:inLength   %length of the one signal
        if(((i-j+1)>0) && (i-j)<inLength)
            A(ind)=A(ind)+(x_post(j)*y_post(i-j+1));
        end
    end
    ind = ind+1;
end

figure();plot(A-xcorr(x,y,100));

%% parabolic Interpolation
[yc, lag] = max(A);
yl = A(lag-1);
yr = A(lag+1);

delta=(yl-yr)./(2.*(yl-2*yc+yr)); 
TapsLag = maxlag+1-lag-delta;

%% Time Delay to Angle trough trigonometric function

constant_val = 1/stdvalues.fs*stdvalues.c/(2*stdvalues.radiusmean);
phi = asin(TapsLag*constant_val)*180/pi

%% Plot Section

% figure();
% sgtitle('Framing a Signal');
% subplot(311);
% plot(Window);
% legend('Hann Window');
% subplot(312);
% plot(X(inLength+1:2*inLength,1)*scale);
% legend('Input Stream');
% subplot(313);
% plot((X(inLength+1:2*inLength,1))*scale.*Window);
% legend('Framed Input Stream');

%% Binary Implementation of AND+OR Convolution
% 
% x_inp = myquant([x zeros(1,maxLag)],LSB);
% y_inp = myquant([y zeros(1,maxLag)],LSB);
% 
% sum_result = zeros(1,2*maxLag+1);
% vecLength = length(x);
% 
% for curLag = -maxLag:+maxLag
%        
%     if curLag > 0
%         shift_x = curLag;
%         shift_y = 0;
%     else   
%         shift_x = 0;
%         shift_y = -curLag;
%     end
%     
%     x_temp = x_inp(1+ shift_x:vecLength+shift_x);
%     y_temp = y_inp(1+ shift_y:vecLength+shift_y);
%     
%     mul = x_temp.*y_temp;
%     
%     sum_result(curLag+maxLag+1) = MultiOrApproxBinary(mul);
%     
% end
% figure();plot(sum_result)
% title('Binary Approximation')
