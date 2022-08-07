function tf = mems_tf_var(alpha,fvec)
% the following lines are just for checking reasons
% Cback = 83.5e-15; Rvent = 91.6e9; Mvent = 6.66e3   % nominal case
% Cback = 97.5e-15; Rvent = 98.4e9; Mvent = 6.97e3   % lower corner case
% Cback = 69.4e-15; Rvent = 85.4e9; Mvent = 6.36e3   % higher corner case


% By choosing the factor alpha between -1 and 1 you get the variation of
% the MEMS highpass cutoff frequency between 16 - 26 Hz

if abs(alpha) > 1
    error('alpha has to be between -1 and 1')
end

dCback = 14e-15;
dRvent = 6e9;
dMvent = 300;


Cback = 83.5e-15 + alpha*dCback;
Rvent = 91.4e9 + alpha*dRvent;
Mvent = 6.66e3 + alpha*dMvent;


            Rport = 2.2*19.8604e6;
            Mport = 3.08124e3;
            Cport = 5.31549e-15; 
           
            Rcore = 0.7*83.47812e6 ;   % fitting value needed because Niccolo uses ROM model !!!
            Rcore = 1*83.47812e6 ;   % fitting value needed because Niccolo uses ROM model !!!
            Mcore = 1.50474e3; 
            Cdia = 0.57*1.50483e-15 ;   % fitting value needed because Niccolo uses ROM model !!!
            Cdia = 2.3*1.50483e-15;   % fitting value needed because Niccolo uses ROM model !!!
            
            
%             Rvent=77.9429e9;
%             Mvent=5.9936e3;
%             Cback=1*98.5915e-15;

% state space representation

A=[-Rport/Mport  -1/Mport    0          0             0                0 
    1/Cport        0         0        -1/Cport       -1/Cport          0
    0              0         0         1/Cback        1/Cback          0
    0             1/Mvent   -1/Mvent  -Rvent/Mvent    0                0
    0             1/Mcore   -1/Mcore    0            -Rcore/Mcore  -1/Mcore
    0              0         0          0             1/Cdia           0 ];

b=[1/Mport  0  0  0  0  0]';
ct=[ 0  0  0  0  0  1];

d=0;


[num,den] = ss2tf(A,b,ct,d);
   % SSt transformation -->calculation of TF--> frequency response
  
   msys = ss(A,b,ct,d);
   c_msys = canon(msys,'modal',Inf);
%    c_msys = canon(msys,'companion');
   
   [A,b,ct,d] = ssdata(c_msys);
   
   [num_mems,den_mems] = ss2tf(A,b,ct,d,1);
  
   tf = freqs(num_mems,den_mems,2*pi*fvec);
   tf = 1.0408*tf;  % compensation to 0dB@1kHz for alpha = 0
       
%    figure
% 
%    subplot(211)
%     plot(fvec,20*log10(abs(tf))); grid on;
%    subplot(212)
%     plot(fvec,angle(tf)); grid on;
%     
%     axis([ 0 10e3 -40 20]);
%     title(' AC response Thor')
%     xlabel('Frequency [Hz]'),
%     ylabel('magnitude response [dB]')
%     
