function [ phi, theta ] = srpPhat( X_Frame, stdvalues )
%% Steered Response Power method using Phat Weighing
% Bandpass filter to filter for all frequencies f>300 and f<6000
% Calculate Delay for each frequency for given search grid and sum over
% them
%% Inputs: 
%%% 1) X_frame is a frame of Data (Window_Length x Microphones)
%%% 2) stdvalues contains information about Array Geometry, search resolution,
%%% sampling rate
%% Outputs:
%%% 1) phi: azimuth
%%% 2) theta: elevation
%% 

%-------------------------------------
% Cauer Bandpass - Passband 300hz - 6000khz
%-------------------------------------
Bandpass = [1.0	1.5091602180533874	1.0	1.0	1.309386604007795	0.9718231494427941
            1.0	-1.9886909525924237	1.0	1.0	-1.9788783710616373	0.9951413232134249
            1.0	1.8538619980330788	1.0	1.0	0.37500030973601706	0.7527775296868646
            1.0	-1.996927779472577	1.0	1.0	-1.8865714937664464	0.9367540846824399
            1.0	1.436761666160152	1.0	1.0	1.4073257957189724	0.996851276478052
            1.0	-1.98675592847237	1.0	1.0	-1.985486639606707	0.9994752214260119];


X = sosfilt(Bandpass,X_Frame); %second-order-section digital filter
Xf_Frame = transpose(fft(X,stdvalues.nfft));



Power_Mat_3D = zeros(length(stdvalues.PhiAngles),length(stdvalues.ThetaAngles),length(linspace(stdvalues.klow,stdvalues.khigh)));
for k = stdvalues.klow:stdvalues.khigh
    
  Xfi = Xf_Frame(:,k);
  S= Xfi*ctranspose(Xfi); 
  
  for thetaind = 1:length(stdvalues.ThetaAngles)
      for phiind = 1:length(stdvalues.PhiAngles)

        Dm = stdvalues.Delay_Matrix(:,:,phiind,thetaind,k);
        Power_Mat_3D(phiind,thetaind,k) = ctranspose(Dm)*S*Dm;

       end
  end
 %Phat Weighting   
  Power_Mat_Weighing(:,:,k) = 1./(abs(ctranspose(Xfi)*Xfi)+eps).*Power_Mat_3D(:,:,k);
end

Pssl_Mat = sum(Power_Mat_Weighing,3);

%-------------------------------------
% Search for Maximum
%-------------------------------------
[phiind, thetaind] = findMaximum (Pssl_Mat);
% change to find Maximum, because other approach gives wrong results with
% negative values

phi = stdvalues.PhiAngles(phiind)*180/pi;
theta = stdvalues.ThetaAngles(thetaind)*180/pi;

end