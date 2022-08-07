load('simulation_data.mat','X','stdvalues');
inLength = 256; % [128, 256, 512] precision vs complexity
scale =300; % magnify the signal. (Will be done Analog)
LSB = 2^-6;
%first microphone
x=myquant(X(1:inLength,1)'.*scale,LSB);
%second microphone
y=myquant(X(1:inLength,2)'.*scale,LSB);

maxLag = 3;

%positive half-wave only
x(x<0) = 0;
y(y<0) = 0;

corr = xcorr(x,y,maxLag);

figure();plot(corr);

store_second = zeros(1,maxLag+1);
store_prim = zeros(1,maxLag+1);

xcorr_result = zeros(1,2*maxLag+1);

for ind = 1:inLength
   % update storag
   cur_second = y(ind);
   cur_prim = x(ind);
   
   store_second = [cur_second, store_second(1:end-1)];
   store_prim = [cur_prim, store_prim(1:end-1)];
   
   %parallel
   for par_ind = 1:maxLag
      xcorr_result(par_ind) = xcorr_result(par_ind)+store_second(par_ind)*store_prim(end); 
   end
   xcorr_result(maxLag+1) = xcorr_result(maxLag+1)+store_second(end)*store_prim(end);
   for par_ind = 1:maxLag
       xcorr_result(maxLag+1+par_ind) = xcorr_result(maxLag+1+par_ind)+store_second(end)*store_prim(end-par_ind);
   end

% bits-to-store
% 2*maxLag*2^Bits als input für die multiplikation
% 2*maxLag+1*2^(2*Bits) als summierungszwischenergebnis
% für maxLag = 3 und Bits = 2^6

% 2*3*2^6 + 2*5*2^12 = 41344 bits
% für Binär ->
% 2*maxLag*Bits + 2*(maxLag+1)*2*Bits = 2*3*6+2*4*12 = 132 bits
% Verhältnis = 1:313
   
%% meaning for parallel
%    xcorr_res(1) = xcorr_res(1)+store_second(1)*store_prim(end);
%    xcorr_res(2) = xcorr_res(2)+store_second(2)*store_prim(end);
%    xcorr_res(3) = xcorr_res(3)+store_second(3)*store_prim(end);
%    xcorr_res(4) = xcorr_res(4)+store_second(end)*store_prim(end);
%    xcorr_res(5) = xcorr_res(5)+store_second(end)*store_prim(3);
%    xcorr_res(6) = xcorr_res(6)+store_second(end)*store_prim(2);
%    xcorr_res(7) = xcorr_res(7)+store_second(end)*store_prim(1);
end

figure();plot(xcorr_result);