function [Corr_result] = My_SequentialCorr(x,y,maxLag,LSB)
%% Corr_fi fixedpoint version of corr
%   Detailed explanation goes here

switch nargin
    case 3
        LSB = eps;
end


% this function uses row vectors
if (iscolumn(x))
   x = transpose(x);
end

if (iscolumn(y))
   y = transpose(y); 
end


% debug start
%q_debug = quantizer([5,4]);   % word length - fraction length in this case unsigned
%xx = num2bin(q_debug,double(bin2num(q_debug,x)))
%yy = num2bin(q_debug,double(bin2num(q_debug,y)))
%q_debug_result = quantizer([8,6]);
%% Cross Correlation
inLength = length(x);
store_second = zeros(1,maxLag+1);
store_prim = zeros(1,maxLag+1);
Corr_result=zeros(1,2*maxLag+1);

for ind = 1:inLength
   % update storag
   cur_second = y(ind);
   cur_prim = x(ind);
   
   store_second = [cur_second, store_second(1:end-1)];
   store_prim = [cur_prim, store_prim(1:end-1)];
   

   
   %parallel
   for par_ind = 1:maxLag
      Corr_result(par_ind) = quant( Corr_result(par_ind)+ double(store_second(par_ind)*store_prim(end)),LSB); 
   end
   Corr_result(maxLag+1) = quant ( Corr_result(maxLag+1)+double( store_second(end)*store_prim(end)),LSB);
   for par_ind = 1:maxLag
       Corr_result(maxLag+1+par_ind) = quant( Corr_result(maxLag+1+par_ind) + double(store_second(end)*store_prim(end-par_ind)),LSB);
   end
   
    % debug start
  % num2bin(q_debug_result,double(Corr_result))
  % second = num2bin(q_debug,double(bin2num(q_debug,store_second)))
   % debug end
   
%% meaning for parallel
%    xcorr_res(1) = xcorr_res(1)+store_second(1)*store_prim(end);
%    xcorr_res(2) = xcorr_res(2)+store_second(2)*store_prim(end);
%    xcorr_res(3) = xcorr_res(3)+store_second(3)*store_prim(end);
%    xcorr_res(4) = xcorr_res(4)+store_second(end)*store_prim(end);
%    xcorr_res(5) = xcorr_res(5)+store_second(end)*store_prim(3);
%    xcorr_res(6) = xcorr_res(6)+store_second(end)*store_prim(2);
%    xcorr_res(7) = xcorr_res(7)+store_second(end)*store_prim(1);
end





end

