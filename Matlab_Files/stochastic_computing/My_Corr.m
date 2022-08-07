function [Corr_result] = My_Corr(x,y,maxLag,LSB)
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

%% Cross Correlation

x_fi = x;
y_fi = fliplr(y);
inLength = length(x_fi);

Corr_result=zeros(1,2*maxLag+1);
ind = 1;
for i=inLength-maxLag:inLength+maxLag   %the taps to be calculated
    Corr_result(ind)=0;
    for j=1:inLength   %length of the one signal
        if(((i-j+1)>0) && (i-j)<inLength)
            Corr_result(ind)= quant( Corr_result(ind)+double((x_fi(j)*y_fi(i-j+1))),LSB);
        end
    end
    

    ind = ind+1;
end

end

