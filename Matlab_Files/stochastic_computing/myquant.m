function [retval] = myquant(X,LSB,method)

if exist('method','var')
    
    if (strcmp(method,'fix') || strcmp(method,'truncate')) %round towards zero
        retval = fix(X./LSB).*LSB;
    elseif strcmp(method,'ceil')
        retval = ceil(X ./ LSB) .*LSB; %towards positive inivity
    elseif strcmp(method,'floor')
        retval = floor(X ./ LSB) .*LSB; %towards negative inivity
    elseif strcmp(method, 'round')
        retval = round(X ./ LSB) .*LSB; %default is rounding
    else
        disp('MyQuant - No method selected - do Roudning');
        retval = round(X ./ LSB) .*LSB; %default is rounding
    end
else
    retval = round(X ./ LSB) .*LSB; %default is rounding
end

end
