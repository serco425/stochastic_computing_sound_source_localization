function [pwm] = pwmVec(duty,vec_length,optional)
%pwmVec Creates one Periode of PWM signal
%   Inputs: duty high time [0 1]
%   length: desired length of the pwm signal

%limit duty cycle to [0,1]

%if ~exist('optional','var') %
%    optional.pwm = 'None';
%end

duty = reshape(duty,length(duty),1); %make to column vector

if strcmp(optional.pwm,'LTC6992')
    
    duty = (duty-0.1).*1.25;
    duty(duty>0.9) = 1;
    duty(duty<0.1) = 0;

    %if duty>0.9
    %    duty = 1;
    %elseif duty<0.1
    %    duty = 0;
    %else
    %    %from 10 to 90 -> 80% now 100% of pwm
    %    duty = (duty-0.1)*1.25;
    %end
else

    duty(duty>1) = 1;
    duty(duty<0) = 0;

    %if duty > 1
    %    duty = 1;
    %elseif duty < 0
    %    duty =0;
    %end
end

%assure right length
length_one = round(duty.*vec_length);
%length_zero = round((1-duty)*vec_length);
%if length_one+length_zero > vec_length %special treatment
%    length_zero = length_zero-1;
%end
pwm = zeros(length(duty),vec_length);
for ind = 1:length(duty)
pwm(ind,:) = [ones(1,length_one(ind)), zeros(1,vec_length-length_one(ind))];
end

end

