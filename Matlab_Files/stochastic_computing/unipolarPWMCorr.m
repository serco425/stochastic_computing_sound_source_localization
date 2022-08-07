function [sum_result] = unipolarPWMCorr(x,y,x_pwm_len,y_pwm_len,maxLag,optional)

if ~exist('optional','var') %
    optional.type = 'None';
end

mul_len = x_pwm_len*y_pwm_len;

x_inp = [x zeros(1,maxLag)];
y_inp = [y zeros(1,maxLag)];

sum_result = zeros(1,2*maxLag+1);
vecLength = length(x);

x_sc = pwmVec(x_inp,x_pwm_len,optional);
y_sc = pwmVec(y_inp,y_pwm_len,optional);


for curLag = -maxLag:+maxLag

    %mul = zeros(vecLength,mul_len);
    
    if curLag > 0
        shift_x = curLag;
        shift_y = 0;
    else   
        shift_x = 0;
        shift_y = -curLag;
    end
    mul = scAndMul(x_sc((1+shift_x):(shift_x+vecLength),:),   y_sc((1+shift_y):(shift_y+vecLength),:),mul_len);
    %for ind = 1:vecLength
    %   mul(ind,:) = scAndMul(x_sc((ind+shift_x),:),y_sc(ind+shift_y,:),mul_len);
    %end

    if strcmp(optional.type,'MultiOrApproxBinary')
        mul_binary = Unary2Binary(mul);
        sum_result(curLag+maxLag+1) = MultiOrApproxBinary(mul_binary);
    else
        sum_sc = scMultiOrAdd(mul,mul_len,optional);
        sum_result(curLag+maxLag+1) = Unary2Binary(sum_sc);        
    end

end

end

