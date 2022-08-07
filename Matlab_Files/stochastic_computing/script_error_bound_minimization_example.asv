% Given n=16
% v_max = 5
% N = 6
% We have two input that is in interval [0,1]
% two are [0,0.5]
% eight that small [0,0.25]
n=16;
N=6;
func = @(c,v,L)((3*c*(c+1))/2+c*(v-1))/(n*(n-1))*L;

for v = 1:15
    L_high = 0;
    L_medium = 0;
    L_low = 0;
    
    c_high = max(n-v,0);
    c_medium = max(n/2-v,0);
    c_low = max(n/4-v,0);
    
    if c_high>0
        L_high = 1;
    end
    if c_medium>0
        L_medium=1;
    end
    
    if c_low>0
        L_low=4;
    end
    
    Err_high = func(c_high,v,L_high);
    
    
    
end