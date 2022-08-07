% 28.07.2020
% This script shows the obvious fact that the maximum induced delay in
% worst case (v=1) is equal to nk-1
% the maximum induced delay decreases with increasing value of v


clear all;
curmax = 0;
config = [0,0, 0];

nmax = 256;
kmax = nmax-1;

for n = 3:nmax
    for v=15:128
        
        Delays = getShifts(n,v);
        
        tmpmax = max(Delays)/(n*(n-1));
        
        if curmax < tmpmax
            curmax = tmpmax;
            config = [n,v, max(Delays), length(Delays)];
        end
    end
end

