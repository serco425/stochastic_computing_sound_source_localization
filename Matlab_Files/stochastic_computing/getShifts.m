function [shift_vec] = getShifts(n,v)
%getShifts returns required relative shifts for 0 error add using or and pwm
%   n = first relative prime
%   v = allowed bits true - if 0 - it is set to 1 to avoid divicion by 0
%   intended to be with matlab function scDoRelativeShift
v = max(v,1); %avoid division by 0

    Large_Shifts = floor((n-2*v)/(v));
    Small_Shifts = floor((n-v)/v);
    N = floor((n^2-2*n*v)/(n*v)+1)*(floor((n-v)/v)+1);
    shift_vec = zeros(1,N);
    ind = 1;
    for large_shift_ind = 0:Large_Shifts
        cur_Shift = large_shift_ind*(n*v);
        shift_vec(ind) = cur_Shift;
        for ind_2 = 1:(Small_Shifts)
           ind = ind+1;
           cur_Shift = cur_Shift+v;
           shift_vec(ind) = cur_Shift;
        end
        ind = ind+1;
    end

end

