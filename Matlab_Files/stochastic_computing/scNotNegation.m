function [not_result] = scNotNegation(a,N_bits)
%scMul takes to logical vectors and ANDs them
%   if bits > then a or b the vectors are repeated
%   a and b are not required to have the same length

rep_a = 1;


%expand a if required
if length(a)<N_bits
   rep_a = ceil(N_bits/length(a));
end



a_long = repmat(a,1,rep_a);

not_result = not(a_long(1:N_bits));
end

