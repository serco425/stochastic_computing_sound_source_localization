function [mul] = scXNOrMul(a,b,N_bits)
%scMul takes to logical vectors and ANDs them
%   if bits > then a or b the vectors are repeated
%   a and b are not required to have the same length


rep_a = 1;
rep_b = 1;
%%expand a if required
if size(a,2)<N_bits
   rep_a = ceil(N_bits/size(a,2));
end

if size(b,2)<N_bits
   rep_b = ceil(N_bits/size(b,2));
end

a_long = repmat(a,1,rep_a);
b_long = repmat(b,1,rep_b);

mul = not(xor(a_long,b_long));
mul = mul(:,1:N_bits);
end

