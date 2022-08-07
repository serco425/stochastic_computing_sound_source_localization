function [multiAdd] = scMuxMultiAdd(mat,N_bits,selvec)
%scMultiAdd Stochastic Mulitplexer (MUX) for more than 2 signals
%   Unweighed add of rows in matrix
%   mat = 3x100 = mux with 3 inputs
%   N_bits if bits > then length of a or b the vectors are repeated

dims = size(mat,1);
rep = 1;
%expand a if required
if size(mat,2)<N_bits
   rep = ceil(N_bits/size(mat,2));
end

if ~exist('selvec','var')
    selvec = randi(dims,1,N_bits);
end

mat_long = repmat(mat,1,rep);
% 2^nextpow2(dims)

multiAdd = zeros(1,N_bits);

linear_index = sub2ind(size(mat_long),selvec,1:N_bits);
multiAdd(1,1:N_bits) = mat_long(linear_index);
end