function [multiAdd] = scSobolMultiAdd(mat,N_bits,selvec)
%scMultiAdd Stochastic Mux for more than 2 signals
%   Unweighed add of rows in matrix
%   mat = 3x100 = mux with 3 inputs
%   N_bits if bits > then length of a or b the vectors are repeated

if ~exist('selvec','var') %create selecten vector with sobol sequence
    % third parameter does not exist, so default it to something
    dims = size(mat,1);   
    Bitwidth = nextpow2(dims);
    
    p = sobolset(Bitwidth,'Skip',1e3);
    sobolseq = (net(p,N_bits))';
    scale = 2.^(0:Bitwidth-1);
    sobolseq_round = round(sobolseq);

    scaled_selvec = sobolseq_round.*transpose(scale);
    selvec = sum(scaled_selvec);
    % selvec = randi(dims,1,N_bits); Wintout sobol sequence
end

rep = 1;

%expand a if required
if size(mat,2)<N_bits
   rep = ceil(N_bits/size(mat,2));
end
mat_long = repmat(mat,1,rep);

multiAdd = zeros(1,N_bits);

linear_index = sub2ind(size(mat_long),1+selvec,1:N_bits);
multiAdd(1,1:N_bits) = mat_long(linear_index);
end