function [vecOrAdd] = scOrAdd(a,b,N_bits)
%scAdd  takes to logical vectors and adds a and b using mux
%  Input:   sel ...if no select is chosen it will be a unweighed uniform distributed one with
%           N_bits if bits > then length of a or b the vectors are repeated

rep_a = 1;
rep_b = 1;

if (iscolumn(a)) % convert to rowvector if required
   a = transpose(a);
end

if (iscolumn(b))
   b = transpose(b); 
end

%expand a if required
if length(a)<N_bits
   rep_a = ceil(N_bits/length(a));
end

if length(b)<N_bits
   rep_b = ceil(N_bits/length(b));
end

a_long = repmat(a,1,rep_a);
b_long = repmat(b,1,rep_b);

vecOrAdd = or(a_long,b_long);
end

