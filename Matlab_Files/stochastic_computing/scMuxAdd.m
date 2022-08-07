function [vecAdd] = scMuxAdd(a,b,N_bits, sel)
%scAdd  takes to logical vectors and adds a and b using mux
%  Input:   sel ...if no select is chosen it will be a unweighed uniform distributed one with
%           N_bits if bits > then length of a or b the vectors are repeated
if ~exist('sel','var')
 % third parameter does not exist, so default it to something
 sel = scVec(0.5,N_bits); %probability = 0.5
end

rep_a = 1;
rep_b = 1;
rep_sel = 1;

%expand a if required
if length(a)<N_bits
   rep_a = ceil(N_bits/length(a));
end

if length(b)<N_bits
   rep_b = ceil(N_bits/length(b));
end

if length(sel)<N_bits
   rep_sel = ceil(N_bits/length(sel));
end

a_long = repmat(a,1,rep_a);
b_long = repmat(b,1,rep_b);
sel_long = repmat(sel,1,rep_sel);

vecAdd = boolean(zeros(1,N_bits)); %init
vecAdd(find(not(sel_long))) =  a_long(find(not(sel_long)));
vecAdd(find(sel_long)) =  b_long(find(sel_long));

end

