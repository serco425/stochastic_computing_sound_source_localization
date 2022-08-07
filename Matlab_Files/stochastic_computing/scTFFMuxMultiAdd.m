function [multiAdd] = scTFFMuxMultiAdd(mat,N_bits, initialization)
%scMultiAdd Stochastic Mulitplexer (MUX) for more than 2 signals
%   Unweighed add of rows in matrix
%   mat = 3x100 = mux with 3 inputs
%   N_bits if bits > then length of a or b the vectors are repeated

if ~exist('initialization','var')
    initialization = 1;
end
    

dims = size(mat,1);
rep = 1;
%expand a if required
if size(mat,2)<N_bits
   rep = ceil(N_bits/size(mat,2));
end


    

mat_long = repmat(mat,1,rep);
multiAdd = zeros(1,N_bits);

% PairWiseReduction
stages = ceil(log2(dims));
%first_stage_select_signals = floor(dims/2);
intermediate_results = zeros(stages+1,dims,N_bits);
tff_signal = zeros(1,N_bits);

%init copy
intermediate_results(1,1:dims,:) = mat_long;
tff_signal(1) = initialization; %initialize first tff to 1

%loop computation
stage_dims = dims;
for ind_stage = 1:stages

    ind_par = 1;
    for ind_dims = 1:floor(stage_dims/2) %could also be used as index
        cur_stage_cur_mux = squeeze(intermediate_results(ind_stage,ind_par:ind_par+1,:)); % init mux with inputs

        xor_signal =  xor(cur_stage_cur_mux(1,:),cur_stage_cur_mux(2,:));
        tff_signal(2:end) = mod(cumsum(xor_signal(1:end-1))+initialization,2) ;
        
        cur_stage_cur_mux(2,:) = tff_signal;
        linear_index = sub2ind(size(cur_stage_cur_mux),xor_signal+1,1:N_bits);
        tmp_intermediate = cur_stage_cur_mux(linear_index);
        intermediate_results(ind_stage+1,ind_dims,:) = tmp_intermediate;
        ind_par = ind_par+2;
    end
        
    if mod(stage_dims,2)~=0 %if not power of two we forward one signal
        intermediate_results(ind_stage+1,ceil(stage_dims/2),:) = intermediate_results(ind_stage,stage_dims,:); % copy forward
    end

    stage_dims = ceil(stage_dims/2);

end
%sum(intermediate_results,3)

multiAdd(1,1:N_bits) = intermediate_results(end,1,:);
end
