function [vecMultiOrAdd] = scMultiOrAdd(mat,N_bits,optional)
%scAdd  takes to logical vectors and adds a and b using or gate
%  Input:   sel ...if no select is chosen it will be a unweighed uniform distributed one with
%           N_bits if bits > then length of size(mat,2) the values are repeated are repeated
%   mat = 3x100 = or gate with 3 inputs

if ~exist('optional','var') %
    delay = 0;
    optional.type = 'none';
else
    if strcmp(optional.type,'RotateSum')
        delay = round(optional.curPrime(1)*optional.curPrime(2)/size(mat,1)); %should depend somehow on 
    elseif strcmp(optional.type, 'RotateFor')                                                   % size of input and length
        delay = optional.rotate;
    elseif strcmp(optional.type, 'ShiftFor')
        delay = optional.rotate;
    elseif strcmp(optional.type, 'CyclicBuffer')
        delay = optional.rotate;
    else
        delay = 0;
        display('scMultiOrAdd_no_choice');
    end
end

dims = size(mat,1);
rep = 1;
%expand a if required
if size(mat,2)<N_bits
   rep = ceil(N_bits/size(mat,2));
end

mat_long = repmat(mat,1,rep);
vecMultiOrAdd = zeros(1,N_bits);

if strcmp(optional.type,'ShiftFor')
    for ind = linspace(dims,1,dims)
        temp_vec = [zeros(1,(ind-1)*delay), mat_long(ind,1:end-(ind-1)*delay)];
        vecMultiOrAdd = temp_vec | vecMultiOrAdd; 
    end
    
elseif strcmp(optional.type, 'CyclicBuffer') % Itroduces an additional zero field
    vecMultiOrAdd = [vecMultiOrAdd zeros(1,delay)];
    
    for ind = 1:dims
        rotate_times = dims-ind;
        temp_vec = [mat_long(ind,:) zeros(1,delay)];
        temp_vec = circshift(temp_vec,delay*rotate_times);
    
        vecMultiOrAdd = temp_vec | vecMultiOrAdd; 
    end
else
    for ind = 1:dims    % if delay is zero - it is just normal OR for all in matrix
       vecMultiOrAdd = [vecMultiOrAdd(1+(delay):end), vecMultiOrAdd(1:delay)];
       vecMultiOrAdd = mat_long(ind,:) | vecMultiOrAdd; 
    end
end


end

