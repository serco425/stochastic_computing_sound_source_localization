function [return_mat] = filterMat(mat,num)
%filterMat Removes all nums from mat, when rows contain different
%numbers of num, then all rows are shrink to one length so it fits into the
%matrix
%   sample input: 1250x5
%   sample output: 890x5 with no num in it
cell_temp = cell(size(mat,2),1);
for i = 1:size(mat,2)
	cell_temp{i} = mat(mat(:,i)~=num,i);
end

min_vec = min(cellfun(@length,cell_temp));
return_mat = zeros(min_vec,size(mat,2));
for i=1:size(mat,2)
   temp = cell2mat(cell_temp(i));
   return_mat(:,i) = temp(1:min_vec);
end
end

