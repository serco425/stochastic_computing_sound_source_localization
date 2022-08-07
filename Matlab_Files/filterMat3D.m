function [return_mat] = filterMat3D(mat,num)
%filterMat Removes all nums from mat, when rows contain different
%numbers of num, then all rows are shrink to one length so it fits into the
%matrix
%   uses filterMat(mat,num)
%   sample input: 1250x5x10
%   sample output: 890x5x10 with no num in it
cell_temp = cell(size(mat,3),1);
for i = 1:size(mat,3)
	cell_temp{i} = filterMat(mat(:,:,i),num);
end

min_mat = min(cellfun(@length,cell_temp));
return_mat = zeros(min_mat,size(mat,2),size(mat,3));
for i=1:size(mat,3)
   temp = cell2mat(cell_temp(i));
   return_mat(:,:,i) = temp(1:min_mat,:);
end

end


