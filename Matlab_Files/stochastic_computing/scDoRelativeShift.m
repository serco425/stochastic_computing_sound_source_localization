function [mat_out] = scDoRelativeShift(mat,shift_vec)
%scDoRelativeShift takes one matirx and a vector, then performs relative
%   shift of rows in mat according to values in shift_vec
%   mat_out:    right shifted version of mat
%               size is increaced to fit the shifted values
%               size(mat_out,2) = size(mat,2)+max(shift_vec)              

mat_out = [mat, zeros(size(mat,1),max(shift_vec))];

for ind = 1:length(shift_vec)
    mat_out(ind,:) = circshift(mat_out(ind,:),shift_vec(ind));
end

end
