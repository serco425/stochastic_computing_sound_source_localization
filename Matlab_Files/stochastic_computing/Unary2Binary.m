function [value] = Unary2Binary(mat)
% Takes an unipolar Bistream and returns the represented value
% takes a matrix or vector
% 3x1000 -> 3x1
if isrow(mat) || iscolumn(mat)
    value = sum(mat)/length(mat);
else
   value = sum(mat,2)/size(mat,2);
end
end

