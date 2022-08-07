function [value] = Bipolar2Binary(mat)
% Takes an bipolar Bistream and returns the represented value
% takes a matrix or vector
% 3x1000 -> 3x1

if isrow(mat) || iscolumn(mat)
    value = 2*sum(mat)/length(mat)-1;
else
   value = 2*sum(mat,2)/size(mat,2)-1;
end

end

