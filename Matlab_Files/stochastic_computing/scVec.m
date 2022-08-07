function [vec] = scVec(probability,length)
%scVec Creates an vector of length n with ones uniformly distributed with proability
% Input probability [0,1]
%       length ... vector length
vec = logical(round(rand(1,length)+(probability-0.5)));

end

