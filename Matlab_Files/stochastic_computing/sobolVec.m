function [vec] = sobolVec(probability,sobol_rand)
%scVec Creates an vector of length n with ones uniformly distributed with proability
% Input probability [0,1]
%       length ... vector length

vec = sobol_rand < probability;
%vec = logical(round(sobol_rand+(probability-0.5)));

end
% Deterministic Methods for Stochastic Computing using Low-Discrepancy Sequences
% Note that, when converting to a bitstream representation,
% a one is generated if the Sobol number is less than the input target
% number.

% February 2021 changed to sobol_rand < probability