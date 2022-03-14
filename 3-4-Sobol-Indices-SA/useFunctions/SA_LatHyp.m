function [sample] = SA_LatHyp(k,bounds,N)
%Creates a random sample size [N x k] using latin hypercube minimizing the
%correlation of the parameters
%   k: number of parameters
%   bounds: [k,2] array containing lower and upper bounds in that order
%   [lower,upper]
%   N: number of parameter vectors (samples)
%   C: correlation, (default = no correlation)

latin = lhsdesign(N,k,'Criterion','correlation','Iterations',20);
sample = (bounds(:,1) + (bounds(:,2)-bounds(:,1)).*latin')';

% rng('shuffle')
% MC = bounds(:,1)+rand([k,N]).*bounds(:,2);
% sample=MC';
end

