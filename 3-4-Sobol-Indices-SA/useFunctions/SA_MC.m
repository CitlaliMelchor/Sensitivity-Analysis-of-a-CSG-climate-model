function [sample] = SA_MC(k,bounds,N)
%Creates a random sample size [N x k] using monte carlo
%   k: number of parameters
%   bounds: [k,2] array containing lower and upper bounds in that order
%   [lower,upper]
%   N: number of parameter vectors (samples)

rng('shuffle')
MC = bounds(:,1)+rand([k,N]).*bounds(:,2);
sample=MC';
end

