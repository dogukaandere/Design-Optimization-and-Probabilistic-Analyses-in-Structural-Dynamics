    function Fx = cdf_GumbelDistribution(x,a,b)
%% cumulative distribution function of the Gumbel distribution
% a (real)
% b (real) > 0 

Fx = exp(-exp(-(x-a)/b));
