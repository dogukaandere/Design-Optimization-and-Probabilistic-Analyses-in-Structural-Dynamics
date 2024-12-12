function Fx = cdf_LogNormalDistribution(x,mue,sigma)
%% cumulative distribution function of the lognormal distribution
%   mue (real)
%   var (real)

if nargin<2
    mue=0; sigma=1;
end

Fx = logncdf(x,mue,sigma);
