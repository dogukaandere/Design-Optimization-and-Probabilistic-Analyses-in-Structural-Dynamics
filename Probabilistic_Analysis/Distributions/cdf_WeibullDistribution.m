function Fx = cdf_WeibullDistribution(x,a,b)
%% cumulative distribution function of the Weibull distribution

Fx = 1-exp(-a.*x^b)