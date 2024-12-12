function [EX, VarX] = EX_VarX_of_WeibullDistribution(a,b)
%% Determines the mean and variance of a Weibull distribution 

g1 = gamma(1/b+1);
g2 = gamma(2/b+1);

EX = a^(-1/b) *g1;

VarX = a^(-2/b) *(g2-g1^2);
