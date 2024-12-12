function [mue,sigma] = EX_VarX_to_LogNormalDistribution(EX, VarX)
%% Determines the parameters of a LogNormal-distribution 
%% from empiric mean and variance 


mue = 2*log(EX) - 1/2*log(VarX+exp(2*log(EX)));
sigma = sqrt(-2*log(EX)+log(VarX+exp(2*log(EX))));


