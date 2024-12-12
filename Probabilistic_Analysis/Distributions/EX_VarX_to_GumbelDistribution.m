function [a,b] = EX_VarX_to_GumbelDistribution(EX, VarX) 
%% Determines the parameters of a Gumbel-distribution 
%% from empiric mean and variance 

% Euler-Mascheroni-constant
gamma = 0.5772156649;

b = sqrt(VarX*6)/pi;
a = EX - b * gamma;


