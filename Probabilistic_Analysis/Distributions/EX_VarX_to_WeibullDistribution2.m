function [a,b] = EX_VarX_to_WeibullDistribution2(EX, VarX, a0,b0)
%% Determines the parameters of a Weibull-distribution from empiric mean and
%% variance
if nargin<3
    a0=1; b0=1;
end
f = @(x) target(x,EX,VarX);
x = fmincon(f,[a0,b0],[],[],[],[],[0,0]);
a=x(1);
b=x(2);

function f = target(x,EX_target,VarX_target)
a=x(1);
b=x(2);
[EX, VarX] = EX_VarX_of_WeibullDistribution(a,b);
f = (EX-EX_target)^2 + (VarX-VarX_target)^2;