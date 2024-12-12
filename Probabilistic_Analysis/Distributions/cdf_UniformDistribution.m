function Fx = cdf_UniformDistribution(x,a,b)
%% cumulative distribution function of the uniform distribution
% a (real)
% b (real)

if x<=a
    Fx = 0;
elseif x>=b
    Fx = 1;
else
    Fx = (x-a)/(b-a);
end
