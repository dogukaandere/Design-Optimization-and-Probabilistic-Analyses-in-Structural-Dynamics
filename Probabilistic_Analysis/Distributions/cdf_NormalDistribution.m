function Fx = cdf_NormalDistribution(x,mue,var)
%% cumulative distribution function of the normal distribution
%   mue: mean value
%   var: variance

if nargin<2
    mue=0; var=1;
end

Fx = normcdf(x,mue,sqrt(var));

% numerisch
% % Fx = numInt_GaussQuad('pdf_NormalDistribution',[mue var],min(x-sqrt(var)*5,mue-sqrt(var)*5),x,20);
% fh = @(x) pdf_NormalDistribution(x,mue,var);
% Fx = numInt_GaussQuad_fh(fh,min(x-sqrt(var)*5,mue-sqrt(var)*5),x,20);
