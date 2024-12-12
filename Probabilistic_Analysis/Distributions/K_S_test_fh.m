function [accept,dmax] = K_S_test_fh(x,fh,sigLevel) 
%% Kolmogorow-Smirnow-Test
% 
% input:        x          Datensatz
%               fh         Verteilungsfunktion
%               sigLevel   Signifikanzniveau 
%                 
% output:       accept     Verteilungsfunktion
%               dmax       absolut größte Differenz




if nargin<3
    sigLevel=0.05;                 
end

%% emperical distribution function
n = length(x);                      
x = sort(x);                       
F_emp = zeros(n+1,1);               
for i=1:n                           
    F_emp(i+1) = F_emp(i)+1/n;     
end

%% assumed distribution
F_ass = zeros(n,1);
for i=1:n
    F_ass(i) = fh(x(i));            
end

%% search for max d and mean d 
dm=zeros(n,1);                      
for i=1:n
    d1 = abs(F_emp(i)-F_ass(i));
    d2 = abs(F_emp(i+1)-F_ass(i));
    dm(i) = max(d1,d2);
end
dmax = max(dm);

%% compare to dmax
dallow = Kolmogorov_Smirnov_table(n,sigLevel);
if dallow>dmax
    accept=1; %akzeptiert
else
    accept=0;
end