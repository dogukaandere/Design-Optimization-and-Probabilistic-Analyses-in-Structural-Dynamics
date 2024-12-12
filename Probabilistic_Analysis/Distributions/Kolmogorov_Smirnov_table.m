function dmax = Kolmogorov_Smirnov_table(n,alpha)
%% Gibt die absolute größte erlaubte Differenz zurück
% 
% input:        n          Anzahl an Daten
%               alpha      Signifikanzniveau 
%                 
% output:       dmax       absolut größte erlaubte Differenz
%               
if nargin<2
    alpha = 0.05;
end

table=[0.9	0.684	0.565	0.494	0.446	0.41	0.381	0.358	0.339	0.322	0.307	0.295	0.284	0.274	0.266	0.258	0.25	0.244	0.237	0.231	0.21	0.19	0.18	1.07;
0.925	0.726	0.597	0.525	0.474	0.436	0.405	0.381	0.36	0.342	0.326	0.313	0.302	0.292	0.283	0.274	0.266	0.259	0.252	0.246	0.22	0.20	0.19	1.14;
0.95	0.776	0.642	0.564	0.51	0.47	0.438	0.411	0.388	0.368	0.352	0.338	0.325	0.314	0.304	0.295	0.286	0.278	0.272	0.264	0.24	0.22	0.21	1.22;
0.975	0.842	0.708	0.624	0.565	0.521	0.486	0.457	0.432	0.41	0.391	0.375	0.361	0.349	0.338	0.328	0.318	0.309	0.301	0.294	0.27	0.24	0.23	1.36;
0.995	0.929	0.828	0.733	0.669	0.618	0.577	0.543	0.514	0.49	0.468	0.45	0.433	0.418	0.404	0.392	0.381	0.371	0.363	0.356	0.32	0.29	0.27	1.63];

if alpha==0.2
    s=1;
elseif alpha==0.15
    s=2;
elseif alpha==0.1
    s=3;
elseif alpha==0.05
    s=4;
elseif alpha==0.01
    s=5;
else
    'not tabled'
end

if n<=20
    dmax=table(s,n);
elseif n<=25
    dmax=table(s,21);
elseif n<=30
    dmax=table(s,22);
elseif n<=35
    dmax=table(s,23);
else
    dmax=table(s,24)/sqrt(n);
end
