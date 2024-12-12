%% Comparison of the explicit derivative with the FDM / Vergleich der expliziten Ableitung mit der FDM 

%% Code Initialization / Code-Initialisierung
clc; 
clear all; 
close all;
%% Parameters for derivative validation / Parameter für die Validierung der Ableitung
q0 = 0.5; %[N/mm]; 
NE = 200; 
s = (-1)^NE;
lengthh = 1500;
lengthhElement = lengthh/NE; 
EModul = 210e3;

%% Profile parameters (reference values b=50mm, d=50mm, t1=1mm, t2=3.57142857mm, constant in all elements)/ 
% / Profilparameter (Referenzwerte b=50mm, d=50mm, t1=1mm, t2=3.57142857mm, in allen Elementen konstant
b00=50;
d00=50;
t100=1;
t200=3.57142857;
b = ones(NE,1)*b00;
d = ones(NE,1)*d00;
t1 = ones(NE,1)*t100;
t2 = ones(NE,1)*t200;

%% FEM calculation / FEM Rechnung
Le = (lengthhElement) * ones(NE,1);     
E(1:NE)= EModul; % Element lengths / Elementlängen
I = (b.*d.^3 - (b-t1).*(d-2*t2).^3)/12; % Moment of inertia I-profile (constant in all elements) / Flächenträgheitsmoment I-Profil (konstant in allen Elementen)
f = 0.5*q0*(lengthhElement+lengthhElement);

[u,Ke,ue,K] = BeamFEM(Le, E, I, f);

Vektor_e = zeros(length(u)-4,1);
Abfrage = (-1)^NE;
if Abfrage == -1;
    Vektor_e(NE-2) = 1; % NE, n is odd / NE,n ist ungerade
else  
    Vektor_e(NE-1) = 1; % NE, n is even / NE,n ist gerade
end

lambda = K\-Vektor_e;
lambda = [0;0;lambda]; 
lambda = [lambda;0;0];

for p = 1:NE
    lambdae(:,p) = lambda(2*p-1:2*p+2);
end

%% Derivative of the maximum deflection at the beam center with respect to b (Adjoint Method) / Ableitung der maximalen Verschiebung in Balkenmitte nach b (Adjungierten-Methode)

for i = 1:NE
    dIidbi = (1/12) * (d(i)^3 - (d(i) - 2 * t2(i))^3); 
    dKidIi = Ke(:,:,i) / I(i);               
    dKidbi = dKidIi * dIidbi;               
    dcdb(i) = lambdae(:,i)' * dKidbi * ue(:,i);  
end
dcdb = transpose(dcdb);
dcdb;

%% Derivative of the maximum deflection at the beam center with respect to d (Adjoint Method) / Ableitung der maximalen Verschiebung in Balkenmitte nach d (Adjungierten-Methode)
for i = 1:NE
    dIiddi = (1/4) * (b(i) * d(i)^2 - (b(i) - t1(i)) * (d(i) - 2 * t2(i))^2);
    dKidIi = Ke(:,:,i) / I(i);               
    dKiddi = dKidIi * dIiddi;               
    dcdd(i) = lambdae(:,i)' * dKiddi * ue(:,i);
end
dcdd = transpose(dcdd);
dcdd;

%% Derivative of the maximum deflection at the beam center with respect to t1 (Adjoint Method) / Ableitung der maximalen Verschiebung in Balkenmitte nach t1 (Adjungierten-Methode)
for i = 1:NE
    dIidt1i = (1/12) * ((d(i) - 2 * t2(i))^3); 
    dKidIi = Ke(:,:,i) / I(i);               
    dKidt1i = dKidIi * dIidt1i;                                                   
    dcdt1(i) = lambdae(:,i)' * dKidt1i * ue(:,i);
end
dcdt1 = transpose(dcdt1);
dcdt1;

%% Derivative of the maximum deflection at the beam center with respect to t2 (Adjoint Method) / Ableitung der maximalen Verschiebung in Balkenmitte nach t2 (Adjungierten-Methode)
for i = 1:NE
    dIidt2i = (1/2) * ((b(i) - t1(i)) * (d(i) - 2 * t2(i))^2);
    dKidIi = Ke(:,:,i) / I(i);               
    dKidt2i = dKidIi * dIidt2i;               
    dcdt2(i) = lambdae(:,i)' * dKidt2i * ue(:,i);
end
dcdt2 = transpose(dcdt2);
dcdt2;

%% Finite difference derivative (Central Differences) of the maximum deflection at the beam center with respect to b
% Finite Differenzen Ableitung (Zentrale Differenzen) der maximalen Verschiebung in Balkenmitte nach b
b = ones(NE,1) * b00;
d = ones(NE,1) * d00;
t1 = ones(NE,1) * t100;
t2 = ones(NE,1) * t200;
for i = 1:NE
    b0 = b;
    % Forward step / Vorwärtsschritt
    b = b0 + 0.0001 * ones(length(b0), 1);
    I(i) = (b(i) * d(i)^3 - (b(i) - t1(i)) * (d(i) - 2 * t2(i))^3) / 12;
    u_p = BeamFEM(Le, E, I, f);

    if s == -1
        c_p = u_p(NE); % NE is odd / NE ist ungerade
    else
        c_p = u_p(NE + 1); % NE is even / NE ist gerade
    end

    % Backward step / Rückwärtsschritt
    b = b0 - 0.0001 * ones(length(b0), 1);
    I(i) = (b(i) * d(i)^3 - (b(i) - t1(i)) * (d(i) - 2 * t2(i))^3) / 12;
    u_m = BeamFEM(Le, E, I, f);

    if s == -1
        c_m = u_m(NE); % NE is odd / NE ist ungerade
    else
        c_m = u_m(NE + 1); % NE is even / NE ist gerade
    end

    dcdb_FDM(i) = (c_p - c_m) / (2 * 0.0001); % Central Differences / Zentrale Differenzen
end
dcdb_FDM;

%% Finite difference derivative (Central Differences) of the maximum deflection at the beam center with respect to d
% Finite Differenzen (Zentrale Differenzen) Ableitung der maximalen Verschiebung in Balkenmitte nach d
b = ones(NE, 1) * b00;
d = ones(NE, 1) * d00;
t1 = ones(NE, 1) * t100;
t2 = ones(NE, 1) * t200;
for i = 1:NE
    d0 = d;
    % Forward step / Vorwärtsschritt
    d = d0 + 0.0001 * ones(length(d0), 1);
    I(i) = (b(i) * d(i)^3 - (b(i) - t1(i)) * (d(i) - 2 * t2(i))^3) / 12;
    u_p = BeamFEM(Le, E, I, f);

    if s == -1
        c_p = u_p(NE); % NE is odd / NE ist ungerade
    else
        c_p = u_p(NE + 1); % NE is even / NE ist gerade
    end

    % Backward step / Rückwärtsschritt
    d = d0 - 0.0001 * ones(length(d0), 1);
    I(i) = (b(i) * d(i)^3 - (b(i) - t1(i)) * (d(i) - 2 * t2(i))^3) / 12;
    u_m = BeamFEM(Le, E, I, f);

    if s == -1
        c_m = u_m(NE); % NE is odd / NE ist ungerade
    else
        c_m = u_m(NE + 1); % NE is even / NE ist gerade
    end

    dcdd_FDM(i) = (c_p - c_m) / (2 * 0.0001); % Central Differences / Zentrale Differenzen
end
dcdd_FDM;

%% Finite difference derivative (Central Differences) of the maximum deflection at the beam center with respect to t1
% Finite Differenzen (Zentrale Differenzen) Ableitung der maximalen Verschiebung in Balkenmitte nach t1
b = ones(NE, 1) * b00;
d = ones(NE, 1) * d00;
t1 = ones(NE, 1) * t100;
t2 = ones(NE, 1) * t200;
for i = 1:NE
    t10 = t1;
    % Forward step / Vorwärtsschritt
    t1 = t10 + 0.0001 * ones(length(t10), 1);
    I(i) = (b(i) * d(i)^3 - (b(i) - t1(i)) * (d(i) - 2 * t2(i))^3) / 12;
    u_p = BeamFEM(Le, E, I, f);

    if s == -1
        c_p = u_p(NE); % NE is odd / NE ist ungerade
    else
        c_p = u_p(NE + 1); % NE is even / NE ist gerade
    end

    % Backward step / Rückwärtsschritt
    t1 = t10 - 0.0001 * ones(length(t10), 1);
    I(i) = (b(i) * d(i)^3 - (b(i) - t1(i)) * (d(i) - 2 * t2(i))^3) / 12;
    u_m = BeamFEM(Le, E, I, f);

    if s == -1
        c_m = u_m(NE); % NE is odd / NE ist ungerade
    else
        c_m = u_m(NE + 1); % NE is even / NE ist gerade
    end

    dcdt1_FDM(i) = (c_p - c_m) / (2 * 0.0001); % Central Differences / Zentrale Differenzen
end
dcdt1_FDM;

%% Finite difference derivative (Central Differences) of the maximum deflection at the beam center with respect to t2
% Finite Differenzen (Zentrale Differenzen) Ableitung der maximalen Verschiebung in Balkenmitte nach t2
b = ones(NE, 1) * b00;
d = ones(NE, 1) * d00;
t1 = ones(NE, 1) * t100;
t2 = ones(NE, 1) * t200;
for i = 1:NE
    t20 = t2; 
    % Forward step / Vorwärtsschritt
    t2 = t20 + 0.0001 * ones(length(t20), 1);
    I(i) = (b(i) * d(i)^3 - (b(i) - t1(i)) * (d(i) - 2 * t2(i))^3) / 12;
    u_p = BeamFEM(Le, E, I, f);

    if s == -1
        c_p = u_p(NE); % NE is odd / NE ist ungerade
    else
        c_p = u_p(NE + 1); % NE is even / NE ist gerade
    end

    % Backward step / Rückwärtsschritt
    t2 = t20 - 0.0001 * ones(length(t20), 1);
    I(i) = (b(i) * d(i)^3 - (b(i) - t1(i)) * (d(i) - 2 * t2(i))^3) / 12;
    u_m = BeamFEM(Le, E, I, f);

    if s == -1
        c_m = u_m(NE); % NE is odd / NE ist ungerade
    else
        c_m = u_m(NE + 1); % NE is even / NE ist gerade
    end

    dcdt2_FDM(i) = (c_p - c_m) / (2 * 0.0001); % Central Differences / Zentrale Differenzen
end
dcdt2_FDM;

%% Plot: Comparison of explicit and FDM derivatives / Vergleich der Ableitungen explizit und FDM
fig1 = figure; % opens a new window / öffnet ein Fenster
movegui(fig1, [300 -50]); % move window to specified screen location / bewegt Fenster an bestimmte Bildschirmposition
sgtitle('Comparison explicit vs. FDM with N=100 / Vergleich explizit mit FDM mit N=100') % Title appears at the top / Titel erscheint oben

% Comparison of dcdb and dcdb_FDM
subplot(2,2,1) % Additional screen appears, top left / Weiterer Screen erscheint, links oben 
x_plot = linspace(1, NE, NE); % generates n points with spacing (x2 - x1) / generiert n Punkte mit Abstand (x2 - x1)
dcdb_FDM_plot = full(dcdb_FDM); % Convert sparse matrix to full storage / Konvertiert Sparse-Matrix zu voller Speicherung
y_plot = [dcdb, dcdb_FDM_plot']; % transpose as full generates row vector / transponiert da full zeilenvektor erzeugt
h = plot(x_plot, y_plot);
title('Derivative with respect to b / Ableitung nach b')
h(1).LineWidth = 2;
h(1).Color = 'red';
h(2).Color = 'blue';
h(2).LineWidth = 1;
legend({'dcdb', 'dcdb_{FDM}'}, 'Location', 'southeast');
xlabel('Number of Elements [-] / Elementanzahl [-]'), ylabel('dcdb');

% Comparison of dcdd and dcdd_FDM
subplot(2,2,2)
x_plot = linspace(1, NE, NE);
dcdd_FDM_plot = full(dcdd_FDM);
y_plot = [dcdd, dcdd_FDM_plot'];
h = plot(x_plot, y_plot);
title('Derivative with respect to d / Ableitung nach d')
h(1).LineWidth = 2;
h(1).Color = 'red';
h(2).Color = 'blue';
h(2).LineWidth = 1;
legend({'dcdd', 'dcdd_{FDM}'}, 'Location', 'southeast');
xlabel('Number of Elements [-] / Elementanzahl [-]'), ylabel('dcdd');

% Comparison of dcdt1 and dcdt1_FDM
subplot(2,2,3)
x_plot = linspace(1, NE, NE);
dcdt1_FDM_plot = full(dcdt1_FDM);
y_plot = [dcdt1, dcdt1_FDM_plot'];
h = plot(x_plot, y_plot);
title('Derivative with respect to t1 / Ableitung nach t1')
h(1).LineWidth = 2;
h(1).Color = 'red';
h(2).Color = 'blue';
h(2).LineWidth = 1;
legend({'dcdt1', 'dcdt1_{FDM}'}, 'Location', 'southeast');
xlabel('Number of Elements [-] / Elementanzahl [-]'), ylabel('dcdt1');

% Comparison of dcdt2 and dcdt2_FDM
subplot(2,2,4)
x_plot = linspace(1, NE, NE);
dcdt2_FDM_plot = full(dcdt2_FDM);
y_plot = [dcdt2, dcdt2_FDM_plot'];
h = plot(x_plot, y_plot);
title('Derivative with respect to t2 / Ableitung nach t2')
h(1).LineWidth = 2;
h(1).Color = 'red';
h(2).Color = 'blue';
h(2).LineWidth = 1;
legend({'dcdt2', 'dcdt2_{FDM}'}, 'Location', 'southeast');
xlabel('Number of Elements [-] / Elementanzahl [-]'), ylabel('dcdt2');

% Relative error between the explicit derivative and the Finite Difference Method / Relative Abweichung zwischen der expliziten Ableitung und der Finite-Differenzen-Methode
e1 = abs((norm(dcdb) - norm(dcdb_FDM)) / norm(dcdb));    % Relative error: dcdb vs dcdb_FDM / Relative Abweichung: dcdb zu dcdb_FDM
e2 = abs((norm(dcdd) - norm(dcdd_FDM)) / norm(dcdd));    % Relative error: dcdd vs dcdd_FDM / Relative Abweichung: dcdd zu dcdd_FDM
e3 = abs((norm(dcdt1) - norm(dcdt1_FDM)) / norm(dcdt1)); % Relative error: dcdt1 vs dcdt1_FDM / Relative Abweichung: dcdt1 zu dcdt1_FDM
e4 = abs((norm(dcdt2) - norm(dcdt2_FDM)) / norm(dcdt2)); % Relative error: dcdt2 vs dcdt2_FDM / Relative Abweichung: dcdt2 zu dcdt2_FDM
error = [e1; e2; e3; e4];                                % Vector with relative errors / Vektor mit den relativen Abweichungen

% Plot: Relative error of the derivative with respect to b (between Adjoint Method and FDM) as a function of the element number / Plot: relative Abweichung der Ableitung nach b (zwischen Adjungierte-Methode und FDM) in Abhängigkeit der Elementanzahl
e1_5 = 4.6725e-06;       % Error between dcdb and dcdb_FDM with 5 elements / Abweichung zwischen dcdb und dcdb_FDM mit 5 Elementen
e1_10 = 1.3616e-05;      % Error between dcdb and dcdb_FDM with 10 elements / Abweichung zwischen dcdb und dcdb_FDM mit 10 Elementen
e1_50 = 8.1499e-05;      % Error between dcdb and dcdb_FDM with 50 elements / Abweichung zwischen dcdb und dcdb_FDM mit 50 Elementen
e1_100 = 2.4151e-04;     % Error between dcdb and dcdb_FDM with 100 elements / Abweichung zwischen dcdb und dcdb_FDM mit 100 Elementen
e1_200 = 8.6503e-04;     % Error between dcdb and dcdb_FDM with 200 elements / Abweichung zwischen dcdb und dcdb_FDM mit 200 Elementen
e1_300 = 0.0100;         % Error between dcdb and dcdb_FDM with 300 elements / Abweichung zwischen dcdb und dcdb_FDM mit 300 Elementen
e1_350 = 0.0220;         % Error between dcdb and dcdb_FDM with 350 elements / Abweichung zwischen dcdb und dcdb_FDM mit 350 Elementen
e1_400 = 0.1322;         % Error between dcdb and dcdb_FDM with 400 elements / Abweichung zwischen dcdb und dcdb_FDM mit 400 Elementen
e1_450 = 0.2134;         % Error between dcdb and dcdb_FDM with 450 elements / Abweichung zwischen dcdb und dcdb_FDM mit 450 Elementen
e1_500 = 0.5614;         % Error between dcdb and dcdb_FDM with 500 elements / Abweichung zwischen dcdb und dcdb_FDM mit 500 Elementen

error_dcdb = 100 * [e1_5; e1_10; e1_50; e1_100; e1_200; e1_300; e1_350; e1_400; e1_450; e1_500];
x = [5; 10; 50; 100; 200; 300; 350; 400; 450; 500];

fig2 = figure;
movegui(fig2, [800 -50]);
plot(x, error_dcdb)
xlabel('Number of Elements [-] / Elementanzahl der Rechnung [-]'), ylabel('Relative error between dcdb and dcdb(FDM) [%] / relative Abweichung zwischen dcdb und dcdb(FDM) [%]')
title('Relative error of derivative with respect to b / relative Abweichung der Ableitung nach b')
