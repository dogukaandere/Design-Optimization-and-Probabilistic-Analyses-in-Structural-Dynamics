function [c, dc] = Objective_Function(x)
q0 = 0.5;
NE = length(x) / 4; 
lengthh = 1500;
lengthhElement = lengthh / NE;
b  = x(1:4:4*NE); 
d  = x(2:4:4*NE); 
t1 = x(3:4:4*NE); 
t2 = x(4:4:4*NE);  
EModul = 210e3; 
Le = (lengthhElement) * ones(NE,1); 
E(1:NE) = EModul; 
I = (b.*d.^3 - (b-t1).*(d-2*t2).^3) / 12; % Moment of inertia / Flächenträgheitsmoment
f = 0.5 * q0 * (lengthhElement + lengthhElement);
[u,Ke,ue,K] = BeamFEM(Le, E, I, f); % Cantilever FEM solution / Kragarm FEM Lösung

Vektor_e = zeros(length(u) - 4, 1);
Abfrage = (-1)^NE;
if Abfrage == -1
    Vektor_e(NE - 2) = 1; % NE, n is odd / NE, n ist ungerade
else  
    Vektor_e(NE - 1) = 1; % NE, n is even / NE, n ist gerade
end
lambda = K \ -Vektor_e;
lambda = [0; 0; lambda]; 
lambda = [lambda; 0; 0];
for p = 1:NE
    lambdae(:, p) = lambda(2*p - 1 : 2*p + 2); % Element-wise lambda / Elementweises Lambda
end

%% Deflection at the beam center / Durchbiegung in der Balkenmitte
s = (-1)^NE;

if s == -1
    c = u(NE); % NE is odd / NE ungerade
else
    c = u(NE + 1); % NE is even / NE gerade
end

%% Derivative of deflection at beam center with respect to design parameters / Ableitung der Durchbiegung in der Balkenmitte nach den Entwurfsparametern

% Derivative of deflection with respect to b (adjoint method) / Ableitung nach b (Adjungierten-Methode)
for i = 1:NE
    dIidbi = (1/12) * (d(i)^3 - (d(i) - 2 * t2(i))^3);    % Derivative of moment of inertia / Ableitung des Flächenträgheitsmoments
    dKidIi = Ke(:, :, i) / I(i);                          % Derivative of element stiffness matrix / Ableitung der Elementsteifigkeitsmatrix
    dKidbi = dKidIi * dIidbi;                             % Derivative with respect to design parameter / Ableitung nach Entwurfsgröße
    dcdb(i) = lambdae(:, i)' * dKidbi * ue(:, i);         % No need for assembled stiffness matrix / Keine Notwendigkeit der assemblierten Steifigkeitsmatrix
end
dcdb = transpose(dcdb);

% Derivative of deflection with respect to d (adjoint method) / Ableitung nach d (Adjungierten-Methode)
for i = 1:NE
    dIiddi = (1/4) * (b(i) * d(i)^2 - (b(i) - t1(i)) * (d(i) - 2 * t2(i))^2); 
    dKidIi = Ke(:, :, i) / I(i);               
    dKiddi = dKidIi * dIiddi;                                             
    dcdd(i) = lambdae(:, i)' * dKiddi * ue(:, i);
end
dcdd = transpose(dcdd);

% Derivative of deflection with respect to t1 (adjoint method) / Ableitung nach t1 (Adjungierten-Methode)
for i = 1:NE
    dIidt1i = (1/12) * ((d(i) - 2 * t2(i))^3); 
    dKidIi = Ke(:, :, i) / I(i);               
    dKidt1i = dKidIi * dIidt1i;                                                 
    dcdt1(i) = lambdae(:, i)' * dKidt1i * ue(:, i);
end
dcdt1 = transpose(dcdt1);

% Derivative of deflection with respect to t2 (adjoint method) / Ableitung nach t2 (Adjungierten-Methode)
for i = 1:NE
    dIidt2i = (1/2) * ((b(i) - t1(i)) * (d(i) - 2 * t2(i))^2); 
    dKidIi = Ke(:, :, i) / I(i);               
    dKidt2i = dKidIi * dIidt2i;                                                                                   
    dcdt2(i) = lambdae(:, i)' * dKidt2i * ue(:, i);
end
dcdt2 = transpose(dcdt2);

%% Derivative vector for deflection at beam center with respect to design parameters / Ableitungsvektor der Durchbiegung in der Balkenmitte
dc(1:4:4*NE,1) = dcdb;                
dc(2:4:4*NE,1) = dcdd;
dc(3:4:4*NE,1) = dcdt1;
dc(4:4:4*NE,1) = dcdt2;

%% Finite Difference Derivatives (Central Differences) of max deflection with respect to b / Zentrale Differenzen der Ableitungen nach b
for i = 1:NE
    b0 = b;
    b = b0 + 0.0001 * ones(length(b0), 1); % Forward step / Vorwärtsschritt
    I(i) = (b(i) * d(i)^3 - (b(i) - t1(i)) * (d(i) - 2 * t2(i))^3) / 12;
    u_p = BeamFEM(Le, E, I, f);

    if s == -1
        c_p = u_p(NE); % NE is odd / NE ist ungerade
    else
        c_p = u_p(NE + 1); % NE is even / NE ist gerade
    end

    b = b0 - 0.0001 * ones(length(b0), 1); % Backward step / Rückwärtsschritt
    I(i) = (b(i) * d(i)^3 - (b(i) - t1(i)) * (d(i) - 2 * t2(i))^3) / 12;
    u_m = BeamFEM(Le, E, I, f);

    if s == -1
        c_m = u_m(NE); % NE is odd / NE ist ungerade
    else
        c_m = u_m(NE + 1); % NE is even / NE ist gerade
    end

    dcdb_FDM(i) = (c_p - c_m) / (2 * 0.0001); % Central differences / Zentrale Differenzen
end
dcdb_FDM;

% Repeat similarly for dcdd_FDM, dcdt1_FDM, and dcdt2_FDM following the same steps / Wiederholen für dcdd_FDM, dcdt1_FDM, und dcdt2_FDM
