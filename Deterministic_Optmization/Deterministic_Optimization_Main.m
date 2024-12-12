%% Deterministic Optimization / Deterministische Optimierung

% Adjusting the code to specific boundary conditions and loads / Anpassung des Codes an spezifische Randbedingungen und Lasten
% Defining a constant cross-section (e.g., I-profile) for the beam / Festlegung eines konstanten Querschnitts (z. B. I-Profil) für den Balken
% Performing FEM simulation for the reference deflection at the beam center / Durchführung der FEM-Simulation für die Referenzdurchbiegung in der Balkenmitte
% Defining the objective function to minimize deflection / Definition der Zielfunktion zur Minimierung der Durchbiegung
% Setting constraints (e.g., max volume, allowable stresses) / Festlegung der Nebenbedingungen (z. B. max. Volumen, zulässige Spannungen)
% Variable profile parameters; profile type remains constant / Variabel gestaltete Profilparameter; Profiltyp bleibt konstant
% Calculating the derivative of the deflection at the beam center with respect to profile parameters / Berechnung der Ableitung der Durchbiegung in der Balkenmitte nach Profilparametern
% Validating derivatives using finite differences / Validierung der Ableitungen mittels finiter Differenzen
% Performing optimization with appropriate algorithms (e.g., gradient methods) / Durchführung der Optimierung mit geeigneten Algorithmen (z. B. Gradientenverfahren)
% Plausibility check of optimized results / Plausibilitätsprüfung der optimierten Ergebnisse
% Comparing the optimized deflection with the reference solution / Vergleich der optimierten Durchbiegung mit der Referenzlösung
%% Code Initialization / Code-Initialisierung
clc;
close all;
clear all;
%% Parameters / Parameter
NE                 = 20;                                        % Number of elements / Anzahl der Elemente
lengthh            = 1500;                                      % Beam length [mm] / Balkenlänge [mm]
EModul             = 210e3;                                     % Young's Modulus [N/mm^2] / E-Modul [N/mm^2]
lengthhElement     = lengthh/NE;                                % Element length [mm] / Elementlängen [mm]
q0                 = 0.5;                                       % Applied distributed load [N/mm] / Aufgebrachte Streckenlast [N/mm]
f                  = 0.5*q0*(lengthhElement+lengthhElement);    % Load on an inner node of the beam due to distributed load q0 [N] / Belastung eines inneren Knotens des Balkens durch Streckenlast q0 [N]
faussen            = 0.5*q0*lengthhElement;                     % Load on an outer node of the beam due to distributed load q0 [N] / Belastung eines äußeren Knotens des Balkens durch Streckenlast q0 [N]
%% Reference Design I-profile with parameters b, d, t1, t2 / Referenzentwurf I-Profil mit den Parametern b, d, t1, t2

b_ref    = 50;  % [mm] b0, d0, t10, t20 constant across all elements / b0, d0, t10, t20 konstant in allen Elementen
d_ref    = 50;  % [mm] 
t1_ref   = 1;   % [mm] 
t2_ref   = 3.571428571; % [mm]

V_zul = 6e5;    % [mm^3] allowable maximum volume of the beam / zulässiges maximales Volumen des Balkens
I_ref=(b_ref*d_ref^3-(b_ref-t1_ref)*(d_ref-2*t2_ref)^3)/12; % Moment of inertia with reference parameters [mm^4] / Flächenträgheitsmoment mit Referenzparametern [mm^4] 
c_ref=(q0*lengthh^4)/(384*EModul*I_ref); % Maximum deflection at beam center with reference parameters [mm] / Maximale Verschiebung in der Balkenmitte mit Referenzparametern [mm]
%% Initial parameters for b, d, t1, t2 / Startparameter für b, d, t1, t2

b0      = b_ref;    % [mm]
d0      = d_ref;    % [mm]
t10     = t1_ref;   % [mm] 
t20     = t2_ref;   % [mm] 
%% Initial vector x0 / Startvektor x0
x0(1:4:4*NE,1)      = b0;    
x0(2:4:4*NE,1)      = d0;
x0(3:4:4*NE,1)      = t10;
x0(4:4:4*NE,1)      = t20;
Le(1:NE)            = lengthh/NE;   % Element lengths / Elementlängen
E(1:NE)             = EModul;

%% fmincon-Solver

c           = @(x) Objective_Function(x);           % Define objective function / Zielfunktion definieren
nonlcon     = @(x) Constraints(x);       % Define constraints / Nebenbedingungen definieren

A_element   = [-1 0 1 0; 0 -1 0 2];           % Inequality constraints for each element / Ungleichheitsbedingungen für jedes Element
Cellarray   = repmat({A_element}, 1, NE);     % Cell array for each element / Zellen-Array für jedes Element
A           = blkdiag(Cellarray{:});          % Block diagonal matrix for all elements / Blockdiagonale Matrix für alle Elemente
b           = zeros(2*NE,1);                  % Inequality vector / Ungleichheitsvektor

Aeq = [];                                     % No equality constraints / Keine Gleichheitsbedingungen
beq = [];

ub           = repmat([50; 50; 25; 25],NE,1); % Upper bounds for design variables / Obere Grenze für die Entwurfsvariablen
lb           = 1*ones(4*NE,1);                % Lower bounds for design variables / Untere Grenze für die Entwurfsvariablen

options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',100000000,'SpecifyObjectiveGradient',true);

% Perform optimization / Optimierung durchführen
[x_opt] = fmincon(c, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

x_opt_A1 = x_opt;

b  = x_opt(1:4:4*NE);   
d  = x_opt(2:4:4*NE);
t1 = x_opt(3:4:4*NE);
t2 = x_opt(4:4:4*NE);

%% Deflection at the beam center with optimized design parameters / Durchbiegung in der Balkenmitte mit optimierten Entwurfsparametern

I = (b.*d.^3 - (b-t1).*(d-2*t2).^3)/12; % Moment of inertia / Flächenträgheitsmoment
[u,Ke,ue,K] = BeamFEM(Le, E, I, f);

s = (-1)^NE;
if s == -1
    c = u(NE); % NE is odd / NE ist ungerade
else
    c = u(NE+1); % NE is even / NE ist gerade
end
c_reduced_by = 1 - c/c_ref;             % Reduction in deflection [-] / Reduktion der Durchbiegung [-]
       

%% Beam volume with optimized design parameters / Balkenvolumen mit optimierten Entwurfsparametern
Ve = Le.*(2.*b.*t2 + (d - 2.*t2).*t1); % Element volume = Le * Ae / Elementvolumen = Le * Ae
BV = sum(Ve);                          % Beam volume / Balkenvolumen
Le = (1500/NE) * ones(NE,1);           % Vector of element lengths / Vektor der Elementlängen

%% Plot: Comparison of design parameters with optimized design parameters / Plot: Vergleich der Entwurfsparameter mit optimierten Entwurfsparametern
Plot_for_Comparison_of_Parameters(b, d, t1, t2, NE);


