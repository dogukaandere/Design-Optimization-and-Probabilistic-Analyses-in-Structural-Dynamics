%% RDO_FOSM_main
% Robustness_Optimization / Robusheitstoptimierung

%% Code Initialization / Code-Initialisierung
clc;
close all;
clear all;

%% Parameters / Parameter
NE               = 4;  
% NE               = 20;                                      % Number of elements / Elementanzahl
% NE                 = 50;                                    % Number of elements (alternative) / Elementanzahl (Alternative)
% NE                 = 100;                                   % Number of elements (alternative) / Elementanzahl (Alternative)
lengthh            = 1500;                                    % Beam length [mm] / Balkenlänge [mm]
EModul             = 210e3;                                   % Young's modulus [N/mm^2] / E-Modul [N/mm^2]
lengthhElement     = lengthh/NE;                              % Element length [mm] / Elementenlängen [mm]
q0                 = 0.5;                                     % Distributed load [N/mm] / Wirkende Streckenlast [N/mm]
f                  = 0.5*q0*(lengthhElement+lengthhElement);  % Load on an internal node of the beam from distributed load q0 [N] / Belastung eines inneren Knotens des Balkens durch Streckenlast q0 [N]
faussen            = 0.5*q0*lengthhElement;                   % Load on an external node of the beam from distributed load q0 [N] / Belastung eines äußeren Knotens des Balkens durch Streckenlast q0 [N]

%% Reference Design I-Profile with parameters b, d, t1, t2 / Referenzentwurf I-Profil mit den Parametern b, d, t1, t2

b_ref    = 50;             % [mm] b0, d0, t10, t20 constant in all elements / b0, d0, t10, t20 in allen Elementen konstant
d_ref    = 50;             % [mm] Height of the beam / Höhe des Balkens
t1_ref   = 1;              % [mm] Thickness of the flange / Dicke des Flansches
t2_ref   = 3.571428571;    % [mm] Thickness of the web / Dicke des Steges

V_zul = 6e5;    % [mm^3] Allowable maximum volume of the beam / Zulässiges maximales Volumen des Balkens
I_ref=(b_ref*d_ref^3-(b_ref-t1_ref)*(d_ref-2*t2_ref)^3)/12; % Area moment of inertia with reference parameters [mm^4] / Flächenträgheitsmoment mit Referenzparametern [mm^4]
c_ref=(q0*lengthh^4)/(384*EModul*I_ref); % Maximum deflection at beam center with reference parameters [mm] / Maximale Verschiebung in Balkenmitte mit Referenzparametern [mm]

% Initial vector x0 / Startvektor x0
x0(1:4:4*NE,1)      = b_ref;        % Assign initial value for width / Zuweisung des Anfangswerts für die Breite
x0(2:4:4*NE,1)      = d_ref;        % Assign initial value for height / Zuweisung des Anfangswerts für die Höhe
x0(3:4:4*NE,1)      = t1_ref;       % Assign initial value for flange thickness / Zuweisung des Anfangswerts für die Flanschdicke
x0(4:4:4*NE,1)      = t2_ref;       % Assign initial value for web thickness / Zuweisung des Anfangswerts für die Stegdicke

Le(1:NE)            = lengthh/NE;   % Element length / Elementlänge
E(1:NE)             = EModul;       % Modulus of elasticity / Elastizitätsmodul

b_ref    = 50;             % [mm] b0, d0, t10, t20 constant in all elements / b0, d0, t10, t20 in allen Elementen konstant
d_ref    = 50;             % [mm] Height of the beam / Höhe des Balkens
t1_ref   = 1;              % [mm] Thickness of the flange / Dicke des Flansches
t2_ref   = 3.571428571;    % [mm] Thickness of the web / Dicke des Steges

V_zul = 6e5;    % [mm^3] Allowable maximum volume of the beam / Zulässiges maximales Volumen des Balkens
I_ref=(b_ref*d_ref^3-(b_ref-t1_ref)*(d_ref-2*t2_ref)^3)/12; % Area moment of inertia with reference parameters [mm^4] / Flächenträgheitsmoment mit Referenzparametern [mm^4]
c_ref=(q0*lengthh^4)/(384*EModul*I_ref); % Maximum deflection at beam center with reference parameters [mm] / Maximale Verschiebung in Balkenmitte mit Referenzparametern [mm]

% Initial vector x0 / Startvektor x0
x0(1:4:4*NE,1)      = b_ref;        % Assign initial value for width / Zuweisung des Anfangswerts für die Breite
x0(2:4:4*NE,1)      = d_ref;        % Assign initial value for height / Zuweisung des Anfangswerts für die Höhe
x0(3:4:4*NE,1)      = t1_ref;       % Assign initial value for flange thickness / Zuweisung des Anfangswerts für die Flanschdicke
x0(4:4:4*NE,1)      = t2_ref;       % Assign initial value for web thickness / Zuweisung des Anfangswerts für die Stegdicke

Le(1:NE)            = lengthh/NE;   % Element length / Elementlänge
E(1:NE)             = EModul;       % Modulus of elasticity / Elastizitätsmodul

%% fmincon Solver / fmincon-Solver 

% Define objective and constraint functions / Definition der Zielfunktion und Nebenbedingungen
c_rdo       = @(x) Objective_Function(x);     % Handle for the objective function / Funktionshandle für die Zielfunktion
nonlcon     = @(x) Constraints(x);            % Handle for the nonlinear constraint function / Funktionshandle für nichtlineare Nebenbedingungen

% Define linear inequality constraints for one element / Definition der linearen Ungleichheitsbedingungen für ein Element
% Constraints: t1_i - b_i <= 0 (flange thickness should not exceed flange width) / t1_i - b_i <= 0 (Flanschdicke darf Flanschbreite nicht überschreiten)
%              2*t2_i - d_i <= 0 (web thickness should not exceed half the beam height) / 2*t2_i - d_i <= 0 (Stegdicke darf die halbe Balkenhöhe nicht überschreiten)
A_element    = [-1 0 1 0; 0 -1 0 2];          % Matrix for constraints / Matrix für Nebenbedingungen
Cellarray = repmat({A_element}, 1, NE);       % Replicate for all elements / Wiederholung für alle Elemente
A  = blkdiag(Cellarray{:});                   % Assemble block diagonal matrix / Blockdiagonalmatrix erstellen
b  = zeros(2*NE,1);                           % Right-hand side for inequality constraints / Rechte Seite der Ungleichungsbedingungen

% No equality constraints for this problem / Keine Gleichheitsbedingungen für dieses Problem
Aeq = [];                                    
beq = [];

% Define bounds for design variables / Definition der Grenzen für die Entwurfsgrößen
ub           = repmat([50; 50; 25; 25],NE,1); % Upper bounds: [b, d, t1, t2] for all elements / Obere Grenzen: [b, d, t1, t2] für alle Elemente
lb           = 1*ones(4*NE,1);                % Lower bounds: Minimum dimensions for safety / Untere Grenzen: Minimale Abmessungen für Sicherheit

% Optimization options / Optionen für den Optimierungsalgorithmus
options = optimoptions('fmincon', ...
    'Display', 'iter', ...                    % Display optimization progress / Anzeige des Optimierungsfortschritts
    'MaxFunctionEvaluations', 1e6, ...        % Maximum number of function evaluations / Maximale Anzahl von Funktionsaufrufen
    'FiniteDifferenceStepSize', 0.01);        % Step size for numerical gradient approximation / Schrittgröße für numerische Gradientenschätzung

% Run the optimization solver / Ausführung des Optimierungssolvers
[x_opt] = fmincon(c_rdo, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

% Extract optimized values for each design variable / Extraktion der optimierten Werte für jede Entwurfsgröße
b  = x_opt(1:4:4*NE);                         % Optimized flange widths / Optimierte Flanschbreiten
d  = x_opt(2:4:4*NE);                         % Optimized beam heights / Optimierte Balkenhöhen
t1 = x_opt(3:4:4*NE);                         % Optimized flange thicknesses / Optimierte Flanschdicken
t2 = x_opt(4:4:4*NE);                         % Optimized web thicknesses / Optimierte Stegdicken

% Reference values for comparison / Referenzwerte zum Vergleich
b_ref  = x0(1:4:4*NE);                        % Reference flange widths / Referenz-Flanschbreiten
d_ref  = x0(2:4:4*NE);                        % Reference beam heights / Referenz-Balkenhöhen
t1_ref = x0(3:4:4*NE);                        % Reference flange thicknesses / Referenz-Flanschdicken
t2_ref = x0(4:4:4*NE);                        % Reference web thicknesses / Referenz-Stegstärken

%% Comparison of compliance between x_opt_RDO_FOSM, x_ref, and x_opt_DO / Vergleich der Nachgiebigkeit zwischen x_opt_RDO_FOSM, x_ref und x_opt_DO

% x_opt_RDO_FOSM     Results from the optimization in Robustness Optimization FOSM / Ergebnisse der Optimierung in Robustness Optimization FOSM
% x_ref              Reference design parameters / Referenzentwurfsparameter
% x_opt_DO           Results from the optimization in Deterministic Optimization / Ergebnisse der Optimierung in Deterministic Optimization

% %% Include optimization vector from Deterministic Optimization / Optimierungsvektor aus Deterministischen Optimierung mitberücksichtigen

% Load results from Deterministic_Optimization_Main.m based on NE / Ergebnisse aus Deterministic_Optimization_Main.m basierend auf NE laden
if NE == 20
    load('x_opt_DO_NE_20.mat');            % Results for NE=20 with index "x_opt_DO" / Ergebnisse für NE=20 mit Index "x_opt_DO"
elseif NE == 50
    load('x_opt_DO_NE_50.mat');            % Results for NE=50 with index "x_opt_DO" / Ergebnisse für NE=50 mit Index "x_opt_DO"
elseif NE == 100
    load('x_opt_DO_NE_100.mat');           % Results for NE=100 with index "x_opt_DO" / Ergebnisse für NE=100 mit Index "x_opt_DO"
else
    error('Unsupported value of NE. Please use NE=20, 50, or 100.');
end 

% Extract optimized design variables / Extraktion der optimierten Entwurfsgrößen
b_opt_DO  = x_opt_DO(1:4:4*NE);          % Optimized flange widths / Optimierte Flanschbreiten
d_opt_DO  = x_opt_DO(2:4:4*NE);          % Optimized beam heights / Optimierte Balkenhöhen
t1_opt_DO = x_opt_DO(3:4:4*NE);          % Optimized flange thicknesses / Optimierte Flanschdicken
t2_opt_DO = x_opt_DO(4:4:4*NE);          % Optimized web thicknesses / Optimierte Stegdicken

% Calculate area moment of inertia / Flächenträgheitsmoment berechnen
I = (b .* d.^3 - (b - t1) .* (d - 2*t2).^3) / 12;                             % Area moment of inertia / Flächenträgheitsmoment
I_opt_DO = (b_opt_DO .* d_opt_DO.^3 - (b_opt_DO - t1_opt_DO) .* (d_opt_DO - 2*t2_opt_DO).^3) / 12; % Optimized area moment of inertia / Optimiertes Flächenträgheitsmoment

% Calculate displacements using FEM for current and Task 1 parameters / Durchbiegung mit FEM berechnen für aktuelle und Aufgabe 1 Parameter
[u] = BeamFEM(Le, E, I, f);                                      % Displacement for current parameters / Durchbiegung für aktuelle Parameter
[u_opt_DO] = BeamFEM((lengthh / NE) * ones(NE, 1), EModul * ones(NE, 1), I_opt_DO, f); % Displacement for Task 1 / Durchbiegung für Aufgabe 1

% Display deflection at beam center for Task 3 / Anzeige der Durchbiegung in der Balkenmitte für Aufgabe 3
disp('Deflection at the beam center: Task 3 / Durchbiegung in der Balkenmitte: Aufgabe 3');
s = (-1)^NE;                                                       % Determine parity of NE / Parität von NE bestimmen
if s == -1
    c_opt_RDO_FOSM = u(NE);                                              % NE is odd / NE ist ungerade
else
    c_opt_RDO_FOSM = u(NE+1);                                            % NE is even / NE ist gerade
end

% Display deflection at beam center for Task 1 / Anzeige der Durchbiegung in der Balkenmitte für Aufgabe 1
disp('Deflection at the beam center: Task 1 / Durchbiegung in der Balkenmitte: Aufgabe 1');
if s == -1
    c_opt_DO = u_opt_DO(NE);                                       % NE is odd / NE ist ungerade
else
    c_opt_DO = u_opt_DO(NE+1);                                     % NE is even / NE ist gerade
end  

% Assign optimized values for Robustness Optimization / Weise optimierte Werte für 
A3.b = b_opt_DO;             % Optimierte Flanschbreiten
A3.d = d_opt_DO;             % Optimierte Balkenhöhen
A3.t1 = t1_opt_DO;           % Optimierte Flanschdicken
A3.t2 = t2_opt_DO;           % Optimierte Stegdicken

% Pass updated A3 to RDO_FOSM / Übergebe die aktualisierten A3 Parameter an RDO_FOSM
disp('Updating parameters in RDO_FOSM...');
output_RDO = RDO_FOSM(A3);

% Display results from RDO_FOSM / Ergebnisse aus RDO_FOSM anzeigen
disp('Results from RDO_FOSM:');
disp(output_RDO);

% Display deflection at beam center for reference design / Anzeige der Durchbiegung in der Balkenmitte für den Referenzentwurf
disp('Deflection at the beam center: Reference design / Durchbiegung in der Balkenmitte: Referenzentwurf');
c_ref = output_RDO.deflection_reference; % Assuming deflection_reference is in the output / Angenommen, Durchbiegung ist Teil von output

%% Vergleich der Standardabweichungen (Varianzen) zwischen x_opt_Aufgabe3, x_ref und x_opt_Aufgabe1
    
[sigma2_c_optA3]  = FOSM(b,d,t1,t2,NE);                            % FOSM mit x_opt_Aufgabe3
[sigma2_c_ref]    = FOSM(b_ref,d_ref,t1_ref,t2_ref,NE);            % FOSM mit x_ref
[sigma2_c_optA1]  = FOSM(b_opt_A1,d_opt_A1,t1_opt_A1,t2_opt_A1,NE);% FOSM mit x_opt_Aufgabe1 

disp('Varianz: Aufgabe3')
sigma2_c_optA3
disp('Varianz: Aufgabe1')
sigma2_c_optA1           %berechnet mit NE=20/NE=50
disp('Varianz: Referenzentwurf')
sigma2_c_ref

disp('Standardabweichung: Aufgabe3')
std_opt_A3 = sqrt(sigma2_c_optA3)
disp('Standardabweichung: Aufgabe1')
std_opt_A1 = sqrt(sigma2_c_optA1)           %berechnet mit NE=20/NE=50
disp('Standardabweichung: Referenzentwurf')
std_ref = sqrt(sigma2_c_ref)

%% Balkenvolumen mit optimiertern Entwurfsparametern
Ve = Le*(2*b.*t2+(d-2*t2).*t1); %Elementvolumen = Le * Ae
BV = sum(Ve);                   %Balkenvolumen
Le = (1500/NE) * ones(NE,1);    %Vektor der Elementlaengen
%% Plot: Vergleich der Entwurfsparameter mit den Referenzgroessen


plot_mit_A1(b,d,t1,t2,b_opt_A1,d_opt_A1,t1_opt_A1,t2_opt_A1,NE);

plot_ohne_A1(b,d,t1,t2,NE);


%% Zielfunktion
function [c_rdo] = Zielfunktion(x)
q0 = 0.5;
NE = length(x)/4; 
lengthh = 1500;
lengthhElement = lengthh/NE;
b  = x(1:4:4*NE); 
d  = x(2:4:4*NE); 
t1 = x(3:4:4*NE); 
t2 = x(4:4:4*NE);  
EModul = 210e3; 
Le = (lengthhElement) * ones(NE,1); 
E(1:NE)= EModul; 
I = (b.*d.^3 - (b-t1).*(d-2*t2).^3)/12; 
f = 0.5*q0*(lengthhElement+lengthhElement);
[u,Ke,ue,K] = KragarmFEM(Le,E,I,f);
Vektor_e = zeros(length(u)-4,1);
    Abfrage= (-1)^NE;
        if Abfrage==-1;
        Vektor_e(NE-2)=1; %NE,n is odd
        else  
        Vektor_e(NE-1)=1; %NE,n is even
        end
lambda=K\-Vektor_e;
lambda=[0;0;lambda]; 
lambda=[lambda;0;0];
for p=1:NE
        lambdae(:,p) = lambda(2*p-1:2*p+2);
    end

[sigma2_c] = FOSM(b,d,t1,t2,NE);     %Bestimmung der Varianz mit FOSM


%% Durchbiegung in der Balkenmitte + Standardabweichung mit einem Wichtungsfaktor w
% Robustheitsoptimierung 
% Optimierung mit Gradienten- basierten Algorithmus 
% Ableitung mit Finiten Differenzen
% If NE is odd:
%    c_rdo=w1*u(NE)+w2*sqrt(sigma2_c);
% If Ne is even:
%    c_rdo=w1*u(NE+1)+w2*sqrt(sigma2_c):

% s = (-1)^NE;
% 
% if s== -1;
%     
%     c_rdo = 1*u(NE) + 10* sqrt(sigma2_c); %NE ungerade (odd) %w1=1 %w2=10
% else
%    
%     c_rdo = 1*u(NE+1) + 10* sqrt(sigma2_c); %NE gerade (even) %w1=1 %w2=10
% end

s = (-1)^NE;

if s== -1;
    
    c_rdo = 0.01*u(NE) + 100* sqrt(sigma2_c); %NE ungerade (odd) %w1=1 %w2=3
else
   
    c_rdo = 0.01*u(NE+1) + 100* sqrt(sigma2_c); %NE gerade (even) %w1=1 %w2=3
end

% s = (-1)^NE;
% 
% if s== -1;
%     
%     c_rdo = 1*u(NE) + 1* sqrt(sigma2_c); %NE ungerade (odd) %w1=1 %w2=1
% else
%    
%     c_rdo = 1*u(NE+1) + 1* sqrt(sigma2_c); %NE gerade (even) %w1=1 %w2=1
% end
end
%% Nebenbedingung
% Die Nebenbedingung ist das maximale Balkenvolumen.
% = Summe der Elementvolumen
function [h, heq] = Nebenbedingungen(x)
NE = length(x)/4;          
b  = x(1:4:4*NE);          
d  = x(2:4:4*NE);          
t1 = x(3:4:4*NE);          
t2 = x(4:4:4*NE);          
V_zul = 6e5; 

Le = (1500/NE) ;                %Element Laenge
Ve = Le*(2*b.*t2+(d-2*t2).*t1); %Elementvolumen = Le * Ae
BV = sum(Ve);                   %Balkenvolumen, Skalarer Wert
h = BV - V_zul;                 %Ungleichheitsbedingung: Maximales Volumen <= 6*10^5 mm^3
heq = [];                       %Keine Gleichheitsbedingungen
end

function [sigma2_c] = FOSM(b,d,t1,t2,NE) 
%% FOSM
lengthh            = 1500;                                      % Balkenlänge           [mm]
lengthhElement     = lengthh/NE;                                % Elementenlängen       [mm]
q0                 = 0.5;                                       % Wirkende Streckenlast [N/mm]        
f                  = 0.5*q0*(lengthhElement+lengthhElement);    % Belastung eines inneren Knotens des Balkens durch Streckenlast q0  [N]
faussen            = 0.5*q0*lengthhElement;                     % Belastung eines äußeren Knotens des Balkens durch Streckenlast q0  [N]
E_mue              = 205.6378e3;                                % Mittelwert [N/mm^2] des E-Moduls bestimmt aus KS-Test (MaterialTestProben.m)      



%% stochastische Eingangsgrößen

%Mittelwertvektor

E_mue= E_mue * ones(NE, 1);
mue_x = [E_mue];

%Varianzen

E_var    = 81.4913e6; %Varianz des E-Moduls bestimmt aus KS-Test (MaterialTestProben.m)

%Kovarianzmatrix
%Elastizitätsmodul wird als streuende Größe betrachtet. 

nx=NE;
lc = 400;                           %Korrelationslänge [mm]
x = [lengthhElement/2:lengthhElement:NE*lengthhElement];                %Schrittweite: Abstand der Mittelpunkte der Elemente
R = ones(nx,nx);                    %Korrelationsmatrix

R = ones(nx,nx);                    %Korrelationsmatrix
for i=1:nx                         
    for j=i+1:nx
        R(i,j) = exp( -(x(i)-x(j))^2 /lc^2);
        R(j,i) = R(i,j);
    end
end

% M = ones(nx,nx);                    %Korrelationsmatrix
% for i=1:nx                         
%     for j=i+1:nx
%         M(i,j) = exp( -(x(j)-x(i))^2 /lc^2);
%         M(j,i) = R(i,j);
%     end
% end

% CovMa_E = Korrelationsmatrix*varianz

CovMa_E = R*E_var; 
CovMa = CovMa_E;

%% Varianz und Standardabweichung der Zielfunktion
[~,dc] = Zielfunktion_FOSM(mue_x);
var_c = 0;
for i=1:length(mue_x)                                            
    for j=1:length(mue_x)
        var_c = var_c + dc(i)*dc(j) * CovMa(i,j);
    end
end
var_c;

sigma2_c = var_c;
standardabweichung=sqrt(sigma2_c);
end
function [c,dc] = Zielfunktion_FOSM(x)
q0 = 0.5;
lengthh = 1500;
NE=length(x);
lengthhElement = lengthh/NE;
b_ref    = 50;  % [mm] b0,d0,t10,t20 in allen Elementen konstant
d_ref    = 50;  % [mm] 
t1_ref   = 1;   % [mm] 
t2_ref   = 3.571428571; % [mm]
b=b_ref*ones(NE, 1);
d=d_ref*ones(NE, 1);
t1=t1_ref*ones(NE, 1);
t2=t2_ref*ones(NE, 1);
Le = (lengthhElement) * ones(NE,1); 
E =x(1:1:NE);
I = (b.*d.^3 - (b-t1).*(d-2*t2).^3)/12; 
f = 0.5*q0*(lengthhElement+lengthhElement);
[u,Ke,ue,K] = KragarmFEM(Le,E,I,f);
Vektor_e = zeros(length(u)-4,1);
    Abfrage= (-1)^NE;
        if Abfrage==-1;
        Vektor_e(NE-2)=1; %NE,n is odd
        else  
        Vektor_e(NE-1)=1; %NE,n is even
        end
lambda=K\-Vektor_e;
lambda=[0;0;lambda]; 
lambda=[lambda;0;0];
for p=1:NE
        lambdae(:,p) = lambda(2*p-1:2*p+2);
    end
%% Durchbiegung in der Balkenmitte
s = (-1)^NE;

if s== -1;
    
    c = u(NE); %NE ungerade (odd)
else
   
    c = u(NE+1); %NE gerade (even)
end
%% Ableitung der Durchbiegung in der Balkenmitte nach dem E-Modul


%Ableitung der Durchbiegung in der Balkenmitte nach E (Adjungierten-Methode)
for i=1:NE;

dKidEi = Ke(:,:,i)/E(i);                     %Nur die Elementsteifigkeitsmatrix durch E(i) geteilt                                                                            
dcdE(i)=lambdae(:,i)'*dKidEi*ue(:,i);
end
dcdE=transpose(dcdE);
dcdE;
%% Ableitungsvektor der Durchbiegung in der Balkenmitte nach dem E-Modul

dc=dcdE;
end


