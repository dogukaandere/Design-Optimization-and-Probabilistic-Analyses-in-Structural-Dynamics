%% First Order Second Moment [FOSM] Main

% First Order Second Moment (FOSM) Analysis / Erstes und zweites Moment Analyse (FOSM)

% Calculate the mean, standard deviation, and probability distribution of deflection at mid-span / Berechnung des Mittelwerts, der Standardabweichung und der Wahrscheinlichkeitsverteilung der Durchbiegung in der Balkenmitte

% Use the FOSM method to approximate these statistics / Verwendung der FOSM-Methode zur Approximation dieser Statistiken

% Assume elasticity modulus as a variable parameter with uncertainty / Elastizitätsmodul als unsichere Variable betrachten

% Determine partial derivatives of deflection with respect to input variables for FOSM / Bestimmung der partiellen Ableitungen der Durchbiegung bezüglich der Eingangsvariablen für FOSM

% Validate the derivatives using finite differences / Validierung der Ableitungen mit finiten Differenzen

%% Code Initialization / Code-Initialisierung
clc;
close all;
clear all;
%% Reference Design with Parameters / Referenzentwurf mit den Parametern

NE                 = 100;                                       % Number of elements / Elementanzahl
lengthh            = 1500;                                      % Beam length [mm] / Balkenlänge [mm]
lengthhElement     = lengthh / NE;                              % Element length [mm] / Elementenlängen [mm]
q0                 = 0.5;                                       % Acting distributed load [N/mm] / Wirkende Streckenlast [N/mm]
f                  = 0.5 * q0 * (lengthhElement + lengthhElement); % Load on an internal node of the beam due to distributed load q0 [N] / Belastung eines inneren Knotens des Balkens durch Streckenlast q0 [N]
faussen            = 0.5 * q0 * lengthhElement;                 % Load on an outer node of the beam due to distributed load q0 [N] / Belastung eines äußeren Knotens des Balkens durch Streckenlast q0 [N]
E_mue              = 205.6378e3;                                % Mean value [N/mm^2] of elasticity modulus determined from KS test (Material_Test_Proben.m) / Mittelwert [N/mm^2] des E-Moduls, bestimmt aus KS-Test (Material_Test_Proben.m)


%% FOSM (First Order Second Moment Method) / FOSM 

% Stochastic input variables / Stochastische Eingangsgrößen

% Mean vector / Mittelwertvektor
E_mue = E_mue * ones(NE, 1);          % Mean value vector of elasticity modulus / Mittelwertvektor des Elastizitätsmoduls
mue_x = [E_mue];                      % Vector of mean values / Vektor der Mittelwerte

% Variances / Varianzen
E_var = 81.4913e6;                    % Variance of elasticity modulus determined from KS test (MaterialTestProben.m) / Varianz des E-Moduls, bestimmt aus KS-Test (MaterialTestProben.m)

% Covariance matrix / Kovarianzmatrix
% Elasticity modulus is considered as a random variable / Elastizitätsmodul wird als streuende Größe betrachtet

nx = NE;                              % Number of elements / Elementanzahl
lc = 400;                             % Correlation length [mm] / Korrelationslänge [mm]
x = [lengthhElement/2:lengthhElement:NE*lengthhElement]; % Step size: Distance between element midpoints / Schrittweite: Abstand der Mittelpunkte der Elemente
R = ones(nx, nx);                     % Correlation matrix / Korrelationsmatrix

% Build the correlation matrix based on spatial correlation / Aufbau der Korrelationsmatrix basierend auf räumlicher Korrelation
for i = 1:nx                         
    for j = i+1:nx
        R(i, j) = exp(- (x(i) - x(j))^2 / lc^2);
        R(j, i) = R(i, j);
    end
end

% Covariance matrix calculation / Berechnung der Kovarianzmatrix
CovMa_E = R * E_var;                  % Covariance matrix of elasticity modulus / Kovarianzmatrix des Elastizitätsmoduls
CovMa = CovMa_E;                      % Final covariance matrix / Endgültige Kovarianzmatrix


%% Mean Value of the Objective Function / Mittelwert der Zielfunktion

% disp('Mean value of deflection at mid-span of the beam') / Mittelwert der Durchbiegung in Balkenmitte
mue_c = Objective_Function(mue_x);             % Mean value of the objective function / Mittelwert der Zielfunktion
disp('Mean value of deflection at mid-span of the beam') % Display message / Anzeige der Nachricht
mue_c = full(mue_c);                           % Convert to full format (if sparse) / Konvertieren in volles Format (falls spärlich)

%% Variance and Standard Deviation of the Objective Function / Varianz und Standardabweichung der Zielfunktion

[~, dc] = Objective_Function(mue_x);         % Gradient of the objective function / Gradient der Zielfunktion
var_c = 0;                                   % Initialize variance / Initialisierung der Varianz
for i = 1:length(mue_x)                                           
    for j = 1:length(mue_x)
        var_c = var_c + dc(i) * dc(j) * CovMa(i, j);  % Calculate variance based on covariance matrix / Berechnung der Varianz basierend auf der Kovarianzmatrix
    end
end
var_c;                                       % Final variance value / Endgültiger Varianzwert

disp('Variance of deflection at mid-span of the beam / Varianz der Durchbiegung in Balkenmitte') % Display message / Anzeige der Nachricht
var_c

disp('Standard deviation of deflection at mid-span of the beam / Standardabweichung der Durchbiegung in Balkenmitte') % Display message / Anzeige der Nachricht
sig_c = sqrt(var_c)                          % Standard deviation as square root of variance / Standardabweichung als Wurzel der Varianz


