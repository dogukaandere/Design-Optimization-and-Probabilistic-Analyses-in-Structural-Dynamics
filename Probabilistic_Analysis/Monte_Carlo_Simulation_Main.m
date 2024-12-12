%% Monte Carlo Methode Main

% Monte Carlo Simulation Analysis / Monte-Carlo-Simulationsanalyse

% Calculate the mean, standard deviation, and probability distribution of deflection at mid-span / Berechnung des Mittelwerts, der Standardabweichung und der Wahrscheinlichkeitsverteilung der Durchbiegung in der Balkenmitte

% Use Monte Carlo simulation to analyze the variability in the deflection response / Verwendung der Monte-Carlo-Simulation zur Analyse der Variabilität der Durchbiegungsantwort

% Assume elasticity modulus as a variable parameter with uncertainty / Elastizitätsmodul als unsichere Variable betrachten

% Run the simulation with the provided sample size to obtain statistical results / Durchführung der Simulation mit der angegebenen Stichprobengröße zur Erzielung statistischer Ergebnisse

%% Code Initialization / Code Initialisierung
clc;
close all;
clear all;

%% Reference Design with Parameters / Referenzentwurf mit Parametern

NE                 = 100;                                       % Number of elements / Elementanzahl
lengthh            = 1500;                                      % Beam length [mm] / Balkenlänge [mm]
lengthhElement     = lengthh / NE;                              % Element length [mm] / Elementenlängen [mm]
q0                 = 0.5;                                       % Acting distributed load [N/mm] / Wirkende Streckenlast [N/mm]
f                  = 0.5 * q0 * (lengthhElement + lengthhElement); % Load on an internal node of the beam due to distributed load q0 [N] / Belastung eines inneren Knotens des Balkens durch Streckenlast q0 [N]
faussen            = 0.5 * q0 * lengthhElement;                 % Load on an outer node of the beam due to distributed load q0 [N] / Belastung eines äußeren Knotens des Balkens durch Streckenlast q0 [N]
E_mue              = 205.6378e3;                                % Mean value [N/mm^2] of elasticity modulus determined from KS test (MaterialTestProben.m) / Mittelwert [N/mm^2] des E-Moduls, bestimmt aus KS-Test (MaterialTestProben.m)
E_var              = 81.4913e6;                                 % Variance of elasticity modulus determined from KS test (MaterialTestProben.m) / Varianz des E-Moduls, bestimmt aus KS-Test (MaterialTestProben.m)

%% Monte Carlo Method / Monte-Carlo-Methode
AnzRea = 10000;                                     % Number of realizations / Anzahl der Realisierungen                                                                                                   

%% Generating Realizations / Erzeugen der Realisierungen
% Elasticity modulus is considered as a random variable / Elastizitätsmodul wird als streuende Größe betrachtet
% Covariance matrix of elasticity modulus / Kovarianzmatrix des Elastizitätsmoduls

nx = NE;                                   % Number of elements / Elementanzahl
lc = 400;                                  % Correlation length [mm] / Korrelationslänge [mm]
x = [lengthhElement/2:lengthhElement:NE*lengthhElement];  % Step size: Distance between element midpoints / Schrittweite: Abstand der Mittelpunkte der Elemente
R = ones(nx, nx);                          % Correlation matrix / Korrelationsmatrix

% Construct the correlation matrix based on spatial correlation / Konstruktion der Korrelationsmatrix basierend auf räumlicher Korrelation
for i = 1:nx                         
    for j = i+1:nx
        R(i, j) = exp(-(x(i) - x(j))^2 / lc^2);
        R(j, i) = R(i, j);
    end
end

% Cov_E = Correlation matrix * variance / Cov_E = Korrelationsmatrix * Varianz
CovMa_E = R * E_var;
B = sqrtm(CovMa_E);                        % Matrix for generating correlated samples / Matrix zur Erzeugung korrelierter Stichproben

% Loop for generating realizations / Schleife zur Erzeugung der Realisierungen
for i = 1:AnzRea
    
    % Generate random values / Zufallszahlen erzeugen
    E_r = B * randn(NE, 1) + E_mue;

    % Evaluate the objective function / Zielfunktion auswerten
    x = [E_r];
    c(i) = Objective_Function(x);          % Store result of objective function / Ergebnis der Zielfunktion speichern
    
end

%% Mittelwert, Varianz und Standardabweichung der Zielfunktion
c=real(c);                          %Realteil der Durchbiegung in Balkenmitte
disp('Mean of deflection at mid-span / Mittelwert der Durchbiegung in Balkenmitte');
mue_c = mean(c)                    %Mittelwert der Zielfunktion

disp('Variance of deflection at mid-span / Varianz der Durchbiegung in Balkenmitte');
c_var = var(c)                     %Varianz der Zielfunktion
disp('Standard deviation of deflection at mid-span / Standardabweichung der Durchbiegung in Balkenmitte');
sig_c = std(c)                     %Standardabweichung der Zielfunktion

%% Visualization / Visualisierung

% Probability density / Wahrscheinlichkeitsdichte
fig1 = figure;
subplot(1, 2, 1);
% histogram(c, 250)                      % Optional histogram with 250 bins / Optionales Histogramm mit 250 Klassen
hist(c, 40);                             % Histogram with 40 bins / Histogramm mit 40 Klassen
title('Probability density of deflection at mid-span NE=100 AnzRea=10000 / Wahrscheinlichkeitsdichte der Durchbiegung in Balkenmitte NE=100 AnzRea=10000');
xlabel('Deflection at mid-span [mm] / Durchbiegung in Balkenmitte [mm]');
ylabel('PDF');

% Cumulative frequency / Summenhäufigkeit
subplot(1, 2, 2);
y_plot = sort(c);                        % Sorted values of deflection / Sortierte Werte der Durchbiegung
plot(sort(c), 0:1/(AnzRea-1):1);         % Plot of sorted values vs. cumulative frequency / Plot der sortierten Werte vs. kumulative Häufigkeit
stairs(sort(c), 0:1/(AnzRea-1):1);       % Stair-step plot for CDF / Treppen-Plot für CDF
title('CDF of deflection at mid-span NE=100 AnzRea=10000 / CDF der Durchbiegung in Balkenmitte NE=100 AnzRea=10000');
xlabel('Deflection at mid-span [mm] / Durchbiegung in Balkenmitte [mm]');
ylabel('CDF');

