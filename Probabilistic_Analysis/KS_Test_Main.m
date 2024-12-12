%% Kolmogrov-Smirnov Test 

% Kolmogorov-Smirnov Test for Elastic Modulus Distribution / Kolmogorov-Smirnov-Test für die Verteilung des Elastizitätsmoduls

% Test the distribution of the elasticity modulus using the Kolmogorov-Smirnov test / Verteilung des Elastizitätsmoduls mithilfe des Kolmogorov-Smirnov-Tests überprüfen

% Evaluate the match with different distributions: Uniform, Normal, Gumbel, Log-normal, and Weibull / Übereinstimmung mit verschiedenen Verteilungen bewerten: Gleich-, Normal-, Gumbel-, Log-Normal- und Weibull-Verteilung

% Select the best-fit distribution based on test results / Auswahl der am besten passenden Verteilung basierend auf den Testergebnissen
%% Code Initialization / Code-Initialisierung
clc; 
clear all; 
close all;                

%% Distributions / Verteilungen
addpath('Distributions')                  % Add path to distributions folder / Pfad zum Distributions-Ordner hinzufügen
[E, ~] = Material_Test_Sample;            % Load material test sample data / Materialtest-Proben-Daten laden
x = E;                                    % Assign sample data to variable x / Probendaten der Variable x zuweisen

EX = mean(x);                             % Mean / Mittelwert
VarX = var(x);                            % Variance / Varianz
stdvX = std(x);                           % Standard deviation / Standardabweichung

% Gumbel Distribution / Gumbel-Verteilung
[a, b] = EX_VarX_to_GumbelDistribution(EX, VarX);   % Parameters of the Gumbel distribution / Parameter der Gumbel-Verteilung
fh_gumb = @(x) cdf_GumbelDistribution(x, a, b);     % CDF function for Gumbel / CDF-Funktion für Gumbel-Verteilung
[accept_gumb, dmax_gumb] = K_S_test_fh(x, fh_gumb); % KS test for Gumbel / KS-Test für Gumbel

% Log-Normal Distribution / Logarithmische Normalverteilung
[mue, sigma] = EX_VarX_to_LogNormalDistribution(EX, VarX);  % Parameters of the LogNormal distribution / Parameter der Log-Normalverteilung
fh_log = @(x) cdf_LogNormalDistribution(x, mue, sigma);     % CDF function for Log-Normal / CDF-Funktion für Log-Normalverteilung
[accept_log, dmax_log] = K_S_test_fh(x, fh_log);            % KS test for Log-Normal / KS-Test für Log-Normalverteilung

% Normal Distribution / Normalverteilung
fh_norm = @(x) cdf_NormalDistribution(x, EX, VarX);         % CDF function for Normal distribution / CDF-Funktion für Normalverteilung
[accept_norm, dmax_norm] = K_S_test_fh(x, fh_norm);         % KS test for Normal distribution / KS-Test für Normalverteilung

% Uniform Distribution / Gleichverteilung
a = min(x); b = max(x);                     % Parameters of Uniform distribution / Parameter der Gleichverteilung
fh_uni = @(x) cdf_UniformDistribution(x, a, b); % CDF function for Uniform / CDF-Funktion für Gleichverteilung
[accept_uni, dmax_uni] = K_S_test_fh(x, fh_uni); % KS test for Uniform distribution / KS-Test für Gleichverteilung

% Weibull Distribution / Weibull-Verteilung
skal = 200;                                 % Scaling factor / Skalierungsfaktor
EX = EX / skal; 
VarX = VarX / skal^2;
[a, b] = EX_VarX_to_WeibullDistribution2(EX, VarX); % Parameters of Weibull distribution / Parameter der Weibull-Verteilung
[EX__, VarX__] = EX_VarX_of_WeibullDistribution(a, b); % Mean and variance of Weibull / Mittelwert und Varianz der Weibull-Verteilung
fh_wb = @(x) cdf_WeibullDistribution(x, a, b);        % CDF function for Weibull / CDF-Funktion für Weibull-Verteilung
[accept_wb, dmax_wb] = K_S_test_fh(x / skal, fh_wb);  % KS test for Weibull / KS-Test für Weibull



%% Plots / Grafiken
hold on;

% Plot CDFs of each distribution / CDFs jeder Verteilung plotten
fplot(fh_gumb, [min(x), max(x)], 'LineWidth', 1);  % Gumbel distribution / Gumbel-Verteilung
fplot(fh_log, [min(x), max(x)], 'LineWidth', 1);   % Log-normal distribution / Log-Normalverteilung
fplot(fh_norm, [min(x), max(x)], 'LineWidth', 1);  % Normal distribution / Normalverteilung
fplot(fh_uni, [min(x), max(x)], 'LineWidth', 1);   % Uniform distribution / Gleichverteilung
fplot(fh_wb, [min(x), max(x)], 'LineWidth', 1);    % Weibull distribution / Weibull-Verteilung

axis([min(x) max(x) -0.2 1.1]);                    % Set axis limits / Achsenbegrenzungen setzen

% CDF: Elasticity Modulus (E modulus) / CDF: Elastizitätsmodul (E-Modul)
n = length(E);
x_plot = sort(E);                                  % Sort data for plotting / Daten zum Plotten sortieren
y_plot = [1/n:1/n:1];                              % Generate CDF values / CDF-Werte generieren
stairs(x_plot, y_plot, 'LineWidth', 2);            % Plot stair-step function / Treppenfunktion plotten

title('Comparison: Distribution function of E modulus with typical distributions / Vergleich: Verteilungsfunktion des E-Moduls mit typischen Verteilungen');  % Plot title / Titel des Plots
legend({'Gumbel', 'Log-Normal', 'Normal', 'Uniform', 'Weibull', 'E-Modul'}, 'Location', 'southeast');  % Legend
xlabel('E-Modulus in GPa'); ylabel('CDF');         % Axis labels / Achsenbeschriftungen

% Display acceptance and deviations / Akzeptanz und Abweichungen anzeigen
disp('All distributions are accepted / Alle Verteilungen werden akzeptiert'); % Display message / Nachricht anzeigen
Akzeptiert = [accept_gumb, accept_log, accept_norm, accept_uni, accept_wb] % Accepted distributions / Akzeptierte Verteilungen

disp('Maximum deviations / Maximale Abweichungen'); % Display message / Nachricht anzeigen
dmax = [dmax_gumb; dmax_log; dmax_norm; dmax_uni; dmax_wb] % Max deviations / Maximale Abweichungen

disp('Smallest maximum deviation: Normal distribution / Kleinste maximale Abweichung: Normalverteilung'); % Display result / Ergebnis anzeigen
dmax_min = min(dmax);                              % Minimum deviation / Kleinste Abweichung
dmax_norm                                          % Display max deviation for Normal distribution / Maximale Abweichung der Normalverteilung anzeigen








                                      