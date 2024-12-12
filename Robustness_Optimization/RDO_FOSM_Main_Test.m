%% Aufgabe3_RDO_FOSM_main

clc;
clear all;
close all;
%% Parameter
NE               = 20;                                        % Elementanzahl
% NE                 = 50;                                        % Elementanzahl
% NE                 = 100;                                        % Elementanzahl
lengthh            = 1500;                                      % Balkenlänge           [mm]
EModul             = 210e3;                                     % E-Modul               [N/mm^2]  
lengthhElement     = lengthh/NE;                                % Elementenlängen       [mm]
q0                 = 0.5;                                       % Wirkende Streckenlast [N/mm]        
f                  = 0.5*q0*(lengthhElement+lengthhElement);    % Belastung eines inneren Knotens des Balkens durch Streckenlast q0  [N]
faussen            = 0.5*q0*lengthhElement;                     % Belastung eines äußeren Knotens des Balkens durch Streckenlast q0  [N]
%% Referenzentwurf I-Profil mit den Parametern b,d,t1,t2

b_ref    = 50;             % [mm] b0,d0,t10,t20 in allen Elementen konstant
d_ref    = 50;             % [mm] 
t1_ref   = 1;              % [mm] 
t2_ref   = 3.571428571;    % [mm] 

V_zul = 6e5;    % [mm^3] zulässiges Maximales Volumen des Balkens
I_ref=(b_ref*d_ref^3-(b_ref-t1_ref)*(d_ref-2*t2_ref)^3)/12; % Flächenträgheitsmoment mit Referenzparametern [mm^4] 
c_ref=(q0*lengthh^4)/(384*EModul*I_ref); % Maximale Verschiebung in Balkenmitte mit Referenzparametern [mm]

%Startvektor x0
x0(1:4:4*NE,1)      = b_ref;    
x0(2:4:4*NE,1)      = d_ref;
x0(3:4:4*NE,1)      = t1_ref;
x0(4:4:4*NE,1)      = t2_ref;

Le(1:NE)            = lengthh/NE;   %Elementlaenge
E(1:NE)             = EModul;

%% fmincon Solver

c_rdo       = @(x) Zielfunktion(x);
nonlcon     = @(x) Nebenbedingungen(x);

A_element    = [-1 0 1 0; 0 -1 0 2];   %Lineare Ungleichheitsbedingungsmatrix für ein Element t1_i - b_i <= 0 und 2*t2_i -d_i <=0
Cellarray = repmat({A_element}, 1, NE);
A  = blkdiag(Cellarray{:});            % für A*x <= b mit x = Entwurfsgrößen
b  = zeros(2*NE,1);                    % für rechte Seite der Nebenbedingungen

Aeq = [];                              %keine Gleichheitsbedingungen
beq = [];

ub           = repmat([50; 50; 25; 25],NE,1); %upper bound für die Entwurfsgrößen
lb           = 1*ones(4*NE,1);                %lower bound für die Entwurfsgrößen

options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',1000000,'FiniteDifferenceStepSize',0.01);

% options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',10000000000000000,'FiniteDifferenceStepSize',0.000001);


[x_opt] = fmincon(c_rdo, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

b  = x_opt(1:4:4*NE);   
d  = x_opt(2:4:4*NE);
t1 = x_opt(3:4:4*NE);
t2 = x_opt(4:4:4*NE);

b_ref  = x0(1:4:4*NE);   
d_ref  = x0(2:4:4*NE);
t1_ref = x0(3:4:4*NE);
t2_ref = x0(4:4:4*NE);

%% Vergleich der Nachgiebigkeit zwischen x_opt_Aufgabe3, x_ref und x_opt_Aufgabe1
% x_opt_Aufgabe3     Ergebnisse der Optimierung in Aufgabe3
% x_opt_ref          Referenzentwurfsparameter
% x_opt_A1           Ergebnisse der Optimierung in Aufgabe1
%% Optimierungsvektor aus Aufgabe 01 mitberücksichtigen

if NE==20
load('x_opt_A1_NE_20.mat')                %Ergebnisse aus Aufgabe1_Main.m mit NE=20 Index "_opt_A1"
end 

if NE==50
load('x_opt_A1_NE_50.mat')                %Ergebnisse aus Aufgabe1_Main.m mit NE=50 Index "_opt_A1"
end 

if NE==100
load('x_opt_A1_NE_100.mat')                %Ergebnisse aus Aufgabe1_Main.m mit NE=100 Index "_opt_A1"
end 

b_opt_A1  = x_opt_A1(1:4:4*NE);   
d_opt_A1  = x_opt_A1(2:4:4*NE);
t1_opt_A1 = x_opt_A1(3:4:4*NE);
t2_opt_A1 = x_opt_A1(4:4:4*NE);

I = (b.*d.^3 - (b-t1).*(d-2*t2).^3)/12; % Flaechentraegheitsmoment
I_opt_A1 = (b_opt_A1.*d_opt_A1.^3 - (b_opt_A1-t1_opt_A1).*(d_opt_A1-2*t2_opt_A1).^3)/12; % Flaechentraegheitsmoment

[u] = KragarmFEM(Le, E, I, f);
[u_opt_A1] = KragarmFEM((lengthh/NE)*ones(NE,1), EModul*ones(NE,1), I_opt_A1, f); %Für Aufgabe1

disp('Durchbiegung in der Balkenmitte: Aufgabe3')
s = (-1)^NE;
if s== -1
    c_opt_A3 = u(NE)%NE ist ungerade
else
    c_opt_A3 = u(NE+1) %NE ist gerade
end
% c_opt_A3= full(c_opt_A3)

disp('Durchbiegung in der Balkenmitte: Aufgabe1')
if s== -1;
    c_opt_A1 = u_opt_A1(NE)%NE ist ungerade
else
    c_opt_A1 = u_opt_A1(NE+1) %NE ist gerade
end  %Nachgiegkeit mit den Optimierten Parametern aus Aufgabe 1 mit NE=20 / NE=50 /NE = 100
% c_opt_A1 = full(c_opt_A1)
disp('Durchbiegung in der Balkenmitte: Referenzentwurf')
c_ref


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
[~,dc] = Zielfunktion_FOSM(mue_x,b,d,t1,t2);
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
function [c,dc] = Zielfunktion_FOSM(x,b,d,t1,t2)
q0 = 0.5;
lengthh = 1500;
NE=length(x);
lengthhElement = lengthh/NE;
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


