clc; clear all; close all;
%Vergleich der expliziten Ableitung mit der FDM
%% Parameter für die Validierung der Ableitung
q0 = 0.5; %[N/mm]; 
NE = 100; 
s = (-1)^NE;
lengthh = 1500;
lengthhElement = lengthh/NE; 
E_mue    = 205.6378e3;


%Profilparameter 
b_ref    = 50;             % [mm]
d_ref    = 50;             % [mm] 
t1_ref   = 1;              % [mm] 
t2_ref   = 3.571428571;    % [mm]
b    = b_ref * ones(NE,1);         % [mm] 
d    = d_ref * ones(NE,1);         % [mm] 
t1   = t1_ref * ones(NE,1);        % [mm] 
t2   = t2_ref * ones(NE,1);        % [mm] 

%% FEM Rechnung
Le = (lengthhElement) * ones(NE,1);     
E    = E_mue * ones(NE,1);         % [mm]
I = (b.*d.^3 - (b-t1).*(d-2*t2).^3)/12; % Flaechentraegheitsmoment I-Profil (konstant in allen Elementen)
f = 0.5*q0*(lengthhElement+lengthhElement);
[u,Ke,ue,K] = KragarmFEM(Le, E, I, f);

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

%% Ableitung der maximalen Verschiebung in Balkenmitte nach E (Adjungierten-Methode)
for i=1:NE;

dKidEi = Ke(:,:,i)/E(i);                     %Nur die Elementsteifigkeitsmatrix durch E(i) geteilt                                                                            
dcdE(i)=lambdae(:,i)'*dKidEi*ue(:,i);
end
dcdE=transpose(dcdE);
dcdE;

%% Finite Differenzen Ableitung der Nachgiebigkeit nach E
E    = E_mue * ones(NE,1);         % [mm]
b    = b_ref * ones(NE,1);         % [mm] 
d    = d_ref * ones(NE,1);         % [mm] 
t1   = t1_ref * ones(NE,1);        % [mm] 
t2   = t2_ref * ones(NE,1);        % [mm]

for i=1:NE;
E0 = E;                        
% E  = E0+0.0001*ones(length(E0),1);
E(i)  = E0(i)+1;
u_p = KragarmFEM( Le,E,I,f );
if s== -1
    c_p = u_p(NE); % NE is odd
else
    c_p = u_p(NE+1); % Ne is even
end

% E  = E0-0.0001*ones(length(E0),1);
E(i)  = E0(i)-1;
u_m = KragarmFEM( Le,E,I,f );
if s== -1
    c_m = u_m(NE); % NE is odd
else
    c_m = u_m(NE+1); %NE is even
end
% dcdE_FDM(i)=(c_p-c_m)/(2*0.0001); 
dcdE_FDM(i)=(c_p-c_m)/(2*1); 
end
dcdE_FDM;

%% Plot: Vergleich der Ableitungen explizit und FDM
fig1 = figure;
% movegui(fig1,[300 -50]);
sgtitle('Vergleich explizit mit FDM für FOSM mit N=100, Ableitung nach E')

%Vergleich dcdE dcdE_FDM
x_plot = linspace(1,NE,NE);
dcdE_FDM_plot = full(dcdE_FDM);
y_plot = [dcdE,dcdE_FDM_plot'];
h = plot(x_plot,y_plot);
h(1).LineWidth = 2;
h(1).Color = 'red';
h(2).Color = 'blue';
h(2).LineWidth = 1;
legend({'dcdE','dcdE_{FDM}'},'Location','southeast');
xlabel('Elementanzahl [-]'), ylabel('dcdE');

