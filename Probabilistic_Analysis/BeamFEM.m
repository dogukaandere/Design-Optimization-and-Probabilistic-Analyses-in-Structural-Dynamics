function [u,Ke,ue,K] = BeamFEM(L,E,I,f)
%% Simple Cantilever Beam FE Program / Einfaches Kragarm-Balken-FE-Programm
% input:    L (vector)   Element lengths / Element-Längen
%           E (vector)   Young's modulus per element / E-Modul pro Element
%           I (vector)   Moment of inertia per element / Flächenträgheitsmoment pro Element
%           f (scalar)   Load at the end of the beam / Last am Balkenende
% output:   u (vector)   All degrees of freedom at all nodes / Alle Freiheitsgrade an allen Knoten
%           w (vector)   Displacements per node / Verschiebungen pro Knoten
%           r (vector)   Rotations per node / Verdrehungen pro Knoten
%           K (4x4xN tensor)  Element stiffness matrices / Elementsteifigkeitsmatrizen
%           ue (4xN matrix)   Degrees of freedom per element / Freiheitsgrade je Element

    n = length(L);
    if length(E) == 1
        E = E * ones(n, 1);
    end
    if length(I) == 1
        I = I * ones(n, 1);
    end
    if length(E) ~= n || length(I) ~= n
        disp('L, E, and I must be the same length'); return; % L, E und I müssen gleich lang sein
    end
    [K, Ke] = TotalBeamStiffness(L, E, I, n);     % Stiffness matrix / Steifigkeitsmatrix
    F = zeros(size(K, 1) - 4, 1);                      % Force vector / Kraftvektor
    F(1:2:end) = f;                                    % Assign load to force vector / Last zuweisen
    RB3 = (2 * (n) + 1); % Boundary condition u_last = 0, r_last = 0 / Randbedingung u_last = 0, r_last = 0
    RB4 = (2 * (n) + 2); % Boundary condition u_last = 0, r_last = 0 / Randbedingung u_last = 0, r_last = 0
    RB1 = 1;             % Boundary condition u_0 = 0, r_0 = 0, left support / Randbedingung u_0 = 0, r_0 = 0, linke Einspannung
    RB2 = 2;             % Boundary condition u_0 = 0, r_0 = 0, right support / Randbedingung u_0 = 0, r_0 = 0, rechte Einspannung
    
    K(RB3:RB4, :) = []; % Apply boundary condition u_last = 0, r_last = 0 / Randbedingung u_last = 0, r_last = 0 anwenden
    K(:, RB3:RB4) = [];
    K(RB1:RB2, :) = []; % Apply boundary condition u_0 = 0, r_0 = 0 / Randbedingung u_0 = 0, r_0 = 0 anwenden
    K(:, RB1:RB2) = [];
    
    u = K \ F;           % Solve the system K*u = f / Gleichungssystem K*u = f lösen
    u = [0; 0; u];       % Add boundary condition at left support / Randbedingung an linker Einspannung ergänzen
    u = [u; 0; 0];       % Add boundary condition at right support / Randbedingung an rechter Einspannung ergänzen
    w = u(1:2:end);      % Extract displacements (u_i) / Verschiebungen extrahieren (u_i)
    r = u(2:2:end);      % Extract rotations (phi_i) / Rotationen extrahieren (phi_i)
    for e = 1:n
        ue(:, e) = u(2 * e - 1 : 2 * e + 2); % Assign degrees of freedom per element / Freiheitsgrade je Element zuweisen
    end
end

function [K, Ke] = TotalBeamStiffness(L, E, I, n)
% Computes the stiffness matrix of a bending beam / Berechnet die Steifigkeitsmatrix eines Biegebalkens
    K = zeros(2 * n + 2);
    Ke = zeros(4, 4, n);
    for i = 1:n
        Ke(:, :, i) = ElementStiffness(L(i), E(i), I(i)); % Element stiffness / Elementsteifigkeit
        K(i * 2 - 1 : i * 2 + 2, i * 2 - 1 : i * 2 + 2) = ...
            K(i * 2 - 1 : i * 2 + 2, i * 2 - 1 : i * 2 + 2) + Ke(:, :, i);
    end
end

function Ke = ElementStiffness(Le, Ee, Ie)
% Element stiffness matrix / Elementsteifigkeitsmatrix
    Ke = Ee * Ie / Le^3 * ...
      [12,   -6 * Le,  -12,  -6 * Le; 
      -6 * Le,  4 * Le^2, 6 * Le, 2 * Le^2;
      -12,    6 * Le,   12,    6 * Le;
      -6 * Le,  2 * Le^2, 6 * Le, 4 * Le^2];
end

