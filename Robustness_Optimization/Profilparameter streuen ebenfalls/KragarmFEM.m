function [u,Ke,ue,K] = KragarmFEM(L,E,I,f)
%% einfaches Kragarm-Balken-FE-Programm
% input:    L (Vektor)  Element-Längen
%           E (Vektor)  E-Modul pro Element
%           I (Vektor)  FTM pro Element
%           f (Skalar)  Last am Balkenende
% output:   u (Vektor)  alle Freieheitsgrade an allen Knoten
%           w (Vektor)  Verschiebungen pro Knoten
%           r (Vektor)  Verdrehungen pro Knoten
%           K (4x4xN-Tensor)  Elementsteifigkeitsmatrizen
%           ue (4xN-Matrix)  Freieheitsgrade je Element

    n=length(L);
    if length(E)==1
        E=E*ones(n,1);
    end
    if length(I)==1
        I=I*ones(n,1);
    end
    if length(E)~= n || length(I)~=n
        disp('L, E und I müssen gleich lang sein'); return;
    end
    [K,Ke]=SteifigkeitGesamtBalken(L,E,I,n);     %Steifigkeitsmatrix
    F=zeros(size(K,1)-4,1);         %Kraftvektor
    F(1:2:end)=f;                     %    "
%     F=sparse(F);                    %    "
    RB3 = (2*(n)+1); %Randbedinung ulast = 0,rlast = 0
    RB4 = (2*(n)+2);   %Randbedinung ulast= 0,rlast = 0
    RB1 = 1; %Randbedinung u0=0,r0=0, linke Einspannung
    RB2 = 2; %Randbedinung u0=0,r0=0, rechte Einspannung   
%     RB3 = (2*(n)-1); %Randbedinung ulast = 0,rlast = 0
%     RB4 = (2*(n));   %Randbedinung ulast= 0,rlast = 0
    K(RB3:RB4,:)=[]; %Randbedinung ulast= 0,rlast = 0
    K(:,RB3:RB4)=[]; %Randbedinung ulast= 0,rlast = 0
    K(RB1:RB2,:)=[]; %Randbedinung u0=0,r0=0
    K(:,RB1:RB2)=[]; %Randbedinung u0=0,r0=0
%     K(RB3:RB4,:)=[]; %Randbedinung ulast= 0,rlast = 0
%     K(:,RB3:RB4)=[]; %Randbedinung ulast= 0,rlast = 0
    u=K\F;           %Loesen des LGS K*u=f
    u=[0;0;u];                      %Ergaenzen um Randbedingung linke Einspannung
    u=[u;0;0];                      %Ergaenzen um Randbedingung rechte Einspannung
    w=u(1:2:end);                   %Auslesen der Verschiebung (u i)
    r=u(2:2:end);                   %Auslesen der Rotation (phi i)
    for e=1:n
        ue(:,e) = u(2*e-1:2*e+2);
    end
end

function [K,Ke] = SteifigkeitGesamtBalken( L,E,I,n )
% Berechnet die Steifigkeitsmatrix eines Biegebalkens
    K=zeros(2*n+2);
    Ke = zeros(4,4,n);
    for i=1:n
        Ke(:,:,i) = ElementSteifigkeit(L(i),E(i),I(i));
        K(i*2-1:i*2+2,i*2-1:i*2+2)=...
            K(i*2-1:i*2+2,i*2-1:i*2+2)+Ke(:,:,i);
    end
%     K=sparse(K);
end

function Ke = ElementSteifigkeit(Le,Ee,Ie)
% Elementsteifigkeitsmatrix
    Ke=Ee*Ie/Le^3*...
      [12,   -6*Le,  -12   -6*Le; 
      -6*Le,  4*Le^2, 6*Le, 2*Le^2;
      -12,    6*Le,   12    6*Le;
      -6*Le,  2*Le^2, 6*Le, 4*Le^2];
end
