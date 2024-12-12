function [h, heq] = Constraints(x)
NE = length(x) / 4;          % Number of elements / Anzahl der Elemente
b  = x(1:4:4*NE);            % Width values / Breitenwerte
d  = x(2:4:4*NE);            % Depth values / Tiefenwerte
t1 = x(3:4:4*NE);            % Thickness t1 / Dicke t1
t2 = x(4:4:4*NE);            % Thickness t2 / Dicke t2
V_zul = 6e5;                 % Allowable maximum volume / Zulässiges maximales Volumen

Le = (1500 / NE);                % Element length / Elementlänge
Ve = Le * (2 * b .* t2 + (d - 2 * t2) .* t1); % Element volume = Le * Ae / Elementvolumen = Le * Ae
BV = sum(Ve);                    % Beam volume, scalar value / Balkenvolumen, Skalarer Wert
h = BV - V_zul;                  % Inequality constraint: Maximum volume <= 6*10^5 mm^3 / Ungleichheitsbedingung: Maximales Volumen <= 6*10^5 mm^3
heq = [];                        % No equality constraints / Keine Gleichheitsbedingungen
end
