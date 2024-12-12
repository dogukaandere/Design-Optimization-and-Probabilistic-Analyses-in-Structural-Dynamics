%% Plot for Comparison of Parameters / Plot für Vergleich von Parametern

function Plot_for_Comparison_of_Parameters(b, d, t1, t2, NE)
fig1 = figure;
sgtitle('Comparison of optimal design parameters (NE=200) with reference values / Vergleich der optimalen Entwurfsparameter (NE=200) mit den Referenzgrößen')

subplot(2,2,1)
% Comparison of parameter b / Vergleich des Parameters b
x_plot = linspace(0, NE, NE); % Generates NE points between 0 and NE / erzeugt NE Punkte zwischen 0 und NE
b_ref = 50 * ones(NE, 1);
y = [b_ref, b];
h = stairs(x_plot, y);
axis([0 NE 0 inf]);
title('Comparison of b with b_{ref} / Vergleich b mit b_{ref}')
h(1).LineWidth = 2;
h(1).Color = 'red';
h(2).Color = 'blue';
h(2).LineWidth = 1;
legend({'b_{ref}', 'b_{opt}'}, 'Location', 'northwest');
xlabel('Element count [-] / Elementanzahl [-]'), ylabel('Parameter b [mm] / Parameter b [mm]');

subplot(2,2,2)
% Comparison of parameter d / Vergleich des Parameters d
x_plot = linspace(0, NE, NE);
d_ref = 50 * ones(NE, 1);
y = [d_ref, d];
h = stairs(x_plot, y);
axis([0 NE 0 inf]);
title('Comparison of d with d_{ref} / Vergleich d mit d_{ref}')
h(1).LineWidth = 2;
h(1).Color = 'red';
h(2).Color = 'blue';
h(2).LineWidth = 1;
legend({'d_{ref}', 'd_{opt}'}, 'Location', 'southwest');
xlabel('Element count [-] / Elementanzahl [-]'), ylabel('Parameter d [mm] / Parameter d [mm]');

subplot(2,2,3)
% Comparison of parameter t1 / Vergleich des Parameters t1
x_plot = linspace(0, NE, NE);
t1_ref = 1 * ones(NE, 1);
y = [t1_ref, t1];
h = stairs(x_plot, y);
axis([0 NE 0 inf]);
title('Comparison of t1 with t1_{ref} / Vergleich t1 mit t1_{ref}')
h(1).LineWidth = 2;
h(1).Color = 'red';
h(2).Color = 'blue';
h(2).LineWidth = 1;
legend({'t1_{ref}', 't1_{opt}'}, 'Location', 'southwest');
xlabel('Element count [-] / Elementanzahl [-]'), ylabel('Parameter t1 [mm] / Parameter t1 [mm]');

subplot(2,2,4)
% Comparison of parameter t2 / Vergleich des Parameters t2
x_plot = linspace(0, NE, NE);
t2_ref = 3.571428571 * ones(NE, 1);
y = [t2_ref, t2];
h = stairs(x_plot, y);
axis([0 NE 0 inf]);
title('Comparison of t2 with t2_{ref} / Vergleich t2 mit t2_{ref}')
h(1).LineWidth = 2;
h(1).Color = 'red';
h(2).Color = 'blue';
h(2).LineWidth = 1;
legend({'t2_{ref}', 't2_{opt}'}, 'Location', 'southwest');
xlabel('Element count [-] / Elementanzahl [-]'), ylabel('Parameter t2 [mm] / Parameter t2 [mm]');

end


