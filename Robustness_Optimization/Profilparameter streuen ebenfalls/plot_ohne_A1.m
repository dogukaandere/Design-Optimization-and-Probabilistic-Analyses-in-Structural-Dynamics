function plot_ohne_A1(b,d,t1,t2,NE)
fig2 = figure;
sgtitle('Vergleich der optimalen Entwurfsparameter (RDO NE=20 w1=1, w2=10) mit den Referenzgroessen')
subplot(2,2,1)
% Vergleich des Parameters b
x_plot = linspace(0,NE,NE);
b_ref=50*ones(NE,1);
y = [b_ref,b];
h = stairs(x_plot,y);
axis([0 NE 0 inf]);
title('Vergleich b mit b_{ref}')
h(1).LineWidth = 2;
h(1).Color = 'red';
h(2).Color = 'blue';
h(2).LineWidth = 1;
legend({'b_{ref}','b_{opt}'},'Location','northwest');
xlabel('Elementanzahl [-]'), ylabel('Parameter b [mm]');

subplot(2,2,2)
% Vergleich des Parameters d
x_plot = linspace(0,NE,NE);
d_ref=50*ones(NE,1);
y = [d_ref,d];
h = stairs(x_plot,y);
axis([0 NE 0 inf]);
title('Vergleich d mit d_{ref}')
h(1).LineWidth = 2;
h(1).Color = 'red';
h(2).Color = 'blue';
h(2).LineWidth = 1;
legend({'d_{ref}','d_{opt}'},'Location','southwest');
xlabel('Elementanzahl [-]'), ylabel('Parameter d [mm]');

subplot(2,2,3)
% Vergleich des Parameters t1
x_plot = linspace(0,NE,NE);
t1_ref=1*ones(NE,1);
y = [t1_ref,t1];
h = stairs(x_plot,y);
title('Vergleich t1 mit t1_{ref}')
ylim([-0.2 1.2]);
h(1).LineWidth = 2;
h(1).Color = 'red';
h(2).Color = 'blue';
h(2).LineWidth = 1;
legend({'t1_{ref}','t1_{opt}'},'Location','southwest');
xlabel('Elementanzahl [-]'), ylabel('Parameter t1 [mm]');

subplot(2,2,4)
% Vergleich des Parameters t2
x_plot = linspace(0,NE,NE);
t2_ref=3.571428571*ones(NE,1);
y = [t2_ref,t2];
h = stairs(x_plot,y);
axis([0 NE 0 inf]);
title('Vergleich t2 mit t2_{ref}')
h(1).LineWidth = 2;
h(1).Color = 'red';
h(2).Color = 'blue';
h(2).LineWidth = 1;
legend({'t2_{ref}','t2_{opt}'},'Location','southwest');
xlabel('Elementanzahl [-]'), ylabel('Parameter t2 [mm]');
end
