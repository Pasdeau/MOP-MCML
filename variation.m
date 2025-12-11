close all; clear all;

T = readtable('variation_d.xlsx');
d = T.whiteAndGrayMatter_d_;
R1 = T.R1;
R3 = T.R3;
T1 = T.T1;
T3 = T.T3;

subplot(2,2,1)
plot(d,R1,'-o')
xlabel('d'); ylabel('R\_MCML');
title('sans\_os : variation de ''d'' de la couche white matter and gray matter')

subplot(2,2,3)
plot(d,R3,'-o')
xlabel('d'); ylabel('R\_SL');
title('sans\_os : variation de ''d'' de la couche white matter and gray matter')

subplot(2,2,2)
% plot(d,T1,'-o')
scatter(d,T1,'o')
hold on
p=polyfit(d,T1,3)
y1=polyval(p,d);
plot(d,y1,'b','LineWidth',2,'DisplayName','T__MCML = -0.0001*d^3+0.0007*d^2-0.0015*d+0.0012');
xlabel('d'); ylabel('T\_MCML');
title('sans\_os : variation de ''d'' de la couche white matter and gray matter')
set(gca,'yscale','log')
legend

subplot(2,2,4)
scatter(d,T3,'o')
hold on
p=polyfit(d,T3,5)
y1=polyval(p,d);
plot(d,y1,'b','LineWidth',2,'DisplayName','T__SL = -0.000004*d^5+0.0000309*d^4-0.000923*d^3+0.001351*d^2-0.001008*d+0.0334');
xlabel('d'); ylabel('T\_SL');
title('sans\_os : variation de ''d'' de la couche white matter and gray matter')
set(gca,'yscale','log')
legend

x=[0.3:0.0001:2.4];
y=exp(-2.1464*x-10.7886);

figure(2)
semilogy(d,T3,'o--','LineWidth',2,'DisplayName','Définition 3');
hold on
plot(x,y,'r','LineWidth',2,'DisplayName','Courbe de régression : T=exp(-10.79-2.15d_v_a_r_i )')
% p=polyfit(d,T3,5)
% y1=polyval(p,d);
% plot(d,y1,'b','LineWidth',2,'DisplayName','T__SL = -0.000004*d^5+0.0000309*d^4-0.000923*d^3+0.001351*d^2-0.001008*d+0.0334');
xlabel('d_v_a_r_i (cm)'); ylabel('T (sans unité)');
%title('sans\_os : variation de ''d'' de la couche white matter and gray matter')
%set(gca,'yscale','log')
legend
box off
set(gca,'FontWeight','Bold','FontSize',20);

% x=[0.3:0.0001:2.4];
% y=exp(-2.1464*x--10.7886);
% subplot(1,2,2)
% scatter(d,log(T3),'o')
% hold on
% p=polyfit(d,log(T3),1)
% y1=polyval(p,d);
% plot(d,y1,'b','LineWidth',2,'DisplayName','Courbe de regression');
% xlabel('Δd'); ylabel('lnT');
% title('Courbe de regression')
% %set(gca,'yscale','log')
% legend