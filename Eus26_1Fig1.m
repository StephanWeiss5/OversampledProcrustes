% Eus26_1Fig1.m
%
% generates Fig 1 of EUSIPCO'26 draft on analytic Procrustes solution

clear all; close all;
FS=12;
Omega = (0:0.01:0.5);
plot(Omega,2*cos(2*pi*Omega),'b');
hold on;
plot(Omega+0.5,2*cos(2*pi*(Omega+0.5)),'b--');
plot(0.25,0,'bo','MarkerFaceColor','b');
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/2$','$\pi$','$3\pi/2$',...
      '$2\pi$','$5\pi/2$','$3\pi$','$7\pi/2$','$4\pi$'},...
    'YTick',(-2:1:2),'YTickLabel',...
     {'$-2$','$-1$','$0$','$1$','$2$'});
ylabel('$\sigma\prime(\Omega)$',...
	'interpreter','latex','fontsize',FS);
axis([0 1 -2.1 2.1]);
grid on;
xlabel('normalised angular frequency $\Omega$','interpreter','latex','fontsize',FS);
legend({'$\sigma^\prime(\Omega)$','$-\sigma^\prime(\Omega-2\pi)$'},...
        'interpreter','latex','location','SouthEast');
set(gcf,'OuterPosition',[230 250 570 250]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Eus26_1Fig1.eps

