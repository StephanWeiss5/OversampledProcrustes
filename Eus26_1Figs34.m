% Generate Figs. 3 and 4 for Eusipco paper
%
%  display switching function and allpass accuracy

close all; clear all;

Nfft = 2^10;
sigma = [-1i/2 0 1/2 0 1/2 0 1i/2];
% zero crossings:
Wo = [0.375 0.75 0.875];
Wo2 = [Wo Wo+1]*pi;
Wo=Wo*2*pi;
W = (0:Nfft-1)'/Nfft;
Sigma = cos(2*pi*W) + sin(6*pi*W);
Sigma_s = cos(pi*W) + sin(3*pi*W);

%------------------------------------------------------------------------------
%  generate and plot switching functions
%------------------------------------------------------------------------------
%  *** higher order design (Markus Lang)
NN = 3*2.^(1:4);
H2 = zeros(Nfft,length(NN));
Angle2 = zeros(Nfft,length(NN));
for n = 1:length(NN),
   [b,a] = SingularValueAllpassSwitchCompact(Wo2,NN(n));
%   [b,a] = AllpassSwitchCompactEusipco(Wo,NN(n));
   b = b(1:2:end); a = a(1:2:end);
   h = impz(b,a,Nfft);
   h2 = circshift(h,-NN(n)/2+2);
   Angle2(:,n) = unwrap(angle(fft(h2,Nfft)))+W*2*pi;
   H2(:,n) = cos(unwrap(angle(fft(h2,Nfft)))+W*pi);
%   plot(w,cos(unwrap(angle(fft(h2,1024)))));
end;


%------------------------------------------------------------------------------
%  Figure 3a
%------------------------------------------------------------------------------
FS=12;
figure(1); clf;
Angle0 = W; % zeros(size(W));
for i =1:3,
    Angle0 = Angle0 - (W>Wo(i)/2/pi);
end;
h = plot(W,Angle0,'-','linewidth',3); set(h(1),'color',[1 1 1]*.75);
hold on;    
plot(W,Angle2(:,1)/pi,'b-'); 
plot(W,Angle2(:,2)/pi,'r--'); 
h = plot(W,Angle2(:,3)/pi,'-.'); set(h(1),'color',[0 0.5 0]);
plot(W,Angle2(:,4)/pi,'k:');
axis([0 1 -2.2 0.7]);
text(0.01,0.3,'(a)','interpreter','latex'); 
ylabel('phase response $\Phi(\Omega)$',...
	'interpreter','latex','fontsize',FS);
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$',...
      '$\pi$','$5\pi/4$','$3\pi/2$','$7\pi/4$','$2\pi$'},...
    'YTick',(-2:1:0),'YTickLabel',...
     {'$-2\pi$','$-\pi$','$0$'});
     grid on;
legend({'ideal','$J=3$','$J=6$','$J=12$','$J=24$'},...
       'interpreter','latex','fontsize',FS-2,...
        'location','SouthWest');
set(gcf,'OuterPosition',[230 250 570 230]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Eus26_1Fig3a.eps
 
%------------------------------------------------------------------------------
%  Figure 3b
%------------------------------------------------------------------------------
figure(2); clf;
h = plot(W,sign(Sigma_s),'-','linewidth',3); set(h(1),'color',[1 1 1]*.75);
hold on;
plot(W,H2(:,1),'b-'); 
plot(W,H2(:,2),'r--'); 
h = plot(W,H2(:,3),'-.'); set(h(1),'color',[0 0.5 0]);
plot(W,H2(:,4),'k:'); 

axis([0 1 -1.1 1.1]);
text(0.01,-.9,'(b)','interpreter','latex'); 
ylabel('freq.~resp.~$f^\prime(\mathrm{e}^{\mathrm{j}\Omega})$',...
	'interpreter','latex','fontsize',FS);
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$',...
      '$\pi$','$5\pi/4$','$3\pi/2$','$7\pi/4$','$2\pi$'},...
    'YTick',(-1:.5:1),'YTickLabel',...
     {'$-1$','$-.5$','$0$','$.5$','$1$'});
     grid on;
legend({'$\mathrm{sgn}(\sigma^\prime(\mathrm{e}^{\mathrm{j}\Omega}))$'},...
       'interpreter','latex','fontsize',FS-2,...
        'location','West');
set(gcf,'OuterPosition',[230 250 570 210]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Eus26_1Fig3b.eps

%------------------------------------------------------------------------------
%  Figure 3c
%------------------------------------------------------------------------------
Ssum = sum(abs(Sigma_s));
Int = zeros(1,size(H2,2));
for n = 1:size(H2,2),
   H2(:,n) = H2(:,n).*Sigma_s;
   Int(n) = (Ssum-sum(H2(:,n)))/Ssum;
end;
figure(3); clf; 
h = plot(W,Sigma_s,'--','linewidth',1); set(h(1),'color',[1 1 1]*.75);
hold on;
h = plot(W,abs(Sigma_s),'-','linewidth',3); set(h(1),'color',[1 1 1]*.75);
plot(W,H2(:,1),'b-');
plot(W,H2(:,2),'r--');
h = plot(W,H2(:,3),'-.'); set(h(1),'color',[0 0.5 0]);
plot(W,H2(:,4),'k:');
for n = 1:3, plot(Wo(n)/(2*pi)*[1 1],[-.5 1.5],'k--'); end;
text(Wo(1)/(2*pi)-.01,1.75,'$\Omega_1$','interpreter','latex');
text(Wo(2)/(2*pi)-.01,1.75,'$\Omega_2$','interpreter','latex');
text(Wo(3)/(2*pi)-.01,1.75,'$\Omega_3$','interpreter','latex');
axis([0 1 -1.2 2]); 
grid on;
text(0.01,-.8,'(c)','interpreter','latex'); 
ylabel('$f^\prime(\mathrm{e}^{\mathrm{j}\Omega}) \sigma^\prime(\mathrm{e}^{\mathrm{j}\Omega})$',...
	'interpreter','latex','fontsize',FS);
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$',...
      '$\pi$','$5\pi/4$','$3\pi/2$','$7\pi/4$','$2\pi$'},...
    'YTick',(-1:.5:2),'YTickLabel',...
     {'$-1$','$-.5$','$0$','$.5$','$1$','$1.5$','$2$'});
     grid on;
xlabel('normalised angular frequency $\Omega$','interpreter','latex','fontsize',FS);
legend({'$\sigma^\prime(\mathrm{e}^{\mathrm{j}\Omega})$',...
       '$|\sigma^\prime(\mathrm{e}^{\mathrm{j}\Omega})|$'},...
       'interpreter','latex','fontsize',FS-2,...
        'location','West');
set(gcf,'OuterPosition',[230 250 570 260]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Eus26_1Fig3c.eps

%------------------------------------------------------------------------------
%  Figures 4a and 4b
%------------------------------------------------------------------------------
[b,a] = SingularValueAllpassSwitchCompact(Wo2,24);
h = impz(b,a,8192);
H = fft(circshift(h,-21),8192);
Angle0 = 0*W;
for i =1:6,
    Angle0 = Angle0 - (W>Wo2(i)/2/pi);
end;
figure(6); clf;
plot(W,Angle0/2,'-','linewidth',3,'color',[1 1 1]*.75);
hold on;
plot((0:8191)/8192,unwrap(angle(H))/2/pi,'b-'); grid on;
plot(Wo2(1)/(2*pi)*[1 1],[.1 -1],'k--'); 
text(Wo2(1)/(2*pi)-.01,-1.2,'$\frac{\Omega_1}{2}$','interpreter','latex');
plot(Wo2(2)/(2*pi)*[1 1],[.1 -1.5],'k--'); 
text(Wo2(2)/(2*pi)-.01,-1.7,'$\frac{\Omega_2}{2}$','interpreter','latex');
plot(Wo2(3)/(2*pi)*[1 1],[.1 -1.75],'k--'); 
text(Wo2(3)/(2*pi)-.01,-1.95,'$\frac{\Omega_3}{2}$','interpreter','latex');
plot(Wo2(4)/(2*pi)*[1 1],[-.75 -3],'k--');
text(Wo2(4)/(2*pi)-.035,-.55,'$\frac{\Omega_1}{2}+\pi$','interpreter','latex');
plot(Wo2(5)/(2*pi)*[1 1],[-1.25 -3],'k--');
text(Wo2(5)/(2*pi)-.035,-1.05,'$\frac{\Omega_2}{2}+\pi$','interpreter','latex');
plot(Wo2(6)/(2*pi)*[1 1],[-1.75 -3],'k--');
text(Wo2(6)/(2*pi)-.035,-1.55,'$\frac{\Omega_3}{2}+\pi$','interpreter','latex');
text(0.01,-.5,'(a)','interpreter','latex'); 
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$',...
      '$\pi$','$5\pi/4$','$3\pi/2$','$7\pi/4$','$2\pi$'},...
    'YTick',(-2:1:0),'YTickLabel',...
     {'$-4\pi$','$-2\pi$','$0$'});
     grid on;
axis([0 1 -3.1 .1]);     
legend({'ideal','$J=24$'},...
       'interpreter','latex','fontsize',FS-2,...
        'location','SouthWest');
xlabel('normalised angular frequency $\Omega$','interpreter','latex','fontsize',FS);
ylabel('$\angle\{P_{\mathrm{OS}}^\prime(\mathrm{e}^{\mathrm{j}\Omega})\}$',...
	'interpreter','latex','fontsize',FS);
set(gcf,'OuterPosition',[230 250 570 270]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Eus26_1Fig4a.eps

figure(7); clf;
stem(-20:40,real([h(2:62)]),'b'); hold on;
plot(-20:40,imag([h(2:62)]),'r*'); grid on;
text(-19,.4,'(b)','interpreter','latex'); 
xlabel('discrete time $\nu$','interpreter','latex','fontsize',FS);
ylabel('$p_{\mathrm{OS}}^\prime[\nu]$',...
	'interpreter','latex','fontsize',FS);
set(gca,'TickLabelInterpreter','latex',...
        'YTick',(-.5:.5:.5),'YTickLabel',...
     {'$-\frac12$','$0$','$\frac12$'});
legend({'real part','imag.~part'},'interpreter','latex',...
          'location','SouthEast');     
set(gcf,'OuterPosition',[230 250 570 230]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Eus26_1Fig4b.eps

