% Generate Fig. 2 for Eusipco paper
%
%  display switching function and allpass accuracy

close all; clear all;

Nfft = 2^10;
W = (0:Nfft-1)/Nfft;
%Sigma = .5 + .5*cos(2*pi*W) + sin(4*pi*W);
Sigma = .5 + .5*sin(2*pi*W) - sin(4*pi*W);
figure(1);
plot(W,Sigma);

% approximate zero crossings:
Wo = [0.0588 0.1575 0.5326 0.7500]*2*pi;

%------------------------------------------------------------------------------
%  generate and plot switching functions
%------------------------------------------------------------------------------
%  *** higher order design (Markus Lang)
NN = 2.^(1:4);
H2 = zeros(Nfft,length(NN));
Angle2 = zeros(Nfft,length(NN));
for n = 1:length(NN),
%   [b,a] = SingularValueAllpassSwitch(Wo,NN(n));
   [b,a] = SingularValueAllpassSwitchCompact(Wo,NN(n));
   h = impz(b,a,Nfft);
   h2 = circshift(h,-NN(n)+2);
   Angle2(:,n) = unwrap(angle(fft(h2,Nfft)));
   H2(:,n) = cos(unwrap(angle(fft(h2,Nfft))));
%   plot(w,cos(unwrap(angle(fft(h2,1024)))));
end;


%------------------------------------------------------------------------------
%  Figure 2a
%------------------------------------------------------------------------------
FS=12;
figure(1); clf;
Angle0 = zeros(size(W));
for i =1:4,
    Angle0 = Angle0 - (W>Wo(i)/2/pi);
end;
h = plot(W,Angle0,'-','linewidth',3); set(h(1),'color',[1 1 1]*.75);
hold on;    
plot(W,Angle2(:,1)/pi,'b-'); 
plot(W,Angle2(:,2)/pi,'r--'); 
h = plot(W,Angle2(:,3)/pi,'-.'); set(h(1),'color',[0 0.5 0]);
plot(W,Angle2(:,4)/pi,'k:');
axis([0 1 -4.2 0.2]);
text(0.01,-3.8,'(a)','interpreter','latex'); 
ylabel('phase response $\Phi(\Omega)$',...
	'interpreter','latex','fontsize',FS);
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$',...
      '$\pi$','$5\pi/4$','$3\pi/2$','$7\pi/4$','$2\pi$'},...
    'YTick',(-4:1:0),'YTickLabel',...
     {'$-4\pi$','$-3\pi$','$-2\pi$','$-\pi$','$0$'});
     grid on;
legend({'ideal','$J=2$','$J=4$','$J=8$','$J=16$'},...
       'interpreter','latex','fontsize',FS-2,...
        'location','NorthEast');
set(gcf,'OuterPosition',[230 250 570 250]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Eus26_1Fig2a.eps
 
%------------------------------------------------------------------------------
%  Figure 2b
%------------------------------------------------------------------------------
figure(2); clf;
h = plot(W,sign(Sigma),'-','linewidth',3); set(h(1),'color',[1 1 1]*.75);
hold on;
plot(W,H2(:,1),'b-'); 
plot(W,H2(:,2),'r--'); 
h = plot(W,H2(:,3),'-.'); set(h(1),'color',[0 0.5 0]);
plot(W,H2(:,4),'k:'); 

axis([0 1 -1.1 1.1]);
text(0.01,-.9,'(b)','interpreter','latex'); 
ylabel('freq. response $f(\mathrm{e}^{\mathrm{j}\Omega})$',...
	'interpreter','latex','fontsize',FS);
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$',...
      '$\pi$','$5\pi/4$','$3\pi/2$','$7\pi/4$','$2\pi$'},...
    'YTick',(-1:.5:1),'YTickLabel',...
     {'$-1$','$-.5$','$0$','$.5$','$1$'});
     grid on;
legend({'$\mathrm{sgn}(\sigma(\mathrm{e}^{\mathrm{j}\Omega})$)'},...
       'interpreter','latex','fontsize',FS-2,...
        'location','East');
set(gcf,'OuterPosition',[230 250 570 230]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Eus26_1Fig2b.eps

%------------------------------------------------------------------------------
%  Figure 2c
%------------------------------------------------------------------------------
Ssum = sum(abs(Sigma));
Int = zeros(1,size(H2,2));
for n = 1:size(H2,2),
   H2(:,n) = H2(:,n).*Sigma.';
   Int(n) = (Ssum-sum(H2(:,n)))/Ssum;
end;
figure(3); clf; 
%subplot(212);
h = plot(W,Sigma,'--','linewidth',1); set(h(1),'color',[1 1 1]*.75);
hold on;
h = plot(W,abs(Sigma),'-','linewidth',3); set(h(1),'color',[1 1 1]*.75);
plot(W,H2(:,1),'b-');
plot(W,H2(:,2),'r--');
h = plot(W,H2(:,3),'-.'); set(h(1),'color',[0 0.5 0]);
plot(W,H2(:,4),'k:');
%plot(W,H2(:,6:end),'--');
for n = 1:4, plot(Wo(n)/(2*pi)*[1 1],[-.5 1.5],'k--'); end;
text(Wo(1)/(2*pi)-.01,1.75,'$\Omega_1$','interpreter','latex');
text(Wo(2)/(2*pi)-.01,1.75,'$\Omega_2$','interpreter','latex');
text(Wo(3)/(2*pi)-.01,1.75,'$\Omega_3$','interpreter','latex');
text(Wo(4)/(2*pi)-.01,1.75,'$\Omega_4$','interpreter','latex');
axis([0 1 -1 2]); 
grid on;
text(0.01,-.8,'(c)','interpreter','latex'); 
ylabel('$f(\mathrm{e}^{\mathrm{j}\Omega}) \sigma(\mathrm{e}^{\mathrm{j}\Omega})$',...
	'interpreter','latex','fontsize',FS);
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$',...
      '$\pi$','$5\pi/4$','$3\pi/2$','$7\pi/4$','$2\pi$'},...
    'YTick',(-1:.5:2),'YTickLabel',...
     {'$-1$','$-.5$','$0$','$.5$','$1$','$1.5$','$2$'});
     grid on;
xlabel('normalised angular frequency $\Omega$','interpreter','latex','fontsize',FS);
legend({'$\sigma(\mathrm{e}^{\mathrm{j}\Omega})$',...
        '$|\sigma(\mathrm{e}^{\mathrm{j}\Omega})|$'},...
       'interpreter','latex','fontsize',FS-2,...
        'location','SouthEast');
set(gcf,'OuterPosition',[230 250 570 280]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Eus26_1Fig2c.eps

return;

%---------------------------------------------------
%  determine zero crossings
%---------------------------------------------------
W0 = [pi/2    2*pi/3;
      7/8*pi  9/8*pi;
      4*pi/3  5*pi/3;
      5*pi/3  2*pi];
Theta = zeros(4,1);
for i = 1:4,
   Wlower = W0(i,1); Wupper = W0(i,2);
   flower = .5 + .5*cos(Wlower) + sin(2*Wlower); fupper =   .5 + .5*cos(Wupper) + sin(2*Wupper);
   fmid = 1;
   while (abs(fmid)>1e-12);
      Wmid = (abs(flower)*Wlower + abs(fupper)*Wupper)/(abs(flower) + abs(fupper));
      fmid =   .5 + .5*cos(Wmid) + sin(2*Wmid);
      if sign(fmid) == sign(flower),
         flower = fmid; Wlower = Wmid;
      else
         fupper = fmid; Wupper = Wmid;
      end;
   end;         
   Theta(i) = Wmid;
end;

%----------------------------------------------------
%  calculate allpasses
%----------------------------------------------------
Nfft = 2^14;
NN = 2.^(1:8);
H2 = zeros(Nfft,length(NN));
Angle2 = zeros(Nfft,length(NN));
for n = 1:length(NN),
    disp(sprintf('allpass design for order %d',NN(n)));
   [b,a] = SingularValueAllpassSwitchCompact(Theta,NN(n));
   h = impz(b,a,Nfft);
   h2 = circshift(h,-NN(n)+2);
   Angle2(:,n) = unwrap(angle(fft(h2,Nfft)));
   H2(:,n) = cos(unwrap(angle(fft(h2,Nfft))));
%   plot(w,cos(unwrap(angle(fft(h2,1024)))));
end;

%----------------------------------------------------
%  numerical integration
%----------------------------------------------------
W = (0:(Nfft-1))'/Nfft*2*pi;
f = .5 + .5*cos(W) + sin(2*W);
F = zeros(length(NN),1);
for n = 1:length(NN),
   F(n) = sum( (abs(H2(:,n).*f)-1).^2 )/Nfft;  
end;
Flimit = sum( (abs(f)-1).^2 )/Nfft;

%---------------------------------------------------
%  Figure 10 --- accuracy of allpass
%---------------------------------------------------
figure(10); clf;
% dummy plots for legend
plot(-1,-1,'b-'); hold on;
h = plot(-1,-1,'-','linewidth',3); set(h(1),'color',[1 1 1]*.75);
% actual plots
h = plot(log2(NN),-1-2*log2(NN),'-','linewidth',3); set(h(1),'color',[1 1 1]*.75);
plot((1:8),log2(F-Flimit),'b-')
axis([1 8 -18 -2]);
grid on;
xlabel('$\log_{2} J$','interpreter','latex','fontsize',FS);
ylabel('$\log_{2} (e_{\mathrm{LS}}-\epsilon_{\mathrm{LS}})$','interpreter','latex','fontsize',FS);
legend({'$e_{\mathrm{LS}}-\epsilon_{\mathrm{LS}}$','power law $1/J^2$'},...
       'interpreter','latex','fontsize',FS-2,...
        'location','SouthWest');
        set(gcf,'OuterPosition',[230 250 570 240]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig10.eps     

  

%-------------------------------------------------------------------------
%    cost function
%-------------------------------------------------------------------------
function costval = CostFunction2(X)
rho=X(1); phi=X(2); theta=X(3); j = sqrt(-1);     % 1st allpass
b1 = [rho*exp(j*phi) -1]*exp(j*theta);
a1 = [1 -rho*exp(-j*phi)]; 
rho=X(4); phi=X(5); theta=X(6);                   % 2nd allpass
b2 = [rho*exp(j*phi) -1]*exp(j*theta);
a2 = [1 -rho*exp(-j*phi)]; 
b = conv(b1,b2); a = conv(a1,a2);
h = impz(b,a,512);
f = zeros(1023,1);
f(512:1023) = h;
f = f + conj(flipud(f));
f2 = zeros(1024,1);
f2(1:512) = f(512:1023);
f2(514:1024) = f(1:511);
F = fft(f2,1024);

W = (0:1023)'/1024;
S = .5 + .5*cos(2*pi*W) + sin(4*pi*W);
costval = 1/sum(real(F.*S));
end
