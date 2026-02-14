% Generate Figs. 4 and 5 for Eusipco'26 paper on analytic Procrustes solution
%
% display result with existing method (TSP'26)

close all; clear all;

%-------------------------------------------------------------------------------
%  parameters
%-------------------------------------------------------------------------------
M = 32; 

if  (exist('Eus26_1Fig5.mat','file')~=2)
  %-------------------------------------------------------------------------------
  %  matrix A(z)
  %-------------------------------------------------------------------------------
  rand('seed',0);
  Sigma = zeros(M,M,3);
  for m=1:M/2,                            % 1/2 of the SVs with 2 zero crossings
    dummy = randn(1,1,3)+1i*rand(1,1,3);
    dummy = dummy + conj(dummy(1,1,3:-1:1));
    Sigma(m,m,:) = dummy/norm(squeeze(dummy));
  end;
  for m = M/2+1:M,                        % 1/2 of the SVs with 1 zero crossing
    theta = rand(1,1);
    if (abs(theta)-.5)<0.05, theta = theta+.25; end;
    Sigma(m,m,1:2) = [exp(-1i*theta*pi) exp(1i*theta*pi)]/sqrt(2);
  end;   
  U = PUPolyMatRand(M,1,0,'complex');                
  V = PUPolyMatRand(M,1,1,'complex');                
  A = PolyMatConv(U,PolyMatConv(Sigma,ParaHerm(V)));   

  A2 = zeros(M,M,9); A2(:,:,1:2:end)=A; 
  B = zeros(M,M,1); B(:,:,1) = eye(M);
  [Q,D,S] = PUProcrustes(A2,B,1024,0,16);

  save Eus26_1Fig5.mat A Sigma Q
else
  load Eus26_1Fig5.mat
end;

%------------------------------------------------------------------------------
%  Figure 5 --- compare real parts
%------------------------------------------------------------------------------
figure(5);
Rows = [1 12 30]; Columns = [9, 15, 26, 31];
AA = zeros(length(Rows),length(Columns),41); AA(:,:,11:2:19)=A(Rows,Columns,:); 
QQ = Q(Rows,Columns,503:543);
t = (-10:30);
for m = 1:length(Rows),
   for n = 1:length(Columns),
      subplot(length(Rows),length(Columns),length(Columns)*(m-1)+n); 
      stem(t,squeeze(real(AA(m,n,:))),'b'); hold on;
      plot(t,squeeze(real(QQ(m,n,:))),'r*'); 
      axis([-5 10 -.15 .2]); grid on;
      if n == 1,
        ylabel(sprintf('$\\ell=%d$',Rows(m)),'interpreter','latex'); 
      end;   
      if m==1,
         title(sprintf('$m=%d $',Columns(n)),'interpreter','latex');
      end;
      if m==length(Rows),
         xlabel('index $n$','interpreter','latex');
      end;
      set(gca,'TickLabelInterpreter','latex',...
      'XTick',(-5:5:10),'XTickLabel',{'$-5$','$0$','$5$','$10$'},...
      'YTick',(-.15:.15:.15),'YTickLabel',...
     {'$-.15$','$.0$','$0.15$'});
   end;
end;
set(gcf,'OuterPosition',[230 250 570 400]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Eus26_1Fig5.eps

%------------------------------------------------------------------------------
%  Figure 6 --- compare imaginary parts
%------------------------------------------------------------------------------
figure(20);
for m = 1:length(Rows),
   for n = 1:length(Columns),
      subplot(length(Rows),length(Columns),length(Columns)*(m-1)+n); 
      stem(t,squeeze(imag(AA(m,n,:))),'b'); hold on;
      plot(t,squeeze(imag(QQ(m,n,:))),'r*'); 
      axis([-5 10 -.15 .2]); grid on;
      if n == 1,
        ylabel(sprintf('$\\ell=%d$',Rows(m)),'interpreter','latex'); 
      end;   
      if m==1,
         title(sprintf('$m=%d $',Columns(n)),'interpreter','latex');
      end;
      if m==length(Rows),
         xlabel('index $n$','interpreter','latex');
      end;
      set(gca,'TickLabelInterpreter','latex',...
      'XTick',(-5:5:10),'XTickLabel',{'$-5$','$0$','$5$','$10$'},...
      'YTick',(-.1:.1:.1),'YTickLabel',...
     {'$-.15$','$.0$','$0.15$'});
   end;
end;
set(gcf,'OuterPosition',[230 250 570 400]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Eus26_1Fig6.eps


