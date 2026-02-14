% Generate Fig. 7 for Eusipco'26 paper on analytic Procrustes solution
%

close all; clear all;

%-------------------------------------------------------------------------------
%  parameters
%-------------------------------------------------------------------------------
filename = 'EusipcoEnsembleResults.txt';

A = dlmread(filename);

Rold = A(find(A(:,1)==1),2:8);
Rnew = A(find(A(:,1)==2),2:8);
M = [2 4 8 16 32];
MedianOld = zeros(5,3);  MedianNew = zeros(5,3);
PercentOld = zeros(5,6); PercentNew = zeros(5,6);
for m = 1:length(M),
   dummy = sort(Rold(find(Rold(:,2)==M(m)),5:7));
   MedianOld(m,:) = dummy(50,:);   
   PercentOld(m,:) = [dummy(6,:) dummy(95,:)];
   dummy = sort(Rnew(find(Rnew(:,2)==M(m)),5:7));
   MedianNew(m,:) = dummy(50,:);   
   PercentNew(m,:) = [dummy(6,:) dummy(95,:)];
end;

% for legend only
leg1 = semilogy([-1 -2],[1 1],'b-'); hold on;
leg2 = semilogy([-1 -2],[1 1],'-','color',[0 0.5 0]);
leg3 = semilogy([-1 -2],[1 1],'r-');
leg4 = semilogy([-1 -2],[1 1],'ko-');
leg5 = semilogy([-1 -2],[1 1],'k*--');

% now for the actual curves
semilogy(log2(M),MedianOld(:,2),'r-');
errorbar(log2(M)+0.025,MedianOld(:,2),(MedianOld(:,2)-PercentOld(:,2)),...
        (PercentOld(:,5)-MedianOld(:,2)),'ro');
semilogy(log2(M),MedianOld(:,3),'-','color',[0 0.5 0]);
errorbar(log2(M)+0.025,MedianOld(:,3),(MedianOld(:,3)-PercentOld(:,3)),...
        (PercentOld(:,6)-MedianOld(:,3)),'o','color',[0 0.5 0]);
%
semilogy(log2(M),MedianNew(:,2),'r--');
errorbar(log2(M)-0.025,MedianNew(:,2),(MedianNew(:,2)-PercentNew(:,2)),...
        (PercentNew(:,5)-MedianNew(:,2)),'r*');
semilogy(log2(M),MedianNew(:,3),'--','color',[0 0.5 0]);
errorbar(log2(M)-0.025,MedianNew(:,3),(MedianNew(:,3)-PercentNew(:,3)),...
        (PercentNew(:,6)-MedianNew(:,3)),'*','color',[0 0.5 0]);
ylabel('LS error $e_{\mathrm{LS}}$, execution time $T$/[s]','interpreter','latex');
%legend({'LS error','exec. time','PU error','oversampled [xx]','proposed'},'interpreter',...
%    'latex','fontsize',10,...
%        'location','West');%%
legend([leg1 leg2 leg3],{'PU error $e_{\mathrm{PU}}$','exec. time $T$','LS error $e_{\mathrm{LS}}$'},...
    'interpreter','latex','fontsize',10,'location','West');%%
axis([0.9 5.1 0.1 100]); 
grid on;     
set(gca,'TickLabelInterpreter','latex','XTick',(1:5),'XTickLabel',...
       {'$2$','$4$','$8$','$16$','$32$'});  
%% right hand side of plot
%%
yyaxis right;
%semilogy(log2(M),MedianOld(:,1),'b-'); 
plot(log2(M),10*log10(MedianOld(:,1)),'b-'); 
%errorbar(log2(M)-.025,10*log10(MedianOld(:,1)),10*log10((MedianOld(:,1)-PercentOld(:,1))),...
%        10*log10((PercentOld(:,4)-MedianOld(:,1))),'bo');
errorbar(log2(M),10*log10(MedianOld(:,1)),10*log10(MedianOld(:,1))-10*log10(PercentOld(:,1)),...
        10*log10(PercentOld(:,4))-10*log10(MedianOld(:,1)),'bo');
%semilogy(log2(M),MedianNew(:,1),'b--');
plot(log2(M),10*log10(MedianNew(:,1)),'b--');
errorbar(log2(M),10*log10(MedianNew(:,1)),10*log10(MedianNew(:,1))-10*log10(PercentNew(:,1)),...
        10*log10(PercentNew(:,4))-10*log10(MedianNew(:,1)),'b*');       
xlabel('spatial dimension $M$','interpreter','latex');
ylabel('paraunitarity error $10\log10\{e_{\mathrm{PU}}\}$','interpreter','latex');
axis([0.9 5.1 -230 -90]); 

ah1 = axes('position',get(gca,'position'),'visible','off');
legend(ah1,[leg4 leg5],{,'oversampled [19]','proposed'},...
'interpreter','latex','fontsize',10, 'Location','SouthWest');

set(gcf,'OuterPosition',[230 250 570 350]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Eus26_1Fig7.eps


