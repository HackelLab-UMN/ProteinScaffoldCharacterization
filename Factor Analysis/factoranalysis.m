%% Input data

factorin=dir('../factor/experiment/*.txt*');
factorname=[];
factor=[];
for i=1:length(factorin)
    factorname= [factorname string(factorin(i).name(1:end-4))];
    expfact=dlmread(['../factor/experiment/' factorin(i).name]);
    otherfact=dlmread(['../factor/other/' factorin(i).name]);
    factor=[factor [expfact;otherfact]];
end

texpname=textscan(fopen('../factor/order/exppdb_order.txt'),'%s','Delimiter','\n');
tothername=textscan(fopen('../factor/order/pdb_order.txt'),'%s','Delimiter','\n');

scafname=[texpname{1};tothername{1}];

torder=textscan(fopen('../matlabOutputs/factor_order.txt'),'%s','Delimiter','\n');
factorname=torder{1};
good=[];
for j=1:size(factor,1)
    if min(factor(j,:))~=-100
        good=[good j];
    end
end

scafname=scafname(good);
factor=factor(good,:);
for i=1:size(factor,2)
    factor(:,i)=(factor(:,i)-mean(factor(:,i)))./std(factor(:,i));
end
binding=dlmread('../matlabOutputs/score.txt');


%% Boxplot Section
% figure
% hold on
% boxplot(factor,'Symbol','','OutlierSize',4,'Colors','k')
% xticks([1:20])
% xticklabels(factorname)
% xtickangle(75)
% shift=linspace(-.2,.2,17);
% key='A':'Q';
% for j=1:17
%     scatter((1:20)+shift(j),factor(j,:),8,'filled','DisplayName',key(j))
% end
% leg=legend('show','Location','eastoutside');
% ylabel('Z-Score')
% title(leg,'Scaffold')
% ylim([-3 3])
% h=gcf;
% set(h,'Units','inches')
% set(h,'Position',[0.1 0.1 3.5 2.5])
% set(gca,'FontName','Arial','fontsize',6)
% print(h,'C:\Users\alexg\Google Drive\Hackel Lab Share\Research Updates\A Golinski\ProtID Final\paper\PNAS\Matlab figures\2_r','-dpng','-r300')
% hold off


%% PCA Section
[coeff,score,latent,tsquared,explained]=pca(factor,'NumComponents',6);
% score=factor;
% fig=figure;
% hold on
% bar(explained,'black')
% ylabel('Percent Variability Explained')
% xlabel('# of Components')
% yyaxis right
% plot(1:20,cumsum(explained),'r-o','MarkerFaceColor','r')
% ax=gca;
% ylim([0 100])
% ax.YColor=('red');
% ylabel('Cumulative Variability Explained')
% hold off
 
%% ICA Section
Mdl=rica(score,6);
score=transform(Mdl,score);
ica=coeff*Mdl.TransformWeights;
% ica=coeff;
% ica=ones(1,20);

% figure
% for k=1:6
%     subplot(6,1,k)
%     hold on
%     bar(1:20,ica(:,k),'k')
%     title(['Independent Component ' num2str(k)]);
%     if k<6
%         xticks([])
%     else
%         xticks(1:20)
%         xticklabels(factorname)
%         xtickangle(45)
%     end
%     set(gca,'FontName','Arial','fontsize',6)
% 
%     hold off
% end
% h=gcf;
% set(h,'Units','inches')
% set(h,'Position',[0.1 0.1 5 7])
% print(h,'C:\Users\alexg\Google Drive\Hackel Lab Share\Research Updates\A Golinski\ProtID Final\paper\PNAS\Matlab figures\ICA','-dpng','-r300')

%% regularizatoin 
expcomp=score(1:17,:);
best=ones(3,1); %rmse,df,alpha
for df=2 %1:size(expcomp,2)
for alpha=.01 %[.01 0.1 .25 .5 .75 1]
    [B,fitinfo]=lassoglm(expcomp,binding,'normal','CV',17,'Alpha',alpha,'DFmax',df);
    %lassoPlot(B,fitinfo,'plottype','CV'); 
    idxLambda1SE = fitinfo.Index1SE;
    min1coefs = find(B(:,idxLambda1SE));
    idxLambdaMinDeviance = fitinfo.IndexMinDeviance;
    mincoefs = find(B(:,idxLambdaMinDeviance));
    if ~isempty(mincoefs)
    lassoin=expcomp(:,mincoefs);
    index=1:17;
    msereal=0;
    rmse=10;
    llopredict=zeros(17,2); %(value,error)
    factorweights=zeros(17,size(lassoin,2)+1);
    for i=1:length(index)
        test = (index == i);
        train = ~test;
        p=regress(binding(train,:),[ones(sum(train==1),1) lassoin(train,:)]);
        ymodel=sum(p'.*[1 lassoin(test,:)]);
        msereal=msereal+(binding(test)-ymodel).^(2); %calulate error
        llopredict(i,1)=ymodel;
        llopredict(i,2)=delta;
        factorweights(i,:)=p;
    end    
    rmse=(msereal/17)^(1/2)
    if rmse<best(1)
        best(1)=rmse;
        best(2)=df;
        best(3)=alpha;
    end
    end
end
end
PFinal=regress(binding(:,:),[ones(length(binding),1) lassoin(:,:)]);
allpredict=sum(PFinal'.*[ones(length(score),1) score(:,mincoefs)],2);
% [lassoin llopredict binding]

%% ROC Curve 
% strgb2=[1,0,0,0,1,0,0,1,0,1,1,0,0,1,1,0,0];
% [X,Y,t,AUC]=perfcurve(strgb2,binding,1);
% figure
% hold on 
% plot(X,Y,'k-o','MarkerFaceColor','k')
% % text(X-0.05,Y-.025,num2str(t,2))
% xlim([-0.05 1.05])
% ylim([-0.05 1.05])
% xlabel('False Positive Rate')
% ylabel('True Positive Rate')
% hold off

%% Distribution of predicted performances
% figure
% hold on 
% histogram(allpredict,10,'Normalization','probability','FaceColor','black')
% xlabel('Predicted Scaffold Binding Performance')
% ylabel('Fraction of Scaffolds')
% hold off

%% example of ICA for figure 1
% figure
% hold on
% bar(20:-1:1,ica(:,1),'black');
% yticks([])
% xticks([])
% %     xlabel('Parameters')
% %     ylabel('Loadings')
% %     title('PC 1','fontsize',6)
% 
% h=gcf;
% set(h,'Units','inches')
% set(h,'Position',[0.1 0.1 1 1])
% set(gca,'FontName','Arial','fontsize',6)
% print(h,'C:\Users\alexg\Google Drive\Hackel Lab Share\Research Updates\A Golinski\ProtID Final\paper\PNAS\Matlab figures\1_a','-dpng','-r300')
% hold off


%% final model parameter coeff
figure
hold on
barh(21:-1:1,[sum(PFinal(2:end)'.*ica(:,mincoefs),2) ; PFinal(1)],'black');
% barh(2:-1:1,[sum(PFinal(2:end)'.*ica(:,mincoefs),2) ; PFinal(1)],'black');
%  yticks([1:2])
%  yticklabels(['Constant' ;flipud(factorname(mincoefs))])

yticks([1:21])
yticklabels(['Constant' ;flipud(factorname)])


xlabel('Predictive Model Coefficents')
h=gcf;
set(h,'Units','inches')
set(h,'Position',[0.1 0.1 2.5 2])
set(gca,'FontName','Arial','fontsize',6)
print(h,'C:\Users\alexg\Google Drive\Hackel Lab Share\Research Updates\A Golinski\ProtID Final\paper\PNAS\Matlab figures\4_a','-dpng','-r300')
hold off
mincoefs

%% LOO prediction figure 
datal={'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q'};
figure
hold on 
strgb2=[1,0,0,0,1,0,0,1,0,1,1,0,0,1,1,0,0];
% scatter(binding(strgb2==1),llopredict((strgb2==1),1),'o','filled','MarkerFaceColor','black')
% scatter(binding(strgb2==0),llopredict((strgb2==0),1),'x','MarkerFaceColor','black')
scatter(binding,llopredict(:,1),6,'ko','filled');
ylabel('LOO Model Prediction')
xlabel('Experimental Binding Performance')
% textlocy=llopredict(:,1);
% textlocy(4)=textlocy(4)+0.01;
% text(binding+0.005,textlocy,datal);
xlim([-.15 .25])
ylim([-.15 .25])
z=refline(0,-0.006);
z.LineStyle='--';
z.Color='black';
line([-0.006 -0.006],[-.15 .25],'Color','black','LineStyle','--')
% legend({"Evolved Strong Binding Variant","Not Functional"})
hold off
h=gcf;
set(h,'Units','inches')
set(h,'Position',[0.1 0.1 2.5 2])
set(gca,'FontName','Arial','fontsize',6)
print(h,'C:\Users\alexg\Google Drive\Hackel Lab Share\Research Updates\A Golinski\ProtID Final\paper\PNAS\Matlab figures\4_bsup','-dpng','-r300')

% 
% 
% 
% 
