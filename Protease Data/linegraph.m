stab=[2.99 2.99 2.99];
med=[2.98 2.72 1.93];
low=[2.90 0.94 0.90];

figure
hold on
scatter(0:2,stab,'filled','o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0])
P1=-0.002;
a=plot((0:2),P1(1)*(0:2)+stab(1),'Color',[1 0 0],'DisplayName','Stable: m = 0.00','LineStyle','--');

scatter(0:2,med,'filled','o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
P2=-0.476;
b=plot((0:2),P2(1)*(0:2)+med(1),'Color',[0 0 1],'DisplayName','Semi-Stable: m = -0.48','LineStyle','--');

scatter(0:2,low,'filled','o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0])
P3=-1.194;
c=plot((0:2),P3(1)*(0:2)+low(1),'Color',[0 1 0],'DisplayName','Non-Stable: m = -1.19','LineStyle','--');

xlim([-0.5 2.5])
ylim([0 3.5])
xticks([0 1 2])
yticks([0 1 2 3])
xticklabels({'No Protease','Low Protease','High Protease'})
ylabel('Protease Resistance')
lg=legend([a b c],'Location','SouthOutside');
lg.Title.String=('Slope \propto Stability');
lg.FontSize=6;
lg.FontName='Arial';
h=gcf;
set(h,'Units','inches')
set(h,'Position',[0.1 0.1 2.5 1])
set(gca,'FontName','Arial','fontsize',6)
print(h,'C:\Users\alexg\Google Drive\Hackel Lab Share\Research Updates\A Golinski\ProtID Final\paper\PNAS\Matlab figures\7_a4n','-dpng','-r300')
hold off