load parsco.mat
scaf=['A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'L';'M';'N';'O';'P';'Q'];
figure 
hold on
scatter(parsco(:,1),parsco(:,2),6,'filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])
parsco(14,2)=parsco(14,2)+0.01;
parsco(9,1)=parsco(6,1)-0.08;
text(parsco(:,1)+.06,parsco(:,2),scaf,'HorizontalAlignment','right','FontSize',6,'FontName','Arial')
% text(-.6,.23,"Spearman Correlation Coefficient: 0.56",'HorizontalAlignment','center')
% text(-.6,.215,'p<0.05','HorizontalAlignment','center')
xlabel('Parental Protease Stability')
ylabel('Binding Performance')

h=gcf;
set(h,'Units','inches')
set(h,'Position',[0.1 0.1 2.5 2.5])
set(gca,'FontName','Arial','fontsize',6)
print(h,'C:\Users\alexg\Google Drive\Hackel Lab Share\Research Updates\A Golinski\ProtID Final\paper\PNAS\Matlab figures\7_b','-dpng','-r300')
hold off