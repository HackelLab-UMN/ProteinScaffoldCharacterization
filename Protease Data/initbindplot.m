cat=["A-Unsort" "A-Bind" "B-Unsort","B-Bind","E-Unsort","E-Bind","H-Unsort","H-Bind","J-Unsort","J-Bind","K-Unsort","K-Bind","N-Unsort","N-Bind","O-Unsort","O-Bind"];
cat2=["Naïve" "Binding"];
load('initbind.mat')
scaf=['A';'B';'E';'H';'J';'K';'N';'O'];
colormap=[ 0 .44 .74 ; .3 .75 .93 ; .49 .18 .56; ...
    .85 .33 .1 ; .64 .08 .18; .47 .67 .19; ...
    .93 .69 .13; 1 0 1];
figure
h=gcf;
set(h,'Units','inches')
set(h,'Position',[0.1 0.1 4.5 2.5])
set(gca,'FontName','Arial','fontsize',6)

j=1:2:16;
for i=1:8
    subplot(1,8,i)
    hold on
    temp=violinplot(initbind(:,j(i):j(i)+1),cat2);
    temp(1).ViolinColor=colormap(i,:);
    temp(2).ViolinColor=colormap(i,:);
    temp(2).MedianColor=[0 0 0];
    xtickangle(45)
    ylim([-1.4 0.5])
    p=ranksum(initbind(:,j(i)),initbind(:,j(i)+1),'tail','right');
    title([scaf(i)])
    set(gca,'FontName','Arial','fontsize',6)
    if i==1
        ylabel('Protease Stability')
    else 
        yticks([])
        ax=gca;
        ax.YAxis.Visible='off';
    end
    hold off
end

print(h,'C:\Users\alexg\Google Drive\Hackel Lab Share\Research Updates\A Golinski\ProtID Final\paper\PNAS\Matlab figures\7_c','-dpng','-r300')
hold off
