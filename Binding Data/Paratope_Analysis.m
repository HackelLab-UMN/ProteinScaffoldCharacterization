clear
nativ=fastaread('native.txt');
nativ=nativ(3:end);
totalreads=0;
uniquereads=0;
par=cell(17,1);

for i=1:size(nativ,1)
    x=strfind(nativ(i).Sequence,'X'); %Find locations of all varied sites
    par{i}=x;  % Residue Number
end
nnk=[2,1,1,3,2,2,2,1,1,1,1,3,2,1,1,1,1,1,1,3]./31;
aaa='PMILVAGFWYCSTNQDEHKR';
binding=zeros(20,7);
for j=1:7
     campaign=zeros(20,1);
     test=fastaread(['./allcorrect/' num2str(j) '.txt']);
     si=size(test,1);
     
     for i=1:si
         x=test(i).Header;
         xx=test(i).Sequence;
         delim=strfind(x,';');
         scaf=str2double(x(1:delim-1));
         freq=str2double(x(delim+1:end));
         totalreads=totalreads+freq;
         uniquereads=uniquereads+1;
         nsites=size(par{scaf},2);
         
         for m=1:nsites
             flag=0;
             n=1;
            while flag==0
                if xx(par{scaf}(m))==aaa(n)
                    campaign(n)=campaign(n)+freq^(1/4);
                    flag=1;
                end
                n=n+1;
                if n>21
                    print('AHHHHH');
                end
            end        
         end
     end
     binding(:,j)=campaign(:)./sum(campaign(:));
end

starting=zeros(20,3);
for j=8:10
     campaign=zeros(20,1);
     test=fastaread(['./allcorrect/' num2str(j) '.txt']);
     si=size(test,1);
     
     for i=1:si
         x=test(i).Header;
         xx=test(i).Sequence;
         delim=strfind(x,';');
         scaf=str2double(x(1:delim-1));
         freq=str2double(x(delim+1:end));
%          totalreads=totalreads+freq;
         nsites=size(par{scaf},2);
         
         for m=1:nsites
             flag=0;
             n=1;
            while flag==0
                if xx(par{scaf}(m))==aaa(n)
                    campaign(n)=campaign(n)+freq^(1/4);
                    flag=1;
                end
                n=n+1;
            end        
         end
     end
     starting(:,j-7)=campaign(:)./sum(campaign(:));
end
     
% figure
% hold on
% scatter((1:20)-0.1,nnk,'k','filled')
% errorbar((1:20)-0.1,mean(starting'),std(starting'),'LineStyle','none','Marker','none','Color','k');
% markers=['o','s','d','p','h','o','s'];
% for i=1:7
%     scatter((1:20)+0.1,binding(:,i),'filled','Marker',markers(i));
% end
% xticklabels(aaa(:));
% xticks(1:20);
% legend({'NNK','Initial','Luciferase', 'CTLA/Avidin', 'PD1/Avidin', 'CTLA/Tris/Carb', 'GFP', 'PE', 'VEGF'})
% xlim([0,21]);
% ylim([0,.4]);
% ylabel('Abundance')
% hold off
% addpath ('./export')
% export_fig('./data/paratopediversity','-tif')

difference=mean(binding')-mean(starting');
figure
hold on
for i=1:20
    y=bar(i,difference(i),'FaceColor','flat');
    if difference(i)>0
       y.FaceColor=[0 0 1];
    else
       y.FaceColor=[1 0 0];
    end
end
% title('Binding Paratope Composition')
xticklabels(aaa(:));
xticks(1:20);
xlim([0,21])
ylabel('Paratope Abundance (Binding-Initial)')
xlabel('Amino Acid')
h=gcf;
set(h,'Units','inches')
set(h,'Position',[0.1 0.1 2.5 2])
set(gca,'FontName','Arial','fontsize',6)
print(h,'C:\Users\alexg\Google Drive\Hackel Lab Share\Research Updates\A Golinski\ProtID Final\paper\PNAS\Matlab figures\5_b','-dpng','-r300')
        
   

