%%% Creates a bubble plot for each binding campagins where the area is
%%% proportional to the square root of the quartic dampened number of reads
%%% in a bin determined by matlab 

close
clear
figure %%% get the figure set up with colors and legend
hold on
col=[153/256 50/256 204/256; 0 1 1; 1 0 1; 0 0 1 ; 0 1 0; 1 0 0; 1 0.65 0];
camp=['1', '2', '3', '4', '5', '6', '7'];
adj=linspace(-.3,.3,7); 

 for j=1:7 %for each campagin 
     test=fastaread(['./allcorrect/' num2str(j) '.txt']);
     si=size(test,1);
     scat=NaN(17,si); %(scaffold, read index)
     
     for i=1:si %for each read
         x=test(i).Header;
         delim=strfind(x,';');
         scat(str2double(x(1:delim-1)),i)=str2double(x(delim+1:end));
     end
     
     scat(:,:)=scat(:,:)./nansum(nansum(scat(:,:))); %normalize to total number of reads per campagin
     for z=1:17 %for each scaffold
         if isnan(max(scat(z,:,1)))==0
            [a,b]=histcounts(((scat(z,:,1))),'BinMethod','auto'); %(number of reads, edges) auto bin each scaffold reads
            c=zeros(length(a),1); %center of bin
            for zz=1:length(b)-1
                 c(zz)=(b(zz)+b(zz+1))/2;
                 if a(zz)==0
                     a(zz)=1e-100; %log scale problems
                 end
             end

            y(j,z)=scatter(ones(1,length(a)).*z+adj(j),c,(a).^(1/2).*16,col(j,:),'filled','DisplayName',camp(j)); %(scaffold, center of bin, size of bubble, color, legend)
        end
     end
     

     


 end
 
 %%% The rest is adjusting the graph to get the lines to show
 
 nshow=17;
 set(gca,'yscale','log')
 xlim([0.5 (nshow+0.5)])

 xticks(1:nshow)
 set(gca,'XTickLabel',{'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q'})
 set(gca,'YTick',[],'ticklength',[0 0])
% set(gca,'fontsize',30)
 leg=legend(y(:,1),'Location','Eastoutside');
 leg.Title.String='Campaign';
 ylm=get(gca,'YLim');

 xlabel('Scaffold')
 ax1=get(gca,'Position');
%  ax1(2)=ax1(2)+0.005;
%  ax1(4)=ax1(4)-0.02;
 ax2 = axes('Position',ax1,'XAxisLocation','bottom','YaxisLocation','left','Color','none');
 set(gca,'yscale','log')
 ylim(ylm)
 %set(gca,'fontsize',30)
 set(gca,'XTick',linspace(0.5,(nshow+0.5),nshow+1),'XtickLabels',[]);
 xlim([0.5 nshow+0.5])
 ylabel({'Campaign Fraction of Sequences within Bin'})

 ax3 = axes('Position',ax1,'XAxisLocation','bottom','YaxisLocation','left','Color','none');
 set(gca,'XTick',linspace(0.5,(nshow+0.5),nshow+1),'XtickLabels',[]);
 xlim([0.5 nshow+0.5])
 set(gca,'yscale','log')
 ylim(ylm)
 set(gca,'YTick',[],'ticklength',[.7 1],'XColor',[.8 .8 .8])
     
 addpath ('./export')
 export_fig('./data/bubble','-tif')