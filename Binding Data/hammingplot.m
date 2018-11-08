clear
addpath ('./export')
warning('off','all')

nnkdist=[0,2,3,4,5,6,8,9,10,11,14,15,16,18,19,22,25,27,29,30,31];
seqcount=5000;
paratopesize=[15,14,12,17,12,11,12,10,8,9,14,9,8,9,11,11,10];
scaffoldsize=[31,47,55,44,45,48,52,40,48,48,49,43,37,33,44,41,43];

aa=0;
bb=0;

nbins=20;
histcenters=linspace(0.05,1.95,20);
d=diff(histcenters/2);
edges=linspace(0,2,21);


nnkpdist=nan(17,nbins);
initseqpdist=nan(17,3,nbins);
finalseqpdist=nan(17,7,nbins);
nnkinitdif=nan(17,1);
finalinitdif=nan(17,1);
initseqpdistscafavg=nan(17,nbins);
finalseqpdistscafavg=nan(17,nbins);

for scaffold=['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p' 'q']
    aa=aa+1;

    concerved=rand(scaffoldsize(aa)-paratopesize(aa),seqcount);
    concerved(concerved<(1/1000))=0;
    concerved(concerved>=(1/1000))=1;
    nnkseqt=[randi(31,seqcount,paratopesize(aa)) concerved'];
    for i=2:21
        nnkseqt(nnkseqt<=nnkdist(i) & nnkseqt>nnkdist(i-1))=nnkdist(i);
    end
    nnkpdist(aa,:)=histcounts(scaffoldsize(aa).*pdist(nnkseqt,'hamming')./paratopesize(aa),edges,'Normalization','probability');
   
        
    bb=0;  
    for zz=[9 10 11]
        bb=bb+1;
        try seqs = fastaread(['./allsorted/' num2str(zz) '-' scaffold '.txt']);
        catch
            continue
        end
        if length(seqs)<3
            continue
        end

        seqs_int=nan(length(seqs),length(seqs(1).Sequence));
        parfor i=1:length(seqs)
            seqs_int(i,:)=aa2int(seqs(i).Sequence)
        end

        initseqpdist(aa,bb,:)=histcounts(scaffoldsize(aa).*pdist(seqs_int,'hamming')./paratopesize(aa),edges,'Normalization','probability');
    end
    
    bb=0;
    bbc=false;
    for zz=[1 2 3 5 6 7 8]
        bb=bb+1;
        try seqs = fastaread(['./allsorted/' num2str(zz) '-' scaffold '.txt']);
        catch
            continue
        end
        if length(seqs)<3
            continue
        end
        bbc=true;
        seqs_int=nan(length(seqs),length(seqs(1).Sequence));
        parfor i=1:length(seqs)
            seqs_int(i,:)=aa2int(seqs(i).Sequence)
        end

        finalseqpdist(aa,bb,:)=histcounts(scaffoldsize(aa).*pdist(seqs_int,'hamming')./paratopesize(aa),edges,'Normalization','probability');
    end
    initseqpdistaverage=mean(initseqpdist(aa,:,:),2);
    nnkinitdif(aa)=kstest2(initseqpdistaverage(:),nnkpdist(aa,:),'Alpha',0.05/17,'Tail','smaller');
    initseqpdistscafavg(aa,:)=initseqpdistaverage(:);
    if bbc
        finalseqpdist2=[];
        for i=1:7
            if ~isnan(finalseqpdist(aa,i,1))
                finalseqpdist2=[finalseqpdist2; finalseqpdist(aa,i,:)];
            end
        end
        finalseqpdistaverage=mean(finalseqpdist2(:,:),1);
        finalinitdif(aa)=kstest2(finalseqpdistaverage(:),initseqpdistaverage(:),'Alpha',0.05/12,'Tail','smaller');
        finalseqpdistscafavg(aa,:)=finalseqpdistaverage(:);
    end
end


figure
hold on
errorbar(histcenters-0.01,mean(nnkpdist),std(nnkpdist),'o','DisplayName','NNK','MarkerSize',2);
errorbar(histcenters,mean(initseqpdistscafavg),std(initseqpdistscafavg),'*','DisplayName','Initial','MarkerSize',2);
errorbar(histcenters+0.01,nanmean(finalseqpdistscafavg,1),nanstd( finalseqpdistscafavg),'s','DisplayName','Binding','MarkerSize',2);
leg=legend('show','Location','Northeast');
leg.FontSize=6;
leg.Box='off';
ylabel('Frequency')
xlabel('Hamming Distance (Mismatched Residues/Paratope Size)')
hold off
h=gcf;
set(h,'Units','inches')
set(h,'Position',[0.1 0.1 2.5 2])
set(gca,'FontName','Arial','fontsize',6)
print(h,'C:\Users\alexg\Google Drive\Hackel Lab Share\Research Updates\A Golinski\ProtID Final\paper\PNAS\Matlab figures\5_a','-dpng','-r300')

