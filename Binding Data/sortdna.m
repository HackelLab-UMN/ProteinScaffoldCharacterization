%%% This script sorts through all of the fasta files and create files for
%%% each campaign and each scaffold within each campagin

%%% Output has Header of "scaffold;count" and Sequence in terms of AA

% The following lines are used to run pandaseq and turn fastq->a
%system('pandaseq -f ../fastq/A1f.fastq -r ../fastq/A1r.fastq -w fasta/0a.fasta -t 0.9 -p TTCCAGACTACGCTCTGCAG -q AGATAGAAATCAGTACTACCGGATCC')
%system('pandaseq -f ../fastq/B1f.fastq -r ../fastq/B1r.fastq -w ../fasta/4.fasta -t 0.9 -p TTCCAGACTACGCTCTGCAG -q AAGTCGATTTTGTTACATCTACACTGT')



clear
warning('off','all')
native=fastaread('native.txt');
for n=1:2
    tag{n}=native(n).Sequence; % 1 is GS linker, 2 is Au5
end


for j= [9 10 11] % for each campagin, j=4 is the gIgG sort, skip for non GP2
    for m=3:length(native)
        out{m-2,1}=native(m).Sequence; %out{:,1}=the ordered sequence
        lnative(m-2)=length(native(m).Sequence); %the number of AA in scaffold
    end
    %Sort through all reads of a campaign, and assign a scaffold
    
    counted=fastaread(['./fasta/' num2str(j) '.fasta']);
    out{1,length(counted)+1}=''; %preallocate 
    tic
    for i=1:length(counted) %for number of reads 
        x=nt2aa(counted(i).Sequence,'ACGTOnly', false, 'Frame', 1, 'AlternativeStartCodons',false);
        if mod(i,5000)==1 %for tracking progress
            disp((i/length(counted))) %percent done
            toc %time elapsed 
            disp(toc/(i/length(counted))-toc) %time left 
        end
        stop=strfind(x,'*');
            if isempty(stop)==0
               x=x(1:stop); %find truncated sequences
            end
        unresolved=strfind(x,'X');
        if length(x)>1 && isempty(unresolved)==1
            y=1; %turn on flag, 1=not yet assigned to scaffold
        end
        if pm(x,tag{1})>=70 && length(x)>=21 %percent match of GS linker
            x=x(20:end);
            if pm(x,tag{2})>=70 %percent match of Au5
                for jj=1:size(out,1)
                    if scafid(x,out{jj,1})==100 && y==1
                        if jj==12 || jj==15  
                            if min(x(22:23)=='MN') %to sort between L(contains MN) and O
                                out{12,i+1}=x(1:lnative(12));
                                y=0;
                            else
                                out{15,i+1}=x(1:lnative(15));
                                y=0;
                            end
                        elseif jj==6 && pm(x,'IAGRFEG')>=80 %to sort between F and I(contains IAGRFEG)
                            out{9,i+1}=x(1:lnative(9));
                            y=0;
                        else
                            out{jj,i+1}=x(1:lnative(jj));
                            y=0;
                        end
                    end   
                end
            end
        end
    end
    
    parfor z=1:size(out,1) %for each scaffold
        tic
        aa=char(out{z,2:end}); %convert cells to char strings 
        [uni, ab, ac]=unique(aa,'rows'); % uni = aa(ab,:) and aa = uni(ac,:).
        ad=hist(ac,length(ab)); %number of counts in bins (each sequence is a unique bin)
        [ae, af]=sort(ad,'descend'); %ae=ad(af), put the counts in order
        uni=uni(af,:); %order the reads in the order of the counts
        id='abcdefghijklmnopq';
        for zz=1:size(uni,1)
            if (contains(uni(zz,:),' ')==0 || contains(uni(zz,:),'*')==1) && isempty(uni(zz,:))==0 %remove all errors
                fastawrite(['./allcorrect/100/' num2str(j) '.txt'],[num2str(z) ';' num2str(ae(zz))] ,uni(zz,:)) 
                %fastawrite(['./allsorted/' num2str(j) '-' id(z) '.txt'],[num2str(z) ';' num2str(ae(zz))] ,uni(zz,:)) 

            end
        end
        toc
    end
    clearvars -except tag native j         
end
