function [ pmatch ] = pm(read, actual)
%%% Percent match of letters in read compared to actual. useful for finding
%%% matches not at the beginning of reads
ncorrect=0;
mcorrect=0;
for i=1:length(read)-length(actual)+1
    ncorrect=0;
    for j=1:length(actual)
        if read(i+j-1)==actual(j)
            ncorrect=ncorrect+1;
        end
    end
    
    if ncorrect>mcorrect
        mcorrect=ncorrect;
    end
    
    if mcorrect==length(actual)
    pmatch=100;
    return 
    end
    
end
pmatch=mcorrect/length(actual)*100.00;
end

