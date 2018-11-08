function [pmatch] = scafid(read,actual)
%%% A stricter version of pm.m, this goes letter by letter matching the
%%% ordered sequence to the read. It skips degenerate sites.
    ncorrect=0;    
    nposs=0;
    pmatch=0;
    if length(read)>=length(actual)
        for i=1:length(actual)
            if actual(i)~='X'
                nposs=nposs+1;
                if read(i)==actual(i)
                    ncorrect=ncorrect+1;
                end
            end
        end
        pmatch=(ncorrect/nposs)*100;
    end
end