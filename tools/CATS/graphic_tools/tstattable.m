function tTable = tstattable(s)
%% T-statistics

tTable = dataset({s.tstat.xhat,'Coef'},{s.tstat.se,'StdErr'},...
    {s.tstat.t,'tStat'},{s.tstat.pval,'pVal'});

end
