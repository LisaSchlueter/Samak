RunList = 'FTpaper';
FT = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','DataEffCor','RunSummary',...
    'ELossFlag','Abdurashitov','exclDataStart',12,'SysBudget',0,'DataType','Real','FSDFlag','SAENZ',...
    'fixPar','1 5 6 7 8 9 10 11','NonPoissonScaleFactor',1);
%%
FT.chi2 = 'chi2CMShape';
FT.ComputeCM;
%%
FT.qUScan('firstPoint',7,'lastPoint',16,'RecomputeFlag','OFF','RelFlag','ON');
