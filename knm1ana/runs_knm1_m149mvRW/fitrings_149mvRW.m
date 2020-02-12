M149 = MultiRunAnalysis('RunList','KNM1_m149mvRW');
M149.exclDataStart=2;
R149 = RingAnalysis('RunAnaObj',M149,'RingList',1:12);
R149.FitRings('SaveResult','ON','RecomputeFlag','OFF');
R149.PlotFits('SavePlot','ON','PlotPar',2);

%%
A149 = MultiRunAnalysis('RunList','KNM1_m149mvRW','RingList',1:10,'AnaFlag','Ring');
A149.FitAllRings;
A149.PlotFit;

%% 