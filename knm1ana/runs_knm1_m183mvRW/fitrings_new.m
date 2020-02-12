 M = MultiRunAnalysis('RunList','KNM1_m183mvRW');
 R = RingAnalysis('RunAnaObj',M,'RingList',1:12);
 %%
 R.FitRings;
 R.PlotFits;