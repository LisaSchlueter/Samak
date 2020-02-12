M4 = MultiRunAnalysis('RunList','KNM1','RingMerge','Full');
M10 = MultiRunAnalysis('RunList','KNM1','RingMerge','Default');
M12 = MultiRunAnalysis('RunList','KNM1','RingMerge','None');
MC4 = MultiRunAnalysis('RunList','KNM1','RingMerge','Full','DataType','Twin');
MC10 = MultiRunAnalysis('RunList','KNM1','RingMerge','Default','DataType','Twin');
MC12 = MultiRunAnalysis('RunList','KNM1','RingMerge','None','DataType','Twin');
MC4.GetPlotColor
MC10.GetPlotColor
MC12.GetPlotColor
%%
R4 = RingAnalysis('RunAnaObj',M4,'RingList',1:4);
R4.FitRings('SaveResult','OFF','RecomputeFlag','ON');

R10 = RingAnalysis('RunAnaObj',M10,'RingList',1:10);
R10.FitRings('SaveResult','OFF','RecomputeFlag','ON');

R12 = RingAnalysis('RunAnaObj',M12,'RingList',1:12);
R12.FitRings('SaveResult','OFF','RecomputeFlag','ON');

RC4 = RingAnalysis('RunAnaObj',MC4,'RingList',1:4);
RC4.FitRings('SaveResult','OFF','RecomputeFlag','ON');

RC10 = RingAnalysis('RunAnaObj',MC10,'RingList',1:10);
RC10.FitRings('SaveResult','OFF','RecomputeFlag','ON');

RC12 = RingAnalysis('RunAnaObj',MC12,'RingList',1:12);
RC12.FitRings('SaveResult','OFF','RecomputeFlag','ON');
%%
R4.PlotFits('SavePlot','ON','PlotPar',3,'Blind','OFF');
RC4.PlotFits('SavePlot','ON','PlotPar',3,'Blind','OFF');
R12.PlotFits('SavePlot','ON','PlotPar',3,'Blind','OFF');
R10.PlotFits('SavePlot','ON','PlotPar',3,'Blind','OFF');
RC10.PlotFits('SavePlot','ON','PlotPar',3,'Blind','OFF');
RC12.PlotFits('SavePlot','ON','PlotPar',3,'Blind','OFF');
R10.PlotFits('SavePlot','ON','PlotPar',2,'Blind','OFF');
R4.PlotFits('SavePlot','ON','PlotPar',2,'Blind','OFF');
R12.PlotFits('SavePlot','ON','PlotPar',2,'Blind','OFF');
RC10.PlotFits('SavePlot','ON','PlotPar',2,'Blind','OFF');
RC4.PlotFits('SavePlot','ON','PlotPar',2,'Blind','OFF');
RC12.PlotFits('SavePlot','ON','PlotPar',2,'Blind','OFF');