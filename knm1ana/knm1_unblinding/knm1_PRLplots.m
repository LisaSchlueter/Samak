D = MultiRunAnalysis('RunList','KNM1','FSDFlag','Sibille0p5eV','fitter',...
'minuit','minuitOpt','min;migrad','exclDataStart',14,'chi2','chi2CMShape');
D.Fit;
%%
D.PlotFit('FitResultsFlag','OFF','yLimRes',[-3,3],'saveplot','pdf','ErrorBarScaling',100);%,'LabelFlag','FinalKNM1');