% fit to ksn1 data
range = 95;
NonPoissonScaleFactor=1.064;
chi2 = 'chi2CMShape';
DataType = 'Real';
freePar = 'E0 Bkg Norm mnu4Sq sin2T4';
RunList = 'KNM1';
RunAnaArg = {'RunList',RunList,...
    'fixPar',freePar,...
    'DataType',DataType,...
    'FSDFlag','Sibille0p5eV',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'chi2',chi2,...
    'ROIFlag','Default',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','OFF',...
    'ISCSFlag','Edep',...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'TwinBias_Q',18573.73};
T = MultiRunAnalysis(RunAnaArg{:});
T.exclDataStart = T.GetexclDataStart(range);
T.Fit;
FitResult_freeSt = T.FitResult; 
%%
T.pullFlag = 11; %9
T.Fit;
FitResult_freeStpull11 = T.FitResult; 


%% test with sterile parameter combi, that has smaller chi2
sin2t4 = 0.0193;
m4=4.7149e+03; 
T.ModelObj.SetFitBiasSterile(m4,sin2t4);   % asign values to steriles
T.fixPar = 'E0 Bkg Norm';                  % fix sterile parameters
T.InitFitPar;
T.pullFlag = [];
T.Fit;
FitResult_fixSt = T.FitResult; 


