% do a systematics breakdown with twin MC

DataType = 'Twin';
exclDataStart = 14;
nSamples = 1000;
savedir = [getenv('SamakPath'),'knm1ana/knm1sensitivity/results/'];
savefile = [savedir,'SystematicsBreakdown_MC'];

try
M = MultiRunAnalysis('RunList','KNM1','DataType',DataType,'fixPar','5 6 7 8 9 10 11',...
    'FSDFlag','Sibille0p5eV','exclDataStart',exclDataStart,'chi2','chi2Stat','fitter','minuit','minuitOpt','min;minos',...
    'TwinBias_mNuSq',-0.96);
M.NonPoissonScaleFactor = 1.0;
M.chi2 = 'chi2Stat';
[FitPar, FitErr, FitChi2min, dof]  = M.FitTwin('nSamples',nSamples);
thisfile = [savefile,'_StatNoNP.mat'];
save(thisfile,'FitPar','FitErr','FitChi2min','dof','M');
catch
end

try
M = MultiRunAnalysis('RunList','KNM1','DataType',DataType,'fixPar','5 6 7 8 9 10 11',...
    'FSDFlag','Sibille0p5eV','exclDataStart',exclDataStart,'chi2','chi2Stat','fitter','minuit','minuitOpt','min;minos',...
    'TwinBias_mNuSq',-0.96);
M.chi2 = 'chi2CMShape';
M.ComputeCM('SysEffect',struct('FSD','ON'),'BkgCM','OFF');
[FitPar, FitErr, FitChi2min, dof]  = M.FitTwin('nSamples',nSamples);
thisfile = [savefile,'_FSD.mat'];
FitCMShape = M.FitCMShape;
save(thisfile,'FitPar','FitErr','FitChi2min','dof','M','FitCMShape');
catch
end

try 
M = MultiRunAnalysis('RunList','KNM1','DataType',DataType,'fixPar','5 6 7 8 9 10 11',...
    'FSDFlag','Sibille0p5eV','exclDataStart',exclDataStart,'chi2','chi2Stat','fitter','minuit','minuitOpt','min;minos',...
    'TwinBias_mNuSq',-0.96);
M.chi2 = 'chi2CMShape';
M.ComputeCM('SysEffect',struct('RF_EL','ON'),'BkgCM','OFF');
[FitPar, FitErr, FitChi2min, dof]  = M.FitTwin('nSamples',nSamples);
thisfile = [savefile,'_EL.mat'];
FitCMShape = M.FitCMShape;
save(thisfile,'FitPar','FitErr','FitChi2min','dof','M','FitCMShape');
catch
end

try 
M = MultiRunAnalysis('RunList','KNM1','DataType',DataType,'fixPar','5 6 7 8 9 10 11',...
    'FSDFlag','Sibille0p5eV','exclDataStart',exclDataStart,'chi2','chi2Stat','fitter','minuit','minuitOpt','min;minos',...
    'TwinBias_mNuSq',-0.96);
M.chi2 = 'chi2CMShape';
M.ComputeCM('SysEffect',struct('RF_BF','ON'),'BkgCM','OFF');
[FitPar, FitErr, FitChi2min, dof]  = M.FitTwin('nSamples',nSamples);
thisfile = [savefile,'_BF.mat'];
FitCMShape = M.FitCMShape;
save(thisfile,'FitPar','FitErr','FitChi2min','dof','M','FitCMShape');
catch
end

try 
M = MultiRunAnalysis('RunList','KNM1','DataType',DataType,'fixPar','5 6 7 8 9 10 11',...
    'FSDFlag','Sibille0p5eV','exclDataStart',exclDataStart,'chi2','chi2Stat','fitter','minuit','minuitOpt','min;minos',...
    'TwinBias_mNuSq',-0.96);
M.chi2 = 'chi2CMShape';
M.ComputeCM('SysEffect',struct('RF_RX','ON'),'BkgCM','OFF');
[FitPar, FitErr, FitChi2min, dof]  = M.FitTwin('nSamples',nSamples);
thisfile = [savefile,'_RX.mat'];
FitCMShape = M.FitCMShape;
save(thisfile,'FitPar','FitErr','FitChi2min','dof','M','FitCMShape');
catch
end

try 
M = MultiRunAnalysis('RunList','KNM1','DataType',DataType,'fixPar','5 6 7 8 9 10 11',...
    'FSDFlag','Sibille0p5eV','exclDataStart',exclDataStart,'chi2','chi2Stat','fitter','minuit','minuitOpt','min;minos',...
    'TwinBias_mNuSq',-0.96);
M.chi2 = 'chi2CMShape';
M.ComputeCM('SysEffect',struct('FSD','OFF'),'BkgCM','ON');
[FitPar, FitErr, FitChi2min, dof]  = M.FitTwin('nSamples',nSamples);
thisfile = [savefile,'_Bkg.mat'];
FitCMShape = M.FitCMShape;
save(thisfile,'FitPar','FitErr','FitChi2min','dof','M','FitCMShape');
catch
end

try 

M = MultiRunAnalysis('RunList','KNM1','DataType',DataType,'fixPar','5 6 7 8 9 10 11',...
    'FSDFlag','Sibille0p5eV','exclDataStart',exclDataStart,'chi2','chi2Stat','fitter','minuit','minuitOpt','min;minos',...
    'TwinBias_mNuSq',-0.96);
M.chi2 = 'chi2CMShape';
M.ComputeCM('SysEffect',struct('Stack','ON'),'BkgCM','OFF');
[FitPar, FitErr, FitChi2min, dof]  = M.FitTwin('nSamples',nSamples);
thisfile = [savefile,'_Stack.mat'];
FitCMShape = M.FitCMShape;
save(thisfile,'FitPar','FitErr','FitChi2min','dof','M','FitCMShape');
catch
end

try 

M = MultiRunAnalysis('RunList','KNM1','DataType',DataType,'fixPar','5 6 7 8 9 10 11',...
    'FSDFlag','Sibille0p5eV','exclDataStart',exclDataStart,'chi2','chi2Stat','fitter','minuit','minuitOpt','min;minos',...
    'TwinBias_mNuSq',-0.96);
M.chi2 = 'chi2CMShape';
M.ComputeCM('SysEffect',struct('FPDeff','ON'),'BkgCM','OFF');
[FitPar, FitErr, FitChi2min, dof]  = M.FitTwin('nSamples',nSamples);
thisfile = [savefile,'_FPD.mat'];
FitCMShape = M.FitCMShape;
save(thisfile,'FitPar','FitErr','FitChi2min','dof','M','FitCMShape');
catch
end

try 

M = MultiRunAnalysis('RunList','KNM1','DataType',DataType,'fixPar','5 6 7 8 9 10 11',...
    'FSDFlag','Sibille0p5eV','exclDataStart',exclDataStart,'chi2','chi2Stat','fitter','minuit','minuitOpt','min;minos',...
    'TwinBias_mNuSq',-0.96);
M.chi2 = 'chi2CMShape';
M.ComputeCM('SysEffect',struct('TC','ON'),'BkgCM','OFF');
[FitPar, FitErr, FitChi2min, dof]  = M.FitTwin('nSamples',nSamples);
thisfile = [savefile,'_TC.mat'];
FitCMShape = M.FitCMShape;
save(thisfile,'FitPar','FitErr','FitChi2min','dof','M','FitCMShape');
catch
end




