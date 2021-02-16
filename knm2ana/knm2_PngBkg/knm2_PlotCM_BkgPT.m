% unblinded fit with penning track background slope
TestFit = 'ON';
BKG_PtSlopeErr = 3*1e-06;
BKG_PtSlope = 3*1e-06;
TwinBias_BKG_PtSlope = 3*1e-06;

%%
RunList = 'KNM2_Prompt';
range     = 40;
freePar   = 'mNu E0 Bkg Norm';
chi2      = 'chi2Stat';%CMShape';%CMShape';
DataType  = 'Twin';
AnaFlag   = 'StackPixel';
RingMerge = 'Full';%'None';
FSDFlag   = 'KNM2';
SigmaSq =  0.0124+0.0025;



RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar',freePar,...
    'RadiativeFlag','ON',...
    'minuitOpt','min ; minos',...
    'FSDFlag',FSDFlag,...
    'ELossFlag','KatrinT2A20',...
    'SysBudget',40,...
    'AnaFlag',AnaFlag,...
    'chi2',chi2,...
    'TwinBias_Q',18573.7,...
    'NonPoissonScaleFactor',1,...
    'FSD_Sigma',sqrt(SigmaSq),...
    'TwinBias_FSDSigma',sqrt(SigmaSq),...
    'RingMerge',RingMerge,...
    'PullFlag',99,...
    'BKG_PtSlope',BKG_PtSlope,...
    'TwinBias_BKG_PtSlope',TwinBias_BKG_PtSlope};%99 = no pull


if ischar(RunList)
        A = MultiRunAnalysis('RunList',RunList,RunAnaArg{:});
else
        A = RunAnalysis('RunNr',RunList,RunAnaArg{:});
end

A.exclDataStart = A.GetexclDataStart(range);
%%
if strcmp(TestFit,'ON')
    savedir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
    savename = sprintf('%sknm2_PlotTestBkgPTCM_%s%s_err%.1fmuCpsPerS.mat',savedir,DataType,RunList,1e6*BKG_PtSlopeErr);
    if exist(savename,'file')
        load(savename,'mNuSqErr','mNuSqErrCM');
        A.chi2 = 'chi2CMShape';
        A.ComputeCM('BkgPTCM','ON','RecomputeFlag','ON','BKG_PtSlopeErr',BKG_PtSlopeErr,'SysEffect',struct('FSD','OFF'),'BkgCM','OFF','nTrials',1e4);
    else
        A.chi2 = 'chi2Stat';
        A.Fit;
        mNuSqErr = 0.5*(A.FitResult.errPos(1)- A.FitResult.errNeg(1));
        %
        A.chi2 = 'chi2CMShape';
        A.ComputeCM('BkgPTCM','ON','RecomputeFlag','ON','BKG_PtSlopeErr',BKG_PtSlopeErr,'SysEffect',struct('FSD','OFF'),'BkgCM','OFF','nTrials',1e4);
        A.Fit;
        mNuSqErrCM = 0.5*(A.FitResult.errPos(1)- A.FitResult.errNeg(1));
        save(savename,'mNuSqErr','mNuSqErrCM','BKG_PtSlopeErr','RunAnaArg');
    end
    fprintf('slope uncertainty (%s, %s,%.1f mucps/s): 1 sigma sys err on m^2 = %.3f eV^2 \n',DataType,RunList,1e6.*BKG_PtSlopeErr,sqrt(mNuSqErrCM^2-mNuSqErr^2))
else
    A.chi2 = 'chi2CMShape';
    A.ComputeCM('BkgPTCM','ON','RecomputeFlag','ON','BKG_PtSlopeErr',BKG_PtSlopeErr,'SysEffect',struct('FSD','OFF'),'BkgCM','OFF','nTrials',1e4);
end
%%

[BkgCMPtShape,BkgCMPtFracShape] = A.FitCM_Obj.DecomposeCM('CovMatFrac',A.FitCM_Obj.CovMatFrac,...
    'exclDataStart',A.exclDataStart,'BkgCM','ON');
%%
A.FitCM_Obj.ComputeCM_BackgroundPT('Display','OFF')
%%
A.FitCM_Obj.SysEffect.BkgPT = 'ON';
A.FitCM_Obj.PlotCM('CovMatInput',A.FitCM_Obj.CovMat,'Mode','CM','qUWindowIndexMax',-135,'savename',sprintf('_Err_%.1fmuCpsPerS_CMtot',A.FitCM_Obj.BKG_PtSlopeErr*1e6),'saveplot','ON');
A.FitCM_Obj.PlotCM('CovMatInput',A.FitCM_Obj.CovMatFracShape,'Mode','Shape','qUWindowIndexMax',-135,'savename',sprintf('_Err_%.1fmuCpsPerS_CMFracShape',A.FitCM_Obj.BKG_PtSlopeErr*1e6),'saveplot','ON');
A.FitCM_Obj.PlotCorr('qUWindowIndexMax',-135,'savename',sprintf('_Err_%.1fmuCpsPerS',A.FitCM_Obj.BKG_PtSlopeErr*1e6),'qUWindowIndexMin',40,'saveplot','ON')
%%
close all
CutIn = 25;
imagesc(BkgCMPtFracShape(1:end-CutIn,1:end-CutIn));
colorbar
PrettyFigureFormat;
title('fractional shape only')