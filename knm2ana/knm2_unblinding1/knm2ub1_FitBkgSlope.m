range   = 40;
freePar = 'mNu E0 Bkg Norm BkgSlope';
chi2    = 'chi2CMShape';%CMShape';
SysBudget = 38;
DataType = 'Real';
AnaFlag = 'StackPixel';%StackPixel';
PullFlag = 99; % 99 == No pull
BkgSlopeSigma = 4.74.*1e-06;

savedir = [getenv('SamakPath'),'knm2ana/knm2_unblinding1/results/BestFit/'];
savename = sprintf('%sknm2ub1_FitBkgSlope_%s_%.0feV_%s_%s_%s.mat',...
    savedir,DataType,range,strrep(freePar,' ',''),chi2,AnaFlag);
if PullFlag~=99
    savename = strrep(savename,'mat',sprintf('_BkgPull%.0f_%.2fmcpskeV.mat',PullFlag,BkgSlopeSigma*1e6));
end

if strcmp(chi2,'chi2Stat+')
    NonPoissonScaleFactor = 1.112;
    chi2 = 'chi2Stat';
else
    NonPoissonScaleFactor = 1;
end

if ~strcmp(chi2,'chi2Stat')
    savename = strrep(savename,'.mat',sprintf('_SysBudget%.0f.mat',SysBudget));
end

if exist(savename,'file')
    load(savename,'FitResult','RunAnaArg','A');
else
    SigmaSq =  0.0124+0.0025;
    
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'chi2','chi2Stat',...
        'DataType',DataType,...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',SysBudget,...
        'AnaFlag',AnaFlag,...
        'chi2','chi2Stat',...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'RingMerge','Full',...
        'PullFlag',PullFlag,...
        'pulls',BkgSlopeSigma};
    A = MultiRunAnalysis(RunAnaArg{:});
    %%
    A.exclDataStart = A.GetexclDataStart(range);
    
    if ~strcmp(chi2,'chi2Stat')
        A.NonPoissonScaleFactor = 1.112;
        A.SetNPfactor;
        
        A.chi2 = chi2;
        A.ComputeCM('BkgCM','OFF');
    end
    
    if strcmp(DataType,'Twin')
        A.ModelObj.RFBinStep = 0.01;
        A.ModelObj.InitializeRF;
    end
    
    A.Fit;
    FitResult = A.FitResult;
    MakeDir(savedir);
    save(savename,'FitResult','RunAnaArg','A')
end
%%
%A.PlotFit;
fprintf('m_nu^2 = %.3f + %.3f %.3f eV^2       , ',FitResult.par(1),FitResult.errPos(1),FitResult.errNeg(1))
fprintf('mean err = %.3f eV^2 \n',(FitResult.errPos(1)-FitResult.errNeg(1))/2)


fprintf('B slope = %.3f +- %.3f mcps/keV --> %.2f sigma \n ',FitResult.par(12)*1e6,FitResult.err(12)*1e6,FitResult.par(12)/FitResult.err(12))

%%
% A.PlotFit('LabelFlag','FinalKNM1',...
%     'saveplot','pdf',...
%     'ErrorBarScaling',50,...
%     'YLimRes',[-2.2,2.9],...
%     'Colors','RGB',...
%     'DisplayStyle','PRL',...
%     'FitResultsFlag','OFF',...
%     'qUDisp','Abs',...
%     'TickDir','Out');