% unblinded fit with penning track background slope
range     = 40;
freePar   = 'mNu E0 Bkg Norm qU';
chi2      = 'chi2Stat';
DataType  = 'Real';
AnaFlag   = 'Ring';
RingMerge = 'Full';%'None';
DopplerEffectFlag = 'FSD';
BKG_PtSlope = 3*1e-06;
TwinBias_BKG_PtSlope = 3*1e-06;
FSDFlag   = 'KNM2';

if strcmp(AnaFlag,'Ring')
    SysBudget = 39;
    if strcmp(RingMerge,'Full')
        AnaStr = AnaFlag;
    else
        AnaStr = sprintf('Ring%s',RingMerge);
    end
else
    SysBudget = 40;
    AnaStr = AnaFlag;
end

savedir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savename = sprintf('%sknm2ubfinal_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s.mat',...
    savedir,BKG_PtSlope*1e6,DataType,range,strrep(freePar,' ',''),chi2,AnaStr,FSDFlag);

if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
elseif strcmp(chi2,'chi2Stat+')
    NonPoissonScaleFactor = 1.112;
    chi2 = 'chi2Stat';
end

if ~strcmp(chi2,'chi2Stat')
    savename = strrep(savename,'.mat',sprintf('_SysBudget%.0f.mat',SysBudget));
end

if strcmp(DataType,'Twin')
    savename = strrep(savename,'.mat',sprintf('_TwinBpng-%.1fmucpsPers.mat',1e6*TwinBias_BKG_PtSlope));
end

if exist(savename,'file')
    load(savename,'FitResult','RunAnaArg','A');
else
    SigmaSq =  0.0124+0.0025;
    
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag',FSDFlag,...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',SysBudget,...
        'AnaFlag',AnaFlag,...
        'chi2',chi2,...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'RingMerge',RingMerge,...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',BKG_PtSlope,...
        'TwinBias_BKG_PtSlope',TwinBias_BKG_PtSlope,...
        'DopplerEffectFlag',DopplerEffectFlag};
    A = MultiRunAnalysis(RunAnaArg{:});
    %%
    A.exclDataStart = A.GetexclDataStart(range);
    
    if strcmp(DataType,'Twin')
        A.ModelObj.RFBinStep = 0.01;
        A.ModelObj.InitializeRF;
    end
    
    A.Fit;
    FitResult = A.FitResult;
    MakeDir(savedir);
    save(savename,'FitResult','RunAnaArg','A','SigmaSq')
end
%%

%A.PlotFit;
fprintf('m_nu^2 = %.3f + %.3f %.3f eV^2       , ',FitResult.par(1),FitResult.errPos(1),FitResult.errNeg(1))
fprintf('mean err = %.3f eV^2 \n',(FitResult.errPos(1)-FitResult.errNeg(1))/2)
fprintf('E_0 = %.3f + %.3f eV  \n',FitResult.par(2)+A.ModelObj.Q_i,FitResult.err(2))
fprintf('chi2 = %.3f (%.0f dof), p = %.3f  \n',FitResult.chi2min,FitResult.dof,1-chi2cdf(FitResult.chi2min,FitResult.dof));

%%

Plot = 'OFF';
if strcmp(Plot,'ON')
    if strcmp(AnaFlag,'StackPixel')
        A.PlotFit('LabelFlag','FinalKNM1',...
            'saveplot','pdf',...
            'ErrorBarScaling',50,...
            'YLimRes',[-2.2,2.9],...
            'Colors','RGB',...
            'DisplayStyle','PRL',...
            'FitResultsFlag','OFF',...
            'qUDisp','Abs',...
            'TickDir','Out');
    else
        % A.PlotResidualsMultiRing('saveplot','ON','YLimRes',[-2.4 2.6])
        A.PlotFit('ring',4,...
            'saveplot','pdf',...
            'ErrorBarScaling',50,...
            'Colors','RGB',...
            'FitResultsFlag','OFF',...
            'qUDisp','Abs',...
            'TickDir','In',...
            'YLimRes',[-2.2,2.9],...
            'DisplayStyle','PRL');
        %A.PlotFitMultiRing('PlotPar','qU','savePlot','ON')
    end
end
%  A.FitCM_Obj.PlotCM('qUWindowIndexMax',10,'qUWindowIndexMin',40,'saveplot',...
%         'ON','Convergence','OFF','CovMatInput',A.FitCMFracShape,'PlotEffect','total',...
%         %         'savename','KNM2_UB1_MultiRingFull');
% A.FitCM_Obj.PlotCorr('qUWindowIndexMax',10,'qUWindowIndexMin',90,'saveplot',...
% 'ON','CovMatInput',A.FitCMFracShape,...
% 'savename',sprintf('KNM2_UB1_%s',AnaFlag));

%%