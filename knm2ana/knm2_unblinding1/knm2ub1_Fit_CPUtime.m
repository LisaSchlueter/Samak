% calculate CPU time for knm-2 fit
range   = 40;
freePar = 'mNu E0 Bkg Norm';
chi2    = 'chi2CMShape';%CMShape';
DataType = 'Real';%Real';
AnaFlag = 'StackPixel';%Ring';
RingMerge = 'Full';%'None';

if strcmp(AnaFlag,'Ring')
    SysBudget = 39;
    if strcmp(RingMerge,'Full')
        AnaStr = AnaFlag;
    else
        AnaStr = sprintf('Ring%s',RingMerge);
    end
else
    SysBudget = 38;
    AnaStr = AnaFlag;
end
savedir = [getenv('SamakPath'),'knm2ana/knm2_unblinding1/results/CPUtime/'];
savename = sprintf('%sknm2ub1_Fit_CPUtime_%s_%.0feV_%s_%s_%s.mat',...
    savedir,DataType,range,strrep(freePar,' ',''),chi2,AnaStr);


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

if exist(savename,'file') 
    load(savename,'tStart','tStop','CPUtimeHours');
else
    tStart = cputime;
    SigmaSq =  0.0124+0.0025;
    
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',SysBudget,...
        'AnaFlag',AnaFlag,...
        'chi2',chi2,...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'RingMerge',RingMerge};
    A = MultiRunAnalysis(RunAnaArg{:});
    %%
    A.exclDataStart = A.GetexclDataStart(range);

    if strcmp(DataType,'Twin')
        A.ModelObj.RFBinStep = 0.01;
        A.ModelObj.InitializeRF;
    end
    
    A.Fit;
    FitResult = A.FitResult;
    tStop = cputime;
    CPUtimeHours = (tStop-tStart)./(60*60);
    MakeDir(savedir);
    save(savename,'tStart','tStop','CPUtimeHours','FitResult','RunAnaArg','A','SigmaSq')
end
%%
%A.PlotFit;
fprintf('CPU time  = %.3f hours \n',CPUtimeHours)