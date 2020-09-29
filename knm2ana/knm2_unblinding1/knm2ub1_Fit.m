range   = 40;
freePar = 'mNu E0 Bkg Norm';
chi2    = 'chi2CMShape';
SysBudget = 38;
DataType = 'Real';

savedir = [getenv('SamakPath'),'knm2ana/knm2_unblinding1/results/'];
savename = sprintf('%sknm2ub1_Fit_%s_%.0feV_%s_%s.mat',...
    savedir,DataType,range,strrep(freePar,' ',''),chi2);

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
        'AnaFlag','StackPixel',...
        'chi2',chi2,...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',1,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq)};
    A = MultiRunAnalysis(RunAnaArg{:});
    %%
    A.exclDataStart = A.GetexclDataStart(range);
   
    if ~strcmp(chi2,'chi2Stat')
        A.NonPoissonScaleFactor = 1.112;
        A.SetNPfactor;
    end
    
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
fprintf('m_nu^2 = %.3f + %.3f - %.3f eV^2 \n',FitResult.par(1),FitResult.errPos(1),FitResult.errNeg(1))