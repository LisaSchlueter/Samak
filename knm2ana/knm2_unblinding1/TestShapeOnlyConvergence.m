range   = 40;
freePar = 'mNu E0 Bkg Norm';
chi2    = 'chi2CMShape';
SysBudget = 38;
DataType = 'Twin';

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
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
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
    
    mNuSqErr2    = zeros(5,1);
    mNuSqErrSym2 = zeros(5,1);
    
    for i=1:5
    A.ComputeCM;
    A.Fit;
    FitResult = A.FitResult;
    mNuSqErr2(i)    = (FitResult.errPos(1)-FitResult.errNeg(1))/2;
    mNuSqErrSym2(i) = FitResult.err(1);
    end
%%

plot(mNuSqErr2-mean(mNuSqErr2))
%hold on;

%plot(mNuSqErrSym2-mean(mNuSqErrSym2),'--')
%%
%A.PlotFit;
fprintf('m_nu^2 = %.3f + %.3f %.3f eV^2       , ',FitResult.par(1),FitResult.errPos(1),FitResult.errNeg(1))
fprintf('mean err = %.3f eV^2 \n',(FitResult.errPos(1)-FitResult.errNeg(1))/2)