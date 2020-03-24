savedir = [getenv('SamakPath'),'knm2ana/knm2_LargeRangeFit/results/'];
savename = sprintf('%sknm2_ApplyLargeRangeFit.mat',savedir);
RecomputeFlag = 'OFF';

if exist(savename,'file')
    load(savename)
else
    E0twins = knm2FS_GetE0Twins('SanityPlot','OFF');
    range = 40;
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'fixPar','mNu E0 Bkg Norm',...         % free Parameter !!
        'DataType','Twin',...
        'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2','chi2Stat',...                    % statistics only
        'ROIFlag','14keV',...
        'SysBudget',33,...
        'TwinBias_Q',E0twins,...
        'NonPoissonScaleFactor',1};
    
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);
    A.ModelObj.RFBinStep = 0.02;
    A.ModelObj.InitializeRF;
    %% broaden FSD --> this time with std from large range fit!
    [E0, E0Err] = knm2_LargeRangeFit();
    Sigma = std(E0);
    
    FSDArg = {'SanityPlot','ON','Sigma',Sigma};
    A.ModelObj.LoadFSD(FSDArg{:});
    A.Fit;
    FitResultStat = A.FitResult;
    
    %% get covariance matrix for uncertainty
    A.ComputeCM('SysEffects',struct('FSD','OFF'),'BkgCM','OFF'); % get covariance matrix object
    A.FitCM_Obj.E0Offsets    = E0';
    A.FitCM_Obj.E0OffsetsErr = E0Err';
    A.FitCM_Obj.nTrials = 1000;
    A.FitCM_Obj.ComputeCM_PlasmaOffsets;
    
    % Compute Statistical Uncertainties,including P/NP fluctuations
    [StatCM, StatCMFrac] = A.ComputeCM_StatPNP;
    [A.FitCMShape,A.FitCMFracShape] = A.FitCM_Obj.DecomposeCM(...
        'CovMatFrac',A.FitCM_Obj.CovMatFracShape,'exclDataStart',A.exclDataStart);
    
    A.FitCMShape = A.FitCMShape + StatCM;
    A.chi2 = 'chi2CMShape';
    A.Fit;
    FitResultCM = A.FitResult;
    
    CMObj = A.FitCM_Obj;
    save(savename,'E0','E0Err','E0twins','RunAnaArg','FitResultStat','FitResultCM','CMObj')
end
%%
fprintf('--------------------------------------\n')
fprintf('mNuSq = %.3f (%.4f +%.4f) eV^2  (stat) \n',FitResultStat.par(1),FitResultStat.errNeg(1),FitResultStat.errPos(1))
fprintf('mNuSq = %.3f (%.4f +%.4f) eV^2  (stat + syst) \n',FitResultCM.par(1),FitResultCM.errNeg(1),FitResultCM.errPos(1))
fprintf('syst err = %.3f eV^2',sqrt(FitResultCM.err(1)^2-FitResultStat.err(1)^2));
fprintf('--------------------------------------\n')
%%
d = importdata(CMObj.CovMatFile);
f1 = figure('Units','normalized','Position',[0.1,0.1,0.45,0.5]);
h1 = histogram(d.FSDSigma_v);
h1.FaceColor = rgb('SkyBlue');
h1.FaceAlpha = 1;
PrettyFigureFormat('FontSize',24);
xlabel(sprintf('\\sigma({\\itE}_0^{fit}) (eV)'));
ylabel('Occurrence');
leg = legend(sprintf('mean = %.0f meV \\sigma = %.0f meV',...
    1e3*mean(d.FSDSigma_v),1e3*std(d.FSDSigma_v)),'EdgeColor',rgb('Silver'));
plotname = strrep(strrep(savename,'results','plots'),'.mat','.pdf');
export_fig(f1,plotname);