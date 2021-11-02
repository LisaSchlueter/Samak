% Fit asimov mc without doppler effect
range     = 40;
freePar   = 'mNu E0 Bkg Norm';

% labelling
savedir = [getenv('SamakPath'),'knm2ana/knm2_FSD/results/'];
savename = sprintf('%sknm2_DifferentFSD_noDoppler_%.0feV_%s_statNP_Uniform_MC.mat',...
    savedir,range,strrep(freePar,' ',''));

if exist(savename,'file')
    load(savename);
else
    TwinBias_BKG_PtSlope = 3*1e-06;
    SigmaSq =  0.0124+0.0025;
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType','Real',...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',1.112,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'RingMerge','None',...
        'PullFlag',99,...;
        'BKG_PtSlope',3*1e-06,...
        'TwinBias_BKG_PtSlope',3*1e-06,...
        'AngularTFFlag','ON'};
    
    A = MultiRunAnalysis(RunAnaArg{:},'FSDFlag','KNM2','DopplerEffectFlag','OFF');
    A.exclDataStart = A.GetexclDataStart(range);
    A.InitModelObj_Norm_BKG;
    
    TBDIS_i = A.ModelObj.TBDIS;
    A.RunData.TBDIS = TBDIS_i;
    A.Fit;
    FitResult_i = A.FitResult;
    
    C = MultiRunAnalysis(RunAnaArg{:},'FSDFlag','SAENZ','DopplerEffectFlag','OFF');
    C.exclDataStart = C.GetexclDataStart(range);
    C.RunData.TBDIS = TBDIS_i;
    C.Fit;
    FitResult_Saenz2000 = C.FitResult;
    
    MakeDir(savedir);
    save(savename,'FitResult_i','FitResult_Saenz2000',...
        'RunAnaArg','A','C','SigmaSq')
    
end
%
mNuSq = zeros(2,1);
mNuSq(1) = FitResult_i.par(1);
mNuSq(2) = FitResult_Saenz2000.par(1);

