% unblinded fit with penning trap background slope
% use saenz as MC default
range     = 40;
freePar   = 'mNu E0 Bkg Norm';

% labelling
savedir = [getenv('SamakPath'),'knm2ana/knm2_FSD/results/'];
savename = sprintf('%sknm2_DifferentFSD_%.0feV_%s_statNP_Uniform_MC2.mat',...
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
    %%
    A = MultiRunAnalysis(RunAnaArg{:},'FSDFlag','SAENZ','DopplerEffectFlag','FSD');
    A.exclDataStart = A.GetexclDataStart(range);
    A.InitModelObj_Norm_BKG;
    
    TBDIS_i = A.ModelObj.TBDIS;
    A.RunData.TBDIS = TBDIS_i;
    A.Fit;
    FitResult_i = A.FitResult;
    
    B = MultiRunAnalysis(RunAnaArg{:},'FSDFlag','SibilleFull','DopplerEffectFlag','OFF');
    B.exclDataStart = B.GetexclDataStart(range);
    B.RunData.TBDIS = TBDIS_i;
    B.Fit;
    FitResult_Knm1FSD = B.FitResult;
    
    C = MultiRunAnalysis(RunAnaArg{:},'FSDFlag','KNM2','DopplerEffectFlag','FSD');
    C.exclDataStart = C.GetexclDataStart(range);
    C.RunData.TBDIS = TBDIS_i;
    C.Fit;
    FitResult_Saenz2000 = C.FitResult;
    
    %%
    MakeDir(savedir);
    save(savename,'FitResult_i','FitResult_Knm1FSD','FitResult_Saenz2000',...
        'RunAnaArg','A','B','C','SigmaSq')
    
end
%
mNuSq = zeros(3,1);
mNuSq(1) = FitResult_i.par(1);
mNuSq(2) = FitResult_Knm1FSD.par(1);
mNuSq(3) = FitResult_Saenz2000.par(1);

    

