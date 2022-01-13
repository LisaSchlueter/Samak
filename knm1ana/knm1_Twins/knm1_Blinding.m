% illustrate model blinding

savedir = [getenv('SamakPath'),'knm1ana/knm1_Twins/results/'];
savefile = sprintf('%sBlindFSD.mat',savedir);

if exist(savefile,'file')
    load(savefile)
else
    RunList = 'KNM1';
    range   = 40;         % 40eV range = 27 subruns
    ModelArg = {'RunList',RunList,...
        'chi2','chi2Stat',...
        'DataType','Real',...
        'fixPar','mNu E0 Norm Bkg',... free parameter
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',1.064,...
        'minuitOpt','min ; minos',...
        'ELossFlag','KatrinT2',...
        'SysBudget',22,...
        'AngularTFFlag','OFF'};
    
    % Init Model Object with original and blind fsd
    Mo = MultiRunAnalysis(ModelArg{:},'FSDFlag','Sibille');
    Mb =  MultiRunAnalysis(ModelArg{:},'FSDFlag','BlindingKNM1','DopplerEffectFlag','FSD');
    Mo.exclDataStart = Mo.GetexclDataStart(range);
    Mb.exclDataStart = Mb.GetexclDataStart(range);
    
    Mo.Fit; FitResult = Mo.FitResult;
    Mb.Fit; FitResult_blind = Mb.FitResult;
    %% get FSD (ground state only)
    exE = Mo.ModelObj.TTexE_G;
    exP = Mo.ModelObj.TTexP_G;
    Var_o = var(exE,exP);
    
    exE_blind = Mb.ModelObj.TTexE_G;
    exP_blind = Mb.ModelObj.TTexP_G;
    Var_blind = var(exE_blind,exP_blind);
    save(savefile,'FitResult','FitResult_blind',...
        'exE','exP','exE_blind','exP_blind','ModelArg');
   
    
end

%%


