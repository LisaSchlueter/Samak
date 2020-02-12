% neutrino mass sensitivity before hard RW fix using MC twins
% 40 eV range
% 2nd Nov 19, Lisa

% settings
RunList = 'KNM2_beforeFix'; 
fixPar = 'mNu E0 Bkg Norm';
AnaFlag = 'Ring';
RecomputeFlag = 'OFF';

% set up model, compute fake run if necessary
RunArg = {'RunList',RunList,...% has no meaning
    'DataType','Real',...
    'FSDFlag','BlindingKNM2',...
    'fixPar',fixPar,...
    'ELossFlag','KatrinT2',...
    'AnaFlag',AnaFlag,...
    'exclDataStart',12,... % 40eV range (28 subruns)
    'chi2','chi2Stat',...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    'pullFlag',0};

if strcmp(AnaFlag,'Ring')
    RunArg = {RunArg{:},'RingMerge','Full'};
end
R = MultiRunAnalysis(RunArg{:});

% label
savedir = [getenv('SamakPath'),'knm2ana/knm2_sensitivity/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('knm2_Sensitivity_%sAsimov_%s%s_%s_%s_pull%.0f.mat',...
    RunList,AnaFlag,R.RingMerge,strrep(fixPar,' ',''),R.chi2,R.pullFlag)];

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
    R.FitResult = FitResults;
else
    R.InitModelObj_Norm_BKG;
    TBDIS_Asimov = R.ModelObj.TBDIS;
    R.RunData.TBDIS = TBDIS_Asimov;
    R.RunData.TBDISE = sqrt(TBDIS_Asimov);
    
    % fit fake run (Asimov) -> sensitivity
    R.Fit;
    FitResults = R.FitResult;
    save(savename,'FitResults','RunArg','-mat');
end
fprintf('KNM2 before fix (symmetric) neutrino mass squared sensitivty: %.3f eV^2 \n',R.FitResult.err(1));

