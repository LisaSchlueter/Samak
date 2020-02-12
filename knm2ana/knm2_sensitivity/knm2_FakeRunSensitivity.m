% Compute neutrino mass sensitivity before data taking with fake MC run
% measurement times: 50 days
% 84 % column density, KNM2 MTD (one random run)
% 40 eV range, stat only
% 2nd Nov 19, Lisa
% init files
InitFile = @ref_FakeRun_KNM2_CD84_50days; %100 column density, 50 days
fixPar = 'mNu E0 Bkg Norm';
AnaFlag = 'StackPixel';

% set up model, compute fake run if necessary
RunArg = {'RunNr',1,...% has no meaning
    'DataType','Fake',...
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
    RunArg = {RunArg,'RingMerge','Full'};
end
R = RunAnalysis(RunArg{:},'FakeInitFile',InitFile);

% label
savedir = [getenv('SamakPath'),'knm2ana/knm2_sensitivity/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('knm2_FakeRunSensitivity_%s%s_%s_%s_pull%.0f.mat',...
    AnaFlag,R.RingMerge,strrep(fixPar,' ',''),R.chi2,R.pullFlag)];

if exist(savename,'file')
    load(savename)
    R.FitResult = FitResults;
else
    % fit fake run (Asimov) -> sensitivity
    R.Fit;
    FitResults = R.FitResult;
    save(savename,'FitResults','RunArg','-mat');
end
fprintf('KNM2 estimated symmetric neutrino mass squared sensitivty: %.3f eV^2 \n',R.FitResult.err(1));



