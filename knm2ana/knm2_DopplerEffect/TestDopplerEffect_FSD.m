% Test doppler effect- implemented over broadening of FSDs
% Asimov Data with doppler effect
% Fit model without doppler effect
% settings
RecomputeFlag = 'OFF';
RunList = 'KNM2_RW1';
range = 40;
RunAnaArg = {'RunList',RunList,...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'DataType','Twin',...
    'fixPar','mNu E0 Bkg Norm',...
    'TwinBias_Q',18573.70};

savedir = [getenv('SamakPath'),'knm2ana/knm2_DopplerEffect/results'];
MakeDir(savedir);
savename = [savedir,sprintf('knm2_DopplerEffect_%s_%.0feV.mat',RunList,range)];

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
    fprintf('Load file %s \n',savename);
else
    % read data and set up model
    M = MultiRunAnalysis(RunAnaArg{:});
    M.exclDataStart = M.GetexclDataStart(range);
    M.InitModelObj_Norm_BKG('RecomputeFlag','OFF'); % get correct background and normalization
    M.ModelObj.BKG_RateSec_i = M.ModelObj.BKG_RateSec;
    
    % calculate Asimov data with doppler effect
    SigmaDE = M.ModelObj.DE_sigma;
    M.ModelObj.LoadFSD('Sigma',SigmaDE); % broadened FSDs
    M.ModelObj.ComputeTBDDS;
    M.ModelObj.ComputeTBDIS;
    AsimovDataDE = M.ModelObj.TBDIS; % use Asimov data to extract Doppler Effect: Data
    qU = M.ModelObj.qU; % save for plot
    
    % Model without Doppler
    M.ModelObj.LoadFSD;
    M.ModelObj.ComputeTBDDS;
    M.ModelObj.ComputeTBDIS;
    AsimovData =  M.ModelObj.TBDIS;
    
    % fit using model without DE
    M.RunData.TBDIS = AsimovDataDE;
    M.Fit;
    
    FitResult = M.FitResult;
    
    save(savename,'FitResult','SigmaDE','AsimovDataDE','AsimovData','RunAnaArg','qU');
    fprintf('save to %s \n',savename);
end

%% result
fprintf(' ------------------------- Test Doppler Effect (DE)  ------------------------- \n')
fprintf('Asimov Data with DE (sigma = %.4f eV), Fit model w/o DE: \n',SigmaDE)
fprintf('Neutrino mass squared = %.4f eV^2 \n',FitResult.par(1));
fprintf(' -----------------------------------------------------------------------------\n')

