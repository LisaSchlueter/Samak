% KNM2 Monte Carlo comparision
% Calculate Differential spectrum

FSDFlag = 'FSD'; % NoFSD
RecomputeFlag = 'ON';
SanityPlot = 'ON';

% label
savedir = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('knm2_ComputeTBDDS_%s.mat',FSDFlag)];

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
else
    RunAnaArg = {'RunNr',56341,...         %
        'fixPar','mNu E0 Bkg Norm',...     % free Parameter !!
        'DataType','Twin',...              % Real, Twin or Fake
        'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2','chi2Stat',...              % statistics only
        'NonPoissonScaleFactor',1,...      % Non-Poiss. background contribution (1 == no NP background)
        'TwinBias_Q',18573.70,...          % MC endpoint (eV)
        'DopplerEffectFlag','OFF',...
        'exclDataStart',1};
    
    R = RunAnalysis(RunAnaArg{:}); % run analysis object
    TBD = R.ModelObj; % tritium beta-decay object
    %%
    TBD.TD = 'DScomparison';
    TBD.SetTBDDSBinning;
    TBD.SetKinVariables;
    TBD.RadiativeFlag = 'OFF'; % set radiative corrections off
    
    if strcmp(FSDFlag,'NoFSD') % no final states
        TBD.TTFSD = 'OFF';
        TBD.HTFSD = 'OFF';
        TBD.DTFSD = 'OFF';
        TBD.ComputePhaseSpace;
    end
    
    TBD.NormFactorTBDDS = 1;
    TBD.ComputeTBDDS; % calculate differential spectrum
    
    TBDDS = TBD.TBDDS;
    Te    = TBD.Te;
    
    % save to txt file for Kafit
    txtStr = [savedir,'knm2_SamakDiffSpec_',FSDFlag];
    Write2Txt('filename',txtStr,...
        'variable',[Te';TBDDS'],...
        'variableName','E Rate',...
        'nCol',2);
    
    % save to h5 file
    h5Str = strrep(savename,'mat','h5');
    if strcmp(RecomputeFlag,'ON')
        system(['rm ',h5Str]);
    end
    
    h5create(h5Str,'/DiffSpec/Rate',[numel(Te),1]);
    h5create(h5Str,'/DiffSpec/Energy',[numel(Te),1]);
    h5write(h5Str,'/DiffSpec/Rate',TBDDS);
    h5write(h5Str,'/DiffSpec/Energy',Te);
    
    % save for Samak
    save(savename,'Te','TBDDS','TBD','FSDFlag');
end

%% sanity plot (optional)
if strcmp(SanityPlot,'ON')
    TBD.PlotTBDDS('type','log','saveplot','OFF');
end