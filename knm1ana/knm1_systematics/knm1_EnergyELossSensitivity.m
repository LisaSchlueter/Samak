% impact of energy loss function on KNM1 neutrino mass sensitivity
%% settings
RunList = 'KNM1';
ELossFlag ='KatrinT2'; 

% Init Model Object and covariance matrix object
if ~exist('A','var')
    A = MultiRunAnalysis('RunList',RunList,...
        'chi2','chi2Stat',...
        'FSDFlag','BlindingKNM1',...
        'ELossFlag',ELossFlag,...
        'DataType','Twin',...
        'fitter','matlab');
end

A.fixPar = '5 6 7 8 9 10 11'; %free fit parameter: neutrino mass, endpoint, normalization, background
A.exclDataStart = 14; % 17==30eVrange (24), 14==40eV range (27 subruns), 2==90eV,

% stat only
A.chi2 = 'chi2Stat';
A.Fit; A.Fit; % for matlab fitter: 1 init fit good to set up parallel computing
mNuSqStat = A.FitResult.err(1);

% stat + FSD sys
A.chi2 = 'chi2CMShape';
A.ComputeCM('SysEffects',struct('RF_EL','ON'),'BkgCM','OFF','nTrials',1000);
A.Fit;
mNuSqCM = A.FitResult.err(1);




