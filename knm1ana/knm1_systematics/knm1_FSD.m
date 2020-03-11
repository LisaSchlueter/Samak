% check FSD variances for different shape uncertainties
% normalization uncertainty here fixed to 0 --> this is an additional effect

FSDNorm_RelErr = 0; % normalization is fixed -> only check for influence of shape
FSDShapeGS_RelErr = 0.04;
FSDShapeES_RelErr = 0.18;

% tritium model
CommongArg = {'RunList','KNM1','chi2','chi2Stat','DataType','Real',...
    'fixPar','5 6 7 8 9 10 11','exclDataStart',1,'chi2','chi2Stat',...
    'NonPoissonScaleFactor',1};
M = MultiRunAnalysis(CommongArg{:});

% covariance matrix
M.ComputeCM('SysEffect',struct('FSD','ON'),'BkgCM','OFF','nTrials',5000,'RecomputeFlag','OFF',...
    'FSDNorm_RelErr',FSDNorm_RelErr,...
    'FSDShapeGS_RelErr',FSDShapeGS_RelErr,...
    'FSDShapeES_RelErr',FSDShapeES_RelErr);
CM = M.FitCM_Obj; CM.RecomputeFlag = 'OFF';
CM.ComputeCM_FSD;
d = importdata(CM.CovMatFile);

%% check the numbers
Isotopologue = 'TT';
FSD_P = d.([Isotopologue,'_P_norm'])';            % exc. probabilities
FSD_E = d.obj.StudyObject.([Isotopologue,'exE']); % exc. energy

% rel. uncertainty on ground state variance
varianeGS = var(FSD_P(:,1:M.ModelObj.TTGSTh),0,2);  % absolute variance
varianeES = var(FSD_P(:,M.ModelObj.TTGSTh:end),0,2);

GSvarRelErr = std(varianeGS)./mean(varianeGS);
ESvarRelErr = std(varianeES)./mean(varianeES);

fprintf('%.2f %% rel. variance uncertainty ground  state \n',GSvarRelErr*100);
fprintf('%.2f %% rel. variance uncertainty excited state \n',ESvarRelErr*100);