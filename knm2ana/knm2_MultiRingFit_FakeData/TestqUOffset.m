% test qUOffset
% 1. data:  No radial dependent qU
% 2. model: Radial dependent qU
% goal: get correct qU-Offset with fit
%% settings
% fake MC run
% measurement time 50 days
% 84 % column density, KNM2 MTD (one random run)
% 40 eV range, stat only
% FPD: 2 or 4 pseudo rings, outermost ring excluded
% 14th Oct 19, Lisa

%% set up model, compute fake run if necessary
% init file
TimeDays = 50;
InitFile = @ref_FakeRun_KNM2_CD84_50days; %84 column density, 50 days

CommonArg = {'RunNr',1,...% has no meaning
    'DataType','Fake',...
    'FSDFlag','Sibille0p5eV',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','Ring',...    % FPD combinations
    'RingMerge','Half',... % Full== combined 12 rings into 4 pseudo rings
    'exclDataStart',12,...  % 12==40eV range (28 subruns)
    'chi2','chi2Stat',...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    'fixPar','mNu E0 Bkg Norm c'}; % free Parameter

R = RunAnalysis(CommonArg{:},'FakeInitFile',InitFile);
qU_i = R.ModelObj.qU;
%% give model a qU-Offset
qUOffset = linspace(0,-0.1,R.nRings);%-0.1,0,0];
qU = qU_i + qUOffset;
R.ModelObj.qU = qU;
R.ModelObj.AdjustRF; % re-initialize model and response function with shifte qU-values   %alternative: SimulateRun('qU',qU) -> gives same result; 
%% fit
R.Fit;
%% result:
Q_i = 18573.7;
E0diff = R.FitResult.par(2)+R.ModelObj.Q_i+R.FitResult.par(2*R.ModelObj.nPixels+9:3*R.ModelObj.nPixels+8)-Q_i;
