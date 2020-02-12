%% set up model, compute fake run if necessary
% init file
TimeDays = 50;
InitFile = @ref_FakeRun_KNM2_CD84_50days; %84 column density, 50 days

CommonArg = {'RunNr',1,...% has no meaning
    'DataType','Fake',...
    'FSDFlag','Sibille0p5eV',...
    'fixPar','11 12 13 14 15 16 21',...% fix FSD parameters and background slope
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...    % FPD combinations
    'RingMerge','Full',... % Full== combined 12 rings into 4 pseudo rings
    'exclDataStart',12,...  % 12==40eV range (28 subruns)
    'chi2','chi2Stat',...
    'fitter','minuit',...
    'minuitOpt','min;migrad'};

R = RunAnalysis(CommonArg{:},'FakeInitFile',InitFile);
% init values
R.ModelObj.ComputeTBDDS;
Te_i = R.ModelObj.Te;
TBDDS_i = R.ModelObj.TBDDS;
%% differential spectrum with larger endpoint
R.ModelObj.ComputeTBDDS('E0_bias',10);
TBDDS_shift = R.ModelObj.TBDDS;
%% Question: is shift of TBDDS equal to shift of endpoint
plot(Te_i,TBDDS_shift,'Color',rgb('IndianRed'));
hold on;
plot(Te_i+10,TBDDS_i,'-k');
hold off;
set(gca,'YScale','log');

