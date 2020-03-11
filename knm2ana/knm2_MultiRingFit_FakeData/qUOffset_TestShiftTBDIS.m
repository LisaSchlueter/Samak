% question: is a shift of the endpoint equivalent as shift of qU-vector in
% the integral spectrum?
% method: 1. compute integral spectrum 2. increase E0 and decrease qU by
% same amount. check difference, pretending that qU didn't change
% 10th Oct, 2019 -> result: it is exactly the same
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
R.ModelObj.ComputeTBDDS; R.ModelObj.ComputeTBDIS;
qU_i = R.ModelObj.qU;
TBDIS_i = R.ModelObj.TBDIS;
%%
qUShift = 10;
R.ModelObj.qU = qU_i - qUShift;
R.ModelObj.AdjustRF;
R.ModelObj.ComputeTBDDS('E0_bias',qUShift);
R.ModelObj.ComputeTBDIS;
TBDIS_shift = R.ModelObj.TBDIS;
%%
plot(qU_i,TBDIS_i./TBDIS_shift);
set(gca,'YScale','log');
xlabel('qU (eV)');
ylabel('ratio TBDIS');
