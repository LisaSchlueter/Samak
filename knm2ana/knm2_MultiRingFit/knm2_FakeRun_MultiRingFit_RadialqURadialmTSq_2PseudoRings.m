% Script to develop and test multi ring fit
% based on KNM2 fake MC run
% October 2019, Lisa

% method: 
% test qU-Offset:
% Data: MC data with radial dependent qU-values. 2 rings with offset (defined in InitFile)
% Model: Test different cases:
%        1) Uniform fit with averaged qU's. Expect a neutrino mass shift
%        2) Multi ring fit with correct ringwise qU-values
%        3) Multi ring fit with uniform qU-values. Try to fit qU-Offset
%% init file
InitFile = @ref_FakeRun_KNM2_CD84_50days_RadialqU2Rings; %100 column density, 50 days
Display = 'OFF';

%% set up model, compute fake run if necessary
CommonArg = {'RunNr',1,...% has no meaning
    'DataType','Fake',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'exclDataStart',11,... % 11==40eV range (28 subruns)
    'chi2','chi2Stat',...
    'RingMerge','Half',...
    'minuitOpt','min;migrad',...
    'NonPoissonScaleFactor',1};

savedir = [getenv('SamakPath'),'knm2ana/knm2_MultiRingFit/plots/'];

%%  multi-ring fit, take qU-values from first ring
M = RunAnalysis(CommonArg{:},'FakeInitFile',InitFile,'AnaFlag','Ring','fixPar','mNu E0 Bkg Norm qU mTSq',...
    'PullFlag',4,'fitter','minuit');
M.ModelObj.ComputeTBDDS('mTSq_bias',[0, -0.2]);
M.ModelObj.ComputeTBDIS;
TBDIS = M.ModelObj.TBDIS;

qUmean = repmat(M.RunData.qU(:,1),[1,2]);%rrepmat(mean(R_MultiRing.RunData.qU,2),[1,4]);
M.SimulateRun('qU',qUmean);

M.RunData.TBDIS = TBDIS;
M.RunData.TBDISE = sqrt(TBDIS);
M.Fit;
if strcmp(Display,'ON')
    [plotHandle, cbHandle]  = PlotRingwiseqUOffset(M);
    savename2 = [savedir,'knm2_FakeRun_MultiRingFit_RadialqUmTSqfit_2Rings.pdf'];
    export_fig(plotHandle,savename2);
end
