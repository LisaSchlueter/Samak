% Script to develop and test multi ring fit
% based on KNM2 fake MC run
% October 2019, Lisa

% method: 
% test qU-Offset:
% Data: MC data with radial dependent qU-values. 4 rings with offset: 0eV, 0.2eV, 0.4eV, 0.6 eV
% Model: Test different cases:
%        1) Uniform fit with averaged qU's. Expect a neutrino mass shift
%        2) Multi ring fit with 4 rings ->
% init file
InitFile = @ref_FakeRun_KNM2_CD84_50days_RadialqU; %100 column density, 50 days
Display = 'ON';

% set up model, compute fake run if necessary
CommonArg = {'RunNr',1,...% has no meaning
    'DataType','Fake',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'exclDataStart',11,... % 11==40eV range (28 subruns)
    'chi2','chi2Stat',...
    'RingMerge','Full',...
    'minuitOpt','min;migrad',...
    'NonPoissonScaleFactor',1};

savedir = [getenv('SamakPath'),'knm2ana/knm2_MultiRingFit/plots/'];

%% Model 1: uniform fit: average (known!) qU-values
R_Uniform   = RunAnalysis(CommonArg{:},'FakeInitFile',InitFile,'AnaFlag','StackPixel','fixPar','mNu E0 Bkg Norm','RingList',1:12);
R_Uniform.Fit;
if strcmp(Display,'ON')
    [plotHandle, cbHandle]  = PlotRingwiseqUOffset(R_Uniform);
    savename1 = [savedir,'knm2_FakeRun_UniformFit_RadialqU_4Rings.pdf'];
    export_fig(plotHandle,savename1);
end
%% multi-ring fit, taken proper (known!) qU-values into account
R_MultiRing = RunAnalysis(CommonArg{:},'FakeInitFile',InitFile,'AnaFlag','Ring','fixPar','mNu E0 Bkg Norm');
R_MultiRing.Fit;
if strcmp(Display,'ON')
    [plotHandle, cbHandle]  = PlotRingwiseqUOffset(R_MultiRing);
    savename2 = [savedir,'knm2_FakeRun_MultiRingFit_RadialqU_4Rings.pdf'];
    export_fig(plotHandle,savename2);
end
%% multi-ring fit 2
R_MultiRing2 = RunAnalysis(CommonArg{:},'FakeInitFile',InitFile,'AnaFlag','Ring','fixPar','mNu E0 Bkg Norm qU',...
    'PullFlag',4);
qUring1 = repmat(R_MultiRing2.RunData.qU(:,1),[1,4]);%repmat(mean(R_MultiRing.RunData.qU,2),[1,4]);
R_MultiRing2.SimulateRun('qU',qUring1);
R_MultiRing2.Fit;
if strcmp(Display,'ON')
    [plotHandle2, cbHandle2]  = PlotRingWiseFitPar(R_MultiRing2,'PlotPar','qU');%,'PlotParRef',[0,0.05,0.1,0.15]);
    %cbHandle2.Limits = cbHandle.Limits;
    savename2 = [savedir,'knm2_FakeRun_MultiRingFit_RadialqUfit_4Rings.pdf'];
    export_fig(plotHandle2,savename2);
end
