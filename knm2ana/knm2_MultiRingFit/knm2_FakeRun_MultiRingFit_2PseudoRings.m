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
Display = 'ON';

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
%% Model 1: uniform fit: average (known!) qU-values
U   = RunAnalysis(CommonArg{:},'FakeInitFile',InitFile,'AnaFlag','StackPixel','fixPar','mNu E0 Bkg Norm','RingList',1:12);
U.Fit;
if strcmp(Display,'ON')
    [plotHandle, cbHandle]  = PlotRingwiseqUOffset(U);
    savename1 = [savedir,'knm2_FakeRun_UniformFit_RadialqU_2Rings.pdf'];
    export_fig(plotHandle,savename1);
end

%% Model 2: multi-ring fit, taken proper (known!) qU-values into account
MR1 = RunAnalysis(CommonArg{:},'FakeInitFile',InitFile,'AnaFlag','Ring','fixPar','mNu E0 Bkg Norm qU');
MR1.Fit;

if strcmp(Display,'ON')
    [plotHandle, cbHandle]  = PlotRingwiseqUOffset(MR1);
    savename2 = [savedir,'knm2_FakeRun_MultiRingFit_RadialqU_2Rings.pdf'];
    export_fig(plotHandle,savename2);
end
%% Model3: multi-ring fit, take qU-values from first ring
MR2 = RunAnalysis(CommonArg{:},'FakeInitFile',InitFile,'AnaFlag','Ring','fixPar','mNu E0 Bkg Norm qU',...
    'PullFlag',4,'fitter','minuit');
qUmean = repmat(MR2.RunData.qU(:,1),[1,2]);%rrepmat(mean(R_MultiRing.RunData.qU,2),[1,4]);
MR2.SimulateRun('qU',qUmean);
MR2.Fit;
if strcmp(Display,'ON')
    [plotHandle, cbHandle]  = PlotRingwiseqUOffset(MR2);
    savename2 = [savedir,'knm2_FakeRun_MultiRingFit_RadialqUfit_2Rings.pdf'];
    export_fig(plotHandle,savename2);
end

% 
% qUOffset = [-0.1:0.02:0.1];
% chi2min = zeros(numel(qUOffset),1);
% FitResults = cell(numel(qUOffset),1);
% 
% for i=1:numel(qUOffset)
%     progressbar(i/numel(qUOffset));
%     MR2.ModelObj.i_qUOffset(2) = qUOffset(i);
%     MR2.ModelObj.SetFitBias(0);
%     fprintf('Before fit: %.3g \n',MR2.ModelObj.qUOffset);
%     MR2.Fit;
%     fprintf('After fit: %.3g \n',MR2.ModelObj.qUOffset);
%     chi2min(i) = MR2.FitResult.chi2min;
%     FitResults{i} = MR2.FitResult;
%     % re-init
%     MR2.ModelObj.i_qUOffset(2) = 0; MR2.ModelObj.SetFitBias(0);
% end

% %%
% plot(qUOffset,chi2min)