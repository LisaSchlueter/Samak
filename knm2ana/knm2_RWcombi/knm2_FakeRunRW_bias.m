% Assume data without plasma time evolution
% Model with plasma time evolution
% calculate neutrino mass squared bias as function of plasma time evolution

%% plasma settings
MultiPosAll  = [-1,0,1].*(0.05:0.025:0.3)';
MultiWeights = 1/3.*ones(3,1);

% label
savedir = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];
savename = [savedir,sprintf('knm2_FakeRunRW_bias_min%.0fmeV_max%.0fmeV_step%.0fmeV.mat',...
    1e3*min(MultiPosAll(:,3)),1e3*max(MultiPosAll(:,3)),1e3*mean(diff(MultiPosAll(:,3))))];

if exist(savename,'file')
else
%% model settings: Knm2-like, (time: 308 x 2 hours)
InitFileCombi =  @ref_FakeRun_KNM2_CD84_308x2hours; %84 column density, 308 runs a 2 hours
CommonArg = {'RunNr',1,...% has no meaning
    'DataType','Fake',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'exclDataStart',11,... % 11==40eV range (28 subruns)
    'chi2','chi2Stat',...
    'RingMerge','Full',...
    'minuitOpt','min;migrad',...
    'NonPoissonScaleFactor',1,...
    'AnaFlag','StackPixel',...
    'fixPar','mNu E0 Bkg Norm',...
    'RingList',1:12};

mNuSq = zeros(size(MultiPosAll,1)+1,1); %init
%% set up RunAnalysis object: MC data + model
D = RunAnalysis(CommonArg{:},'FakeInitFile',InitFileCombi);
TBDIS_Asimov = D.RunData.TBDIS;
D.Fit;
mNuSq(1) = D.FitResult.par(1);

%%
for i=1:size(MultiPosAll,1)
    progressbar(i/size(MultiPosAll,1));
    FSDArg = {'Dist','Gauss','Sigma',0.001,'MultiPos',MultiPosAll(i,:),'MultiWeights',MultiWeights};
    D.ModelObj.LoadFSD(FSDArg{:});
    D.Fit;
    mNuSq(i+1) = D.FitResult.par(1);
end
save(savename,'mNuSq','MultiPosAll','MultiWeights','TBDIS_Asimov');
end

exclIndex1 = mNuSq<=mNuSq(end);  % don't display fits, which didn't converge properly
exclIndex2 = mNuSq>=mNuSq(1);
exclIndex = logical(exclIndex1.*exclIndex2);

MultiPos = [0;MultiPosAll(:,3)];
x = MultiPos(exclIndex).*1e3;
y = (mNuSq(exclIndex)).*1e3;

f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
p1 = plot(x,y,'LineWidth',2);
PrettyFigureFormat('FontSize',24);
xlabel(sprintf('\\Delta plasma potentials (mV)'));
ylabel(sprintf('\\Delta m_\\nu^2 (meV^2)'))






