% test qU-Offset with randomized MC data
% similar script to knm2_FakeRun_MultiRingFit_2PseudoRings , but with stat. fluct

%% set up model, compute fake run if necessary
%84 column density, 50 days, radial qU dependence, 2 rings
InitFile = @ref_FakeRun_KNM2_CD84_50days_RadialqU2Rings;
RunAnaArg = {'RunNr',1,...% has no meaning
    'DataType','Fake',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'exclDataStart',11,... % 11==40eV range (28 subruns)
    'chi2','chi2Stat',...
    'RingMerge','Half',...
    'minuitOpt','min;migrad',...
    'NonPoissonScaleFactor',1,...
    'fitter','minuit',...
    'pullFlag',4,...
    'FakeInitFile',InitFile,...
    'AnaFlag','Ring',...
    'fixPar','mNu E0 Bkg Norm qU mTSq'}; % free parameter!

nSamples = 1000;
%% Model: multi-ring fit, take qU-values from first ring
MR2 = RunAnalysis(RunAnaArg{:});
qUmean = repmat(MR2.RunData.qU(:,1),[1,2]);
MR2.SimulateRun('qU',qUmean);
TBDIS_i = MR2.RunData.TBDIS;

%% load if possible
if contains(RunAnaArg{end},'mTSq')
    extra_str = '_mTSq';
elseif MR2.pullFlag == 0
   extra_str = '_NoPull';
else
   extra_str = ''; 
end

savedir = [getenv('SamakPath'),'knm2ana/knm2_MultiRingFit/results/'];
savefile = [savedir,...
    sprintf('FitResults_knm2_FakeRun_MultiRingFit_2PseudoRings_StatFluct%s_%.0fSamples.mat',extra_str,nSamples)];
if exist(savefile,'file')
    d = importdata(savefile);
    FitResults = d.FitResults;
else
    %% Add stat fluctuations
    TBDIS = TBDIS_i+(randn([size(MR2.RunData.TBDIS),nSamples]).*MR2.RunData.TBDISE);
    FitResults = cell(nSamples,1);
    
    for i=1:nSamples
        progressbar(i/nSamples)
        MR2.RunData.TBDIS = TBDIS(:,:,i);
        MR2.Fit;
        FitResults{i} = MR2.FitResult;
    end
    
    DataAsimov = MR2.RunData;
    
    MakeDir(savedir);
    save(savefile,'TBDIS','FitResults','DataAsimov','RunAnaArg','InitFile','-mat');
end

%% plot and save histograms

SaveAs = strrep(strrep(savefile,'results','plots'),'.mat','');
PlotHist_MCrandomFit('FitResults',FitResults,'PlotPar','E0','SaveAs',SaveAs);
PlotHist_MCrandomFit('FitResults',FitResults,'PlotPar','mNu','SaveAs',SaveAs);
PlotHist_MCrandomFit('FitResults',FitResults,'PlotPar','qU','SaveAs',SaveAs,'Ring',2);
PlotHist_MCrandomFit('FitResults',FitResults,'PlotPar','mTSq','SaveAs',SaveAs,'Ring',2);

%% plot and save correlations
mypar  = cell2mat(cellfun(@(x) x.par,FitResults,'UniformOutput',false))';
CorrHistPlot('qU','N',mypar,'SaveAs',SaveAs,'Ring',2)
CorrHistPlot('E0','qU',mypar,'SaveAs',SaveAs,'Ring',2,'PlotCorMat','OFF');
CorrHistPlot('mNu','qU',mypar,'SaveAs',SaveAs,'Ring',2,'PlotCorMat','OFF');

