% test qU-Offset with randomized MC data
% similar script to knm2_FakeRun_MultiRingFit_2PseudoRings , but with stat. fluct

%% set up model, compute fake run if necessary
%84 column density, 50 days, radial qU dependence, 4 rings
InitFile = @ref_FakeRun_KNM2_CD84_50days_radialqU;
RunAnaArg = {'RunNr',1,...% has no meaning
    'DataType','Fake',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'exclDataStart',11,... % 11==40eV range (28 subruns)
    'chi2','chi2Stat',...
    'RingMerge','Full',...
    'minuitOpt','min;migrad',...
    'NonPoissonScaleFactor',1,...
    'fitter','minuit',...
    'pullFlag',4,...
    'FakeInitFile',InitFile,...
    'AnaFlag','Ring',...
    'fixPar','mNu E0 Bkg Norm qU'}; % free parameter!

nSamples = 1000;
% Model: multi-ring fit, take qU-values from first ring
MR2 = RunAnalysis(RunAnaArg{:});
qUmean = repmat(MR2.RunData.qU(:,1),[1,4]);
MR2.SimulateRun('qU',qUmean);
TBDIS_i = MR2.RunData.TBDIS;
DataAsimov = MR2.RunData;

% load if possible
if contains(RunAnaArg{end},'mTSq')
    extra_str = '_mTSq';
    
    MR2.ModelObj.ComputeTBDDS('mTSq_bias',[0, 0.1 0.2 0.3]);
    MR2.ModelObj.ComputeTBDIS;
    TBDIS_i = MR2.ModelObj.TBDIS;
    DataAsimov = MR2.RunData;
    
elseif MR2.pullFlag == 0
    extra_str = '_NoPull';
else
    extra_str = '';
end

savedir = [getenv('SamakPath'),'knm2ana/knm2_RingFit/results/'];
savefile = [savedir,...
    sprintf('FitResults_knm2_FakeRun_MultiRingFit_4PseudoRings_StatFluct%s_%.0fSamples.mat',extra_str,nSamples)];
if exist(savefile,'file')
    d = importdata(savefile);
    FitResults = d.FitResults;
else
    
    % Add stat fluctuations
    TBDIS = TBDIS_i+(randn([size(MR2.RunData.TBDIS),nSamples]).*MR2.RunData.TBDISE);
    FitResults = cell(nSamples,1);
    
    for i=1:nSamples
        progressbar(i/nSamples)
        MR2.RunData.TBDIS = TBDIS(:,:,i);
        try 
            MR2.Fit;
            FitResults{i} = MR2.FitResult;
        catch
            fprintf('Fit failed for some reason \n')
        end
        if i==50  % save every 50 fits
            save(savefile,'TBDIS','FitResults','DataAsimov','RunAnaArg','InitFile','-mat');
        elseif mod(i,50)==0
            save(savefile,'TBDIS','FitResults','DataAsimov','RunAnaArg','InitFile','-mat','-append');
        end
    end
    
    MakeDir(savedir);
    save(savefile,'TBDIS','FitResults','DataAsimov','RunAnaArg','InitFile','-mat','-append');
end

%% plot and save histograms
SaveAs = strrep(strrep(savefile,'results','plots'),'.mat','');
PlotHist_MCrandomFit('FitResults',FitResults,'PlotPar','E0','SaveAs',SaveAs);
PlotHist_MCrandomFit('FitResults',FitResults,'PlotPar','mNu','SaveAs',SaveAs);
PlotHist_MCrandomFit('FitResults',FitResults,'PlotPar','qU','SaveAs',SaveAs,'Ring',2);
PlotHist_MCrandomFit('FitResults',FitResults,'PlotPar','qU','SaveAs',SaveAs,'Ring',3);
PlotHist_MCrandomFit('FitResults',FitResults,'PlotPar','qU','SaveAs',SaveAs,'Ring',4);

%% plot Correlations
% if get rid of fits that failed
FailIndex = cell2mat(cellfun(@(x) isempty(x),FitResults,'UniformOutput',false))';
FitResults = FitResults(~FailIndex);
mypar  = cell2mat(cellfun(@(x) x.par,FitResults,'UniformOutput',false))';

CorrHistPlot('qU','N',mypar,'SaveAs',SaveAs,'Ring',2)
CorrHistPlot('E0','qU',mypar,'SaveAs',SaveAs,'Ring',2,'PlotCorMat','OFF');
CorrHistPlot('E0','qU',mypar,'SaveAs',SaveAs,'Ring',3,'PlotCorMat','OFF');
CorrHistPlot('E0','qU',mypar,'SaveAs',SaveAs,'Ring',4,'PlotCorMat','OFF');
CorrHistPlot('qU','mNu',mypar,'SaveAs',SaveAs,'Ring',2,'PlotCorMat','OFF');