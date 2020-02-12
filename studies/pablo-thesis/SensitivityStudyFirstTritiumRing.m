clear;
% Screening corrections does not exist in new TBD (commented at the end, taken from TritiumSpectrumeV)
addpath(genpath('../../../Samak2.0'));

% Apply any pulls you wish (recommended for FT: 2 on mnu, none on the others)
pulls    = [Inf,Inf,Inf,Inf,Inf,Inf];

% Choose type of fit (chi2)
chi2name = 'chi2Stat';

% Choose data range to analyze
dataStart = 9; % 9 = -200 eV; 1 = -1600 eV

% Choose fixed parameters (optional)
fixPar   = '5 6';

% Choose fitter
fitter = 'minuit';

% Number of Monte Carlo Simulations
Ntrials = [10000];

DataEffCorr = 'OFF';

Nring = 11;

RunList = 'StackCD100all';


RunInfo = load('FAKEFTALLring');



for ring = 1:Nring
    
        MRA = MultiRunAnalysis('AnaFlag','Ring',...
    'chi2','chi2Stat','RunList',RunList,...
    'fitter',fitter,'exclDataStart',dataStart,...
    'fixPar',fixPar,...
    'DataEffcorr','OFF', 'ring',ring); %'ROI+PileUp
    
    RunModel = ref_RunAnalysis('FAKEFTALLring','','',...
        'FPD_Segmentation','RING','FPD_Ring',ring,...
        'nTeBinningFactor',20,'recomputeRF','ON','useParallelRF','OFF');
    RunModel.ComputeTBDDS();
    RunModel.ComputeTBDIS();
    
    for trials = 1:length(Ntrials)
        N = Ntrials(trials);
        ResultsTableRing = zeros(N,Nring,10);
        for mc = 1:N
            %qUdata = mean(RunInfo.qU(:,RunModel.ring{ring}),2);
            qUdata = RunInfo.qU;
            TBDISdata = RunInfo.TBDIS{mc,ring};
            Data = [qUdata,TBDISdata,sqrt(TBDISdata)];
            
            F = FITC('SO',RunModel,'DATA',Data,'fitter',fitter,...
                'chi2name',chi2name,...
                'pulls',pulls,...
                'fixPar','5 6',...
                'exclDataStart',dataStart,'COVMAT',MRA.FitCM);
            
            mnu = F.RESULTS{1}(1)+RunModel.mnuSq_i;
            e0 = F.RESULTS{1}(2)+RunModel.Q_i;
            bck = F.RESULTS{1}(3)+RunModel.BKG_RateAllFPDSec;
            norm = F.RESULTS{1}(4)+1;
            ResultsTableRing(mc,ring,:) = [mnu,e0,bck,norm,F.RESULTS{2}(1:4),F.RESULTS{3},F.RESULTS{5}];
            
            fid = fopen(['countingDir/f' num2str(mc)],'w');
            fclose(fid);
            progress = length(dir('countingDir'))-2;
            disp(progress)
            
        end
        save(['FT_MC_RING',num2string(ring),'MC',num2str(N),'.mat'],'ResultsTableRing')
    end
    
    
    
end
