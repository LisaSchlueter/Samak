%clear;
% Screening corrections does not exist in new TBD (commented at the end, taken from TritiumSpectrumeV)
addpath(genpath('../../../Samak2.0'));

% Choose type of fit (chi2)
chi2name = 'chi2CM';

% Choose data range to analyze
dataStart = 1; % 9 = -200 eV; 1 = -1600 eV


% Choose fitter
fitter = 'minuit';

% Number of Monte Carlo Simulations
Ntrials = [10000];

DataEffCorr = 'OFF';

RunInfo = load('FAKEFTALLmpix');

PixelList = 1;

RunList = 'StackCD100all';

% Choose fixed parameters (optional)

fixpar = num2str([2+2*length(PixelList)+1 2+2*length(PixelList)+2]);
    
RunModel = ref_RunAnalysis('FAKEFTALLex2','','',...
    'FPD_Segmentation','OFF',...
    'nTeBinningFactor',50,'HTFSD','OFF','TTFSD','OFF');

RunModel.ComputeTBDDS();
RunModel.ComputeTBDIS();

% Choose fixed parameters (optional)
for ring = 1:11
    MRA = MultiRunAnalysis('AnaFlag','StackPixel',...
    'chi2','chi2CM','RunList',RunList,...
    'fitter',fitter,'exclDataStart',dataStart,...
    'fixPar',fixpar,...
    'DataEffcorr','OFF'); %'ROI+PileUp

for trials = 1:length(Ntrials)
    N = Ntrials(trials);
    ResultsTableMpix = zeros(N,2+2+2+4*length(PixelList));
    for mc = 1:N
        qUdata = RunInfo.qU;
        TBDISdata = RunInfo.TBDIS{mc,1};
        Data = [qUdata,TBDISdata,sqrt(TBDISdata)];
        
        F = FITC('SO',RunModel,'DATA',Data,'fitter',fitter,...
            'chi2name',chi2name,...
            'fixPar',fixpar,...
            'exclDataStart',dataStart,...
            'COVMAT',MRA.FitCM);
        
        mnu = F.RESULTS{1}(1)+RunModel.mnuSq_i;
        e0 = F.RESULTS{1}(2)+RunModel.Q_i;
        bck = F.RESULTS{1}(3:3+length(PixelList)-1)+RunModel.BKG_RateSec;
        bck_model = F.RESULTS{1}(3:3+length(PixelList)-1);
        norm = F.RESULTS{1}(3+length(PixelList):end-2)+1;
        
        ResultsTableMpix(mc,:) = [mnu,e0,bck,norm,F.RESULTS{2}(1:end-2),F.RESULTS{3},F.RESULTS{5}];
       disp(mc) 
    end
    save(['FT_MC_OFF_SYS','MC',num2str(mc),'.mat'],'ResultsTableMpix')
    
    
    
end
end
