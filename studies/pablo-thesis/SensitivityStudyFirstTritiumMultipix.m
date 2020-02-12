%clear;
% Screening corrections does not exist in new TBD (commented at the end, taken from TritiumSpectrumeV)
addpath(genpath('../../../Samak2.0'));

% Choose type of fit (chi2)
chi2name = 'chi2P';

% Choose data range to analyze
dataStart = 9; % 9 = -200 eV; 1 = -1600 eV


% Choose fitter
fitter = 'matlab';

% Number of Monte Carlo Simulations
Ntrials = [1000];

DataEffCorr = 'OFF';

%RunInfo = load('FAKEFTALLmpix');


PixelList = [1:112,114:124];

%PixelList = 1:10;

RunModel = ref_RunAnalysis('FAKEFTALLmpix','','',...
    'FPD_Segmentation','MULTIPIXEL','FPD_Pixel',PixelList,...
    'nTeBinningFactor',5,'HTFSD','OFF','TTFSD','OFF');

RunModel.ComputeTBDDS();
RunModel.ComputeTBDIS();

% Choose fixed parameters (optional)

%fixpar = num2str([1,2,3+length(PixelList):2+2*length(PixelList),2+2*length(PixelList)+1,2+2*length(PixelList)+2]);
%fixpar2 = num2str([3:3+length(PixelList)-1,2+2*length(PixelList)+1,2+2*length(PixelList)+2]);
fixpar = num2str([2+2*length(PixelList)+1 2+2*length(PixelList)+2]);

for trials = 1:length(Ntrials)
    N = Ntrials(trials);
    ResultsTableMpix = zeros(N,2+2+2+4*length(PixelList));
    for mc = 1:N
        qUdata = RunInfo.qU;
        TBDISdata = RunInfo.TBDIS{mc,1};
        Data = [qUdata,TBDISdata,sqrt(TBDISdata)];
        i_mnu = 0;
        i_Q = 0;
        
        F = FITC('SO',RunModel,'DATA',Data,'fitter',fitter,...
            'chi2name',chi2name,...
            'fixPar',fixpar,...
            'exclDataStart',dataStart,...
            'i_mnu',i_mnu,'i_Q',i_Q);
        
        mnu = F.RESULTS{1}(1)+RunModel.mnuSq_i;
        e0 = F.RESULTS{1}(2)+RunModel.Q_i;
        bck = F.RESULTS{1}(3:3+length(PixelList)-1)+RunModel.BKG_RateSec;
        bck_model = F.RESULTS{1}(3:3+length(PixelList)-1);
        norm = F.RESULTS{1}(3+length(PixelList):end-2)+1;
        norm_model = norm - 1;
%         
%         RunModel.ComputeTBDDS();
%         RunModel.ComputeTBDIS();
%         
%         F = FITC('SO',RunModel,'DATA',Data,'fitter',fitter,...
%             'chi2name',chi2name,...
%             'fixPar',fixpar2,...
%             'exclDataStart',dataStart,...
%             'i_mnu',i_mnu,'i_Q',i_Q,'i_B',bck_model,'i_N',norm_model);
        
        ResultsTableMpix(mc,:) = [mnu,e0,bck,norm,F.RESULTS{2}(1:end-2),F.RESULTS{3},F.RESULTS{5}];
       disp(mc) 
    end
    save(['FT_MC_MPIX','MC',num2str(mc),'.mat'],'ResultsTableMpix')
    
    
    
end
