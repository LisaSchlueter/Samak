clear;
% Screening corrections does not exist in new TBD (commented at the end, taken from TritiumSpectrumeV)
addpath(genpath('../../../Samak2.0'));

% Load RunList
RunList = load('RunListFTFullCD.mat');
RunList = RunList.RunListFTFullCD;
RunList((RunList == 40610) | (RunList == 40612) | (RunList == 40539)) = [];
% RunList = 1;
RunList = 'StackCD100all';
% Apply any pulls you wish (recommended for FT: 2 on mnu, none on the others)
%pulls    = [Inf,Inf,Inf,Inf,Inf,Inf];

% Choose type of fit (chi2)
chi2name = 'chi2CM';

% Choose data range to analyze
dataStart = 9; % 9 = -200 eV; 1 = -1600 eV

% Choose fixed parameters (optional)
fixPar   = '5 6';

% Choose fitter
fitter = 'minuit';

% Number of Monte Carlo Simulations
Ntrials = [1];

DataEffCorr = 'OFF';

for trials = 1:length(Ntrials)
    
    N = Ntrials(trials);
    
    ResultsTable = zeros(N,10);
    
    for mc = 1:N
        MRA_OFF = MultiRunAnalysis('AnaFlag','StackPixel',...
            'chi2',chi2name,'RunList',RunList,...
            'fitter',fitter,'exclDataStart',dataStart,...
            'fixPar',fixPar,'pulls',pulls,...
            'DataEffcorr',DataEffCorr); %'ROI+PileUp'
        
        %     MRA_OFF.ModelObj.AddStatFluctTBDIS();
        %     MRA_OFF.RunData.TBDIS = MRA_OFF.ModelObj.TBDIS();
        %     MRA_OFF.ModelObj.ComputeTBDIS();
        
        MRA_OFF.Fit();
        
        mnu = MRA_OFF.FitResult.par(1);
        e0 = MRA_OFF.FitResult.par(2);
        bck = MRA_OFF.FitResult.par(3);
        norm = MRA_OFF.FitResult.par(4);
        mnu_e = MRA_OFF.FitResult.err(1);
        e0_e = MRA_OFF.FitResult.err(2);
        bck_e = MRA_OFF.FitResult.err(3);
        norm_e = MRA_OFF.FitResult.err(4);
        chi2 = MRA_OFF.FitResult.chi2min;
        dof = MRA_OFF.FitResult.dof;
        ResultsTable(mc,:) = [mnu,e0,bck,norm,mnu_e,e0_e,bck_e,norm_e,chi2,dof];
        disp(mc); 
    end
    
    %save(['FT_MC_OFF',num2str(N),'.mat'],'ResultsTable')
    
    Nring = 11;
    ResultsTableRing = zeros(N,Nring,10);
    
    for mc = 1:0
        for ring = 1:Nring
            MRA_RING = MultiRunAnalysis('AnaFlag','Ring',...
                'chi2',chi2name,'RunList',RunList,...
                'ring',ring,...
                'fitter',fitter,'exclDataStart',dataStart,...
                'fixPar',fixPar,'pulls',pulls,...
                'DataEffcorr',DataEffCorr,'ReadFake',mc); %'ROI+PileUp'
            
            %     MRA_OFF.ModelObj.AddStatFluctTBDIS();
            %     MRA_OFF.RunData.TBDIS = MRA_OFF.ModelObj.TBDIS();
            %     MRA_OFF.ModelObj.ComputeTBDIS();
            
            MRA_RING.Fit();
            
            mnu = MRA_RING.FitResult.par(1);
            e0 = MRA_RING.FitResult.par(2);
            bck = MRA_RING.FitResult.par(3);
            norm = MRA_RING.FitResult.par(4);
            mnu_e = MRA_RING.FitResult.err(1);
            e0_e = MRA_RING.FitResult.err(2);
            bck_e = MRA_RING.FitResult.err(3);
            norm_e = MRA_RING.FitResult.err(4);
            chi2 = MRA_RING.FitResult.chi2min;
            dof = MRA_RING.FitResult.dof;
            ResultsTableRing(mc,ring,:) = [mnu,e0,bck,norm,mnu_e,e0_e,bck_e,norm_e,chi2,dof];
            disp(mc); disp(ring);
        end
    end
    
    %save(['FT_MC_RING',num2str(N),'.mat'],'ResultsTableRing')
    pixellist = [1:112,114:124];
    %pixellist = 1:4;
    ResultsTableMpix = zeros(N,2+2+2+4*length(pixellist));
    fixpar = num2str([2+2*length(pixellist)+1 2+2*length(pixellist)+2]);
    for mc = 1:N
        MRA_MPIX = MultiRunAnalysis('AnaFlag','MultiPixel',...
            'chi2','chi2Stat','RunList',RunList,...
            'PixList',pixellist,...
            'fitter','matlab','exclDataStart',dataStart,...
            'fixPar',fixpar,...
            'DataEffcorr',DataEffCorr); %'ROI+PileUp'
        
        %     MRA_OFF.ModelObj.AddStatFluctTBDIS();
        %     MRA_OFF.RunData.TBDIS = MRA_OFF.ModelObj.TBDIS();
        %     MRA_OFF.ModelObj.ComputeTBDIS();
        
        MRA_MPIX.Fit();
        
        mnu = MRA_MPIX.FitResult.par(1);
        e0 = MRA_MPIX.FitResult.par(2);
        bck = MRA_MPIX.FitResult.par(3:3+length(pixellist));
        norm = MRA_MPIX.FitResult.par(3+length(pixellist)+1:end-2);
        mnu_e = MRA_MPIX.FitResult.err(1);
        e0_e = MRA_MPIX.FitResult.err(2);
        bck_e = MRA_MPIX.FitResult.err(3:3+length(pixellist));
        norm_e = MRA_MPIX.FitResult.err(3+length(pixellist)+1:end-2);
        chi2 = MRA_MPIX.FitResult.chi2min;
        dof = MRA_MPIX.FitResult.dof;
        ResultsTableMpix(mc,:) = [mnu,e0,bck,norm,mnu_e,e0_e,bck_e,norm_e,chi2,dof];
        disp(mc);
    end
    
    save(['FT_MC_MPIX',num2str(N),'.mat'],'ResultsTableMpix')
    
    
end