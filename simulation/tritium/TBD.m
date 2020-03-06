% ----------------------------------------------------------------------- %
%
% Class computing the electron energy spectrum in the
% beta decay of tritium, next to the end point energy
%
% ----------------------------------------------------------------------- %
%
% Theoretical Corrections
%  - the conventional Fermi function
%  - Effect of the finite nuclear size on the solution of the Dirac equation for the electron
%  - Screening of the latter by the orbital electron
%  - Exchange between the beta particle and the orbital electron;
%  - Convolution of the lepton and nucleonic wavefunctions through the nuclear volume;
%  - Recoil of the 3He ion in its generation of the Coulomb field;
%  - Effect of the finite nuclear mass in the energy dependence of the phase space;
%  - V-A interference;
%  - Weak magnetism;
%  - Radiative corrections
%  - Excited States
%  - Doppler Effect
%
% Detector Corrections
%  - Energy resolution
%  - Flat background
%  - MAC-E Filter Transmission Function
%
% Option:
%  - Pixel-wise computation
%  - Handling of 148-FPD pixels simultaneously
%
% Th. Lasserre, 2014-18
% P. Morales,   2017-18
% CEA - Saclay
% TUM/IAS and MPP
%
% ----------------------------------------------------------------------- %

classdef TBD < handle & WGTSMACE %!dont change superclass without modifying parsing rules!
    
    properties (Constant = true, Hidden = true)
        
        % Physics Constants
        fsC = 0.007297352569816;    % Fine Structure Cte
        me = 510.998918e3;          % Electron Mass eV
        Gf = 8.963859529398992e-50; % Fermi Cte Mev.m^3 PDG2010
        Mn = 939565.360e3;          % Neutron mass eV
        Mp = 938272.029e3;          % Proton mass eV
        gA = 1.247; gV = 1.000;     % Weak Interaction Coupling constants
        kb = 1.38064852e-23;        % [m2 kg s^-2 K^-1] Boltzmann Constant
        kbeV = 8.617343e-5;         % [eV K^-1]
        M = 2*3.0160492e-3/6.022140857e23; % [kg] mass of tritium molecule
        c = 299792458;              % [m/s] speed of light
        
        % Tritium Beta Decay
        MH3 = 2.809431917483801e9;  % Tritium mass eV, PRC77,055502(2008)
        MHe3 = 2.809413380751987e9; % Helium-3 mass eV, PRC77,055502(2008)
        TdecayC = 1.78283e-09;      % s^-1
        Z = 2;                      % 3He
        
        % Properties of the broadening kernel of the Doppler effect
        E_cms = 18573.7;              % eV Energy center of mass
        e_vel = 8.0833265115e+07;  % [m/s] velocity of the electron at E_0
        %e_vel = sqrt(18573.7*2/(510998.918/299792458^2));
        
        % ----------------------------------------------------------
        % Nuclear radius in unit of me (Elton Formula)
        % Ref: Phys.Rev.C84:024617,2011
        % Radius = 0.0029*3^(1/3) + 0.0063*3^-(1/3) - 0.017*3^(-1);
        % ----------------------------------------------------------
        Rad = 0.002884033115634;
    end
    
    properties (Access=public, Hidden = true)
        
        % Kinematics - Wilkinson Nomenclature
        dhWV0; dhWgamma; dhWw; dhWw0; dhWwbar;
        dhWp; dhWpbar; dhWy; dhWybar; V0bias;
        
        % Pixel-wise qU vector
        qUpixel;
    end
    
    properties (Access=public)
        
        % Kinematics (Electron)
        Q;Q_i;  %18575;% Tritium Q-value, eV (EPT=185833)
        Ee;     % total energy of the electron
        Ee0;    % end point TTexPressed in total neutrino energy
        Te;     % kinetic energy of the electron
        TeMin;  % minimum e kinetic energy considered
        TeMax;  % maximum e kinetic energy considered
        pe;     % electron momentum
        beta;   % electron momentum / electron total energy
        nTe;    % kinetic energy Bins
        TeStep; % energy bin steps width
        
        % Kinematics (Neutrino)
        Enu;    % total energy of the electron neutrino
        Tnu;    % kinetic energy of the electron neutrino
        pnu;    % neutrino momentum
        
        % Neutrino mass/mixing properties
        mnuSq;mnuSq_i;    % keV neutrino 1, 2, 3 masses
        mnu4Sq;mnu4Sq_i;  % keV sterile neutrino mass squared
        sin2T4;sin2T4_i;  % sin2Theta14 mixing (forget the 2x2theta)
        mTSq; mTSq_i       %--> tachyonic neutrino mass sigma^2 (1 per ring, innermost ring is fixed to 0)
        
        % Phase Space Correction for Negative Mass Squared
        PS_Wein93;
        PS_Wein93Coeff;
        
        % Beta Spectrum - Differential
        PhaseSpace;      % Primitive Beta Spectrum
        TBDDS;           % Differential Beta Spectrum
        TBDDSE;          % Error on Differential Beta Spectrum
        NormFactorTBDDS; % Normalization Factor
        mate;            % Migration Matrix For E Fast Resolution
        
        % Beta Spectrum - Integral
        TBDIS;           % Integral Beta Spectrum
        TBDISE;          % Error on Integral Beta Spectrum
        qUOffset; i_qUOffset;
        
        % Screening Correction (electron)
        ScreeningFlag; SCorr; SCorrBias = 0;
        % Finite Extension of Nucleus Charge, L0
        FiniteExtChargeFlag; RadBias; L0Corr; L0CorrBias = 0;
        % Weak interaction finite size correction, C
        WintFiniteSizeFlag; CCorr; CCorrBias = 0;
        % Electron-Electron Exchange Term
        EEexchangeFlag; ECorr; ECorrBias = 0;
        % Recoil Coulomb Effect
        RecoilCoulombFlag; QCorr; QCorrBias = 0;
        % 3He Recoil Corrections:
        RecoilWmVmAFlag; RecCorr; RecCorrBias = 0;
        % Radiative Corrections
        RadiativeFlag; RadCorr; RadCorrBias = 0; RadType = 1;
        
        % Doppler Effect
        DopplerEffectFlag; % turns on and off the Doppler Effect
        Eminus;            % Extra bins added to the end of the spectrum to account for convolution edge effect
        Eplus;             % Extra bins added to the beginning of the spectrum to account for conv. edge effect
        e_velParal;        % [m/s] electron vel. component parallel to emission direction
        cosO;              % cosine of the emission angle
        recomputeKernel;   % flag to recompute kernel
        DE_sigma_e_vel;    % [eV] standard deviation for the kernel in terms of the electron velocity
        DE_sigma;          % [eV] standard deviation for the kernel in terms of the electron energy (usually around 130 meV)
        DE_Bias;
        
        % Broadening kernel g parameters (also Doppler Effect)
        BulkVelT2;        % [m/s] bulk velocity of tritium (mean = 13)
        GaussianKernelDE; % broadening kernel
        ProbKernDE;       % probability under the curve of the kernel (should = 1)
        
        % Normalization
        CumFrac;               % Fraction of the decays that happen in the tail of the spectrum, 
        normFit; normFit_i; 
        normFitallPixels; normFitallPixels_i;
        NormFactorTBDDSPixels; % Normalization factor for all pixels, Vector
        Norm_BiasPixels;       % Normalization bias for all pixels, Vector
        
        % FSD T-T
        TTNormGS;TTNormGS_i;
        TTNormES;TTNormES_i;
        TTexE_G; TTexP_G;
        TTexE_E; TTexP_E;
        
        % FSD D-T
        DTNormGS;DTNormGS_i;
        DTNormES;DTNormES_i;
        DTexE_G; DTexP_G;
        DTexE_E; DTexP_E;
        
        % FSD H-T
        HTNormGS;HTNormGS_i;
        HTNormES;HTNormES_i;
        HTexE_G; HTexP_G;
        HTexE_E; HTexP_E;
        
        % FSD T minus ion
        TmNormGS;TmNormGS_i;
        TmNormES;TmNormES_i;
        TmexE_G; TmexP_G;
        TmexE_E; TmexP_E;
        
        % FSD broadening
        FSD_Sigma         % std of gaussian or width of rectangle    
        FSD_MultiPos      % rel. shifts 
        FSD_MultiWeights  % weighting
        FSD_MultiSigma    % std of gaussian or width of rectangle for each peak
        FSD_Dist          % Gauss or Rect
        
        % Binning
        nTeBinningFactor; % Binning for the differential spectrum
        Binning
        BinningSteps
        
        % Fit: Bias
        mnuSq_Bias;
        E0_bias;
        Bkg_Bias;
        BkgSlope_Bias;
        Norm_Bias;
        mnu4Sq_Bias;
        sin2T4_Bias;
        TTES_bias;
        TTGS_bias;
        DTES_bias;
        DTGS_bias;
        HTGS_bias;
        HTES_bias;
        qUOffset_bias;    % 1 value per spectrum (1 per ring, innermost ring is fixed to 0)
        mTSq_bias; % --> sigma^2 (1 per ring, innermost ring is fixed to 0)
        FracTm_bias; % fraction of T- ions
        
        % Integration method
        IStype; % SIMPSFAST is the only one that works right now
        
        % Response Function
        RF;
        UseParallelRF; % Use the parallelization included in Matlab to calculate the RF
        onlyBKGchange = false;
        prevBKG;
    end
    
    methods
        
        function obj    = TBD(varargin)
            
            p = inputParser;
            % p.addParameter('nPixels',148,@(x)isfloat(x) && x>0);
            
            p.addParameter('PS_Wein93','OFF',@(x)ismember(x,{'OFF','ON'}));
            p.addParameter('PS_Wein93Coeff',0.72,@(x)isfloat(x) && x>0);
            
            % TBD: Parameters
            p.addParameter('mnuSq_i',0,@(x)isfloat(x));
            p.addParameter('mTSq_i','',@(x)isfloat(x));
            p.addParameter('Q_i',18575,@(x)isfloat(x));
            p.addParameter('mnu4Sq_i',0,@(x)isfloat(x));
            p.addParameter('sin2T4_i',0,@(x)isfloat(x));%Statflu
            p.addParameter('normFit_i',0,@(x)isfloat(x));
            p.addParameter('normFitallPixels_i',0,@(x)isfloat(x));
            p.addParameter('i_qUOffset','',@(x)isfloat(x));
            p.addParameter('Binning','Default',@(x)ismember(x,{'Default','NoB'}));
            p.addParameter('BinningSteps',0.1,@(x)isfloat(x));

            % TBD: Flag for Theoretical Corrections
            p.addParameter('ScreeningFlag','OFF',@(x)ismember(x,{'OFF','ON'}));
            p.addParameter('FiniteExtChargeFlag','OFF',@(x)ismember(x,{'OFF','ON'}));
            p.addParameter('WintFiniteSizeFlag','OFF',@(x)ismember(x,{'OFF','ON'}));
            p.addParameter('EEexchangeFlag','OFF',@(x)ismember(x,{'OFF','ON'}));
            p.addParameter('RecoilCoulombFlag','OFF',@(x)ismember(x,{'OFF','ON'}));
            p.addParameter('RadiativeFlag','OFF',@(x)ismember(x,{'OFF','ON'}));
            p.addParameter('RecoilWmVmAFlag','OFF',@(x)ismember(x,{'OFF','ON'}));
            
            % Doppler Effect Flag
            p.addParameter('DopplerEffectFlag','OFF',...
                @(x)ismember(x,{'OFF','matConv','numConv','FSD_Knm1','FSD'}))
            p.addParameter('Eminus',5,@(x)isfloat(x) && x>0);
            p.addParameter('Eplus',5,@(x)isfloat(x) && x>0);
            p.addParameter('recomputeKernel','ON',@(x)strcmp(x,'ON') || strcmp(x,'OFF'));
            
            % Broadening kernel parameters
            p.addParameter('BulkVelT2',0,@(x)isfloat(x));
            % p.addParameter('DE_sigma',0.095056255902670,@(x)isfloat(x) && x>0); %0.134429846285963
            
            % FSD broadening
            p.addParameter('FSD_Sigma','',@(x)isfloat(x) || isempty(x)); % std of gaussian or width of rectangle
            p.addParameter('FSD_MultiPos','',@(x) isfloat(x) || isempty(x));     % rel. shifts
            p.addParameter('FSD_MultiWeights','',@(x) isfloat(x) || isempty(x));  % weighting
            p.addParameter('FSD_MultiSigma','',@(x) isfloat(x) || isempty(x));  % weighting
            p.addParameter('FSD_Dist','Gauss',@(x)ismember(x,{'Gauss','Rect'}));       % Gauss or Rect
            
            % Binning
            p.addParameter('nTeBinningFactor',5,@(x)isfloat(x) && x>0);
            
            % Integration Method
            p.addParameter('IStype','SIMPFAST',@(x)ismember(x,{'GaussKronrod','TRAPZFAST','SUM','SIMPFAST'}));
            
            % Parallelization
            p.addParameter('UseParallelRF','OFF',@(x)ismember(x,{'OFF','ON'}));
            
            % Parse unmatched parameters to WGTSMACE.m
            p.KeepUnmatched=1;
            p.parse(varargin{:});
            if( isempty(fieldnames(p.Unmatched)) ); Unmatched={}; else
                Unmatched = reshape(...
                    [fieldnames(p.Unmatched),struct2cell(p.Unmatched)]',...
                    [1,length(fieldnames(p.Unmatched))*2]);
            end
            obj=obj@WGTSMACE(Unmatched{:}); %Parse to superclass WGTSMACE.m
            
            % KATRIN: GENERAL SETTINGS
            obj.nTeBinningFactor    = p.Results.nTeBinningFactor;
            obj.PS_Wein93           = p.Results.PS_Wein93;
            obj.PS_Wein93Coeff      = p.Results.PS_Wein93Coeff;
            obj.Binning             = p.Results.Binning;
            obj.BinningSteps        = p.Results.BinningSteps;
            
            % TBD: Parameters
            obj.sin2T4_i            = p.Results.sin2T4_i;
            obj.mnuSq_i             = p.Results.mnuSq_i;
            obj.mTSq_i              = p.Results.mTSq_i;
            obj.mnu4Sq_i            = p.Results.mnu4Sq_i;
            obj.normFit_i           = p.Results.normFit_i;
            obj.normFitallPixels_i  = p.Results.normFitallPixels_i;
            obj.Q_i                 = p.Results.Q_i;
            obj.i_qUOffset          = p.Results.i_qUOffset;
            
            % TBD: Flag Theoretical Corrections
            obj.ScreeningFlag       = p.Results.ScreeningFlag;
            obj.FiniteExtChargeFlag = p.Results.FiniteExtChargeFlag;
            obj.WintFiniteSizeFlag  = p.Results.WintFiniteSizeFlag;
            obj.EEexchangeFlag      = p.Results.EEexchangeFlag;
            obj.RecoilCoulombFlag   = p.Results.RecoilCoulombFlag;
            obj.RadiativeFlag       = p.Results.RadiativeFlag;
            obj.RecoilWmVmAFlag     = p.Results.RecoilWmVmAFlag;
            
            % Doppler Effect Flag
            obj.DopplerEffectFlag   = p.Results.DopplerEffectFlag;
            obj.Eminus              = p.Results.Eminus;
            obj.Eplus               = p.Results.Eplus;
            obj.recomputeKernel     = p.Results.recomputeKernel;
            
            % Broadening kernel parameters
            obj.BulkVelT2           = p.Results.BulkVelT2;
            % obj.DE_sigma            = p.Results.DE_sigma;
            
            % FSD broadening
            obj.FSD_Sigma              = p.Results.FSD_Sigma;        % std of gaussian or width of rectangle
            obj.FSD_MultiPos           = p.Results.FSD_MultiPos;     % rel. shifts
            obj.FSD_MultiWeights       = p.Results.FSD_MultiWeights; % weighting
            obj.FSD_MultiSigma         = p.Results.FSD_MultiSigma;
            obj.FSD_Dist               = p.Results.FSD_Dist;             % Gauss or Rect
            
            % Integration Method
            obj.IStype              = p.Results.IStype;
            
            % Parallelization 
            obj.UseParallelRF       = p.Results.UseParallelRF;
            
            obj.computeKernel;
            
            if isempty(obj.i_qUOffset)
                obj.i_qUOffset = zeros(1,obj.nPixels);
            end
            
           if isempty(obj.mTSq_i)
               obj.mTSq_i = zeros(1,obj.nPixels);
           end
           
            % Initialize TD 
            InitializeTD(obj);
            
            % Binning - To get from the TD's
            SetTBDDSBinning(obj);
            
            % Initialize column density: only relevant for multipixel
            % InitializeColumnDensity(obj);
            
            % Initialize Energy Loss Function Parameters
             InitializeELossFunction(obj);
            
             % Apply pixelwise Ba correction
             if ~strcmp(obj.MACE_Ba_Setting,'Data')
                 SetPixel_MACE_BaEa(obj)
             end
             
            % Load FSD Properties
            ComputeTritiumPurity(obj);
            ComputeIsotropologActivityWeight(obj);
            LoadFSD(obj);
            
            % Backgrounds
            InitializeDetectorEfficiency(obj);
            InitializeBackground(obj);
                      
            % Init
            SetFitBias(obj,0);
            
            % Kinematics - WARNING - MOVED - HERE AFTER
            SetKinVariables(obj);
            
            % Initialize Response Function (RF)
            InitializeRF(obj);
            
            % Compute Normalization Factor
            ComputeFracTBDtail(obj);         % <-- calculated live
            ComputeNormFactorTBDDS(obj);

            % Mutlipixelsegmentation for Bias and Normalization
            if strcmp(obj.FPD_Segmentation,'MULTIPIXEL')
                SetFitBiasallPixels(obj,0);
            end
            
            if(strcmp(obj.quiet,'OFF')); fprintf(2,'Finished TBD.m constructor.\n'); end
            
        end % constructor
        
    end % methods
    
    methods
        
        % Response Function
        function AssignRF(obj)
            % TODO: put response function function here
            obj.RF = 0;    
        end
        
        % FSD
        function          weightedFSD(obj)
            % Theory based on: "Spectrum, Response and Statistical Model for Neutrino
            % Mass Measurements with the KATRIN Experiment"
            
            IsTomax = 3; % number of isotopologues
            Jmax    = 6; % quantum number of the molecular angular momentum
            
            %% Initialization
            EnerLev_T2 = [0., 0.004967179, 0.014886286, 0.029727069, 0.049444275, 0.073978639];
            EnerLev_DT = [0., 0.00619, 0.01855, 0.03703, 0.06156, 0.09205];
            EnerLev_HT = [0., 0.00985, 0.02949, 0.05880, 0.09761, 0.14569];
            % precalculated factor gs*gj to make calculations faster
            gsgj = [0.25 2.25 1.25 5.25 2.25 8.25 1 3 5 7 9 11 1 3 5 7 9 11];
            
            %% weights
            PJ(1:6) = gsgj(1:6).*exp(-(EnerLev_T2)./(obj.kbeV*obj.WGTS_Temp));
            PJ(1:6) = PJ(1:6)/sum(PJ(1:6)); % normalization
            PJ(7:12) = gsgj(7:12).*exp(-(EnerLev_DT)./(obj.kbeV*obj.WGTS_Temp));
            PJ(7:12) = PJ(7:12)/sum(PJ(7:12)); % normalization
            PJ(13:18) = gsgj(13:18).*exp(-(EnerLev_HT)./(obj.kbeV*obj.WGTS_Temp));
            PJ(13:18) = PJ(13:18)/sum(PJ(13:18)); % normalization
            
            %% Reading of FSD per angular momentum
            
            CommonName  = [getenv('SamakPath'),'/inputs/FSD/FSD_Roll_'];
            FSDbinningtemp = importdata([CommonName,'T2_J0.txt']);
            FSDbinning = FSDbinningtemp(:,1)';
            clear FSDbinningtemp;
            FSDJ = zeros(Jmax*IsTomax,length(FSDbinning));
            for ii = 1:Jmax*IsTomax
                if ii < 7
                    isot = 'T2';
                elseif ii < 13
                    %isot = 'DT';
                    isot = 'T2';
                else
                    %isot = 'HT';
                    isot = 'T2';
                end
                
                if mod(ii,6) == 1; J = 0; end
                
                FSDFileName = [CommonName,isot,'_J',num2str(J),'.txt'];
                FSDJtemp = importdata(FSDFileName);
                FSDJ(ii,:) = FSDJtemp(:,2)';
                J = J + 1;
            end
            clear FSDJtemp
            
            %% weighted FSD
            FSDweighted = zeros(1,length(FSDbinning));
            PJweigth = zeros(1,18);
            for ii = 1:Jmax*IsTomax
                if ii < 7
                    PJweigth(ii) = obj.WGTS_MolFrac_TT*PJ(ii);
                elseif ii < 13
                    PJweigth(ii) = obj.WGTS_MolFrac_DT*PJ(ii);
                else
                    PJweigth(ii) = obj.WGTS_MolFrac_HT*PJ(ii);
                end
                FSDweighted = FSDweighted + PJweigth(ii)*FSDJ(ii,:);
            end
            
            FSDmine = [FSDbinning',FSDweighted'];
            save([getenv('SamakPath'),'/inputs/FSD/FSD_ROLL_',num2str(obj.WGTS_Temp),...
                'K_TT',num2str(round(obj.WGTS_MolFrac_TT*100)),'DT',num2str(round(obj.WGTS_MolFrac_DT*100)),...
                'HT',num2str(round(obj.WGTS_MolFrac_HT*100)),'.mat'],'FSDmine')
            
        end
        
        function  LoadFSD(obj,varargin)
                 % Load FSD Excitation Energies and Probabilities
            % If requested by the user
            % For: T-T / D-T / H-T / T^-
            p=inputParser;
            p.addParameter('Sigma',obj.FSD_Sigma,@(x)isfloat(x) || isempty(x));                % broadening of FSD
            p.addParameter('MultiPos',obj.FSD_MultiPos,@(x) isfloat(x) || isempty(x));         %3 gaussians instead of using 1 gaussian per energy (for 3 RW settings): relative position
            p.addParameter('MultiWeights',obj.FSD_MultiWeights,@(x) isfloat(x) || isempty(x)); %3 gaussians instead of using 1 gaussian per energy: relative weight
            p.addParameter('BinningFactor',2,@(x) isfloat(x) || isempty(x));                   % enhance binning: twice, 3 times,... as much bins
            p.addParameter('SanityPlot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('ZoomPlot','OFF',@(x)ismember(x,{'ON','OFF'})); % save also zoom to 1 final state
            p.addParameter('Dist',obj.FSD_Dist,@(x)ismember(x,{'Gauss','Rect'}));
            p.parse(varargin{:});
            
            Sigma              = p.Results.Sigma;
            MultiPos           = p.Results.MultiPos;
            MultiWeights       = p.Results.MultiWeights;
            BinningFactor      = p.Results.BinningFactor;
            SanityPlot         = p.Results.SanityPlot;
            ZoomPlot           = p.Results.ZoomPlot;
            Dist               = p.Results.Dist; 
            
            nPseudoRings = numel(obj.MACE_Ba_T);
            
            if ~isempty(MultiPos) % divide FSD into three peaks or not
                nPeaks = numel(MultiPos)./nPseudoRings;
            else
                nPeaks = 1;
            end
            
            if ~isempty(Sigma)
                if numel(Sigma)==1
                    Sigma = squeeze(repmat(Sigma,[3,nPeaks,nPseudoRings])); % 3 isotopologues x nPeaks x nRings
                elseif numel(Sigma)==nPeaks % sigma is not ringwise
                    Sigma = squeeze(repmat(Sigma,[3,1])); 
                elseif numel(Sigma)==nPeaks*nPseudoRings
                    Sigma = squeeze(repmat(Sigma,[1,1,3])); 
                    Sigma = permute(Sigma,[3,1,2]);
                end
            end

            if strcmp(obj.DopplerEffectFlag,'FSD')
                % include Doppler Effect by Modification of FSD
                 if isempty(Sigma)
                     Sigma = obj.DE_sigma.*squeeze(ones(3,nPeaks,nPseudoRings));
                 else
                     Sigma = sqrt(Sigma.^2 + obj.DE_sigma.^2);
                 end
            end
            
            FSDConvArg = {'BinningFactor',BinningFactor,'Dist',Dist,...
                         'MultiWeights',MultiWeights,'MultiPos',MultiPos,...
                         }; %arguments for FSD convolution (optional)
            %% T-T FSD
            FSDdir = [getenv('SamakPath'),'/inputs/FSD/'];
            switch obj.TTFSD
                case {'DOSS','DOSSNOEE'}
                    ttfsdfilename = [FSDdir,'FSD_DOSS_T2_rebinned.dat'];
                case {'SAENZ','SAENZNOEE'}
                    ttfsdfilename = [FSDdir,'FSD_Saenz_T2mod.dat']; %./data/FSD/FSD_Saenz_T2mod_6.dat'
                    % ttfsdfile = importdata([getenv('SamakPath'),'/inputs/FSD/FSD_Saenz_T2_doppler.txt']); 
                case 'ROLL'
                    ttfsdfilename = [FSDdir,'FSD_ROLL_',num2str(obj.WGTS_Temp),...
                        'K_TT',num2str(round(obj.WGTS_MolFrac_TT*100)),'DT',num2str(round(obj.WGTS_MolFrac_DT*100)),...
                        'HT',num2str(round(obj.WGTS_MolFrac_HT*100)),'.mat'];
                    if ~(exist(ttfsdfilename, 'file') == 2)
                        obj.weightedFSD();
                    end
                case 'HT' %take fsd from HT instead
                    ttfsdfilename = [FSDdir,'FSD_Saenz_HTmod.txt'];
                case {'BlindingKNM1'}
                    if ~strcmp(obj.DopplerEffectFlag,'FSD_Knm1')
                        ttfsdfilename = [FSDdir,'FSD_Blinding_KNM1_T2.txt'];
                    elseif strcmp(obj.DopplerEffectFlag,'FSD_Knm1')
                        ttfsdfilename = [FSDdir,'FSD_KNM1_DOPPLER_T2.txt'];
                    end
                case {'BlindingKNM2'}
                    ttfsdfilename = [FSDdir,'FSD_KNM2_T2_Blinding.txt'];
                case {'WGTS100K'} % Test WGTS@100K
                    ttfsdfilename = [FSDdir,'FSD_Saenz_100K_doppler1_T2.txt'];
                case {'Sibille'}
                    ttfsdfilename = [FSDdir,'FSD_KNM1_T2_Doppler.txt'];  
                case {'SibilleFull'}
                    ttfsdfilename = [FSDdir,'FSD_KNM1_T2_Doppler0p5eV_FullRange.txt'];    
                case {'Sibille0p5eV'}
                    ttfsdfilename = [FSDdir,'FSD_KNM1_T2_Doppler0p5eV.txt'];   
            end
            ttfsdfile = importdata(ttfsdfilename);
            [obj.TTexE, TTexE_index] = sort(ttfsdfile(:,1)); %sort from small to large excitation energies           
            obj.TTexE = obj.TTexE';
            obj.TTexP = (ttfsdfile(TTexE_index,2))';
            if ~isempty(Sigma)   %broaden FSDs
                [obj.TTexE,obj.TTexP] = FSD_Convfun(obj.TTexE,obj.TTexP,...
                    squeeze(Sigma(1,:,:)),...
                    FSDConvArg{:},'SanityPlot',SanityPlot,'ZoomPlot',ZoomPlot,...
                    'filename',ttfsdfilename);
            end
            obj.TTGSTh = obj.GetFSDTh(obj.TTexE);          % Limit ground / excited states
            obj.TTNormGS_i = sum(obj.TTexP(:,1:obj.TTGSTh),2);
            obj.TTNormES_i = sum(obj.TTexP(:,obj.TTGSTh:end),2);
            
            %% D-T FSD
            switch obj.DTFSD
                case 'DOSS'
                    dtfsdfilename = [FSDdir,'FSD_DOSS_DT_rebinned.txt'];
                case 'HTFSD' % take HT FSD instead if DT (for testing)
                    dtfsdfilename = [FSDdir,'FSD_Saenz_HTmod.txt'];
                 %    dtfsdfile = importdata([getenv('SamakPath'),'/inputs/FSD/FSD_Saenz_T2_doppler.txt']);
                case 'TTFSD' % take HT FSD instead if TT (for testing)
                    dtfsdfilename = [FSDdir,'FSD_Saenz_T2mod.dat'];
                case 'BlindingKNM1' % Blinding Test KNM1
                    if ~strcmp(obj.DopplerEffectFlag,'FSD_Knm1')
                        dtfsdfilename = [FSDdir,'FSD_Blinding_KNM1_DT.txt'];
                    elseif strcmp(obj.DopplerEffectFlag,'FSD_Knm1')
                        dtfsdfilename = [FSDdir,'FSD_KNM1_DOPPLER_DT.txt'];
                    end
                case {'BlindingKNM2'}
                    dtfsdfilename = [FSDdir,'FSD_KNM2_DT_Blinding.txt'];
                case 'WGTS100K' % Test WGTS@100K
                    dtfsdfilename = [FSDdir,'FSD_Saenz_100K_doppler1_HT.txt'];
                case {'Sibille'}
                    dtfsdfilename = [FSDdir,'FSD_KNM1_DT_Doppler.txt'];
                case {'SibilleFull'}
                    dtfsdfilename = [FSDdir,'inputs/FSD/FSD_KNM1_DT_Doppler0p5eV_FullRange.txt'];
                case {'Sibille0p5eV'}
                    dtfsdfilename = [FSDdir,'inputs/FSD/FSD_KNM1_DT_Doppler0p5eV.txt'];
            end
            dtfsdfile = importdata(dtfsdfilename);
            [obj.DTexE, DTexE_index] = sort(dtfsdfile(:,1));
            obj.DTexE = obj.DTexE';
            obj.DTexP = (dtfsdfile(DTexE_index,2))';
            if ~isempty(Sigma)  %broaden FSDs
                [obj.DTexE,obj.DTexP] = FSD_Convfun(obj.DTexE,obj.DTexP,squeeze(Sigma(2,:,:)),...
                                                  FSDConvArg{:},'filename',dtfsdfilename);
            end
            obj.DTGSTh = obj.GetFSDTh(obj.DTexE); % Limit ground / excited states
            obj.DTNormGS_i = sum(obj.DTexP(:,1:obj.DTGSTh),2);
            obj.DTNormES_i = sum(obj.DTexP(:,obj.DTGSTh:end),2);
            
            %% H-T FSD
            switch obj.HTFSD
                case 'SAENZ'
                    htfsdfilename = [FSDdir,'FSD_Saenz_HTmod.txt'];
                   % htfsdfile = importdata([getenv('SamakPath'),'/inputs/FSD/FSD_Saenz_HT_doppler.txt']);
                case 'BlindingKNM1'
                    if ~strcmp(obj.DopplerEffectFlag,'FSD_Knm1')
                        htfsdfilename = [FSDdir,'FSD_Blinding_KNM1_HT.txt'];
                    elseif strcmp(obj.DopplerEffectFlag,'FSD_Knm1')
                        htfsdfilename = [FSDdir,'FSD_KNM1_DOPPLER_HT.txt'];
                    end
                case {'BlindingKNM2'}
                    htfsdfilename = [FSDdir,'FSD_KNM2_HT_Blinding.txt'];
                case 'WGTS100K' % Test WGTS@100K
                    htfsdfilename = [FSDdir,'FSD_Saenz_100K_doppler1_HT.txt'];
                case {'Sibille'}
                    htfsdfilename = [FSDdir,'FSD_KNM1_HT_Doppler.txt'];
                case {'SibilleFull'}
                    htfsdfilename = [FSDdir,'FSD_KNM1_HT_Doppler0p5eV_FullRange.txt'];
                case {'Sibille0p5eV'}
                    htfsdfilename = [FSDdir,'FSD_KNM1_HT_Doppler0p5eV.txt'];
            end
            htfsdfile = importdata(htfsdfilename);
            [obj.HTexE, HTexE_index] = sort(htfsdfile(:,1));
            obj.HTexE = obj.HTexE';
            obj.HTexP = (htfsdfile(HTexE_index,2))';
            if ~isempty(Sigma)  %broaden FSDs
                [obj.HTexE,obj.HTexP] = FSD_Convfun(obj.HTexE,obj.HTexP,squeeze(Sigma(3,:,:)),...
                                       FSDConvArg{:},'filename',htfsdfilename);
            end
            obj.HTGSTh = obj.GetFSDTh(obj.HTexE);   % Limit ground / excited states
            obj.HTNormGS_i = sum(obj.HTexP(:,1:obj.HTGSTh),2);
            obj.HTNormES_i = sum(obj.HTexP(:,obj.HTGSTh:end),2);
           %% T minus ion
            switch obj.TmFSD
                case 'OFF'
                case 'SAENZ'
                    tmfsdfile = importdata([FSDdir,'FSD_Tminus.mat']);
                    [obj.TmexE, TmexE_index] = sort(tmfsdfile(:,1));
                    obj.TmexE = obj.TmexE';
                    obj.TmexP = (tmfsdfile(TmexE_index,2))';
                    obj.TmGSTh = 15;% % Index, Limit ground / excited states
                    obj.TmNormGS_i = sum(obj.TmexP(1:obj.TmGSTh));
                    obj.TmNormES_i = sum(obj.TmexP(obj.TmGSTh+1:end));
                    obj.TmNormGS =  obj.TmNormGS_i;
                    obj.TmNormES =  obj.TmNormES_i;
            end
        end
        function Threshold = GetFSDTh(obj,TexE)
            % Compute FSD threshold
            energyTh = 10; %threshold energy (eV), which divides GS and ES
            if size(TexE,1)>1
                Threshold = zeros(size(TexE,1),1);
                for i=1:size(TexE,1)
                    Threshold(i) = find(TexE(i,:)<energyTh,1,'last');
                end
            else
                indexTH = find(TexE<energyTh);
                Threshold = indexTH(end);
            end
            Threshold = ceil(mean(Threshold));
        end
        function [GES,TexE,TexP] = ComputeFSD_GSES(obj,TexE,TexP,TGSTh,state)
            % Normalization Factors Ground/Excited States
            switch state
                case 'ground'
                    TexE  = TexE(:,1:TGSTh);        %TGSTh = limit ground state-excited states
                    TNorm_i = sum(TexP(:,1:TGSTh),2); % initial Normalization
                    TexP  = TexP(:,1:TGSTh)./TNorm_i;  %normalize to GS to 1
                case 'excited'
                    TexE  = TexE(:,TGSTh+1:end);
                    TNorm_i = sum(TexP(:,TGSTh+1:end),2); % initial Normalization
                    TexP  = TexP(:,TGSTh+1:end)./TNorm_i; %normalize ES to 1
            end
            
            mNuSq_local = obj.mnuSq-2.*obj.mTSq;
            Q_local = obj.Q;
            
            if sum(obj.mTSq)==0
                mNuSq_local = obj.mnuSq;
                %   TexP_M = repmat(TexP,obj.nTe,1); % probability
                %   TexE_M = repmat(TexE,obj.nTe,1); % energy
                TexP_M = squeeze(permute(repmat(TexP,1,1,obj.nTe),[3,2,1])); % probability
                TexE_M = squeeze(permute(repmat(TexE,1,1,obj.nTe),[3,2,1])); % energy
                Te_M = squeeze(repmat(obj.Te,[1,size(TexP')])); % energies of electron
                pe_M = squeeze(repmat(obj.pe,[1,size(TexP')])); % momentum of electron
                Q_M = squeeze(Q_local'.*ones([obj.nTe,size(TexP')]));
                me_M = squeeze(obj.me.*ones([obj.nTe,size(TexP')]));
                sin2T4_M = squeeze(obj.sin2T4.*ones([obj.nTe,size(TexP')]));
                mnuSq_M = squeeze(mNuSq_local'.*ones([obj.nTe,size(TexP')]));
                mnu4Sq_M = squeeze(obj.mnu4Sq.*ones([obj.nTe,size(TexP')]));
            else
                % Ground/Excited State
                TexP_M = repmat(TexP,obj.nTe,1); % probability
                TexE_M = repmat(TexE,obj.nTe,1); % energy
                
                Te_M = repmat(obj.Te,1,numel(TexP),obj.nPixels); % energies of electron
                pe_M = repmat(obj.pe,1,numel(TexP),obj.nPixels); % momentum of electron
                Q_M = permute(Q_local'.*ones(obj.nPixels,obj.nTe,numel(TexP)),[2 3 1]);
                me_M = obj.me*ones(obj.nTe,numel(TexP));
                sin2T4_M = obj.sin2T4*ones(obj.nTe,numel(TexP));
                mnuSq_M = permute(mNuSq_local'.*ones(obj.nPixels,obj.nTe,numel(TexP)),[2 3 1]);
                mnu4Sq_M = obj.mnu4Sq*ones(obj.nTe,numel(TexP));
            end
            % normal phase space formula, but
            % actual energy of the electron =: initial energy of electron - excitation energy of daughter molecule
            % 1 differential spectrum per excitation energy of the daughter molecule (1 phase space per FSD state)
            % Combination: weight by the exciation probability and sum (weighted mean)
            %              GES = squeeze(real(0 + sum(...
            %                  ((Q_M-Te_M-TexE_M)>=0)...% if (energy of electron - excitation energy of daughter molecule) < Endpoint (Q)
            %                  .*pe_M.*(Te_M+me_M).*(Q_M-TexE_M-Te_M).*(((ones(obj.nTe,numel(TexP))-sin2T4_M)....
            %                  .*(((Q_M-Te_M-TexE_M).^2-mnuSq_M)>=0)...
            %                  .*((Q_M-TexE_M-Te_M).^2-mnuSq_M).^.5) ...
            %                  + sin2T4_M.*(((Q_M-Te_M-TexE_M).^2-mnu4Sq_M)>=0)...
            %                  .*((Q_M-TexE_M-Te_M).^2-mnu4Sq_M).^.5).*TexP_M,2)));
            GES = squeeze(real(0 + sum(...
                ((Q_M-Te_M-TexE_M)>=0)...% if (energy of electron - excitation energy of daughter molecule) < Endpoint (Q)
                .*pe_M.*(Te_M+me_M).*(Q_M-TexE_M-Te_M).*(((ones([obj.nTe,size(TexP')])-sin2T4_M)....
                .*(((Q_M-Te_M-TexE_M).^2-mnuSq_M)>=0)...
                .*((Q_M-TexE_M-Te_M).^2-mnuSq_M).^.5) ...
                + sin2T4_M.*(((Q_M-Te_M-TexE_M).^2-mnu4Sq_M)>=0)...
                .*((Q_M-TexE_M-Te_M).^2-mnu4Sq_M).^.5).*TexP_M,2)));
        end
        
        % Kinematics / Binning
        function          SetKinVariables(obj)
            %
            % Set the kinematical parameters
            %
            % T. Lasserre
            % Updated: Feb. 9 2018
            %
            
            % Kinematical Variables - Classical Nomenclature
            obj.Ee     = obj.me+obj.Te;             % total electron energy
            obj.pe     = (obj.Ee.^2-obj.me.^2).^.5; % electron momentum
            obj.beta   = obj.pe./obj.Ee;            % beta (relativistic)
            
            % Kinematical Variables - Wilkinson Nomenclature
            obj.V0bias   = 0; % Biasing Factor of the screening potential
            obj.dhWV0    = 1.45*obj.fsC^2.*(1+obj.V0bias);
            obj.dhWw     = obj.Ee/obj.me;
            obj.dhWw0    = (obj.Q + obj.me)/obj.me;
            obj.dhWwbar  = obj.dhWw - obj.dhWV0;
            obj.dhWp     = sqrt(obj.dhWw.^2-1);
            obj.dhWpbar  = sqrt(obj.dhWwbar.^2-1);
            obj.dhWy     = obj.fsC * obj.Z * obj.dhWw ./ obj.dhWp;
            obj.dhWybar  = obj.fsC * obj.Z * obj.dhWwbar ./ obj.dhWpbar;
            obj.dhWgamma = sqrt(1-(obj.fsC * obj.Z).^2);
        end
        function          SetTBDDSBinning(obj,varargin)
            p = inputParser;
            p.addParameter('TeStep',0.1,@(x)isfloat(x));
            p.parse(varargin{:});
            TeStep_local = p.Results.TeStep;
            
            qUtmp      = mean(obj.qU,2);    % in case of segmented FPD: mean qU
            qUtmp(1)   = min(min(obj.qU));  % extreme qU values at end and beginning
            qUtmp(end) = max(max(obj.qU));
            
            if ismember(obj.TD,{'FTpaper','StackCD100all','StackCD100_3hours'}) || contains(obj.TD,'FTpaper') % first tritium binning!
                obj.Te = interp1(linspace(1,obj.nqU,obj.nqU)',...
                    qUtmp,linspace(1,obj.nqU,obj.nqU*obj.nTeBinningFactor-(obj.nTeBinningFactor-1))');              
            elseif strcmp(obj.TD,'DScomparison')
                obj.Te = ((18573.70-90):0.1:(18573.70+135))';
            elseif strcmp(obj.TD,'RFcomparison')
                obj.Te = (18540:0.01:18635)'; 
            else
                % HARDCODED - 14/05/2018 - Thierry
                nPoints = ceil((qUtmp(end)-qUtmp(1))/TeStep_local); % 0.1eV Binning necessary for doppler effect
                obj.Te = linspace(qUtmp(1),qUtmp(end),nPoints)'; 
            end

            % Truncating Te values > qUmax
            obj.TeMin = min(obj.Te);%min(min(obj.qU));
            obj.TeMax = max(obj.Te);%max(max(obj.qU));
%             if obj.nTeBinningFactor<=4
%                 obj.Te = obj.Te(obj.Te<=obj.TeMax);
%                 obj.Te = reshape(obj.Te,[obj.nTe obj.nPixels]);
%             end
            obj.nTe               = size(obj.Te,1);
            obj.TeStep            = obj.Te(2,1)-obj.Te(1,1);
        end
        
        function          FillTableqUpixel(obj)
            % Build qU vector for all pixel
            for p=1:1:obj.nPixels
                obj.qUpixel(p,:) =    [];
            end
        end
        function InitializeRF(obj,varargin)   
            p=inputParser;
            p.addParameter('saveRF','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('RebinMode','Fast',@(x)ismember(x,{'Fast','Integral'})); % for RF broadening
            
            p.parse(varargin{:});
            saveRF    = p.Results.saveRF;
            RebinMode = p.Results.RebinMode;
            
            % load or calculate response function according to settings
            %% load response function
            if strcmp(obj.recomputeRF,'OFF') %try loading
                loadSuccess = LoadOrSaveRF(obj,'load'); % loadSuccess==1 when loading was successful
            elseif strcmp(obj.recomputeRF,'ON') %do not attempt loading
                loadSuccess = 0;
            end
            %% calculate response function
            if loadSuccess==0
                tfloc = obj.GetRF(); %function handle of response function
                switch obj.FPD_Segmentation
                    case {'SINGLEPIXEL','MULTIPIXEL'}
                        switch obj.UseParallelRF
                            case 'OFF'
                                for p = 1:obj.nPixels
                                    for ii = 1:obj.nqU
                                        tmpRF = tfloc(obj.Te(:,p),obj.qU(ii,p),'pixel',obj.FPD_PixList(p));
                                        obj.RF(:,ii,p) = tmpRF;
                                    end
                                end
                            case 'ON'
                                %parTe_all = obj.Te;
                                parTe = obj.Te;
                                parqU = obj.qU;
                                parnqU = obj.nqU;
                                parPixel_all = obj.FPD_PixList;
                                parRF = zeros(obj.nTe,obj.nqU,obj.nPixels);
                                for p = 1:obj.nPixels
                                    %parTe = parTe_all(:,p);
                                    parPixel = parPixel_all(p);
                                    parfor ii = 1:parnqU
                                        parRF(:,ii,p) = tfloc(parTe,parqU(ii,p),'pixel',parPixel);
                                    end
                                    
                                end
                                obj.RF = parRF;
                        end
                    case 'OFF'
                        switch obj.UseParallelRF
                            case 'OFF'
                                for ii = 1:obj.nqU
                                    obj.RF(:,ii) = tfloc(obj.Te,obj.qU(ii));
                                end
                            case 'ON'
                                parTe = obj.Te;
                                parqU = obj.qU;
                                parRF = zeros(obj.nTe,obj.nqU);
                                 %                               parfor ii = 1:obj.nqU
                                for ii = 1:obj.nqU
                                    parRF(:,ii) = tfloc(parTe,parqU(ii));
                                end
                                obj.RF = parRF;
                        end
                    case 'RING'
                        switch obj.UseParallelRF
                            case 'OFF'
                                for ri = 1:obj.nRings
                                    for ii = 1:obj.nqU
                                        if strcmp(obj.KTFFlag,'WGTSMACE')
                                            tmpRF = tfloc(obj.Te(:,ri),obj.qU(ii,ri),'pixel',ri);
                                            obj.RF(:,ii,ri) = tmpRF;
                                        else
                                            tmpRF = tfloc(obj.Te(:,ri),obj.qU(ii,ri),'pixel',ri);
                                            obj.RF(:,ii,ri) = tmpRF;
                                        end
                                    end
                                end
                            case 'ON'
                                parqU = obj.qU;
                                parnqU = obj.nqU;
                                parRF = zeros(obj.nTe,obj.nqU,obj.nRings);
                                for ri = 1:obj.nRings
                                    parTe = obj.Te;
                                    parfor ii = 1:parnqU
                                        parRF(:,ii,ri) = tfloc(parTe,parqU(ii,ri),'pixel',ri);
                                    end
                                    
                                end
                                obj.RF = parRF;
                        end
                end
                
                obj.RF(obj.RF<0)=0;
                if strcmp(saveRF,'ON')
                    LoadOrSaveRF(obj,'save')
                end
            end % if loadSuccess
            
            if obj.MACE_Sigma>0
                
                if numel(obj.MACE_Sigma)==1
                    obj.MACE_Sigma = repmat(obj.MACE_Sigma,[obj.nqU,1]); % subrunwise broadening
                end
                
                out = zeros(size(obj.RF));
                for i=1:obj.nqU
                    fprintf('Broadening response function qU(%.0f) - std = %.0f meV \n',i,1e3*obj.MACE_Sigma(i));
                    out(:,i) = TF_Convfun(obj.RF(:,i)',obj.qU(i),obj.Te,obj.MACE_Sigma(i),'RebinMode',RebinMode);
                end
                obj.RF = out;
            end
            
        end
        % Fit Bias
        function          SetFitBias(obj,flag)
            %
            % Set the fitting bias parameters
            % (single pixel fit)
            %
            % T. Lasserre
            % Updated: Feb. 9 2018
            %
            if flag == 0
                obj.Q           = obj.Q_i;
                obj.qUOffset    = obj.i_qUOffset;
                obj.BKG_RateSecallPixels    = obj.BKG_RateSecallPixels_i;             
                obj.mnuSq_Bias  = 0;
                obj.E0_bias     = 0;
                obj.Bkg_Bias    = zeros(1,obj.nPixels);
                obj.BkgSlope_Bias  = obj.BKG_Slope_i;
                obj.Norm_Bias   = zeros(1,obj.nPixels);
                obj.mnu4Sq_Bias = 0;
                obj.sin2T4_Bias = 0;
                obj.TTGS_bias   = 0;
                obj.TTES_bias   = 0;
                obj.DTGS_bias   = 0;
                obj.DTES_bias   = 0;
                obj.HTGS_bias   = 0;
                obj.HTES_bias   = 0;  
                obj.qUOffset_bias     = obj.i_qUOffset;  
                obj.mTSq_bias = obj.mTSq_i;
                obj.FracTm_bias = 0;
                obj.WGTS_MolFrac_Tm =  obj.WGTS_MolFrac_Tm_i;
            end
            if flag == 1
                % Add Bias
                obj.Q           = obj.Q_i       + obj.E0_bias;
                obj.mnuSq       = obj.mnuSq_i   + obj.mnuSq_Bias;
                obj.normFit     = obj.normFit_i + obj.Norm_Bias;
                obj.BKG_RateSec = obj.BKG_RateSec_i + obj.Bkg_Bias;
                obj.BKG_Slope   = obj.BKG_Slope_i + obj.BkgSlope_Bias;
                obj.mnu4Sq      = obj.mnu4Sq_i  + obj.mnu4Sq_Bias;
                obj.sin2T4      = obj.sin2T4_i  + obj.sin2T4_Bias;
                obj.TTNormGS    = obj.TTNormGS_i + obj.TTGS_bias;
                obj.TTNormES    = obj.TTNormES_i + obj.TTES_bias;
                obj.DTNormGS    = obj.DTNormGS_i + obj.DTGS_bias;
                obj.DTNormES    = obj.DTNormES_i + obj.DTES_bias;
                obj.HTNormGS    = obj.HTNormGS_i + obj.HTGS_bias;
                obj.HTNormES    = obj.HTNormES_i + obj.HTES_bias;
                obj.qUOffset    = obj.i_qUOffset + obj.qUOffset_bias;             
                obj.mTSq = obj.mTSq_i + obj.mTSq_bias;
                obj.WGTS_MolFrac_Tm =  obj.WGTS_MolFrac_Tm_i + obj.FracTm_bias;
                % if ( (strcmp(obj.TTFSD,'SAENZNOEE')) || (strcmp(obj.TTFSD,'DOSSNOEE')) )
                %    obj.TTNormES = 0;
                % end
            end
        end
        
        function          SetFitBiasallPixels(obj,flag)
            %
            % Set the fitting bias parameters
            % (multi pixel fit)
            %
            % T. Lasserre
            % Updated: Feb. 9 2018
            %
            switch obj.FPD_Segmentation
                case {'MULTIPIXEL'}
                    dim = obj.nPixels;
                case 'RING'
                    dim = obj.nRings;
            end
            if flag == 0
                obj.Norm_BiasPixels = zeros(1,dim);
            end
            
            if flag == 1
                obj.normFitallPixels = obj.normFitallPixels_i + obj.Norm_BiasPixels;
                obj.mnu4Sq    = obj.mnu4Sq_i   + obj.mnu4Sq_Bias;
                obj.sin2T4    = obj.sin2T4_i   + obj.sin2T4_Bias;
            end
        end
        
        % Spectrum Computation - 1 pixel Equivalent
        function          ComputePhaseSpace(obj)
            % ----------------------------------------------------------
            % Beta Decay Differential Spectrum: Phase Space Factor
            % With or without Final State Distributions for
            % T-T - fTT purity
            %     - TTFSD = 1) Doss: J.Phys.,B41,125701
            %     - TTFSD = 2) Saenz: Phys.Rev.Lett. 84 (2000) no.2, 242
            % D-T - fDT purity
            % H-T - fHT purity
            % Spectrum = fTT (T-T) + fDT (D-T)/2 + fHT (H-T)/2
            % ----------------------------------------------------------
            
            % Initialization
            switch obj.FPD_Segmentation
                case {'OFF'}%,'RING'}
                    dim = {obj.nTe,1};
                 case 'RING'
                     dim = {obj.nTe,obj.nRings}; 
            end
            obj.PhaseSpace = zeros(dim{:});
            TTps = zeros(dim{:});
            DTps = zeros(dim{:});
            HTps = zeros(dim{:});
            
            mNuSq_local = obj.mnuSq-2.*obj.mTSq;
            Q_local = obj.Q; %+ obj.qUOffset;
            % NO FSD
            if ((strcmp(obj.TTFSD,'OFF')) ...
                    && (strcmp(obj.DTFSD,'OFF')) ...
                    && (strcmp(obj.HTFSD,'OFF')))
                
                obj.PhaseSpace = real(0 + ((Q_local-obj.Te)>=0).*...
                    (obj.pe.*(obj.Te+obj.me).*(Q_local-obj.Te).*((1-obj.sin2T4)...
                    .*(((Q_local-obj.Te).^2-mNuSq_local) >= 0)...
                    .*((Q_local-obj.Te).^2-mNuSq_local).^0.5)...
                    + obj.sin2T4.*(((Q_local-obj.Te).^2-obj.mnu4Sq)>=0)... % sterile neutrino
                    .*((Q_local-obj.Te).^2-obj.mnu4Sq).^0.5));
                obj.PhaseSpace(isnan(obj.PhaseSpace)) = 0;
                
            else
                
                if obj.WGTS_MolFrac_TT>0
                    switch obj.TTFSD
                        case {'DOSS','SAENZ','SAENZNOEE','ROLL','BlindingKNM1','WGTS100K','Sibille','Sibille0p5eV','SibilleFull','BlindingKNM2'}
                            % Ground State
                            [GS,obj.TTexE_G,obj.TTexP_G] = ...
                                obj.ComputeFSD_GSES(obj.TTexE,obj.TTexP,obj.TTGSTh,'ground');
                            
                            % Exited States
                            [ES,obj.TTexE_E,obj.TTexP_E] = ...
                                obj.ComputeFSD_GSES(obj.TTexE,obj.TTexP,obj.TTGSTh,'excited');
                            
                            % Build Ground + Exited States
                            if ( (strcmp(obj.TTFSD,'SAENZNOEE')) || (strcmp(obj.TTFSD,'DOSSNOEE')) )
                                obj.TTNormES = 0;
                            end
                            TTps = (obj.TTNormGS)'.*GS + (obj.TTNormES)'.*ES;
                            if strcmp(obj.TTFSD,'ROLL')
                                % FSD with ROLL already includes also DT and HT
                                obj.PhaseSpace = obj.PhaseSpace + TTps;
                            else
                                obj.PhaseSpace = obj.PhaseSpace + obj.WGTS_MolFrac_TT.*TTps;
                            end
                    end
                end
                
                if obj.WGTS_MolFrac_DT>0
                    switch obj.DTFSD
                        case {'DOSS','HTFSD','TTFSD','BlindingKNM1','WGTS100K','Sibille','Sibille0p5eV','SibilleFull','BlindingKNM2'}
                            % Ground State
                            [GS,obj.DTexE_G,obj.DTexP_G] = ...
                                obj.ComputeFSD_GSES(obj.DTexE,obj.DTexP,obj.DTGSTh,'ground');
                            
                            % Exited States
                            [ES,obj.DTexE_E,obj.DTexP_E] = ...
                                obj.ComputeFSD_GSES(obj.DTexE,obj.DTexP,obj.DTGSTh ,'excited');
                            
                            % Build Ground + Exited States
                            %obj.DTNormES = 0;
                            DTps = obj.DTNormGS' .* GS + obj.DTNormES' .* ES;
                            obj.PhaseSpace = obj.PhaseSpace + 0.5*obj.WGTS_MolFrac_DT.*DTps;   
                    end
                end
                
                if obj.WGTS_MolFrac_HT>0
                    switch obj.HTFSD
                        case {'SAENZ','BlindingKNM1','WGTS100K','Sibille','Sibille0p5eV','SibilleFull','BlindingKNM2'}
                            
                            % Ground State
                            [GS,obj.HTexE_G,obj.HTexP_G] = ...
                                obj.ComputeFSD_GSES(obj.HTexE,obj.HTexP,obj.HTGSTh,'ground');
                            
                            % Exited States
                            [ES,obj.HTexE_E,obj.HTexP_E] = ...
                                obj.ComputeFSD_GSES(obj.HTexE,obj.HTexP,obj.HTGSTh,'excited');
                            
                            % Build Ground + Exited States
                            HTps = obj.HTNormGS' .* GS + obj.HTNormES' .* ES;
                            obj.PhaseSpace = obj.PhaseSpace + 0.5*obj.WGTS_MolFrac_HT.*HTps;
                    end
                    
                end
            end
            
            switch obj.TmFSD
                case 'OFF'
                    % no contribution
                    % warning: different definition than for other isotopologues
                case 'SAENZ'
                    if obj.WGTS_MolFrac_Tm~=0
                        % Ground State
                        [GS,obj.TmexE_G,obj.TmexP_G] = ...
                            obj.ComputeFSD_GSES(obj.TmexE,obj.TmexP,obj.TmGSTh,'ground');
                        
                        % Exited States
                        [ES,obj.TmexE_E,obj.TmexP_E] = ...
                            obj.ComputeFSD_GSES(obj.TmexE,obj.TmexP,obj.TmGSTh,'excited');
                        
                        % Build Ground + Exited States
                        Tmps = obj.TmNormGS .* GS + obj.TmNormES .* ES;
                        obj.PhaseSpace = obj.PhaseSpace + obj.WGTS_MolFrac_Tm.*Tmps;
                    end
            end
        end
        
        function          ComputePhaseSpaceWei93(obj)
            % ----------------------------------------------------------
            % Beta Decay Differential Spectrum: Phase Space Factor
            %
            % Consider the case of negative m^2 value according to
            % C. Weinheimer et al.
            % PLB Volume 300, Issue 3, 11 February 1993, Pages 210-216
            %
            % With or without Final State Distributions for
            % T-T - fTT purity
            %     - TTFSD = 1) Doss: J.Phys.,B41,125701
            %     - TTFSD = 2) Saenz: Phys.Rev.Lett. 84 (2000) no.2, 242
            % D-T - fDT purity
            % H-T - fHT purity
            % Spectrum = fTT (T-T) + fDT (D-T) + fHT (H-T)
            % ----------------------------------------------------------
            
            % Init - 0.716
            m    = (obj.mnuSq>=0).*0 + (obj.mnuSq<0) .* (obj.PS_Wein93Coeff .* sqrt(-obj.mnuSq)); % Wei93
            obj.PhaseSpace = zeros(obj.nTe,1);
            
            % NO FSD
            if ((isequal(obj.TTFSD,'OFF')) ...
                    && (isequal(obj.DTFSD,'OFF')) ...
                    && (isequal(obj.HTFSD,'OFF')))
                
                if obj.mnuSq>=0
                    obj.PhaseSpace = ((obj.Q-obj.Te)>=0).*...
                        (obj.pe.*(obj.Te+obj.me).*(obj.Q-obj.Te).*((1-obj.sin2T4)...
                        .*(((obj.Q-obj.Te).^2-obj.mnuSq) >= 0)...
                        .*((obj.Q-obj.Te).^2-obj.mnuSq).^.5...
                        + obj.sin2T4.*(((obj.Q-obj.Te).^2-obj.mnu4Sq)>=0)...
                        .*(obj.Q-obj.Te).^2-obj.mnu4Sq).^.5);
                else
                    obj.PhaseSpace = ((obj.Q-obj.Te)>=0).*...
                        (obj.pe.*(obj.Te+obj.me).*(obj.Q - obj.Te + m .* exp(-(obj.Q-obj.Te)./m-1)).*((1-obj.sin2T4)...
                        .*(((obj.Q-obj.Te).^2-obj.mnuSq) >= 0)...
                        .*((obj.Q-obj.Te).^2-obj.mnuSq).^.5...
                        + obj.sin2T4.*(((obj.Q-obj.Te).^2-obj.mnu4Sq)>=0)...
                        .*(obj.Q-obj.Te).^2-obj.mnu4Sq).^.5);
                end
                obj.PhaseSpace(isnan(obj.PhaseSpace)) = 0;
                
                % good
                %                 obj.PhaseSpace = ((obj.Q-obj.Te)>=0).*...
                %                     ((obj.mnuSq>=0).*...
                %                     (obj.pe.*(obj.Te+obj.me).*(obj.Q-obj.Te).*((1-obj.sin2T4)...
                %                     .*(((obj.Q-obj.Te).^2-obj.mnuSq) >= 0)...
                %                     .*((obj.Q-obj.Te).^2-obj.mnuSq).^.5...
                %                     + obj.sin2T4.*(((obj.Q-obj.Te).^2-obj.mnu4Sq)>=0)...
                %                     .*(obj.Q-obj.Te).^2-obj.mnu4Sq).^.5)...
                %                     + (obj.mnuSq<0).*...
                %                     (obj.pe.*(obj.Te+obj.me).*(obj.Q - obj.Te + m .* exp(-(obj.Q-obj.Te)./m-1)).*((1-obj.sin2T4)...
                %                     .*(((obj.Q-obj.Te).^2-obj.mnuSq) >= 0)...
                %                     .*((obj.Q-obj.Te).^2-obj.mnuSq).^.5...
                %                     + obj.sin2T4.*(((obj.Q-obj.Te).^2-obj.mnu4Sq)>=0)...
                %                     .*(obj.Q-obj.Te).^2-obj.mnu4Sq).^.5));
                
            else
                switch obj.TTFSD
                    case {'DOSS','DOSSNOEE','SAENZ','SAENZNOEE','Sibille','Sibille0p5eV','SibilleFull'}
                        if ( (strcmp(obj.TTFSD,'SAENZNOEE')) || (strcmp(obj.TTFSD,'DOSSNOEE')) )
                            obj.TTNormES = 0;
                        end
                        
                        % Normalization Factors Ground/Excited States
                        %obj.TTexE_G  = obj.TTexE(1:obj.TTGSTh);
                        %obj.TTexP_G  = obj.TTexP(1:obj.TTGSTh)./sum(obj.TTexP(1:obj.TTGSTh)).*obj.TTNormGS;
                        
                        %obj.TTexE_E  = obj.TTexE(obj.TTGSTh:end);
                        %obj.TTexP_E  = obj.TTexP(obj.TTGSTh:end)./sum(obj.TTexP(obj.TTGSTh:end)).*obj.TTNormES;
                        
                        obj.TTexE_G  = obj.TTexE(1:obj.TTGSTh);
                        obj.TTNormGS = sum(obj.TTexP(1:obj.TTGSTh));
                        obj.TTexP_G  = obj.TTexP(1:obj.TTGSTh)/obj.TTNormGS;
                        
                        obj.TTexE_E  = obj.TTexE(obj.TTGSTh:end);
                        obj.TTNormES = sum(obj.TTexP(obj.TTGSTh:end));
                        obj.TTexP_E  = obj.TTexP(obj.TTGSTh:end)/obj.TTNormES;
                        
                        
                        
                        % Ground State
                        clear TTexP_G_M; TTexP_G_M = repmat(obj.TTexP_G,obj.nTe,1);
                        clear TTexE_G_M; TTexE_G_M = repmat(obj.TTexE_G,obj.nTe,1);
                        clear Q_M; Q_M  = repmat(repmat(obj.Q,obj.nTe,1),1,numel(obj.TTexP_G));
                        clear Te_M; Te_M = repmat(obj.Te,1,numel(obj.TTexP_G));
                        clear pe_M; pe_M = repmat(obj.pe,1,numel(obj.TTexP_G));
                        clear me_M; me_M = repmat(repmat(obj.me,obj.nTe,1),1,numel(obj.TTexP_G));
                        clear sin2T4_M; sin2T4_M = repmat(repmat(obj.sin2T4,obj.nTe,1),1,numel(obj.TTexP_G));
                        clear mnuSq_M; mnuSq_M  = repmat(repmat(obj.mnuSq,obj.nTe,1),1,numel(obj.TTexP_G));
                        clear mnu4Sq_M; mnu4Sq_M = repmat(repmat(obj.mnu4Sq,obj.nTe,1),1,numel(obj.TTexP_G));
                        clear m_M; m_M = repmat(repmat(m,obj.nTe,1),1,numel(obj.TTexP_G));
                        
                        clear GS;
                        if obj.mnuSq>=0
                            GS = 0 + sum(...
                                ((Q_M-Te_M-TTexE_G_M)>=0)...
                                .*pe_M.*(Te_M+me_M) ...
                                .*(Q_M-TTexE_G_M-Te_M)...
                                .*(((ones(obj.nTe,numel(obj.TTexP_G))-sin2T4_M)....
                                .*(((Q_M-Te_M-TTexE_G_M).^2-mnuSq_M)>=0).*((Q_M-TTexE_G_M-Te_M).^2-mnuSq_M).^.5)...
                                + sin2T4_M.*(((Q_M-Te_M-TTexE_G_M).^2-mnu4Sq_M)>=0)...
                                .*((Q_M-TTexE_G_M-Te_M).^2-mnu4Sq_M).^.5)...
                                .*TTexP_G_M,2);
                        else
                            GS = 0 + sum(...
                                ((Q_M-Te_M-TTexE_G_M)>=0)...
                                .*pe_M.*(Te_M+me_M) ...
                                .*(Q_M-TTexE_G_M-Te_M + m_M .* exp(-(Q_M-Te_M)./m_M - ones(obj.nTe,numel(obj.TTexP_G)) ))...
                                .*(((ones(obj.nTe,numel(obj.TTexP_G))-sin2T4_M)....
                                .*(((Q_M-Te_M-TTexE_G_M).^2-mnuSq_M)>=0).*((Q_M-TTexE_G_M-Te_M).^2-mnuSq_M).^.5)...
                                + sin2T4_M.*(((Q_M-Te_M-TTexE_G_M).^2-mnu4Sq_M)>=0)...
                                .*((Q_M-TTexE_G_M-Te_M).^2-mnu4Sq_M).^.5)...
                                .*TTexP_G_M,2);
                        end
                        GS(isnan(GS)) = 0;
                        %-TTexE_G_M
                        
                        % Exited States
                        clear TTexP_E_M; TTexP_E_M = repmat(obj.TTexP_E,obj.nTe,1);
                        clear TTexE_E_M; TTexE_E_M = repmat(obj.TTexE_E,obj.nTe,1);
                        clear Q_M; Q_M  = repmat(repmat(obj.Q,obj.nTe,1),1,numel(obj.TTexP_E));
                        clear Te_M; Te_M = repmat(obj.Te,1,numel(obj.TTexP_E));
                        clear pe_M; pe_M = repmat(obj.pe,1,numel(obj.TTexP_E));
                        clear me_M; me_M = repmat(repmat(obj.me,obj.nTe,1),1,numel(obj.TTexP_E));
                        clear sin2T4_M; sin2T4_M = repmat(repmat(obj.sin2T4,obj.nTe,1),1,numel(obj.TTexP_E));
                        clear mnuSq_M; mnuSq_M  = repmat(repmat(obj.mnuSq,obj.nTe,1),1,numel(obj.TTexP_E));
                        clear mnu4Sq_M; mnu4Sq_M = repmat(repmat(obj.mnu4Sq,obj.nTe,1),1,numel(obj.TTexP_E));
                        clear m_M; m_M = repmat(repmat(m,obj.nTe,1),1,numel(obj.TTexP_E));
                        
                        %                         clear ES; ES = 0 + sum(...
                        %                             ((Q_M-Te_M-TTexE_E_M)>=0)...
                        %                             .*pe_M.*(Te_M+me_M) ...
                        %                             .*(Q_M-TTexE_E_M-Te_M)...
                        %                             .*(((ones(obj.nTe,numel(obj.TTexP_E))-sin2T4_M)....
                        %                             .*(((Q_M-Te_M-TTexE_E_M).^2-mnuSq_M)>=0).*((Q_M-TTexE_E_M-Te_M).^2-mnuSq_M).^.5)...
                        %                             + sin2T4_M.*(((Q_M-Te_M-TTexE_E_M).^2-mnu4Sq_M)>=0)...
                        %                             .*((Q_M-TTexE_E_M-Te_M).^2-mnu4Sq_M).^.5)...
                        %                             .*TTexP_E_M,2);
                        
                        clear ES;
                        if obj.mnuSq>=0
                            ES = 0 + sum(...
                                ((Q_M-Te_M-TTexE_E_M)>=0)...
                                .*pe_M.*(Te_M+me_M) ...
                                .*(Q_M-TTexE_E_M-Te_M)...
                                .*(((ones(obj.nTe,numel(obj.TTexP_E))-sin2T4_M)....
                                .*(((Q_M-Te_M-TTexE_E_M).^2-mnuSq_M)>=0).*((Q_M-TTexE_E_M-Te_M).^2-mnuSq_M).^.5)...
                                + sin2T4_M.*(((Q_M-Te_M-TTexE_E_M).^2-mnu4Sq_M)>=0)...
                                .*((Q_M-TTexE_E_M-Te_M).^2-mnu4Sq_M).^.5)...
                                .*TTexP_E_M,2);
                        else
                            ES = 0 + sum(...
                                ((Q_M-Te_M-TTexE_E_M)>=0)...
                                .*pe_M.*(Te_M+me_M) ...
                                .*(Q_M-TTexE_E_M-Te_M + m_M .* exp(-(Q_M-Te_M)./m_M - ones(obj.nTe,numel(obj.TTexP_E)) ))...
                                .*(((ones(obj.nTe,numel(obj.TTexP_E))-sin2T4_M)....
                                .*(((Q_M-Te_M-TTexE_E_M).^2-mnuSq_M)>=0).*((Q_M-TTexE_E_M-Te_M).^2-mnuSq_M).^.5)...
                                + sin2T4_M.*(((Q_M-Te_M-TTexE_E_M).^2-mnu4Sq_M)>=0)...
                                .*((Q_M-TTexE_E_M-Te_M).^2-mnu4Sq_M).^.5)...
                                .*TTexP_E_M,2);
                        end
                        ES(isnan(ES)) = 0;
                        %-TTexE_E_M
                        % Build Ground + Exited States
                        TTps = (obj.TTNormGS ).* GS  + (obj.TTNormES ).* ES;
                        obj.PhaseSpace = obj.PhaseSpace + obj.WGTS_MolFrac_TT.*TTps;
                end
                
                switch obj.DTFSD
                    case {'DOSS','HTFSD','TTFSD'}
                        % Ground State
                        [GS,obj.TTexE_G,obj.TTexP_G] = ...
                            obj.ComputeFSD_GSES(obj.TTexE,obj.TTexP,obj.TTGSTh,'ground');
                        
                        % Exited States
                        [ES,obj.TTexE_E,obj.TTexP_E] = ...
                            obj.ComputeFSD_GSES(obj.TTexE,obj.TTexP,obj.TTGSTh,'excited');
                        
                        % Build Ground + Exited States
                        DTps = obj.DTNormGS .* GS + obj.DTNormES .* ES;
                        obj.PhaseSpace = obj.PhaseSpace + obj.WGTS_MolFrac_DT.*DTps;
                end
                
                switch obj.HTFSD
                    case {'SAENZ','Sibille','Sibille0p5eV','SibilleFull'}
                        % Ground State
                        [GS,obj.HTexE_G,obj.HTexP_G] = ...
                            obj.ComputeFSD_GSES(obj.HTexE,obj.HTexP,obj.HTGSTh,'ground');
                        
                        % Exited States
                        [ES,obj.HTexE_E,obj.HTexP_E] = ...
                            obj.ComputeFSD_GSES(obj.HTexE,obj.HTexP,obj.HTGSTh,'excited');
                        
                        % Build Ground + Exited States
                        HTps = obj.HTNormGS .* GS + obj.HTNormES .* ES;
                        obj.PhaseSpace = obj.PhaseSpace + obj.WGTS_MolFrac_HT.*HTps;
                end
                
            end
        end
                
        function f      = ComputeFermiCorr(obj)
            % ----------------------------------------------------------
            %  Fermi Function
            %  Relativistic / Nuclear Physics A526 (1991) 131-140)
            % ----------------------------------------------------------
            
            f = 4.*(2.*obj.dhWp.*obj.Rad).^(-2*(1-obj.dhWgamma)) ...
                .* abs(gammac(obj.dhWgamma + 1i.*obj.dhWy)).^2 ...
                .* gammac(2.*obj.dhWgamma + 1).^-2 ...
                .* exp(pi.*obj.dhWy);         
        end
        
        function l0     = ComputeFiniteExtChargeCorr(obj)
            % ----------------------------------------------------------
            % Finite Extention of Nucleus Charge, L0
            % Ref: Nuclear Physics A526 (1991) 131-140
            % Ref: Phys.Rev.C84:024617,2011
            % ----------------------------------------------------------
            l0 = (13/60.*(obj.fsC*obj.Z)^2 - ...
                obj.dhWw.*obj.Rad.*obj.fsC.*obj.Z.*(41-26.*obj.dhWgamma) ./ ...
                (15*(2*obj.dhWgamma-1)) - ...
                obj.fsC.*obj.Z.*obj.Rad.*obj.dhWgamma*(17-2.*obj.dhWgamma) ./ ...
                (30.*obj.dhWw.*(2.*obj.dhWgamma-1)));
            l0 = 1 + (1+obj.L0CorrBias) .* l0;
        end
        
        function rc     = ComputeRadiativeCorr(obj)
            % ----------------------------------------------------------
            % Radiative Corrections
            % Case 0 : None
            % Case 1 : Repko Approximation, Phys. Rev. C28 (1983) 2433
            % Case 2 : Repko Generic, Phys. Rev. C28 (1983) 2433
            % Case 3 : Generic, Phys. Rev. 164, 1767 (1967)
            % ----------------------------------------------------------
            switch obj.RadType
                case 1 % Approximation W.W. Repko and C. Wu, Phys. Rev. C28 (1983) 2433
                    Delta = obj.me + obj.Q;
                    tbeta =  0.5./obj.beta .* log((1+obj.beta)./(1-obj.beta))-1;
                    rc =  2*obj.fsC/pi * ...
                        ( tbeta.*(log(2)-1.5+ Delta./obj.Ee -1 ) + ...
                        0.25*(tbeta+1).*(2*(1+obj.beta.^2) + 2*log(1-obj.beta) + 1/6*(Delta./obj.Ee-1).^2) ...
                        - 2 + 0.5*obj.beta -17/36*(obj.beta).^2 + 5/6*(obj.beta).^3);
                case 2 % General Formula  W.W. Repko and C. Wu, Phys. Rev. C28 (1983) 2433
                    Delta = obj.me + obj.Q;
                    tbeta =  0.5./obj.beta .* log((1+obj.beta)./(1-obj.beta))-1;
                    Lbeta1 = zeros(obj.nTe,1); Lbeta2 = zeros(obj.nTe,1);
                    Lbeta3 = zeros(obj.nTe,1); Lbeta4 = zeros(obj.nTe,1);
                    Lbeta5 = zeros(obj.nTe,1);
                    for i=1:1:obj.nTe
                        Lbeta1(i) = integral(@(t)log(1-t)./t,0,obj.beta(i));% L(beta)
                        Lbeta2(i) = integral(@(t)log(1-t)./t,0,-obj.beta(i));% L(-beta)
                        Lbeta3(i) = integral(@(t)log(1-t)./t,0,2*obj.beta(i)/(1+obj.beta(i)));% L(2beta/(1+beta))
                        Lbeta4(i) = integral(@(t)log(1-t)./t,0,(1-obj.beta(i))/2);% L((1-beta)/2)
                        Lbeta5(i) = integral(@(t)log(1-t)./t,0,(1+obj.beta(i))/2);% L((1+beta)/2)
                    end
                    rc =  2*obj.fsC/pi * ...
                        ( tbeta.*(log(2)-1.5+ Delta./obj.Ee -1 ) + ...
                        0.25*(tbeta+1).*(2*(1+obj.beta.^2) - 2*log(2./(1-obj.beta)) + 1/6*(Delta./obj.Ee-1).^2) ...
                        + 0.5./obj.beta .*(Lbeta1 - Lbeta2 + Lbeta3 + 0.5*Lbeta4 -0.5*Lbeta5));
                case 3 % A. Sirlin, Phys. Rev. 164, 1767 (Dec 1967) --> complex number at endpoint...
                    obj.Ee0 = obj.Q + obj.me; Lbeta3 = zeros(obj.nTe,1);
                    for i=1:1:obj.nTe
                        Lbeta3(i) = integral(@(t)log(1-t)./t,0,2*obj.beta(i)/(1+obj.beta(i)));% L(2beta/(1+beta))
                    end
                    rc =  obj.fsC./(2*pi) * (3*log(obj.Mp/obj.me) -3/4 + ...
                        4.*((atanh(obj.beta))./obj.beta - 1) .* ...
                        ((obj.Ee0-obj.Ee)./(3*obj.Ee) -1.5 + log(2*(obj.Ee0-obj.Ee)./obj.me)) ...
                        + 4./obj.beta .* Lbeta3 + ...
                        atanh(obj.beta)./obj.beta .* (2.*(1+obj.beta.^2) + (obj.Ee0-obj.Ee).^2./(6.*obj.Ee.^2) - 4.*atanh(obj.beta)));
            end
            switch obj.RadType
                case {1,2}
                    rc =  ...
                        (Delta./obj.me - obj.Ee./obj.me).^(2*obj.fsC/pi.*tbeta) .* (1 + (1+obj.RadCorrBias) * rc);
                case 3
                    rc = (1 + (1+obj.RadCorrBias) * rc);
            end
        end
        
        function wfs      = ComputeWintFiniteSizeCorr(obj)
            % ----------------------------------------------------------
            % Weak interaction finite size correction, C
            % Ref: Nuclear Physics A526 (1991) 131-140
            % Ref: Phys.Rev.C84:024617,2011
            % ----------------------------------------------------------
            C0  = -233/630*(obj.fsC*obj.Z)^2 - 1/5*(obj.dhWw0*obj.Rad)^2 + 2/35*obj.dhWw0*obj.Rad*obj.fsC*obj.Z;
            C1  = -21/35*obj.Rad*obj.fsC*obj.Z + 4/9*obj.dhWw0*obj.Rad^2;
            C2  = -4/9*obj.Rad^2;
            wfs = (C0 + C1.*obj.dhWw + C2.*obj.dhWw.^2);
            wfs = 1 + (1+obj.CCorrBias).* wfs;
        end
        
        function ee     = ComputeEEexchangeCorr(obj,type)
            % ----------------------------------------------------------
            % Electron-Electron Exchange Term
            % Ref: Nuclear Physics A526 (1991) 131-140
            % ----------------------------------------------------------
            a0 = 1/obj.fsC; %  x2/2, since 1/2 bohr radius & Z=2
            eta = -2 ./ (obj.dhWp .*a0);
            alphaeta = eta.^4 .* exp(2*eta.*atan(-2./eta))./(1+0.25*eta.^2).^2;
            switch type
                case 1 % Wilkinson NPA 526 (1991)
                    ee = (729/256*alphaeta.^2 + 27/16*alphaeta);
                case 2 % New version from Thierry included the first 10 exited states of He+, From Hax85
                    ee = (2.462*alphaeta.^2 + 0.905*alphaeta);
            end
            ee = 1 + (1+obj.ECorrBias) .* ee;
        end
        
        function rc     = ComputeRecoilCoulombCorr(obj)
            % ----------------------------------------------------------
            % Recoil Coulomb Effect (electron moving in a moving field
            % due to the displacement of the recoil daughter nucleus)
            % Ref: Nuclear Physics A526 (1991) 131-140
            % ----------------------------------------------------------
            lambdaT = 1.265; % Physics of Atomic Nuclei, Vol. 65, No. 10, 2002, pp. 1795-1797
            B = (1-lambdaT^2)/(1+3*lambdaT^2);
            rc  = (pi*obj.fsC*obj.Z)/(obj.MHe3/obj.me)...
                ./obj.dhWp.*(1+B.*(obj.dhWw0-obj.dhWw)./(3*obj.dhWw));
            rc = 1 - (1+obj.QCorrBias).* rc;
        end
        
        function he3rec = ComputeRecoilWmVmACorr(obj)
            % ----------------------------------------------------------
            % 3He Recoil Corrections:
            % - Modification of the phase space factor 2body->3body
            % - Weak Magnetism
            % - V-A currents interference
            % Ref: Nuclear Physics A526 (1991) 131-140
            % ----------------------------------------------------------
            lambdaT = 1.265; % Physics of Atomic Nuclei, Vol. 65, No. 10, 2002, pp. 1795-1797
            mu = 5.106588;   % NIST
            aR = 2*(5*lambdaT^2 + 2*lambdaT*mu+1)/(obj.MHe3/obj.me);
            bR = 2*lambdaT*(lambdaT+mu)/(obj.MHe3/obj.me);
            cR = 1 + 3*lambdaT^2 - bR*obj.dhWw0;
            he3rec = (aR.*obj.dhWw - bR./obj.dhWw)/cR;
            he3rec = 1+ (1+obj.RecCorrBias).* he3rec;
        end

        function          ComputeTBDDS(obj,varargin)
            % ----------------------------------------------------------
            % Compute Tritium Decay Differential Spectrum
            % Phase Space
            % Fermi
            % Isotropologs T-T / D-T / H-T
            % FSD (ro-vib states, electronic excitation)
            % KATRIN-like Normalization
            % ----------------------------------------------------------
            
            % Inputs
            % Fit Parameters are parsed here as biases
            % In SetFitBias(obj,1) these biases are added to parameter's inital values
            p = inputParser;
            p.addParameter('mSq_bias',0.,@(x)isfloat(x));
            p.addParameter('E0_bias',0.,@(x)isfloat(x));
            p.addParameter('N_bias',0.,@(x)isfloat(x));
            p.addParameter('DE_bias',0.,@(x)isfloat(x));
            p.addParameter('B_bias',0.,@(x)isfloat(x));
            p.addParameter('BSlope_bias',0,@(x)isfloat(x));
            p.addParameter('TTGS_bias',0.,@(x)isfloat(x));
            p.addParameter('TTES_bias',0.,@(x)isfloat(x));
            p.addParameter('DTGS_bias',0.,@(x)isfloat(x));
            p.addParameter('DTES_bias',0.,@(x)isfloat(x));
            p.addParameter('HTGS_bias',0.,@(x)isfloat(x));
            p.addParameter('HTES_bias',0.,@(x)isfloat(x));
            p.addParameter('mnu4Sq_Bias',0.,@(x)isfloat(x));
            p.addParameter('sin2T4_Bias',0.,@(x)isfloat(x));
            p.addParameter('qUOffset_bias',zeros(1,obj.nPixels),@(x)isfloat(x));  
            p.addParameter('mTSq_bias',zeros(1,obj.nPixels),@(x)isfloat(x));  
            p.addParameter('NormFlag','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('FracTm_bias',0,@(x)isfloat(x)); %fraction of T- ions
            p.parse(varargin{:});
            obj.mnuSq_Bias       = p.Results.mSq_bias;
            obj.E0_bias          = p.Results.E0_bias;
            obj.Norm_Bias        = p.Results.N_bias;
            obj.DE_Bias          = p.Results.DE_bias;
            obj.Bkg_Bias         = p.Results.B_bias;
            obj.BkgSlope_Bias    = p.Results.BSlope_bias;
            obj.TTGS_bias        = p.Results.TTGS_bias;
            obj.TTES_bias        = p.Results.TTES_bias;
            obj.DTGS_bias        = p.Results.DTGS_bias;
            obj.DTES_bias        = p.Results.DTES_bias;
            obj.HTGS_bias        = p.Results.HTGS_bias;
            obj.HTES_bias        = p.Results.HTES_bias;            
            obj.mnu4Sq_Bias      = p.Results.mnu4Sq_Bias;
            obj.sin2T4_Bias      = p.Results.sin2T4_Bias;
            obj.qUOffset_bias    = p.Results.qUOffset_bias;
            obj.FracTm_bias      = p.Results.FracTm_bias;
            obj.mTSq_bias = p.Results.mTSq_bias;
            NormFlag             = p.Results.NormFlag;
                        
            % Save previous values of parameters
            if isempty(obj.mnuSq) 
                % if initialization
                mnuSq_prev = obj.mnuSq_i;
                E0_prev = obj.Q_i;
                TTGS_prev = obj.TTNormGS_i;
                TTES_prev = obj.TTNormES_i;
                DTGS_prev = obj.DTNormGS_i;
                DTES_prev = obj.DTNormES_i;
                HTGS_prev = obj.HTNormGS_i;
                HTES_prev = obj.HTNormES_i;
                mnu4Sq_prev = obj.mnu4Sq_i;
                sin2T4_prev = obj.sin2T4_i;
                mTSq_prev = obj.mTSq_i;
                FracTm_prev = obj.WGTS_MolFrac_Tm_i;
            else             
                mnuSq_prev = obj.mnuSq;
                E0_prev = obj.Q;
                TTGS_prev = obj.TTNormGS;
                TTES_prev = obj.TTNormES;
                DTGS_prev = obj.DTNormGS;
                DTES_prev = obj.DTNormES;
                HTGS_prev = obj.HTNormGS;
                HTES_prev = obj.HTNormES;
                mnu4Sq_prev = obj.mnu4Sq;
                sin2T4_prev = obj.sin2T4;
                mTSq_prev = obj.mTSq;
                FracTm_prev = obj.WGTS_MolFrac_Tm;
            end
            
            % Initialization (add bias)
            SetFitBias(obj,1)
            
            % If all biases are 0, it means it is the inizializtion so
            % do the simulation of the spectrum (calculate TBDDS)
            Init_condition = (obj.mnuSq_Bias == 0) && ...
                (obj.E0_bias == 0) && all(obj.Norm_Bias == 0) && ...
                all(obj.Bkg_Bias == 0) && all(obj.qUOffset_bias == 0) && ...
                all(obj.BkgSlope_Bias == 0) && all(obj.mTSq_bias == 0) && (obj.FracTm_bias==0);
            
            % Check if parameters changed, which require recalculation of phace space
            PhaseSpaceChange_condition = ...
                (mnuSq_prev ~= obj.mnuSq) || ... % neutrino mass
                (E0_prev ~= obj.Q) || ...        % endpoint
                any(mTSq_prev ~= obj.mTSq) || ... % neutrino mass offset (TBDDS energy smearing)
                ((~isempty(obj.TTNormGS)) && any(TTGS_prev ~= obj.TTNormGS)) || ... %FSDs
                ((~isempty(obj.TTNormES)) && any(TTES_prev ~= obj.TTNormES)) || ...
                ((~isempty(obj.DTNormGS)) && any(DTGS_prev ~= obj.DTNormGS)) || ...
                ((~isempty(obj.DTNormES)) && any(DTES_prev ~= obj.DTNormES)) || ...
                ((~isempty(obj.HTNormGS)) && any(HTGS_prev ~= obj.DTNormGS)) || ...
                ((~isempty(obj.HTNormES)) && any(HTES_prev ~= obj.DTNormES)) || ...
                (mnu4Sq_prev ~= obj.mnu4Sq) || (sin2T4_prev ~= obj.sin2T4) || ...; % steriles
                (FracTm_prev ~= obj.WGTS_MolFrac_Tm);
             
            if Init_condition || PhaseSpaceChange_condition
                % (re)calculate phase space
                
                switch obj.DopplerEffectFlag
                    case {'matConv','numConv'}
                        % Enlargement to avoid convolution edge effects.
                        obj.enlargeSpectrum();
                        % Already precalculated in class, use only if parameters change.
                        if strcmp(obj.recomputeKernel,'ON') || strcmp(obj.DopplerEffectFlag,'numConv')
                            obj.computeKernel();
                        end
                end
                
                if FracTm_prev ~= obj.WGTS_MolFrac_Tm % if molecular fraction of T- changed
                    obj.ComputeTritiumPurity;    %(re-)calculate: tritium purity
                    obj.ComputeNormFactorTBDDS;  %(re-)calculate: normalization factor
                end
              
                % Phase Space - Negative Mass Squared....
                switch obj.PS_Wein93
                    case 'ON'
                        obj.ComputePhaseSpaceWei93();
                    case 'OFF'
                        obj.ComputePhaseSpace();
                end
                
                % Apply Fermi Function
                obj.TBDDS =  obj.PhaseSpace .* obj.ComputeFermiCorr();
                
                % Corrections to the allowed diff. beta spectrum
                switch obj.RadiativeFlag
                    case 'ON'
                        obj.TBDDS = obj.TBDDS .* real(obj.ComputeRadiativeCorr());
                end
                switch obj.RecoilWmVmAFlag
                    case 'ON'
                        obj.TBDDS = obj.TBDDS .* obj.ComputeRecoilWmVmACorr() ;
                end
                switch obj.FiniteExtChargeFlag
                    case 'ON'
                        obj.TBDDS = obj.TBDDS .* obj.ComputeFiniteExtChargeCorr() ;
                end
                switch obj.EEexchangeFlag
                    case 'ON'
                        obj.TBDDS = obj.TBDDS .* obj.ComputeEEexchangeCorr(2) ;
                end
                switch obj.ScreeningFlag
                    case 'ON'
                        obj.TBDDS = obj.TBDDS .* obj.ComputeScreeningCorr();
                end
                switch obj.WintFiniteSizeFlag
                    case 'ON'
                        obj.TBDDS = obj.TBDDS .* obj.ComputeWintFiniteSizeCorr() ;
                end
                switch obj.RecoilCoulombFlag
                    case 'ON'
                        obj.TBDDS = obj.TBDDS .* obj.ComputeRecoilCoulombCorr() ;
                end
                
                switch obj.DopplerEffectFlag
                    case {'matConv','numConv'}
                        obj.TBDDS = obj.KernelSpectrumConv(obj.TBDDS);
                        obj.TBDDS = obj.restoreSpectrum(obj.TBDDS);
                end
            end
            
            if strcmp(NormFlag,'ON')
                obj.TBDDS  = (1+obj.normFit).*(obj.TBDDS./simpsons(obj.Te,obj.TBDDS)).*obj.NormFactorTBDDS;
            end
            
            if any(obj.qUOffset~=0) 
                % shift differential spectrum in case qUOffset changed in fit
                TBDDStmp =  obj.TBDDS; % init
                for r=1:obj.nRings
                    if obj.qUOffset(r)~=0
                        %fprintf('Ring %.0f: qUOffset = %.2g eV \n',r,obj.qUOffset(r))
                        TBDDStmp(:,r) = interp1(obj.Te,obj.TBDDS(:,r),obj.Te+obj.qUOffset(r),'spline','extrap'); 
                    end
                end
                obj.TBDDS = TBDDStmp;
            end
            
            obj.TBDDSE = obj.TBDDS.^0.5;
            obj.TBDDSE(obj.TBDDSE == 0) = 1;
        end
        
        function          ComputeTBDIS(obj,varargin)
                switch obj.IStype
                    case 'GaussKronrod'
                        for qUi = 1:1:obj.nqU
                            TBDDSf = @(e) interp1(obj.Te,obj.TBDDS.*obj.KTF(obj.Te,obj.qU(qUi)),e);
                            obj.TBDIS(qUi) = ...
                                (quadgk(TBDDSf,obj.qU(qUi),obj.qUmax)+obj.BKG_RateSec).*obj.TimeSec.*...
                                obj.qUfrac(qUi);
                        end
                    case 'TRAPZFAST'
                        obj.TBDIS = ((trapz(obj.Te,obj.TBDDS.*obj.RF(:,:)',1)+obj.BKG_RateSec).*obj.TimeSec.*obj.qUfrac')';
                    case 'SIMPFAST' %Same as TRAPZFAST but with simpsons instead of trapz
                        if strcmp(obj.FPD_Segmentation,'RING')
                            TBDDSandRF = permute(repmat(obj.TBDDS,1,1,obj.nqU),[1,3,2]).*obj.RF;
                        else
                            TBDDSandRF = repmat(obj.TBDDS,1,obj.nqU).*obj.RF;
                        end
                        %                    if strcmp(obj.FPD_Segmentation,'OFF')
                        TBDISwithoutBCK = squeeze(simpsons(obj.Te,TBDDSandRF));
                        if isrow(TBDISwithoutBCK); TBDISwithoutBCK = TBDISwithoutBCK'; end
                        BKG =  (obj.BKG_Slope.*(obj.qU-18574))+repmat(abs(obj.BKG_RateSec),obj.nqU,1);
                        obj.TBDIS = (TBDISwithoutBCK + BKG).*obj.TimeSec.*obj.qUfrac;
                        
                        % CAUTION: Absolute background above! Beware!
                        
                        %                     else
                        %                         TBDISwithoutBCK = (obj.TBDIS - repmat(obj.prevBKG,obj.nqU,1)).*obj.TimeSec.*obj.qUfrac;
                        %                         obj.TBDIS = (TBDISwithoutBCK + repmat(obj.BKG_RateSec,obj.nqU,1)).*obj.TimeSec.*obj.qUfrac;
                        %                         obj.onlyBKGchange = false;
                        %                    end
                    case 'SUM' % Slightly faster than TRAPZFAST, but only works for equally-spaced Te values
                        obj.TBDIS = ((obj.TeStep*(sum(obj.TBDDS.*obj.RF(:,:)'))+obj.BKG_RateSec).*obj.TimeSec.*obj.qUfrac')';
                end
            
        if strcmp(obj.FPD_PileUpEff,'ON') && strcmp(obj.FPD_ROIEff,'ON')
                switch obj.FPD_Segmentation
                    case 'OFF'
                        RatePerPixel = obj.TBDIS./obj.qUfrac./obj.TimeSec/numel(obj.FPD_PixList);
                    case {'SINGLEPIXEL','MULTIPIXEL'}
                        RatePerPixel = obj.TBDIS./obj.qUfrac./obj.TimeSec;
                    case 'RING'
                        if obj.FPD_Ring == 1
                            RatePerPixel = obj.TBDIS./obj.qUfrac./obj.TimeSec./4;
                        else
                            RatePerPixel = obj.TBDIS./obj.qUfrac./obj.TimeSec./12;
                        end
                end
                RoiCorr = obj.qUEfficiencyCorrectionFactor(obj.qU);
                PileUpCorr = obj.PileUpEfficiencyCorrectionFactor(RatePerPixel);
                obj.TBDIS = obj.TBDIS.*RoiCorr.*PileUpCorr;
        end
        
%         %  apply qU Offset (1 per FPD segmentation e.g. 1 per ring)
%         if any(obj.qUOffset~=0)
%             TBDIStmp = zeros(obj.nqU,obj.nPixels);
%             if obj.nPixels==1
%                 TBDIStmp = interp1(obj.qU,obj.TBDIS,obj.qU+obj.qUOffset,'spline','extrap');
%             else
%                 for p=1:obj.nPixels
%                     TBDIStmp(:,p) = interp1(obj.qU(:,p),obj.TBDIS(:,p),obj.qU(:,p)+obj.qUOffset(p),'spline','extrap');
%                 end
%             end
%             obj.TBDIS = TBDIStmp;
%         end
        
        obj.TBDISE = sqrt(obj.TBDIS);
        end
        
        function f      = ComputeTBDISf(obj,qu,mb,qb,bb,nb,gsb,esb)
            % Function for Matlab Fit
            
            ComputeTBDDS(obj,...
                'mSq_bias',mb,...
                'E0_bias',qb,...
                'N_bias',nb,...
                'B_bias',bb,...
                'TTGS_bias',gsb,...
                'TTES_bias',esb);
            
            obj.TBDDS = (1+nb).*obj.TBDDS;
            ComputeTBDIS(obj);
            
            qUieq = []; qUieq = (find(abs(obj.qU-qu)<1e-3));
            if qUieq>0
                f = obj.TBDIS(qUieq);
            else
                f = -1;
                return;
            end
        end
        
        function f      = ComputeTBDIS4f(obj,qu,mnuSqb,qb,bb,nb,gsb,esb,mnu4Sqb,sin2T4b)
            % Function for Matlab Fit
            ComputeTBDDS(obj,...
                'mSq_bias',mnuSqb,...
                'E0_bias',qb,...
                'N_bias',nb,...
                'B_bias',bb,...
                'TTGS_bias',gsb,...
                'TTES_bias',esb,...
                'mnu4Sq_Bias',mnu4Sqb,...
                'sin2T4_Bias',sin2T4b );
            
            obj.TBDDS = (1+nb).*obj.TBDDS;
            f = ComputeTBDIS(obj);
            
            qUieq = []; qUieq = (find(abs(obj.qU-qu)<1e-3));
            if qUieq>0
                f = obj.TBDIS(qUieq);
            else
                f = -1;
                return;
            end
        end
                
        function          ComputeFracTBDtail(obj)
            % Compute Fraction of decay in TBDDS Tail
            % Check if FSD are included in the spectrum to decide on the
            % filename.
            if contains(obj.TTFSD,'Blinding')
               FSDlabel = 'Blinding';
            else
               FSDlabel = '';
            end
            
            % Check if integral of full spectrum under these conditions
            % already exists, if yes, load it
            FullDSIntDir = [getenv('SamakPath'),'inputs/CumFrac/'];
            MakeDir( FullDSIntDir);
            FullDSIntName = sprintf('%sTT%0.5g_DT%0.5g_HT%0.5g_DE%s_Rad%s_RhoD%0.5g%s.mat',...
                FullDSIntDir,...
                obj.WGTS_MolFrac_TT,obj.WGTS_MolFrac_DT,obj.WGTS_MolFrac_HT,...
                obj.DopplerEffectFlag,obj.RadiativeFlag,obj.WGTS_CD_MolPerCm2,FSDlabel);
            
            if obj.WGTS_MolFrac_Tm~=0
                FullDSIntName = strrep(FullDSIntName,sprintf('HT%0.5g',obj.WGTS_MolFrac_HT),...
                    sprintf('HT%0.5g_Tm%0.5g',obj.WGTS_MolFrac_HT,obj.WGTS_MolFrac_Tm));
            end
            
            if exist(FullDSIntName,'file') == 2
                TempStructFullDS = load(FullDSIntName);
                FullDSInt = TempStructFullDS.FullDSInt;
            else % if not, calculate it
                if obj.Q - obj.qUmin < 50
                    stepsizeCF = 0.1;
                elseif obj.Q - obj.qUmin < 101
                    stepsizeCF = 0.2;
                elseif obj.Q - obj.qUmin < 1001
                    stepsizeCF = 0.5;
                elseif obj.Q - obj.qUmin < 2001
                    stepsizeCF = 1;
                else
                    stepsizeCF = 2;
                end
                
                saveTe = obj.Te;
                
                obj.Te = (0.01:stepsizeCF:obj.qUmax)';
                obj.nTe = length(obj.Te);
                obj.SetKinVariables();
                obj.ComputeTBDDS('NormFlag','OFF');
                
                FullDSInt = simpsons(obj.Te,obj.TBDDS);
                if strcmp(obj.FPD_Segmentation,'RING')
                    FullDSInt = mean(FullDSInt,2); % save FPD average
                end
                save(FullDSIntName,'FullDSInt','-v7.3','-nocompression');
                
                obj.Te = saveTe;
                obj.nTe = length(saveTe);
                obj.SetKinVariables();
            end

            obj.ComputeTBDDS('NormFlag','OFF');
            
            measuredDSInt = simpsons(obj.Te,obj.TBDDS);
            
            obj.CumFrac = measuredDSInt/FullDSInt;
                  
        end
        
        function          ComputeNormFactorTBDDS(obj)
            % Normalize Spectrum according to KATRIN-settings
            switch obj.FPD_Segmentation
                case 'OFF'
%                     P(1)=8.273976836242512e-18;
%                     P(2)=0.141530047843744;
%                     CorrectionRhoD = (P(1)*obj.WGTS_CD_MolPerCm2+P(2));
%                     CorrectionRhoD = CorrectionRhoD./CorrectionRhoD;
                    CorrectionRhoD=1;
                    obj.NormFactorTBDDS = obj.TdecayC ...
                        .*(2*pi*obj.WGTS_FTR_cm^2*obj.WGTS_CD_MolPerCm2.*CorrectionRhoD) ...% mol tritium atoms
                        .*0.5*(1-cos(asin(sqrt(obj.WGTS_B_T./obj.MACE_Bmax_T)))) ...
                        .*(obj.FPD_MeanEff*obj.FPD_Coverage)...
                        .*obj.CumFrac*obj.WGTS_epsT .* numel(obj.FPD_PixList)/148;
%                        .*obj.qUEfficiencyCorrectionFactor(obj.Te)...
%fprintf(2,'\n \n RhoD Fit %.3f \n\n',obj.WGTS_CD_MolPerCm2);
                case 'RING'
                    nPix = cell2mat(cellfun(@(x) numel(x),obj.FPD_RingPixList,'UniformOutput',false)');
                    if strcmp(obj.FPD_RingMerge,'None')
                        nPix(~ismember(obj.FPD_RingPixList))=[];
                    end
                        obj.NormFactorTBDDS = obj.TdecayC ...
                        .*(2*pi*obj.WGTS_FTR_cm^2*obj.WGTS_CD_MolPerCm2) ...% mol tritium atoms
                        .*0.5*(1-cos(asin(sqrt(obj.WGTS_B_T./obj.MACE_Bmax_T)))) ...
                        .*(obj.FPD_MeanEff*obj.FPD_Coverage)...
                        .*obj.CumFrac*obj.WGTS_epsT .* nPix/148;
%                        .*obj.qUEfficiencyCorrectionFactor(obj.Te)...
                case {'SINGLEPIXEL','MULTIPIXEL'}
                    obj.NormFactorTBDDS = (obj.TdecayC ...
                        .*(2*pi*obj.WGTS_FTR_cm^2*obj.WGTS_CD_MolPerCm2allPixels') ...
                        .*0.5.*(1-cos(asin(sqrt(obj.WGTS_B_T./obj.MACE_Bmax_T)))) ...
                        .*(obj.FPD_Eff_pix'.*obj.FPD_Coverage) ...
                        .*obj.CumFrac.*obj.WGTS_epsT);
%                        .*obj.qUEfficiencyCorrectionFactor(obj.Te)...
            end
        end
        
        function          AddStatFluctTBDDS(obj)
            obj.TBDDS   = obj.TBDDS + sqrt(obj.TBDDS).*randn(obj.nTe,1);
            obj.TBDDSE  = obj.TBDDS.^0.5;
        end
        
        function          AddStatFluctTBDIS(obj)
            if obj.TimeSec > 2*60*60*24 
                obj.TBDIS = obj.TBDIS + obj.TBDISE.*randn(obj.nqU,obj.nPixels);
            else 
                obj.TBDIS = poissrnd(obj.TBDIS);
            end
            
            if min(min(obj.TBDIS)) < 0
                error('NEGATIVE COUNTS! STOP!');
            elseif min(min(obj.TBDIS)) < 3
                warning('LESS THAN 3 COUNTS IN AT LEAST ONE BIN! STOP!');
            elseif min(min(obj.TBDIS)) < 50
                warning('LESS THAN 50 COUNT IN AT LEAST ONE BIN! WARNING!');
            end
            
            obj.TBDISE  = sqrt(obj.TBDIS);
        end
        
        function          AddStatSystFluctTBDIS(obj,varargin)
            % Parser
            p = inputParser;
            p.addParameter('CovMat','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SystCorr','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('CM',[]);
            p.parse(varargin{:});
            CovMat     =    p.Results.CovMat;
            SystCorr   =    p.Results.SystCorr;
            CM = p.Results.CM;
            
            switch CovMat
                case 'OFF'
                    obj.TBDIS   = obj.TBDIS + sqrt(obj.TBDIS).*randn(obj.nqU,1);
                case 'ON'
                    
%                     nCMTrials = 10000;
%                     fName = sprintf('../../inputs/CovMat/WGTSMACE_CovMat_%gTrials_%s_Sec_%g.mat',...
%                         nCMTrials,obj.TD,obj.TimeSec);
%                     if exist(fName, 'file') == 2
%                         CM=importdata(fName);
%                     else
%                         A.ComputeTBDISCovarianceMatrix('nTrials',nCMTrials);
%                         CM=importdata(fName);
%                     end
                    
                    switch SystCorr
                        case 'OFF'
                            obj.TBDIS   = obj.TBDIS + ...
                                sqrt(obj.TBDIS + diag(CM.WGTSMACE_CovMat)) ...
                                .*randn(obj.nqU,1);
                        case 'ON'
                            %CNum=chol(CM.WGTSMACE_CovMat);
                            CNum=chol(CM);
%                             obj.TBDIS   = obj.TBDIS + ...
%                                 sqrt(obj.TBDIS).*randn(obj.nqU,1) + ...
%                                 sqrt(abs(CNum' * obj.TBDISE)) * randn(1);
                            
                            obj.TBDIS   = obj.TBDIS + sqrt(obj.TBDIS).*randn(obj.nqU,1)...
                                + (CNum' * ones(1,obj.nqU)')*randn(1);
                            
                    end
                    obj.TBDISE  = obj.TBDIS.^0.5;
            end
        end
        
        % Plots / Display
        function          PlotTBDDS(obj,varargin)
            
            p = inputParser;
            p.addParameter('fign',999,@(x)isfloat(x) && x>0);
            p.addParameter('type','lin',@(x)ismember(x,{'lin','log'}));
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            fign  = p.Results.fign;
            type  = p.Results.type;
            saveplot   = p.Results.saveplot;
            
            if size(obj.TBDDS,2)>1 %multipixel/multiring
              te = repmat(obj.Te,1,size(obj.TBDDS,2))-obj.Q;
            else
                te = obj.Te-obj.Q;
            end
            tbdds = obj.TBDDS;
            
            f1 = figure(fign);
            set(f1,'Units','normalized','Position',[0.1,0.1,0.5,0.5]);
             h = plot(te,tbdds,'-','LineWidth',3,'Color',rgb('DodgerBlue'));
            switch type
                case 'lin'
                    set(gca, 'YScale', 'lin');
                case 'log'
                    set(gca, 'YScale', 'log');
            end
            %errorbar_tick(h,1000);
            xlabel(sprintf('Energy - {\\itE}_0 (eV)'));
            ylabel('Rate per energy bin');
            %title('Tritium Beta Decay','FontSize',14)
            set(gca,'FontSize',12);
            if size(te,2)>1
                a = legend(h,arrayfun(@(x) sprintf('Differential spectrum ring %.0f',x),1:size(te,2),'UniformOutput',0));% legend(a,'boxoff');
                legend boxoff;
            else
                a = legend(h,'Differential spectrum');% legend(a,'boxoff');
            end
             PrettyFigureFormat('FontSize',22);
            a.EdgeColor = rgb('Silver');
            a.FontSize = get(gca,'FontSize');
            if strcmp(saveplot,'ON')
                export_fig(fign,'./plots/TBD_phasespace.pdf')
            end
        end
        
        function          PlotTBDIS(obj,varargin)
            
            p = inputParser;
            p.addParameter('fign',999,@(x)isfloat(x) && x>0);
            p.addParameter('type','lin',@(x)ismember(x,{'lin','log'}));
            p.addParameter('pub',0,@(x)isfloat(x) && x>=0);
            p.addParameter('scalingErr',1,@(x)isfloat(x) && x>=0);
            p.addParameter('cps','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('CovMat','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('qUrange',310,@(x)isfloat(x));
            p.parse(varargin{:});
            fign         = p.Results.fign;
            type         = p.Results.type;
            scalingErr   = p.Results.scalingErr;
            pub          = p.Results.pub;
            cps          = p.Results.cps;
            CovMat       = p.Results.CovMat;
            qUrange      = p.Results.qUrange;
            
            figure(fign)
            e = obj.qU(obj.qU(:,1)>(obj.Q-qUrange),:)-obj.Q;
            tmpis = obj.TBDIS(obj.qU(:,1)>(obj.Q-qUrange),:);
            tmpise = obj.TBDISE(obj.qU(:,1)>(obj.Q-qUrange),:)*scalingErr;
            switch CovMat
                case 'ON'
                    nCMTrials = 5000;
                    fName = [getenv('SamakPath'),sprintf('/inputs/CovMat/WGTSMACE_CovMat_%gTrials_%s_Sec_%g.mat',...
                        nCMTrials,obj.TD,obj.TimeSec)];
                    if exist(fName, 'file') == 2
                        CM=importdata(fName);
                    else
                        obj.ComputeTBDISCovarianceMatrix('nTrials',nCMTrials);
                        CM=importdata(fName);
                    end
                    tmpise_statsyst = sqrt((diag(CM.WGTSMACE_CovMat)) + (obj.TBDIS));
                    tmpise_statsyst = tmpise(obj.TBDIS>0)*scalingErr;
            end
            switch cps
                case 'OFF'
                    h1 = errorbar(e,tmpis,tmpise,'-s','MarkerSize',3,'MarkerFaceColor',.8*[1 1 1],'LineWidth',1);
                    % errorbar_tick(h1,200);
                    switch CovMat
                        case 'ON'
                            hold on
                            h2 = errorbar(e,tmpis,tmpise_statsyst,'ks','MarkerSize',3,'MarkerFaceColor',.4*[1 1 1],'Color','Red','LineWidth',1);
                            %errorbar_tick(h2,200);
                            hold off
                    end
                    
                    ylabel('Counts','FontSize',14);
                case 'ON'
                    tmpis = tmpis./obj.qUfrac(obj.qU(:,1)>(obj.Q-qUrange),:)./obj.TimeSec;
                    h1 = errorbar(e,tmpis,1e-99*tmpis,'-s','MarkerSize',3,'MarkerFaceColor',.8*[1 1 1],'LineWidth',1);
                    ylabel('Counts per second','FontSize',14);
            end
            switch type
                case 'lin'
                    set(gca, 'YScale', 'lin');
                case 'log'
                    set(gca, 'YScale', 'log')  ;
            end
            axis([e(1) e(end) 0.9*min(min(tmpis)) 1.2*max(max(tmpis))]);
            grid on
            xlabel('E-E_0 (eV)','FontSize',14);
            title('Tritium Beta Decay -  Integral Spectrum','FontSize',14)
            set(gca,'FontSize',12);
            switch CovMat
                case 'OFF'
                    a = legend(sprintf('1\\sigma uncertainty \\times %.0f',scalingErr)); legend('boxoff');
                case 'ON'
                    a = legend(sprintf('1\\sigma uncertainty \\times %.0f',scalingErr),'Stat. + Syst. Uncertainties x 10'); legend('boxoff'); 
            end
            if size(tmpis,2)>0 
                 a = legend(h1,arrayfun(@(x) sprintf('ring %.0f: 1\\sigma uncertainty \\times %.0f',x,scalingErr),...
                     1:size(tmpis,2),'UniformOutput',0)); 
                 legend('boxoff');
             
            end
            PrettyFigureFormat();
            if pub>0
                publish_figure(fign,'./figs/TBD_phasespace.eps')
            end
            
            
        end
        
        function          PlotTBDISallPixels(obj,varargin)
            
            p = inputParser;
            p.addParameter('fign',999,@(x)isfloat(x) && x>0);
            p.addParameter('type','lin',@(x)ismember(x,{'lin','log'}));
            p.addParameter('pub',0,@(x)isfloat(x) && x>=0);
            p.addParameter('exclDataStart',1,@(x)isfloat(x));
            p.parse(varargin{:});
            fign         = p.Results.fign;
            type         = p.Results.type;
            pub          = p.Results.pub;
            exclDataStart = p.Results.exclDataStart;
            
            figure(fign)
            [x, y] = meshgrid(obj.qU(exclDataStart:end)-obj.Q,obj.nPixels:1:1);
            ribbon((obj.TBDIS(exclDataStart:end,:)));colormap('summer')
            ylabel('qU')
            xlabel('Pixel / Ring')
            zlabel('Integral Spectrum All Pixels / Rings')
            
            switch type
                case 'lin'
                    set(gca, 'YScale', 'lin');
                case 'log'
                    set(gca, 'YScale', 'log')  ;
            end
            
            title('Tritium Beta Decay - KATRIN Integral Spectrum all Pixels simultaneously')
            set(gca,'FontSize',14);
            
            if pub>0
                publish_figure(fign,'./figs/TBDISallPixels.eps')
            end
            
        end
        
        function DisplayTDBInfo(obj,varargin)
            p=inputParser;
            p.addParameter('Output','screen',@(x)ismember(x,{'file','screen'}));
            p.addParameter('filename','',@(x)ischar(x));
            p.parse(varargin{:});
            Output = p.Results.Output;
            filename = p.Results.filename;
            
            switch Output
                case 'file'
                    try
                        fileID = fopen(filename,'w+');
                    catch
                        fprintf('file doesnt exist. enter valid filename!\n');
                        return
                    end
                case 'screen'
                    fileID = 1;
            end
            fprintf(fileID,'--------------------------------------------------------------\n');
            fprintf(fileID,'Tritium Beta Decay Spectrum with following properties:\n');
            fprintf(fileID,'--------------------------------------------------------------\n');
            fprintf(fileID,'   - EndPoint: %g eV \n',obj.Q);
            fprintf(fileID,'   - (Kinetic) Energy Range: [%g,%g] eV - %g bins \n',obj.TeMin,obj.TeMax,obj.nTe);
            fprintf(fileID,'   - 2 body-decay phase space + relativistic Fermi function \n');
            fprintf(fileID,'   - Screening correction: %s \n',obj.ScreeningFlag);
            fprintf(fileID,'   - Finite Extention of Nucleus Charge: %s \n',obj.FiniteExtChargeFlag);
            fprintf(fileID,'   - Weak interaction finite size correction: %s \n',obj.WintFiniteSizeFlag);
            fprintf(fileID,'   - Electron-Electron Exchange correction: %s \n',obj.EEexchangeFlag);
            fprintf(fileID,'   - Recoil Coulomb correction: %s \n',obj.RecoilCoulombFlag);
            fprintf(fileID,'   - Radiative corrections: %s (type %g)\n',obj.RadiativeFlag,obj.RadType);
            fprintf(fileID,'   - Weinheimer93 Phase Space: %s\n',obj.PS_Wein93);
            fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n');
            fprintf(fileID,'- - - - - Neutrinos \n');
            fprintf(fileID,'   - Effective Neutrino mass: %g eV \n',sqrt(obj.mnuSq_i));
            fprintf(fileID,'   - Sterile Neutrino -  mass: %g eV - mixing: %g\n',sqrt(obj.mnu4Sq_i),obj.sin2T4_i);
            fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n');
            fprintf(fileID,'- - - - - WGTS \n');
            fprintf(fileID,'   - T-T - FSD: %s (activity fraction = %g)\n',obj.TTFSD,obj.WGTS_MolFrac_TT);
            fprintf(fileID,'   - D-T - FSD: %s (activity fraction = %g)\n',obj.DTFSD,obj.WGTS_MolFrac_DT);
            fprintf(fileID,'   - H-T - FSD: %s (activity fraction = %g)\n',obj.HTFSD,obj.WGTS_MolFrac_HT);
            fprintf(fileID,'   - T-Decay Rate out of WGTS: %g Counts in %g sec \n',obj.NormFactorTBDDS./obj.CumFrac,obj.TimeSec);
            fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n');
            fprintf(fileID,'- - - - - MACE \n');
            fprintf(fileID,'   - Integral Spectrum - Time Distribution: %s \n',obj.TD);
            fprintf(fileID,'   - Transmission Function: %s \n',obj.KTFFlag);
            fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n');
            fprintf(fileID,'- - - - - Backgrounds \n');
            fprintf(fileID,'   - Flag: %s  - Type: %s\n',obj.BKG_Flag,obj.BKG_Type);
            switch obj.FPD_Segmentation
                case 'OFF'
                    fprintf(fileID,'   - Rate: %g mcps per qU (renormalized for 148 pixel)\n',obj.BKG_RateSec*1e3./numel(obj.FPD_PixList).*148);
                case 'RING'
                case 'PIXEL'
                    fprintf(fileID,'   - Rate: %g: %g mcps per qU\n',obj.FPD_PixList,obj.BKG_RateSec_i*1e3);
            end
            fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n');
            fprintf(fileID,'- - - - - FPD \n');
            fprintf(fileID,'   - Segmentation: %s  \n',obj.FPD_Segmentation);
            switch obj.FPD_Segmentation
                case 'OFF'
                    fprintf(fileID,'   - Mean Efficiency: %g  \n',obj.FPD_MeanEff);
                case 'RING'
                case 'PIXEL'
            end
            if ~isempty(obj.TBDIS)
                fprintf(fileID,'   - Tritium e- FPD: %g cps at (qU>%g eV)\n',obj.TBDIS(1)/(obj.TimeSec.*obj.qUfrac(1)),obj.TeMin);
            end
            fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n');
            fprintf(fileID,'--------------------------------------------------------------\n');
            switch Output
                case 'file'
                    fclose(fileID);
            end
            
        end
        
        function          PlotFSD(obj,varargin)
            
            p = inputParser;
            p.addParameter('fign',999,@(x)isfloat(x) && x>0);
            p.addParameter('type','lin',@(x)ismember(x,{'lin','log'}));
            p.addParameter('pub','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('holdf','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('TT','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('DT','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('HT','ON',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            fign  = p.Results.fign;
            type  = p.Results.type;
            pub   = p.Results.pub;
            holdf = p.Results.holdf;
            TT = p.Results.TT;
            DT = p.Results.DT;
            HT = p.Results.HT;
            
            switch TT
                case 'ON'
                    f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.4]);
                    h1 = plot(obj.TTexE,obj.TTexP*100,'Color',rgb('DodgerBlue'),'LineWidth',2);
                    switch type
                        case 'lin'
                            set(gca, 'YScale', 'lin');
                        case 'log'
                            set(gca, 'YScale', 'log')  ;
                    end
                    % grid on
                    %axis([0.1 250 1e-3 50]);
                    %ylim([1e-3 50]);
                    xlim([min(obj.TTexE) max(obj.TTexE)]);
                    xlabel('Excitation energy (eV)');
                    ylabel('Probability (%)');
                    %fsdtitle=sprintf('TT Ground and Excited States - %s (\\sigma^2=%.1f eV^2)',obj.TTFSD,var([obj.TTexE_G obj.TTexE_E],[obj.TTexP_G*obj.TTNormGS obj.TTexP_E*obj.TTNormES]));
                    %title(fsdtitle,'FontSize',14);                
                    strG = sprintf(' Ground States: (P=%.1f%% \\sigma^2=%.3f eV^2)',obj.TTNormGS*100,var(obj.TTexE(1:obj.TTGSTh),obj.TTexP(1:obj.TTGSTh)));
                    strE = sprintf('\n Excited States: (P=%.1f%% \\sigma^2=%.1f eV^2)',obj.TTNormES*100,var(obj.TTexE(obj.TTGSTh+1:end),obj.TTexP(obj.TTGSTh+1:end)));
                    leg = legend([strG,strE],'Location','NorthEast'); %legend(a,'boxoff','FontSize',12);
                    legend boxoff
                    leg.Title.String = sprintf('FSD %s T_2',obj.TTFSD);
                    PrettyFigureFormat('FontSize',20);
                    if pub>0
                        figname = sprintf('./plots/TBD_FSD_TT_%s.pdf',obj.TTFSD);
                        publish_figurePDF(fign,figname)
                    end
            end
            
            switch DT
                case 'ON'
                    figure(fign+1)
                    hold on
                    h1 = stairs(-obj.DTexE_G+1.7,obj.DTexP_G*100,'Color','Black','LineWidth',2);
                    hold on
                    h2 = stairs(-obj.DTexE_E+1.7,obj.DTexP_E*100,'Color','Red','LineWidth',2);
                    hold off
                            set(gca, 'XScale', 'lin');
                    switch type
                        case 'lin'
                            set(gca, 'YScale', 'lin');
                        case 'log'
                            set(gca, 'YScale', 'log')  ;
                    end
                    grid on
                    ylim([1e-3 50]);
                    xlim([-250 10]);
                    %axis([0.1 250 1e-3 10]);
                    xlabel('Binding Energy (eV)','FontSize',14);
                    ylabel('Probability (%)','FontSize',14);
                    fsdtitle=sprintf('DT Ground and Excited States - %s (\\sigma^2=%.1f eV^2)',obj.DTFSD,var([obj.DTexE_G obj.DTexE_E],[obj.DTexP_G*obj.DTNormGS obj.DTexP_E*obj.DTNormES]));
                    title(fsdtitle,'FontSize',14);
                    set(gca,'FontSize',12);
                    strG = sprintf('Ground States: (P=%.1f%% \\sigma^2=%.3f eV^2)',obj.DTNormGS*100,var(obj.DTexE_G,obj.DTexP_G));
                    strE = sprintf('Excited States: (P=%.1f%% \\sigma^2=%.1f eV^2)',obj.DTNormES*100,var(obj.DTexE_E,obj.DTexP_E));
                    leg = legend([h1 h2],strG,strE,'Location','NorthWest');
                    %legend(a,'boxoff','FontSize',12);
                    PrettyFigureFormat;
                    if pub>0
                        figname = sprintf('./plots/TBD_FSD_DT_%s.pdf',obj.DTFSD);
                        publish_figurePDF(fign+1,figname)
                    end
            end
            
            switch HT
                case 'ON'
                    figure(fign+2)
                    hold on
                    h1 = stairs(-obj.HTexE_G+1.7,obj.HTexP_G*100,'Color','Black','LineWidth',2);
                    hold on
                    h2 = stairs(-obj.HTexE_E+1.7,obj.HTexP_E*100,'Color','Red','LineWidth',2);
                    hold off
                            set(gca, 'XScale', 'lin');
                    switch type
                        case 'lin'
                            set(gca, 'YScale', 'lin');
                        case 'log'
                            set(gca, 'YScale', 'log')  ;
                    end
                    grid on
                    %axis([0.1 250 1e-3 10]);
                     ylim([1e-3 50]);
                    xlim([-250 10]);
                    xlabel('Binding Energy (eV)','FontSize',14);
                    ylabel('Probability (%)','FontSize',14);
                    fsdtitle=sprintf('HT Ground and Excited States - %s (\\sigma^2=%.1f eV^2)',obj.HTFSD,var([obj.HTexE_G obj.HTexE_E],[obj.HTexP_G*obj.HTNormGS obj.HTexP_E*obj.HTNormES]));
                    title(fsdtitle,'FontSize',14);
                    set(gca,'FontSize',12);
                    strG = sprintf('Ground States: (P=%.1f%% \\sigma^2=%.3f eV^2)',obj.HTNormGS*100,var(obj.HTexE_G,obj.HTexP_G));
                    strE = sprintf('Excited States: (P=%.1f%% \\sigma^2=%.1f eV^2)',obj.HTNormES*100,var(obj.HTexE_E,obj.HTexP_E));
                    leg = legend([h1 h2],strG,strE,'Location','NorthWest');
                    PrettyFigureFormat;
                    if pub>0
                        figname = sprintf('./plots/TBD_FSD_HT_%s.pdf',obj.DTFSD);
                        publish_figurePDF(fign+2,figname)
                    end
            end
            
        end
        function PlotNuMassImprint(obj,varargin)
            % plot neutrino mass imprint
            % use chi2 display: (TBDIS1eV-TBDIS0eV)^2/sqrt(TBDIS0eV)
            % advantage: effect of increased background is taken into account
            % when taking just the ratio: background is eliminated
            p=inputParser;
            p.addParameter('qUrange',40,@(x)isfloat(x));
            p.parse(varargin{:});
            qUrange = p.Results.qUrange;
            
            qUIndex = logical(obj.qU>=(18574-qUrange));
            qU = obj.qU(qUIndex);%linspace(min(obj.qU(qUIndex)),max(obj.qU(qUIndex)),500);
            time = obj.qUfrac(qUIndex).*obj.TimeSec;
            
            % current background
            Bkg_prev = obj.BKG_RateSec;
            obj.ComputeTBDDS('mSq_bias',0); obj.ComputeTBDIS;
            TBDISmNuSq0eV = obj.TBDIS(qUIndex)./time;
            obj.ComputeTBDDS('mSq_bias',1); obj.ComputeTBDIS;
            TBDISmNuSq1eV = obj.TBDIS(qUIndex)./time;
            
            % 10mcps background
            obj.BKG_RateSec_i=0.01; obj.SetFitBias(1);
            obj.ComputeTBDDS('mSq_bias',0); obj.ComputeTBDIS;
            TBDISmNuSq0eV_lowBkg = obj.TBDIS(qUIndex)./time;%interp1(obj.qU(qUIndex),obj.TBDIS(qUIndex)./time,qU,'lin');
            obj.ComputeTBDDS('mSq_bias',1); obj.ComputeTBDIS;
            TBDISmNuSq1eV_lowBkg = obj.TBDIS(qUIndex)./time;%interp1(obj.qU(qUIndex),obj.TBDIS(qUIndex)./time,qU,'lin');

             % 1000mcps background
            obj.BKG_RateSec_i=1; obj.SetFitBias(1);
            obj.ComputeTBDDS('mSq_bias',0); obj.ComputeTBDIS;
            TBDISmNuSq0eV_highBkg = obj.TBDIS(qUIndex)./time;%interp1(obj.qU(qUIndex),obj.TBDIS(qUIndex)./time,qU,'lin');
            obj.ComputeTBDDS('mSq_bias',1); obj.ComputeTBDIS;
            TBDISmNuSq1eV_highBkg = obj.TBDIS(qUIndex)./time;%interp1(obj.qU(qUIndex),obj.TBDIS(qUIndex)./time,qU,'lin');
            obj.BKG_RateSec_i = Bkg_prev;obj.SetFitBias(1);
            
            Chi2_lowBkg_tmp  = (TBDISmNuSq1eV_lowBkg-TBDISmNuSq0eV_lowBkg).^2./TBDISmNuSq0eV_lowBkg;
            Norm = max(Chi2_lowBkg_tmp);
            Chi2_lowBkg = Chi2_lowBkg_tmp./Norm;
            Chi2_highBkg = ((TBDISmNuSq1eV_highBkg-TBDISmNuSq0eV_highBkg).^2./TBDISmNuSq0eV_highBkg)./Norm;
            Chi2         = ((TBDISmNuSq1eV-TBDISmNuSq0eV).^2./TBDISmNuSq0eV)./Norm;
  
            fig1= figure('Units','normalized','Position',[0.1, 0.1, 0.5,0.5]);
            p1 = plot(qU-18574,Chi2,'x-','LineWidth',3);
            hold on;
            p2 = plot(qU-18574,Chi2_lowBkg,'x-','LineWidth',3);
            p3 = plot(qU-18574,Chi2_highBkg,'x-','LineWidth',3);
            PrettyFigureFormat('FontSize',24);
            xlabel('qU -18574 (eV)');
            ylabel(sprintf('\\chi^2_{m_\\nu signal} (arb. units)'));
            grid on;
            title(sprintf('\\chi^2_{m_\\nu signal} = [IS(m_\\nu^2=0 eV^2)-IS(m_\\nu^2=1 eV^2)]^2/IS(m_\\nu^2=0 eV^2)'))
            leg = legend([p2,p1,p3],sprintf(' %.0f mcps',0.01*1e3),sprintf(' %.0f mcps',Bkg_prev*1e3),sprintf(' %.0f mcps',1*1e3));
            leg.Title.String = 'Background'; legend boxoff;
        end
        % Doppler Effect
        function     computeKernel(obj) %REMEMBER: include E_cms in doppler calls
            % The kernel is calculated according to DOI: 10.1140/epjc/s10052-019-6686-7
            
            % The standard deviations of the kernel funtion depending on
            % velocity (DE_sigma_e_vel) and energy (DE_sigma).
            obj.DE_sigma_e_vel = sqrt(obj.kb*obj.WGTS_Temp/obj.M);              % std of velocity distribution
            obj.DE_sigma = obj.DE_sigma_e_vel*sqrt(2*obj.E_cms*obj.me/obj.c^2); % std of energy distribution
            
%             % range for the kernel
%             obj.e_velParal = (obj.Te - obj.Q)/((obj.me/obj.c^2)*obj.e_vel);
%             obj.cosO = linspace(cos(asin(sqrt(obj.WGTS_B_T./obj.MACE_Bmax_T))),1,obj.nTe);
%             [e_velParal_M,cosO_M] = meshgrid(obj.e_velParal,obj.cosO);
%             prefac = (1-cos(asin(sqrt(obj.WGTS_B_T./obj.MACE_Bmax_T)))).^(-1);
%             
%             % integration of the emission angle
%             g_vm = prefac./(sqrt(2*pi)*obj.DE_sigma_e_vel)*...
%                 simpsons(obj.cosO,exp(-0.5*((e_velParal_M-cosO_M*obj.BulkVelT2)/obj.DE_sigma_e_vel).^2));
%             
%             obj.GaussianKernelDE = 1/((obj.me/obj.c^2)*obj.e_vel)*g_vm;
%             obj.ProbKernDE = simpsons(obj.Te - obj.E_cms,obj.GaussianKernelDE);
%             g_fwhm = obj.GaussianKernelDE - 0.5*max(obj.GaussianKernelDE);
%             g_fwhm_index = find(0 < g_fwhm);
%             x1 = g_fwhm_index(1); xf = g_fwhm_index(end);
%             x0_left = obj.Te(x1) - g_fwhm(x1)*(obj.Te(x1-1)-obj.Te(x1))/(g_fwhm(x1-1)-g_fwhm(x1));
%             x0_right = obj.Te(xf) - g_fwhm(xf)*(obj.Te(xf+1)-obj.Te(xf))/(g_fwhm(xf+1)-g_fwhm(xf));
%             FWHM = x0_right - x0_left;
%             obj.DE_sigma = FWHM/(2*sqrt(2*log(2)));
        end % computeKernel
        
        function          plotKernel(obj)
            figure(1);
            hold on
            obj.enlargeSpectrum;
            plot(obj.Te-obj.E_cms,obj.GaussianKernelDE,'-x');
            leg_cell = {['u = ',num2str(obj.BulkVelT2),' m/s']};
            hold off
            PrettyFigureFormat;
            xlabel('E_{lab}-E_{cms} (eV)');
            ylabel('probability');
            xlim([-0.5,0.5])
            legend(leg_cell);
            obj.restoreSpectrum(obj.TBDDS);
        end % plotKernel
        
        function          enlargeSpectrum(obj)
            % add neccessary bins left and right
            nTeStepsMinus = obj.nTeBinningFactor*obj.Eminus*obj.TeStep;
            nTeStepsPlus = obj.nTeBinningFactor*obj.Eplus*obj.TeStep;
            obj.TeMin = obj.qUmin - nTeStepsMinus;
            obj.TeMax = obj.qUmax + nTeStepsPlus;
            obj.Te = [(obj.TeMin:obj.TeStep:obj.qUmin-obj.TeStep)'...
                ;obj.Te;(obj.qUmax+obj.TeStep:obj.TeStep:obj.TeMax)'];
            obj.nTe = length(obj.Te);
            
            % Set properties again
            obj.SetKinVariables()
            % Memory Allocation
            obj.TBDDS = zeros(obj.nTe,1);
        end % enlargeSpectrum
        
        function convolutedspectrum = KernelSpectrumConv(obj,spec)
            switch obj.DopplerEffectFlag
                case 'numConv' % [Work in progress, don't use]
                    convolutedspectrum = obj.TeStep*conv(spec,obj.GaussianKernelDE');
                    nKernel = length(obj.GaussianKernelDE);
                    nconvspec = length(convolutedspectrum);
                    quarterconv = round(nconvspec/4);
                    middlekernel = round(nKernel/2);
                    indexkernel = (1:nKernel)';
                    indexEcms = indexkernel(obj.Te == obj.Q_i);
                    differenceposition = abs(middlekernel - indexEcms);
                    begin = quarterconv + differenceposition + 1;
                    endspec = begin + nKernel - 1;
                    convolutedspectrum = convolutedspectrum(begin:endspec);
                case 'matConv' % convolution using matrix
                    obj.mate = obj.eres(obj.Te,obj.TeMin,obj.TeMax,obj.nTe,0.95056255902670);
                    %obj.mate = obj.eres(obj.Te,obj.TeMin,obj.TeMax,obj.nTe,obj.DE_sigma);
                    %figure(5)
                    %imagesc(obj.mate);
%                     figure(6)
%                     stairs(obj.Te);
                    %eres(E,Emin,Emax,Nbins,sigma)
                    %E = obj.Te; Emin = obj.TeMin; Emax = obj.TeMax; Nbins = obj.nTe; sigma = obj.DE_sigma;
                    convolutedspectrum = obj.mate*spec;
            end
        end % KernelSpectrumConv
        
        function restoredspectrum = restoreSpectrum(obj,spec)
            % Spectrum is restored to the original size
            restoredspectrum = spec(obj.Te>=obj.qUmin & obj.Te<=obj.qUmax);
            %if obj.nPixels == 1
            obj.Te = obj.Te(obj.Te>=obj.qUmin & obj.Te<=obj.qUmax);
            obj.nTe = length(obj.Te);
            %end
        end
        
        % Miscellaneous
        function r      = W(obj,sigma_E)
            % ----------------------------------------------------------
            % Apply Energy Resolution Function
            % ----------------------------------------------------------
            r = 0.5*(erf((obj.TeMax-obj.Te)/sqrt(2)/sigma_E./sqrt(obj.Te))...
                -erf((obj.Te-obj.TeMin)/sqrt(2)/sigma_E./sqrt(obj.Te)));
        end      
        
        function AdjustMolFrac(obj)
            %Call this function after changing isotopologues concentration 
            ComputeTritiumPurity(obj);
            ComputeFracTBDtail(obj);  
            ComputeNormFactorTBDDS(obj);  
        end
        function AdjustRF(obj)
            %Call this function after changing response function parameter
            %(column density, Bs, Bmax, Ba)
            InitializeRF(obj);
            ComputeFracTBDtail(obj);
            ComputeNormFactorTBDDS(obj);
        end
    end %methods
    
    methods(Static)
        
        function r = eres(E,Emin,Emax,Nbins,sigma)
            % Produces Energy resolution matrix from Ebins to Nbins:
            % r is [nB nE].
            Estep = (Emax-Emin)/Nbins;
            dn = repmat(reshape(linspace(Emin,Emax-Estep,Nbins),[Nbins 1]),[1 numel(E)]);
            up = repmat(reshape(linspace(Emin+Estep,Emax,Nbins),[Nbins 1]),[1 numel(E)]);
            f = @(x,a,b,sigma).5*double((erf((b-x)./(sqrt(2)*sigma))-erf((a-x)./(sqrt(2)*sigma))));
            r = f(repmat(reshape(E,[1 numel(E)]),[Nbins 1]),dn,up,sigma);
        end % eres
    end
    
end % class
