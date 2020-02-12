% ----------------------------------------------------------------------- %
%
% Class computing the electron energy spectrum in the
% Krypton 83m Decays
% ----------------------------------------------------------------------- %
%
% Detector Corrections
%  - Energy resolution
%  - Flat background
%  - MAC-E Filter Transmission Function
%
% Th. Lasserre, 2017
% CEA - Saclay
% TUM - IAS
%
% ----------------------------------------------------------------------- %



classdef Kr_Lisa < handle & FPD & WGTSMACE

    
    properties (Constant = true, Hidden = true)
        
        % Physics Constants
        fsC = 0.007297352569816;    % Fine Structure Cte
        me = 510.998918e3;          % Electron Mass eV
        Gf = 8.963859529398992e-50; % Fermi Cte Mev.m^3 PDG2010
        Mn = 939565.360e3;          % Neutron mass eV
        Mp = 938272.029e3;          % Proton mass eV
        gA = 1.247; gV = 1.000;     % Weak Interaction Coupling constants
        kb = 1.38064852e-23;        % [m2 kg s^-2 K^-1] Boltzmann Constant
        M = 82.914136099e-3/6.022140857e23; % [kg] mass of atomic krypton 83
        c = 299792458;              % [m/s] speed of light
        
        % Properties of the broadening kernel of the Doppler effect
        %         E_cms = 18575;              % eV Energy center of mass
        %         v_e = 8.083326511496116e+07;% [m/s] velocity of the electron at E_0
        E_cms = 30474;              % eV Energy center of mass
        v_e   = 9.9223e+07;% [m/s] velocity of the electron at 30474 keV...
        
    end
    
    properties (Access=public, Hidden = true)
        % 83mKr Parameters
        Z=36;                     % Krypton
        
        % L3_32 - single pixel
        %-->Line specification: transmission 32keV & Shell 3     

        L3_32_E_i     = 30472.5; L3_32_E;     % Mean Energy, eV
        L3_32_W_i     = 1.17;    L3_32_W;     % Width, eV
        L3_32_Phi0_i  = 100;     L3_32_Phi0;  % Intensity
        L3_32_B_i     = 0.823;   L3_32_B;     % Branching Ratio

        % L3_32 - all pixels
        L3_32_Phi0allPixels_i  = 100*ones(1,148);  % Intensity All Pixels
        L3_32_Phi0allPixels;                       % Intensity All Pixels
        
        % K32 - single pixel
        K32_E_i      = 17824.2;  K32_E;       % Mean Energy, eV
        K32_W_i      = 2.71;     K32_W;       % Width, eV
        K32_Phi0_i   = 100;      K32_Phi0;    % Intensity
        K32_B_i      = 0.795;    K32_B;       % Branching Ratio

        %values from report Martin Slezak gave me 'KATRIN krypton mode'
        p_main=1-0.193;            %probability of main line->no shake up/off effect  
        p_sat3= 0.193;
        
        L3_32_E_s3;                  %Energy Difference: data-fit
        L3_32_E_s3_i = 30472.5-23;  %Mean Energy of Satellite Line 3
        L3_32_Phi0_s3_i = (0.193/(1-0.193))*100 %(p_sat3/p_main).*L3_32_Phi0_i;   %Intensity satellite
        L3_32_Phi0_s3;

        % K32 - all pixels
        K32_Phi0allPixels_i  = 100*ones(1,148);  % Intensity All Pixels
        K32_Phi0allPixels;                       % Intensity All Pixels
        
        % all lines stored
        tmpline;       %tmp
        tmplineS1;     %tmp line satellites
        tmplineS2;
        tmplineS3;
        KrLine;
        
    end
    
    properties (Access=public)
        
        % Kinematics (Electron)
        Te;             % kinetic energy of the electron
        TeMin;          % minimum e kinetic energy considered
        TeMax;          % maximum e kinetic energy considered
        nTe;            % kinetic energy Bins
        TeStep;         % energy bin steps width
        qUTeMatch;      % matches between Te and qU arrays
        
        % Krypton Spectrum - Differential
        KrDS;           % Differential Krypton Spectrum

        KrDSE;          % Error on Differential Krypton Spectrum
        KrDS_s3;	% Differential Spectrum satellite
	KrDSE_s3;
	NormFactorKrDS; % Normalization Factor
        mate;           % Migration Matrix For E Fast Resolution
        
        % Krypton Spectrum - Integral
        KrIS;           % Integral Krypton Spectrum
        KrISE;          % Error on Krypton Spectrum
        
        % Krypton Spectrum - Integral - ALL PIXELS
        KrDSallPixels;  % Differential Krypton Spectrum - Vector[obj.nPixels]
        
        KrDSEallPixels; % Error on Differential Krypton Spectrum
        KrISallPixels;  % Integral Krypton Spectrum - Vector[obj.nPixels]
        KrISEallPixels; % Error on Krypton Spectrum - Vector[obj.nPixels]
        nPixels;
        
        % [ADD]
        % Doppler Effect
        DopplerEffectFlag; % turns on and off the Doppler Effect
        ConvFlag;       % enables to choose between Gaussian Convolution and Voigt Function for the Kr Spectrum.
        qUminTD;        % Store the original TD qUmin during the convolution
        qUmaxTD;        % Store the original TD qUmax during the convolution
        Eminus;         % Extra bins added to the end of the spectrum to account for convolution edge effect
        Eplus;          % Extra bins added to the beginning of the spectrum to account for conv. edge effect
        Omax;          % [rad] max emission angle
        v_m;            % [m/s] electron vel. component parallel to emission direction
        cosO;           % cosine of the emission angle
        setlog;         % instruction to set logarithmic scale on comparison of the spectra with and without Doppler Effect
        recomputeKernel;% switch to recompute kernel
        realConv;       % calculate the convolution with the usual method
        o_v;            % [eV] standard deviation for the kernel in terms of the electron velocity
        StandDev;            % [eV] standard deviation for the kernel in terms of the electron energy (usually around 130 meV)
        StandDev_num;        % [eV] standard deviation for the kernel in terms of the electron energy calculated numerically
        DE_Bias;
        
        % Broadening kernel g parameters
        E_lab;          % [eV] electron energy in the laboratory frame
        E_range;        % [eV] energy around the endpoint to taken into account for the kernel of the Doppler Effect
        T;              % [K] temperature if the tritium
        u;              % [m/s] bulk velocity of tritium (mean = 13)
        g;              % broadening kernel
        P;              % probability under the curve of the kernel (should = 1)
        % [ADD] END
        
	%Mutiple Peaks
        MultiPeaksFlag;

        % Binning
        nTeBinningFactor;
        
        % Fit: Bias - Single Pxiel
        E_Bias;
        W_Bias;
        B_Bias;
        Phi0_Bias;
        E_Bias_s3; %satellite
        Phi0_s3_Bias;
        % Miscellaneous
        CPS;            % Output of IS in counts per second
    end
    
    methods

        function obj = Kr_Lisa(varargin)

            fprintf(2,'Processing Kr Constructor ...\n');
            
            p = inputParser;
           
            %: GENERAL SETTINGS

            p.addParameter('Normalization','NOMINAL',@(x)ismember(x,{'NOMINAL','COUNTS'}));
            p.addParameter('CPS','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('TimeYear',3,@(x)isfloat(x) && x>0);
            p.addParameter('TD','KrK32',@(x)ismember(x,{'DR20','DR30','DR40','DR50','Kle15','Kle15Ext','Kle15Dyn','Flat20','Flat30','Flat50','Flat60','Flat100','KrK32','KrL3_32','KrL3_32_HS','MoSFitterJuly2017', 'KrL3_32_Satellites'}));
            
            % WGTS Parameters
            p.addParameter('WGTS_CosMaxAAngle',0.6324,@(x)isfloat(x) && x>0);
            p.addParameter('WGTS_Tp',0.95,@(x)isfloat(x) && x>0);
            p.addParameter('WGTS_DTHTr',0.1,@(x)isfloat(x) && x>0);
            p.addParameter('WGTS_FTR_cm',4.1,@(x)isfloat(x) && x>0);
            p.addParameter('WGTS_CD_MolPerCm2',5e17,@(x)isfloat(x) && x>0);
            p.addParameter('WGTS_B_T',3.6,@(x)isfloat(x) && x>0);
            
            % MACE Parameters
            p.addParameter('MACE_Bmax_T',6,@(x)isfloat(x) && x>0);
            p.addParameter('MACE_Ba_T',3e-4,@(x)isfloat(x) && x>0);
            p.addParameter('MACE_R_eV',0.93,@(x)isfloat(x) && x>0);
            p.addParameter('HVRipples','OFF',@(x)ismember(x,{'OFF','ON'}));
            p.addParameter('HVRipplesP2PV',520e-3,@(x)isfloat(x) && x>0);
            
            % Both WGTS & MACE
            p.addParameter('KTFFlag','Kle',@(x)ismember(x,{'OFF','MACE','MACE+WGTSIS','TL','SG','Kle','MACE+RIPPLE'}));
            
            % Krypton
            p.addParameter('L3_32_Phi0_i',100,@(x)isfloat(x) && x>0);
            p.addParameter('K32_Phi0_i',100,@(x)isfloat(x) && x>0);
            p.addParameter('nPixels',1,@(x)isfloat(x) && x>0);
            p.addParameter('L3_32_Phi0allPixels',100*ones(1,148));
            p.addParameter('K32_Phi0allPixels',100*ones(1,148));
            p.addParameter('BKG_RateSecallPixels',1*ones(1,148));
            
            % [ADD]
            % Doppler Effect Flag
            p.addParameter('DopplerEffectFlag','OFF',@(x)ismember(x,{'OFF','ON'}));
            p.addParameter('ConvFlag','DopplerEffect', @(x)ismember(x,{'DopplerEffect','Voigt','OFF'}))
            p.addParameter('Eminus',6,@(x)isfloat(x) && x>0);
            p.addParameter('Eplus',6,@(x)isfloat(x) && x>0);
            p.addParameter('Omax',50.8*pi/180,@(x)isfloat(x) && x>0 && x<2*pi)
            p.addParameter('setlog','OFF',@(x)strcmp(x,'ON') || strcmp(x,'OFF'));
            p.addParameter('recomputeKernel','ON',@(x)strcmp(x,'ON') || strcmp(x,'OFF'));
            p.addParameter('realConv','OFF',@(x)strcmp(x,'ON') || strcmp(x,'OFF'));

            % Broadening kernel parameters
            p.addParameter('E_range',1,@(x)isfloat(x) && x>0);
            p.addParameter('T',100,@(x)isfloat(x) && x>0);
            p.addParameter('u',0,@(x)isfloat(x));
            p.addParameter('StandDev',0.058,@(x)isfloat(x) && x>0); % 0.134429846285963
            
            % [ADD] END
            
 	    %Multiple Peaks
            p.addParameter('MultiPeaksFlag','OFF',@(x)ismember(x,{'OFF','ON'}));

            % FPD
            p.addParameter('FPD_Segmentation','OFF',@(x)ismember(x,{'OFF','RING','PIXEL'}));
            p.addParameter('FPD_Pixel',1,@(x)isfloat(x) && x>-1); % 0=All Pixels
            p.addParameter('FPD_Ring',1,@(x)isfloat(x) && x>-1);  % 0=All Rings
            
            % Background Parameters
            p.addParameter('BKG_Flag','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('BKG_Type','FLAT',@(x)ismember(x,{'FLAT','SLOPE'}));
            p.addParameter('BKG_RateAllFPDSec',1e-3,@(x)isfloat(x) && x>0);
            p.addParameter('BKG_RateRingSec',1e-3,@(x)isfloat(x) && x>0);
            p.addParameter('BKG_RatePixelSec',1e-3,@(x)isfloat(x) && x>0);
            
            % Binning
            p.addParameter('nTeBinningFactor',50,@(x)isfloat(x) && x>0);
            
            p.parse(varargin{:});

            % KATRIN: GENERAL SETTINGS
            obj.Normalization       = p.Results.Normalization;
            obj.CPS                 = p.Results.CPS;
            obj.TimeYear            = p.Results.TimeYear;
            obj.TD                  = p.Results.TD;
            obj.nTeBinningFactor    = p.Results.nTeBinningFactor;
            % WGTS Parameters
            obj.WGTS_CosMaxAAngle   = p.Results.WGTS_CosMaxAAngle;
            obj.WGTS_Tp             = p.Results.WGTS_Tp;
            obj.WGTS_DTHTr          = p.Results.WGTS_DTHTr;
            obj.WGTS_FTR_cm         = p.Results.WGTS_FTR_cm;
            obj.WGTS_CD_MolPerCm2   = p.Results.WGTS_CD_MolPerCm2;
            obj.WGTS_B_T            = p.Results.WGTS_B_T;
            % MACE Parameters
            obj.MACE_Bmax_T         = p.Results.MACE_Bmax_T;
            obj.MACE_Ba_T           = p.Results.MACE_Ba_T;
            obj.MACE_R_eV           = p.Results.MACE_R_eV;
            % WGTSMACE: GENERAL SETTINGS
            obj.KTFFlag             = p.Results.KTFFlag;
            
            % [ADD]
            % Doppler Effect Flag
            obj.DopplerEffectFlag   = p.Results.DopplerEffectFlag;
            obj.ConvFlag            = p.Results.ConvFlag;
            obj.Eminus              = p.Results.Eminus;
            obj.Eplus               = p.Results.Eplus;
            obj.Omax                = p.Results.Omax;
            obj.DopplerEffectFlag   = p.Results.DopplerEffectFlag;
            obj.setlog              = p.Results.setlog;
            obj.recomputeKernel     = p.Results.recomputeKernel;
            obj.realConv            = p.Results.realConv;
            % [ADD] END
            
            % Krypton
            obj.K32_Phi0_i          = p.Results.K32_Phi0_i;
            obj.L3_32_Phi0_i        = p.Results.L3_32_Phi0_i;
            % Krypton all Pixels
            obj.nPixels             = p.Results.nPixels;
            obj.L3_32_Phi0allPixels = p.Results.L3_32_Phi0allPixels;
            obj.K32_Phi0allPixels   = p.Results.K32_Phi0allPixels;
            obj.BKG_RateSecallPixels= p.Results.BKG_RateSecallPixels;
            
            % Broadening kernel parameters
            obj.E_range     = p.Results.E_range;
            obj.T           = p.Results.T;
            obj.u           = p.Results.u;
            obj.StandDev    = p.Results.StandDev;
            
	    %Mutiple Peaks option
            obj.MultiPeaksFlag = p.Results.MultiPeaksFlag;

            % FPD
            obj.FPD_Segmentation    = p.Results.FPD_Segmentation;
            obj.FPD_Pixel           = p.Results.FPD_Pixel;
            obj.FPD_Ring            = p.Results.FPD_Ring;
            % Background Parameters
            obj.BKG_Flag            = p.Results.BKG_Flag;
            obj.BKG_Type            = p.Results.BKG_Type;
            obj.BKG_RateAllFPDSec   = p.Results.BKG_RateAllFPDSec;
            obj.BKG_RateRingSec     = p.Results.BKG_RateRingSec;
            obj.BKG_RatePixelSec    = p.Results.BKG_RatePixelSec;
            
            % Init
            SetFitBias(obj,0);
            
            % WGTS:  Initialization
            
            % MACE:  Initialization
            
            % FPD:   Initialization
            InitializeBackground(obj);
            switch obj.FPD_Segmentation
                case 'OFF'
                    obj.KTF = obj.GetMeanTFAllFPD();
                case 'RING'
                    
                case 'PIXEL'
                    SetPixel_MACE_BaEa(obj,'ReadConfigFile','OFF')
                    if obj.FPD_Pixel == 0
                        fprintf(2,'Processing All Pixels\n');
                        %fprintf(2,'List of Bad Pixels (not processed yet...)\n');

                        obj.KTFallPixels = obj.GetTFPixelWiseFPD_All(); %BUG

                        obj.L3_32_Phi0allPixels  = obj.L3_32_Phi0allPixels_i;
                        obj.K32_Phi0allPixels    = obj.K32_Phi0allPixels_i;
                        obj.BKG_RateSecallPixels = obj.BKG_RateSecallPixels_i;
                    else
                        obj.KTF = obj.GetTFPixelWiseFPD(obj.FPD_Pixel);
                        fprintf(2,'Processing Pixels %g\n',obj.FPD_Pixel);
                    end
            end
            
            % KATRIN Initialization
            InitializeTD(obj,'TD',obj.TD);
            
            % Binning - To get from the TD's
            SetKrDSBinning(obj,obj.qUmin,obj.qUmax,obj.nTeBinningFactor);
            
            % qU / Te match
            obj.qUTeMatch = obj.FindqUTeMatch(1.01e-2);
            if numel(obj.qUTeMatch) < obj.nqU
                fprintf('Anomaleous qU versus Te match  (%.0f / %0.f) --> Aborting... \n',...
                    numel(obj.qUTeMatch),obj.nqU); return;
            elseif numel(obj.qUTeMatch) > obj.nqU
                fprintf('Anomaleous qU versus Te match  (%.0f / %0.f) --> Aborting... \n',...
                    numel(obj.qUTeMatch),obj.nqU); return;
            end
            
            % Compute Normalization Factor
            ComputeNormFactorKrDS(obj);
            
           addpath(genpath('../../../Samak2.0')); 
           % addpath('fminuit'); %digits(6); 
           % addpath('tools')
            
        end % constructor
        
    end % methods
    
    methods
        
        % Initializations
        function SetTeBinning(obj,temin,temax,ten)
            obj.Te     = linspace(temin,temax,ten)';
            obj.TeStep = obj.Te(2)-obj.Te(1);
            % Memory Allocation
            obj.KrDS = zeros(obj.nTe,1);
        end
        
        % Kinematics
        function SetKrDSBinning(obj,TeMin,TeMax,nTeBinningFactor)
            obj.TeMin             = TeMin;
            obj.TeMax             = TeMax;
            obj.nTeBinningFactor  = nTeBinningFactor;
            obj.nTe               = floor(TeMax-TeMin)*nTeBinningFactor+1;
            obj.Te                = linspace(TeMin,TeMax,obj.nTe)';
            obj.TeStep            = obj.Te(2)-obj.Te(1);
            % Memory Allocation
            obj.KrDS              = zeros(obj.nTe,1);
        end
        
        function c = FindqUTeMatch(obj,tolerance)
            % Find match between Te and qU arrays
            % With a tolerance (1e-2 recommanded, 1 meV);
            % Return:
            %  - r: qU index array
            %  - c: qU matches in Te array
            out = abs(bsxfun(@minus,obj.qU(:),obj.Te(:).')) < tolerance;
            [r c] = find(out([1:1:obj.nqU],:,1)==1);
        end
        
        % Compute Kr DS/IS Spectra
        function SetFitBias(obj,flag)
            if flag == 0
                obj.E_Bias     = 0;

                obj.E_Bias_s3  = 0;
                obj.W_Bias     = 0;

                obj.L3_32_B    = obj.L3_32_B_i;
                obj.K32_B      = obj.K32_B_i;
                % Single Pixel
                if (obj.FPD_Pixel>0)
                    obj.B_Bias     = 0;
                    obj.Phi0_Bias  = 0;

                    % All Pixels
                end
            end
            
            if flag == 1
                
                % Add Bias
                if (obj.FPD_Pixel>0)
                    obj.BKG_RateSec = obj.BKG_RateSec_i + obj.B_Bias;
                end
                
                % L3_32
                obj.L3_32_E     = obj.L3_32_E_i    + obj.E_Bias;
                obj.L3_32_W     = obj.L3_32_W_i    + obj.W_Bias;
                obj.L3_32_E_s3  = obj.L3_32_E_s3_i + obj.E_Bias_s3;
                
                
                if (obj.FPD_Pixel>0)
                    obj.L3_32_Phi0  = obj.L3_32_Phi0_i + obj.Phi0_Bias;
                    obj.L3_32_Phi0_s3 = obj.L3_32_Phi0_s3_i + obj.Phi0_s3_Bias;
                end
                
                % K32
                obj.K32_E       = obj.K32_E_i      + obj.E_Bias;
                obj.K32_W       = obj.K32_W_i      + obj.W_Bias;
                if (obj.FPD_Pixel>0)
                    obj.K32_Phi0    = obj.K32_Phi0_i   + obj.Phi0_Bias;
		    

                end
                
                
            end
            
            % To use for the creation of different Spectra
            if flag == 2
                obj.L3_32_E = obj.E_Bias;
                obj.L3_32_W = obj.W_Bias;
                obj.L3_32_Phi0 = obj.Phi0_Bias;
                obj.BKG_RateSec = obj.B_Bias;
            end
        end
        
        function SetFitBiasallPixels(obj,flag)
            if flag == 0
                obj.E_Bias     = 0;
                obj.W_Bias     = 0;
            end
            
            if flag == 1
                obj.L3_32_E     = obj.L3_32_E_i    + obj.E_Bias;
                obj.L3_32_W     = obj.L3_32_W_i    + obj.W_Bias;
                obj.K32_E       = obj.K32_E_i      + obj.E_Bias;
                obj.K32_W       = obj.K32_W_i      + obj.W_Bias;
            end
        end
        
        %function          SetLineInitParallPixels(obj,phi0All,bgAll)
        function SetLineInitParallPixels(obj)
            % Initialize Line Amplitude / Background
            % Based on Pixel-by-Pixel Fit
            switch obj.TD
                case 'KrK32'
                    %obj.K32_Phi0allPixels     = phi0All;
                    obj.K32_Phi0allPixels      = [5.82628718       6.4286757      6.71639758      6.56315701      6.33290149      6.41701042      5.97436507      6.14602028      6.56148191      6.19344321      5.93257474      5.65381996      6.21320543      5.89449769      5.63022818      6.04370412      5.85948031      6.14821642      6.39563041      6.50371692      6.14831694      6.36712887      5.82546359      5.87035486      6.32731804      5.88735241      6.71222882      5.79079133      6.58247217      6.71621052      6.40714188      6.02533709      5.77897097      6.58000529      6.89317295      6.59911641      5.58028314      6.35409377       6.3497417      6.13370344      6.77962071      6.05710631      5.89266257      5.66067427      6.53469541      6.14776446      6.62915575      6.60711978       6.6949361       6.5645906      6.09649441       5.7204478      6.33493424       6.1695564      6.70080694      6.38812609      6.89622266      6.07385701      5.92775761      6.62880562       6.4312363      6.42296968       6.3113004      6.06871301      6.69419998      6.45056315      6.35933422      6.29066952      7.67969842      6.27551314      6.84176215      6.84659048      6.51523667      6.68165404      6.05769158      6.13066627      6.13099858      6.29601724      6.70877844      6.24240345      6.81332952      6.68289131      5.85069547      6.72549655      5.74476265      5.95824086      6.43193681      6.57775774      5.82993973      6.57886822      6.24471608      6.04261633       6.5553628      6.00479354      6.31783028      6.41248604       6.3513146      5.58809486      5.68926138      3.49046438      2.28453493       6.8061936      6.07413302      6.76719097        6.383045      6.26842968      6.63760354      6.15903123      5.88327695      6.90334992      6.68328787      4.76537155      1.99912696      3.72258561       5.2925069      5.96119904      6.42041988      6.60083894      5.98967546      6.06505053      6.17798006      6.22080614      6.81628629      1.48557244      2.21929603      21416.2961      52.7562104      3.09423532      1.87792368     0.124088358      4.19329348      6.04236964      6.06575283      6.20300717      5.58936436      52204.9111      4.68946374    0.0316168668      96650.4952      237.322488     0.141541299      3.50394603      2.98778905      5.95850886      6.15748575       3.8220212      1.40957399     -1.52188163];
                    obj.BKG_RateSecallPixels   = [3.4444      2.9013      2.9933      3.1836      3.1718      3.2189      3.2634      3.1204      2.8651      3.0362      3.0765      3.3151      3.0144      3.2387      3.2319      3.1543      3.1448      2.9433      3.0609      2.8784      2.9491      3.0975      3.1683      3.2471      3.1327      3.1001      2.6788      3.1097      2.7522       2.867       2.866      2.7495      3.1562      2.8297       2.746      3.0666      3.2601       3.019      2.8848      3.0356      2.7423       2.665      2.7002      2.6488      2.4879      2.7116      2.7438      2.8309      2.6771      2.8165      2.6841       2.846      2.3384      2.3102       2.083      2.1538      1.9535      2.2066      2.4962       2.538      2.7512      2.9306      2.5931      2.5801      1.8248      1.9219      1.9362      1.9718      1.7219      1.9709      2.0599      1.9818      2.4131      2.1378      2.2094      2.1218      1.8689      1.7798      1.7079      1.6499      1.5365      1.7775      1.8687      1.9534      2.0986      2.2066      2.0232      1.9007       1.699      1.4346      1.3636      1.3031      1.2524      1.6416      1.7008      1.7925      1.8105      2.1218      2.1589       1.425     0.97931      1.1211      1.1338      1.0101        1.04      1.1334      1.3544      1.5361      1.8423      1.6541      1.5722      1.4886     0.69856     0.50351     0.60041     0.74103      0.6409     0.77577      1.0915      1.3723      1.4853       1.634      1.2468     0.54325    0.073323     0.12734     0.26398      0.8208     0.63641     0.25674     0.52326     0.97061      1.2783      1.2178     0.88438     0.14931      0.4282     0.43534     0.33116     0.34383     0.61636     0.87393     0.48542      0.6858     0.78734     0.66582     0.23379     0.23758 ];
                case 'KrL3_32' % WRONG INIT
                    %obj.L3_32_Phi0allPixels   = phi0All;
                    obj.L3_32_Phi0allPixels    = [5.82628718       6.4286757      6.71639758      6.56315701      6.33290149      6.41701042      5.97436507      6.14602028      6.56148191      6.19344321      5.93257474      5.65381996      6.21320543      5.89449769      5.63022818      6.04370412      5.85948031      6.14821642      6.39563041      6.50371692      6.14831694      6.36712887      5.82546359      5.87035486      6.32731804      5.88735241      6.71222882      5.79079133      6.58247217      6.71621052      6.40714188      6.02533709      5.77897097      6.58000529      6.89317295      6.59911641      5.58028314      6.35409377       6.3497417      6.13370344      6.77962071      6.05710631      5.89266257      5.66067427      6.53469541      6.14776446      6.62915575      6.60711978       6.6949361       6.5645906      6.09649441       5.7204478      6.33493424       6.1695564      6.70080694      6.38812609      6.89622266      6.07385701      5.92775761      6.62880562       6.4312363      6.42296968       6.3113004      6.06871301      6.69419998      6.45056315      6.35933422      6.29066952      7.67969842      6.27551314      6.84176215      6.84659048      6.51523667      6.68165404      6.05769158      6.13066627      6.13099858      6.29601724      6.70877844      6.24240345      6.81332952      6.68289131      5.85069547      6.72549655      5.74476265      5.95824086      6.43193681      6.57775774      5.82993973      6.57886822      6.24471608      6.04261633       6.5553628      6.00479354      6.31783028      6.41248604       6.3513146      5.58809486      5.68926138      3.49046438      2.28453493       6.8061936      6.07413302      6.76719097        6.383045      6.26842968      6.63760354      6.15903123      5.88327695      6.90334992      6.68328787      4.76537155      1.99912696      3.72258561       5.2925069      5.96119904      6.42041988      6.60083894      5.98967546      6.06505053      6.17798006      6.22080614      6.81628629      1.48557244      2.21929603      21416.2961      52.7562104      3.09423532      1.87792368     0.124088358      4.19329348      6.04236964      6.06575283      6.20300717      5.58936436      52204.9111      4.68946374    0.0316168668      96650.4952      237.322488     0.141541299      3.50394603      2.98778905      5.95850886      6.15748575       3.8220212      1.40957399     -1.52188163];
                    obj.BKG_RateSecallPixels   = [3.4444      2.9013      2.9933      3.1836      3.1718      3.2189      3.2634      3.1204      2.8651      3.0362      3.0765      3.3151      3.0144      3.2387      3.2319      3.1543      3.1448      2.9433      3.0609      2.8784      2.9491      3.0975      3.1683      3.2471      3.1327      3.1001      2.6788      3.1097      2.7522       2.867       2.866      2.7495      3.1562      2.8297       2.746      3.0666      3.2601       3.019      2.8848      3.0356      2.7423       2.665      2.7002      2.6488      2.4879      2.7116      2.7438      2.8309      2.6771      2.8165      2.6841       2.846      2.3384      2.3102       2.083      2.1538      1.9535      2.2066      2.4962       2.538      2.7512      2.9306      2.5931      2.5801      1.8248      1.9219      1.9362      1.9718      1.7219      1.9709      2.0599      1.9818      2.4131      2.1378      2.2094      2.1218      1.8689      1.7798      1.7079      1.6499      1.5365      1.7775      1.8687      1.9534      2.0986      2.2066      2.0232      1.9007       1.699      1.4346      1.3636      1.3031      1.2524      1.6416      1.7008      1.7925      1.8105      2.1218      2.1589       1.425     0.97931      1.1211      1.1338      1.0101        1.04      1.1334      1.3544      1.5361      1.8423      1.6541      1.5722      1.4886     0.69856     0.50351     0.60041     0.74103      0.6409     0.77577      1.0915      1.3723      1.4853       1.634      1.2468     0.54325    0.073323     0.12734     0.26398      0.8208     0.63641     0.25674     0.52326     0.97061      1.2783      1.2178     0.88438     0.14931      0.4282     0.43534     0.33116     0.34383     0.61636     0.87393     0.48542      0.6858     0.78734     0.66582     0.23379     0.23758 ];
                case 'KrL3_32_HS'
                    %obj.L3_32_Phi0allPixels   = phi0All;
                    obj.L3_32_Phi0allPixels    = [56.5360514      57.3480234      56.8081639      56.6129319      56.6651322      56.7521378      56.4331566      56.7765422       57.139682      58.3203923      56.8035697       57.509378      57.3134507      56.3585074      57.1910465      56.4915403      56.8331837      57.7819882      57.5332718       57.099742      57.2539468      57.4333427      56.3805174      57.2877082      57.2493351      57.1763282      56.8984839      57.3429887      56.6084765      56.8360138      57.2879331       56.243831      56.5566012       56.433578      57.7451086      57.2440605      56.2422358      57.1043767      58.2373336      56.2438112      57.3310334      57.2233094      57.9683267      57.2562738      56.6051708      56.5921302      57.4050199      56.7143223       57.088038      57.0066145      56.8674904      56.9054365      56.9243977       56.853341      57.2926679      56.6583768      57.8361421      55.8483627      57.0842264      57.7685739      56.2571647      56.8828266      57.5296765      56.5998278       56.468245      56.9482438      57.4223088      56.9444115      57.0761593      56.7717449      57.4750167      58.3111842      56.8557925      58.5809177      57.6444866       57.098027      56.0523479      57.4647075      56.3255289      57.0641549      56.9155181      56.4277599      56.9948574      57.6738179      57.0001469      57.4781666      56.9106279      57.5659522      56.8692609       57.220538      56.5720869      57.2962624      58.7179734      57.7112876      57.6941043      57.7661375      57.2406036      57.3168744      58.7801178      33.7032266      21.2751451      57.5375493      57.8107991      57.2103688       57.766228      56.1421808      57.6496601      58.0498333      57.9669253      57.5510846       57.945322      47.7537967      20.6921401      31.3327892       47.956548      57.2190382      56.4510299      56.4505047      57.5335541      57.0622262      57.4642475      57.9871565      57.4842943      13.4053608      7.63893862     0.622593967      1.95900786      27.3466029      17.7221708      2.74630054      37.1536122      56.8468545      57.0864788      57.4134843      49.1118462      10.6281815      -900.90334    0.0820601151      1408.90234   -0.0783445208      346.301721      2.95014592      25.0924102      57.5671276        55.41964      38.2596145       13.169156      16315.5142];
                    obj.BKG_RateSecallPixels   = [27.028       27.035      26.7847      26.6815      26.7719      26.7425      27.0082       27.049      26.6303      26.4654      26.8427      26.5389      26.6843      27.0508      26.8212      26.8827      26.9802       26.814      26.5839       26.917      26.9161      26.8336      27.0479      26.8629      26.6802      26.6738      26.7857      26.7144      26.7131      26.6946      26.2719       26.503      26.7376       26.909      26.3921      26.6299      26.9731      26.8272      26.6449      26.8906      26.0948      25.9395      25.5633      25.9772      26.3128      26.3866      26.5622      26.8294      26.6298      26.8163      26.5302      26.7241      25.5653      25.6414      25.0528      25.1815      25.1884      25.4453      25.7088      25.7313      26.5097      26.3861      26.1544       26.138      24.6754      23.9798      23.4669      23.6876      23.8851      24.6901      24.9391       25.154      25.5291      25.1862      25.3653      25.1967      23.6484      22.4949      21.9366      21.5055      21.6197      22.3821      23.2923      23.8047      24.8789      24.5595      24.9828      23.9907      20.6212      19.3194      18.5695      18.7807      18.9962      20.3384      21.6166      22.8306      23.2777       23.323      22.6612      16.5671      14.1386      16.3435      14.9948      14.4406      14.6788      16.3339      18.3081      19.8449      20.8193      21.4299       20.779      13.9804      11.2238      9.59267      8.70241      8.17437      9.10707       11.868      14.9024      17.7692      18.9467      18.6027      16.5375      8.33971       5.0646      2.55393      2.25583      2.88224      2.47342      2.44123      6.48736      12.0344      14.7303      15.4808      12.9101      7.89288      1.40932      1.16804     0.716186     0.936746      1.66727      3.37436      3.89498      8.21662      10.1593      8.52242      6.07942      2.06325];            
            end
                %obj.BKG_RateSecallPixels   = bgAll;
        end
        
        function ComputeKrDS(obj,varargin)
            %
            % Compute Differential Spectrum of a Krypton Line
            %
            
            % Inputs
            p = inputParser;
            p.addParameter('E_Bias',0,@(x)isfloat(x));
            p.addParameter('W_Bias',0,@(x)isfloat(x));
            p.addParameter('B_Bias',0,@(x)isfloat(x));
            p.addParameter('Phi0_Bias',0,@(x)isfloat(x));

            p.addParameter('E_Bias_s3',0,@(x)isfloat(x));%satellite
            p.addParameter('Phi0_s3_Bias',0,@(x)isfloat(x));    

            p.parse(varargin{:});
            obj.E_Bias                 = p.Results.E_Bias;
            obj.E_Bias_s3              = p.Results.E_Bias_s3;

            obj.W_Bias                 = p.Results.W_Bias;
            obj.B_Bias                 = p.Results.B_Bias;
            obj.Phi0_Bias              = p.Results.Phi0_Bias;
            obj.Phi0_s3_Bias           = p.Results.Phi0_s3_Bias;
            
            % Initialization
            SetFitBias(obj,1);  
            
            % Kinemtics

            obj.KrDS = zeros(obj.nTe,1); obj.KrDSE = obj.KrDS;
            
            if strcmp(obj.ConvFlag,'DopplerEffect')


                    % The usual value for the standard deviation is already set in the class,
                    % but if the parameters change, the lines below recalculate it.
                    if strcmp(obj.recomputeKernel,'ON')
                        obj.computeKernel();
                    end
                    % Enlargement of the Spectrum to avoid convolution edge effects.
                    obj.enlargeSpectrum();
             end    
        
            %Defining Lorentzfunction for later use
            obj.KrLine = @(B,Phi0,GammaH,Eh,Ec) Phi0*...
                (GammaH/(2*pi))./((Ec-Eh).^2+GammaH.^2/4);
     

            switch obj.TD
                case 'KrK32'
                    E_temp = obj.K32_E; W_temp = obj.K32_W;
                    Phi0_temp = obj.K32_Phi0; B_temp = obj.K32_B;
%                     if strcmp(obj.DopplerEffectFlag,'OFF')
%                         z = (obj.Te-obj.K32_E + 1i.*obj.K32_W)/(obj.StandDev*1.414213562373095); %2^0.5
%                         obj.tmpline = obj.K32_Phi0*real(fadf(z))/(obj.StandDev*2.506628274631000); %(2*pi)^0.5
%                     elseif strcmp(obj.DopplerEffectFlag,'ON')
%                         obj.tmpline = (obj.KrLine(obj.K32_B,obj.K32_Phi0,obj.K32_W,obj.K32_E,obj.Te));
%                     end
                case {'KrL3_32','KrL3_32_HS','MoSFitterJuly2017','KrL3_32_Satellites'}
                    E_temp = obj.L3_32_E; W_temp = obj.L3_32_W;
		            E_temp_s3= obj.L3_32_E_s3; %energy of satellite Line 1 (+Bias) 
                    Phi0_temp = obj.L3_32_Phi0; B_temp = obj.L3_32_B;
                    Phi0_temp_s3= obj.L3_32_Phi0_s3;
                    
%                     if strcmp(obj.ConvFlag,'Voigt')
%                         z = (obj.Te-obj.L3_32_E + 1i.*obj.L3_32_W)/(obj.StandDev*1.414213562373095); %2^0.5
%                         obj.tmpline = obj.L3_32_Phi0*real(fadf(z))/(obj.StandDev*2.506628274631000); %(2*pi)^0.5
%                     elseif strcmp(obj.ConvFlag,'DopplerEffect')
%                         obj.tmpline = (obj.KrLine(obj.L3_32_B,obj.L3_32_Phi0,obj.L3_32_W,obj.L3_32_E,obj.Te));
%                     end
            end
            
            switch obj.ConvFlag
		case {'OFF','DopplerEffect'}

                  if strcmp(obj.MultiPeaksFlag,'OFF')  %create only 1 line  
                      
                    obj.tmpline = (obj.KrLine(B_temp,Phi0_temp,W_temp,E_temp,obj.Te));
                    
                  else %create multiple lines 
                     obj.tmpline = (obj.KrLine(B_temp,Phi0_temp,W_temp,E_temp,obj.Te));
                     obj.tmplineS3 = (obj.KrLine(B_temp,Phi0_temp_s3,W_temp,E_temp_s3,obj.Te));
                 end
   		case 'Voigt'
                    z = (obj.Te-E_temp + 1i.*W_temp)/(obj.StandDev*1.414213562373095*pi); %2^0.5
                    obj.tmpline = Phi0_temp*real(fadf(z))/(obj.StandDev*2.506628274631000); %(pi)^0.5
            end

            obj.KrDS   = obj.tmpline;
	        obj.KrDS_s3 = obj.tmplineS3;
            %             %obj.KrDS  = obj.tmpline .* obj.NormFactorKrDS;
            %             %obj.KrDS  = (obj.KrDS./trapz(obj.Te,obj.KrDS)) .* obj.NormFactorKrDS;
            obj.KrDSE  = obj.KrDS.^0.5;
            obj.KrDSE_s3= obj.KrDS_s3.^0.5;

            if strcmp(obj.ConvFlag,'DopplerEffect')
                    % Do the convolution
                    obj.KernelSpectrumConv();
                    % Restore Spectrum to the original size.
      		    obj.restoreSpectrum();
            end 
  
	    switch obj.MultiPeaksFlag
                case 'ON'
            	       obj.KrDS= (Phi0_temp*obj.KrDS + Phi0_temp_s3* obj.KrDS_s3 );%+...
                       %obj.p_sat1 * obj.KrDS_s1 + obj.p_sat2 * obj.KrDS_s2 );
                end
            end

        function ComputeKrDSallPixels(obj,varargin)
            %
            % Compute Differential Spectrum of a Krypton Line
            % Simultaneously for all pixels
            % (does not work with doppler yet)
            %
            
            % Inputs
            p = inputParser;
            p.addParameter('E_Bias',0,@(x)isfloat(x));
            p.addParameter('W_Bias',0,@(x)isfloat(x));
            %p.addParameter('L3_32_Phi0allPixels',100*ones(1,obj.nPixels));
            %p.addParameter('K32_Phi0allPixels',100*ones(1,obj.nPixels));
            
            p.parse(varargin{:});
            obj.E_Bias                 = p.Results.E_Bias;
            obj.W_Bias                 = p.Results.W_Bias;
            %obj.L3_32_Phi0allPixels    = p.Results.L3_32_Phi0allPixels;
            %obj.K32_Phi0allPixels      = p.Results.K32_Phi0allPixels;
            
            % Initialization
            SetFitBiasallPixels(obj,1);
            SetLineInitParallPixels(obj);
            
            switch obj.DopplerEffectFlag
                case 'ON'
                    % In case one wants to recalcute the standard deviation
                    % of the Gaussian
                    if strcmp(obj.recomputeKernel,'ON')
                        obj.computeKernel();
                    end
                    % Enlargement of the Spectrum to avoid convolution edge
                    % effects.
                    obj.enlargeSpectrum();
            end
            
            % Line Definition
            obj.KrLine = @(B,Phi0,GammaH,Eh,Ec) Phi0*(GammaH/(2*pi))./((Ec-Eh).^2+GammaH.^2/4);
            
            % Loop on All Pixels
            switch obj.TD
                case 'KrK32'
                    obj.KrDSallPixels = arrayfun(@(x) obj.KrLine(obj.K32_B,x,obj.K32_W,obj.K32_E,obj.Te),obj.K32_Phi0allPixels,'UniformOutput',false);
                    %                     for i=1:1:obj.nPixels
                    %                         obj.tmpline = (obj.KrLine(obj.K32_B,obj.K32_Phi0allPixels(i),obj.K32_W,obj.K32_E,obj.Te));
                    %                         obj.KrDSallPixels{i}      = obj.tmpline;
                    %                     end
                case {'KrL3_32','KrL3_32_HS','MoSFitterJuly2017'}
                    % disp(obj.L3_32_Phi0allPixels);
                    obj.KrDSallPixels = arrayfun(@(x) obj.KrLine(obj.L3_32_B,x,obj.L3_32_W,obj.L3_32_E,obj.Te),obj.L3_32_Phi0allPixels,'UniformOutput',false);
                    %                     for i=1:1:obj.nPixels
                    %                         obj.tmpline = (obj.KrLine(obj.L3_32_B,obj.L3_32_Phi0allPixels(i),obj.L3_32_W,obj.L3_32_E,obj.Te));
                    %                         obj.KrDSallPixels{i}       = obj.tmpline;
                    %                     end
            end
            %obj.KrDSallPixels = cell2mat(obj.KrDSallPixels);
            obj.KrDSEallPixels = cell2mat(obj.KrDSallPixels).^.5;
            
            switch obj.DopplerEffectFlag
                case 'ON'
                    % Do the convolution
                    obj.KernelSpectrumConvAllPixels();
                    % Restore Spectrum to the original size.
                    obj.restoreSpectrum();
            end
        end
        
        function ComputeKrIS(obj,varargin)
            % Inputs
            p = inputParser;
            p.addParameter('IStype','TRAPEZ',@(x)ismember(x,{'TRAPEZ','ACCURATE','TRAPEZRipples'}));
            p.parse(varargin{:});
            IStype  = p.Results.IStype;
            
            obj.KrIS = zeros(obj.nqU,1); obj.KrISE = obj.KrIS;
            switch IStype
                case 'ACCURATE' % Real integral --> slow but accurate
                    for qUi = 1:1:obj.nqU
                        KrDSf = @(e) interp1(obj.Te,obj.KrDS.*obj.KTF(obj.Te,obj.qU(qUi),obj.MACE_R_eV),e);
                        obj.KrIS(qUi) = ...
                            (integral(KrDSf,obj.qU(qUi),obj.qUmax)+obj.BKG_RateSec)...
                            *obj.TimeSec.*obj.qUfrac(qUi);
                        obj.KrISE(qUi) = obj.KrIS(qUi).^.5;
                    end
                case 'TRAPEZ'  % Trapezoidal integral --> faster, less accurate
                    for qUi = 1:1:obj.nqU
                        %qUieq = min(find(abs(obj.Te-obj.qU(qUi))<1e-2, 1));
                        qUieq = obj.qUTeMatch(qUi);
                        if qUieq < obj.nTe
                            obj.KrIS(qUi) = ...
                                (trapz(obj.Te(qUieq:end),obj.KrDS(qUieq:end).*obj.KTF(obj.Te(qUieq:end),...
                                obj.qU(qUi),obj.MACE_R_eV))).*obj.TimeSec.*obj.qUfrac(qUi)+obj.BKG_RateSec;
                            obj.KrISE(qUi) = obj.KrIS(qUi).^.5;
                        else
                            obj.KrIS(qUi)  = obj.BKG_RateSec*obj.TimeSec.*obj.qUfrac(qUi);
                            obj.KrISE(qUi) = obj.KrIS(qUi).^.5;
                        end
                    end
                case 'NOLOOP'
                    
                case 'TRAPEZRipples'  % Trapezoidal integral + HV Ripples
                    x = -pi/2:pi/10:pi/2; %
                    for qUi = 1:1:obj.nqU
                        qUieq = obj.qUTeMatch(qUi);
                        if qUieq <= obj.nTe
                            tmp = (trapz(obj.Te,obj.KrDS.*obj.KTF(obj.Te,obj.qU(qUi)+obj.HVRipplesP2PV/2.*sin(x),obj.MACE_R_eV))+obj.BKG_RateSec).*...
                                obj.TimeSec.*obj.qUfrac(qUi);
                            %fprintf(2,'XXX %g %g \n', obj.qU(qUi),mean(obj.qU(qUi)+obj.HVRipplesP2PV/2.*sin(x)));
                            obj.KrIS(qUi)  = mean(tmp);
                            %                            obj.KrIS(qUi)  = trapz(x,tmp)/(2);
                            obj.KrISE(qUi) = obj.KrIS(qUi).^.5;
                        else
                            obj.KrIS(qUi)  = obj.BKG_RateSec*obj.TimeSec.*obj.qUfrac(qUi);
                            obj.KrISE(qUi) = obj.KrIS(qUi).^.5;
                        end
                    end
            end
            
            switch obj.CPS
                case 'ON'
                    obj.KrISE  = obj.KrISE  ./ ((obj.TimeSec*obj.qUfrac)); %BEFORE KrIS!!!
                    obj.KrIS   = obj.KrIS   ./ (obj.TimeSec*obj.qUfrac );
            end
            
        end

        function f      = ComputeKrISf(obj,eb,wb,pb,bb,esb)

            % Function for Matlab Fit
            
            ComputeKrDS(obj,...
                'E_bias',eb,...
                'W_bias',wb,...
                'Phi0_bias',pb,...
                'B_bias',bb,...
		'E_bias_s3',esb);


            ComputeKrIS(obj,'IStype','TRAPEZ');
            f = obj.KrIS;
            
        end
        
        function ComputeKrISallPixels(obj,varargin)
            %
            % Compute Krypton Integral Spectra
            % Simultaneously for all pixels
            %
            
            % Inputs
            p = inputParser;
            p.addParameter('IStype','TRAPEZRipples',@(x)ismember(x,{'TRAPEZ','ACCURATE','TRAPEZRipples'}));
            %p.addParameter('BKG_RateSecallPixels',1*ones(1,obj.nPixels));
            
            p.parse(varargin{:});
            IStype  = p.Results.IStype;
            %obj.BKG_RateSecallPixels   = p.Results.BKG_RateSecallPixels;
            SetLineInitParallPixels(obj);
            
            switch IStype
                case 'TRAPEZ'  % Trapezoidal integral --> faster, less accurate
                    for i=1:1:obj.nPixels
                        localKTF  = obj.KTFallPixels{i};
                        localKrDS = obj.KrDSallPixels{i};
                        for qUi = 1:1:obj.nqU
                            qUieq = obj.qUTeMatch(qUi);
                            
                            if qUieq < obj.nTe
                                obj.KrISallPixels{qUi,i} = ...
                                    (trapz(obj.Te(qUieq:end),localKrDS(qUieq:end).*localKTF(obj.Te(qUieq:end),obj.qU(qUi),obj.MACE_R_eV))+obj.BKG_RateSecallPixels(i)).*...
                                    obj.TimeSec.*obj.qUfrac(qUi);
                                obj.KrISEallPixels{qUi,i} = obj.KrISallPixels{qUi,i}.^.5;
                            else
                                obj.KrISallPixels{qUi,i}  = obj.BKG_RateSecallPixels(i)*obj.TimeSec.*obj.qUfrac(qUi);
                                obj.KrISEallPixels{qUi,i} = obj.KrISallPixels{qUi,i}.^.5;
                            end
                            
                            switch obj.CPS
                                case 'ON'
                                    obj.KrISEallPixels{qUi,i}  = cell2mat(obj.KrISEallPixels(qUi,i))  ./ ((obj.TimeSec*obj.qUfrac(qUi))); %BEFORE KrIS!!!
                                    obj.KrISallPixels{qUi,i}   = cell2mat(obj.KrISallPixels(qUi,i))   ./ (obj.TimeSec*obj.qUfrac(qUi));
                                    %fprintf('CPS - %.0f - %.0f - %.3f \n',qUi,i,(obj.KrISallPixels{qUi,i}));
                            end
                        end
                        
                        %                 switch obj.CPS
                        %                     case 'ON'
                        %                       obj.KrISEallPixels(:,i)  = (cellfun(@(x,y) x./y,{cell2mat(obj.KrISEallPixels(:,i))},{obj.TimeSec*obj.qUfrac},'UniformOutput',false));
                        %                       obj.KrISallPixels(:,i)   = (cellfun(@(x,y) x./y,{cell2mat(obj.KrISallPixels(:,i))},{obj.TimeSec*obj.qUfrac},'UniformOutput',false));
                        %                 end
                        
                    end % Pixel Loop
                case 'TRAPEZRipples'  % Trapezoidal integral + HV Ripples
                    for i=1:1:obj.nPixels
                        localKTF = obj.KTFallPixels{i};
                        localKDS = obj.KrDSallPixels{i};
                        x = -pi/2:pi/10:pi/2; %
                        for qUi = 1:1:obj.nqU
                            qUieq = obj.qUTeMatch(qUi);
                            if qUieq <= obj.nTe
                                tmp = (trapz(obj.Te,localKDS.*localKTF(obj.Te,obj.qU(qUi)+obj.HVRipplesP2PV/2.*sin(x),obj.MACE_R_eV))+obj.BKG_RateSecallPixels(i)).*...
                                    obj.TimeSec.*obj.qUfrac(qUi);
                                obj.KrISallPixels{qUi,i}  = mean(tmp);
                                obj.KrISEallPixels{qUi,i} = obj.KrISallPixels{qUi,i}.^.5;
                            else
                                obj.KrISallPixels{qUi,i}  = obj.BKG_RateSec*obj.TimeSec.*obj.qUfrac(qUi);
                                obj.KrISEallPixels{qUi,i} = obj.KrISallPixels{qUi,i}.^.5;
                            end
                            switch obj.CPS
                                case 'ON'
                                    obj.KrISEallPixels{qUi,i}  = cell2mat(obj.KrISEallPixels(qUi,i))  ./ ((obj.TimeSec*obj.qUfrac(qUi))); %BEFORE KrIS!!!
                                    obj.KrISallPixels{qUi,i}   = cell2mat(obj.KrISallPixels(qUi,i))   ./ (obj.TimeSec*obj.qUfrac(qUi));
                            end
                        end
                    end
            end
        end
        
        %         function f      = ComputeTBDISf(obj,qu,mb,qb,bb,nb,gsb,esb)
        %             % Function for Matlab Fit
        %
        %             ComputeTBDDS(obj,...
        %                 'm_bias',mb,...
        %                 'Q_bias',qb,...
        %                 'N_bias',nb,...
        %                 'TTGS_bias',gsb,...
        %                 'TTES_bias',esb);
        %
        %             ComputeTBDIS(obj,'IStype','TRAPEZ');
        %
        %             qUieq = []; qUieq = (find(abs(obj.qU-qu)<1e-3));
        %             if qUieq>0
        %                 f = (1+nb).*obj.TBDIS(qUieq);
        %             else
        %                 f = -1;
        %                 return;
        %             end
        %         end
        %         function f      = ComputeTBDIS4f(obj,qu,mnuSqb,qb,bb,nb,gsb,esb,mnu4Sqb,sin2T4b)
        %             % Function for Matlab Fit
        %
        %             ComputeTBDDS(obj,...
        %                 'm_bias',mnuSqb,...
        %                 'Q_bias',qb,...
        %                 'N_bias',nb,...
        %                 'B_bias',bb,...
        %                 'TTGS_bias',gsb,...
        %                 'TTES_bias',esb,...
        %                 'mnu4Sq_Bias',mnu4Sqb,...
        %                 'sin2T4_Bias',sin2T4b );
        %
        %             ComputeTBDIS(obj,'IStype','TRAPEZ');
        %
        %             qUieq = []; qUieq = (find(abs(obj.qU-qu)<1e-3));
        %             if qUieq>0
        %                 f = (1+nb).*obj.TBDIS(qUieq);
        %             else
        %                 f = -1;
        %                 return;
        %             end
        %         end
        %
        
        % Normalize Spectrum
        function    ComputeNormFactorKrDS(obj)
            % Normalize Spectrum according to KATRIN-settings
            obj.TimeSec = obj.TimeYear*obj.Year2Sec;
            
            %             switch obj.Normalization
            %                 case {'NOMINAL','CPS'}
            %                     obj.NormFactorKrDS = obj.TimeSec;
            %             end
        end
        
        % Add Statistical Fluctuation
        function          AddStatFluctKrDS(obj)
            obj.KrDS   = obj.KrDS + sqrt(obj.KrDS).*randn(obj.nTe,1);
            %obj.KrDSE  = obj.KrDS.^0.5;
        end
        
        function          AddStatFluctKrIS(obj)
            obj.KrIS   = obj.KrIS + sqrt(obj.KrIS).*randn(obj.nqU,1);
            %obj.KrISE  = (obj.KrIS.^0.5);
        end
        
        function          AddStatFluctKrISallPixels(obj)
            obj.KrISallPixels   =  (cell2mat((obj.KrISallPixels)) + cell2mat(obj.KrISEallPixels).*randn(obj.nqU,obj.nPixels));
            obj.KrISallPixels   = mat2cell(obj.KrISallPixels,obj.nqU,obj.nPixels);
            %obj.KrISE  = (obj.KrIS.^0.5);
        end
        
        % Plots / Display
        function          PlotKrDS(obj,varargin)
            
            p = inputParser;
            p.addParameter('fign',999,@(x)isfloat(x) && x>0);
            p.addParameter('type','lin',@(x)ismember(x,{'lin','log'}));
            p.addParameter('pub',0,@(x)isfloat(x) && x>=0);
            p.parse(varargin{:});
            fign  = p.Results.fign;
            type  = p.Results.type;
            pub   = p.Results.pub;
            
            figure(fign)
            switch type
                case 'lin'
                    h = errorbar(obj.Te,obj.KrDS,obj.KrDSE,...
                        'ks','MarkerSize',3,'MarkerFaceColor',.8*[1 1 1],'LineWidth',1);
                    set(gca, 'YScale', 'lin');
                case 'log'
                    h = errorbar(obj.Te(obj.KrDS>0)-obj.Q,obj.KrDS(obj.KrDS>0),obj.KrDSE(obj.KrDS>0),...
                        'ks','MarkerSize',3,'MarkerFaceColor',.8*[1 1 1],'LineWidth',1);
                    set(gca, 'YScale', 'log')  ;
            end
            %errorbar_tick(h,1000);
            grid on
            xlabel('E (eV)','FontSize',14);
            ylabel('Counts per Energy Bin','FontSize',14);
            title('Krypton Gaseous Source','FontSize',14)
            set(gca,'FontSize',12);
            a = legend(h,'Differential Spectrum');% legend(a,'boxoff');
            PrettyFigureFormat;
            if pub>0
                publish_figure(fign,'./figs/TBD_phasespace.eps')
            end
            
        end
        
        function          PlotKrIS(obj,varargin)
            
            p = inputParser;
            p.addParameter('fign',999,@(x)isfloat(x) && x>0);
            p.addParameter('type','lin',@(x)ismember(x,{'lin','log'}));
            p.addParameter('pub',0,@(x)isfloat(x) && x>=0);
            p.addParameter('scalingErr',1,@(x)isfloat(x) && x>=0);
            p.addParameter('cps','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('CovMat','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            fign         = p.Results.fign;
            type         = p.Results.type;
            scalingErr   = p.Results.scalingErr;
            pub          = p.Results.pub;
            cps          = p.Results.cps;
            CovMat       = p.Results.CovMat;
            
            figure(fign)
            e = obj.qU(obj.KrIS>0);
            tmpis = obj.KrIS(obj.KrIS>0);
            tmpise = obj.KrISE(obj.KrIS>0)*scalingErr;
            switch CovMat
                case 'ON'
                    nCMTrials = 5000;
                    fName = sprintf('./data/CovMat/WGTSMACE_CovMat_%gTrials_%s_Year_%g.mat',...
                        nCMTrials,obj.TD,obj.TimeYear);
                    if exist(fName, 'file') == 2
                        CM=importdata(fName);
                    else
                        %obj.ComputeKrISCovarianceMatrix('nTrials',nCMTrials);
                        CM=importdata(fName);
                    end
                    tmpise_statsyst = sqrt((diag(CM.WGTSMACE_CovMat)) + (obj.KrIS));
                    tmpise_statsyst = tmpise(obj.KrIS>0)*scalingErr;
            end
            switch cps
                case 'OFF'
                    h1 = errorbar(e,tmpis,tmpise,'ks','MarkerSize',3,'MarkerFaceColor',.8*[1 1 1],'LineWidth',1);
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
                    tmpis = tmpis./obj.qUfrac(obj.KrIS>0)./obj.TimeSec;
                    h1 = errorbar(e,tmpis,1e-99*tmpis,'ks','MarkerSize',3,'MarkerFaceColor',.8*[1 1 1],'LineWidth',1);
                    ylabel('Counts per second','FontSize',14);
            end
            switch type
                case 'lin'
                    set(gca, 'YScale', 'lin');
                case 'log'
                    set(gca, 'YScale', 'log')  ;
            end
            axis([e(1) e(end) 0.5*min(tmpis) 1.2*max(tmpis)]);
            grid on
            xlabel('qU (eV)','FontSize',14);
            title('Krypton - KATRIN Integral Spectrum','FontSize',14)
            set(gca,'FontSize',12);
            switch CovMat
                case 'OFF'
                    a = legend('Stat. Uncertainty  '); legend('boxoff');
                case 'ON'
                    a = legend('Stat. Uncertainty  ','Stat. + Syst. Uncertainties '); legend('boxoff');
                    
            end
            PrettyFigureFormat();
            if pub>0
                publish_figure(fign,'./figs/Kr_Line.eps')
            end
            
            
        end
        
        function          DisplayKrInfo(obj)
            fprintf(2,'--------------------------------------------------------------\n');
            fprintf(2,'Krypton Spectrum:\n');
            fprintf(2,'--------------------------------------------------------------\n');
            switch obj.TD
                case {'KrL3_32','KrL3_32_HS','MoSFitterJuly2017'}
                    fprintf('   - Line Position: %g eV \n',obj.L3_32_E_i);
                    fprintf('   - Line Width: %g eV \n',obj.L3_32_W_i);
                    fprintf('   - Line Amplitude: %g cps \n',obj.L3_32_Phi0_i);
                    fprintf('   - Line BR: %g \n',obj.L3_32_B_i);
                case 'KrK32'
                    fprintf('   - Line Position: %g eV \n',obj.K32_E_i);
                    fprintf('   - Line Width: %g eV \n',obj.K32_W_i);
                    fprintf('   - Line Amplitude: %g cps\n',obj.K32_Phi0_i);
                    fprintf('   - Line BR: %g eV \n',obj.K32_B_i);
            end
            fprintf(2,'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n');
            fprintf('   - (Kinetic) Energy Range: [%g,%g] eV - %g bins \n',obj.TeMin,obj.TeMax,obj.nTe);
            fprintf('   - Te Binning Factor: %g\n',obj.nTeBinningFactor);
            fprintf(2,'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n');
            fprintf(2,'- - - - - WGTS \n');
            fprintf('   - cos(theta)    = %g\n',obj.WGTS_CosMaxAAngle);
            fprintf('   - Radius        = %g cm\n',obj.WGTS_FTR_cm);
            fprintf('   - B field       = %g T\n',obj.WGTS_B_T);
            switch obj.DopplerEffectFlag
                case 'ON'
                    fprintf('   - Doppler Effect: ON                     \n');
                    fprintf('    - Doppler Effect: T = %g K               \n',obj.T);
                    fprintf('    - Doppler Effect: Real Convolution  = %s \n',obj.realConv);
                    fprintf('    - Doppler Effect: Energy Broadening = %g eV\n',obj.StandDev);
                    fprintf('    - Doppler Effect: recomputeKernel   = %s   \n',obj.recomputeKernel);
                    fprintf('    - Doppler Effect: DeltaE = [-%g;+%g] bins\n',obj.Eminus,obj.Eplus);
                case 'OFF'
                    fprintf(2,'Doppler Effect: OFF                    \n');
            end
            fprintf(2,'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n');
            fprintf(2,'- - - - - MACE \n');
            fprintf('   - Bmax           = %g T\n',obj.MACE_Bmax_T);
            fprintf('   - Ba             = %g T\n',obj.MACE_Ba_T);
            fprintf('   - Resolution     = %g eV\n',obj.MACE_R_eV);
            fprintf(2,'- - - - -  \n');
            fprintf('   - Integral Spectrum - Time Distribution: %s \n',obj.TD);
            fprintf('   - Transmission Function: %s \n',obj.KTFFlag);
            fprintf('   - HV Ripples      : %s \n',obj.HVRipples);
            switch  obj.HVRipples
                case 'ON'
                    fprintf('    - HV Ripples Value: %g V (Peak-to-Peak, Sinusoid)\n',obj.HVRipplesP2PV);
            end
            fprintf(2,'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n');
            fprintf(2,'- - - - - Backgrounds \n');
            fprintf('   - Flag: %s  - Type: %s\n',obj.BKG_Flag,obj.BKG_Type);
            switch obj.FPD_Segmentation
                case 'OFF'
                    fprintf('   - Rate over all FPD: %g mcps per qU\n',obj.BKG_RateSec_i*1e3);
                case 'RING'
                case 'PIXEL'
                    fprintf('   - Rate per FPD pixel: %g mcps per qU\n',obj.BKG_RateSec_i*1e3);
            end
            fprintf(2,'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n');
            fprintf(2,'- - - - - FPD \n');
            fprintf('   - Segmentation: %s  \n',obj.FPD_Segmentation);
            switch obj.FPD_Segmentation
                case 'OFF'
                    fprintf('   - Mean Efficiency: %g  \n',obj.FPD_MeanEff);
                case 'RING'
                case 'PIXEL'
                    fprintf('   - Pixel: %g  (0=All)\n',obj.FPD_Pixel);
            end
            if ~isempty(obj.KrIS)
                fprintf('   - Krypton e- FPD: %g cps at (qU>%g eV)\n',obj.KrIS(1)/(obj.TimeSec.*obj.qUfrac(1)),obj.TeMin);
            end
            fprintf(2,'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n');
            fprintf(2,'- - - - - KATRIN \n');
            fprintf('   - Normalization   = %s \n',obj.Normalization);
            fprintf('   - CPS Mode        = %s \n',obj.CPS);
            fprintf('   - Time (in years) = %g \n',obj.TimeYear);
            
            fprintf(2,'--------------------------------------------------------------\n');
            disp(' ');
            
        end
        
        
        % [ADD] Doppler Effect
        function computeKernel(obj)
            % The kernel is calculate according to [citation]
            % Here A does not need the spectrum yet.
            
            % The standard deviations of the kernel funtion depending on
            % velocity (o_v) and energy (StandDev).
            obj.o_v = sqrt(obj.kb*obj.T/obj.M);
            obj.StandDev = obj.o_v*sqrt(2*obj.E_cms*obj.me/obj.c^2);
            
            % range for the kernel
            obj.E_lab = linspace(obj.E_cms - obj.E_range,...
                obj.E_cms + obj.E_range,obj.nTe);
            %             obj.E_lab = linspace(min(A.Te),...
            %                 max(A.Te),A.nTe);
            
            obj.v_m = (obj.E_lab - obj.E_cms)/((obj.me/obj.c^2)*obj.v_e);
            obj.cosO = linspace(cos(obj.Omax),1,obj.nTe);
            [v_M,cosO_M] = meshgrid(obj.v_m,obj.cosO);
            prefac = (1-cos(obj.Omax)).^(-1);
            
            % integration of the emission angle
            g_vm = prefac./(sqrt(2*pi)*obj.o_v)*...
                trapz(obj.cosO,exp(-0.5*((v_M-cosO_M*obj.u)/obj.o_v).^2));
            %             g_vm = zeros(1,length(obj.v_m));
            %             for vv = 1:length(obj.v_m)
            %                 g_vm(vv) = prefac.*0.5*1/obj.u*sum(erf(linspace((obj.v_m(vv)-cos(obj.Omax)*obj.u)/(sqrt(2)*obj.o_v),...
            %                 (obj.v_m(vv)-obj.u)/(sqrt(2)*obj.o_v),A.nTe)));
            %             end
            
            obj.g = 1/((obj.me/obj.c^2)*obj.v_e)*g_vm;
            obj.P = trapz(obj.E_lab - obj.E_cms,obj.g);
            g_fwhm = obj.g - 0.5*max(obj.g);
            g_fwhm_index = find(0 < g_fwhm);
            x1 = g_fwhm_index(1); xf = g_fwhm_index(end);
            x0_left = obj.E_lab(x1) - g_fwhm(x1)*(obj.E_lab(x1-1)-obj.E_lab(x1))/(g_fwhm(x1-1)-g_fwhm(x1));
            x0_right = obj.E_lab(xf) - g_fwhm(xf)*(obj.E_lab(xf+1)-obj.E_lab(xf))/(g_fwhm(xf+1)-g_fwhm(xf));
            FWHM = x0_right - x0_left;
            %FWHM = -(obj.E_lab(g_fwhm(1)) - obj.E_lab(g_fwhm(end)));
            %obj.StandDev_num = FWHM/(2*sqrt(2*log(2)));
            obj.StandDev = FWHM/(2*sqrt(2*log(2)));
        end % computeKernel
        
        function plotKernel(obj)
            figure(1);
            hold on
            plot(obj.E_lab-obj.E_cms,obj.g);
            leg_cell = {['u = ',num2str(obj.u),' m/s']};
            hold off
            xlabel('E_{lab}-E_{cms} [eV]');
            ylabel('probability');
            legend(leg_cell);
        end % plotKernel
        
        function enlargeSpectrum(obj)
            % add neccessary bins left and right
            nTeStepsMinus = obj.nTeBinningFactor*obj.Eminus*obj.TeStep;
            nTeStepsPlus = obj.nTeBinningFactor*obj.Eplus*obj.TeStep;
            obj.TeMin = obj.qUmin - nTeStepsMinus;
            obj.TeMax = obj.qUmax + nTeStepsPlus;
            obj.Te = [(obj.TeMin:obj.TeStep:obj.qUmin-obj.TeStep)'...
                ;obj.Te;(obj.qUmax+obj.TeStep:obj.TeStep:obj.TeMax)'];
            obj.nTe = length(obj.Te);
            % Memory Allocation
            obj.KrDS = zeros(obj.nTe,1);
            ComputeNormFactorKrDS(obj);
            InitializeBackground(obj);
            
        end % enlargeSpectrum
        
        function KernelSpectrumConv(obj)
            if strcmp(obj.realConv,'ON') % [Work in progress, don't use]
                %A.TBDDS = trapz(A.TBDDS.*);
                obj.KrDS = conv(obj.g,obj.KrDS);
            else % convolution using matrix
                obj.mate = obj.eres(obj.Te,obj.TeMin,obj.TeMax,obj.nTe,obj.StandDev);
                obj.KrDS = obj.mate*obj.KrDS;

	switch obj.MultiPeaksFlag
                    case 'ON'    
                obj.KrDS_s3 = obj.mate*obj.KrDS_s3;
                end 
            end

        end % KernelSpectrumConv
        
        function KernelSpectrumConvAllPixels(obj)
            if strcmp(obj.realConv,'ON') % ![Work in progress, don't use]
                %A.TBDDS = trapz(A.TBDDS.*);
                obj.KrDS = conv(obj.g,obj.KrDS);
            else % convolution using matrix
                obj.mate = obj.eres(obj.Te,obj.TeMin,obj.TeMax,obj.nTe,obj.StandDev);
                obj.KrDS = obj.mate*obj.KrDS;
            end
        end % KernelSpectrumConv
        
        function restoreSpectrum(obj)
            % Spectrum and error are restores to the original size
            obj.KrDS = obj.KrDS(obj.Te>=obj.qUmin & obj.Te<=obj.qUmax);
            obj.KrDSE = obj.KrDSE(obj.Te>=obj.qUmin & obj.Te<=obj.qUmax);
            obj.Te = obj.Te(obj.Te>=obj.qUmin & obj.Te<=obj.qUmax);
            
	switch obj.MultiPeaksFlag
         case 'ON'
             obj.KrDS_s3 = obj.KrDS_s3(obj.Te>=obj.qUmin & obj.Te<=obj.qUmax);
             obj.KrDSE_s3 = obj.KrDSE_s3(obj.Te>=obj.qUmin & obj.Te<=obj.qUmax);
     		end

            % properties are recalculated
            obj.nTe = length(obj.Te);
            obj.InitializeBackground();
            ComputeNormFactorKrDS(obj);
        end
        
        % [ADD] END
        
        
        
        % Miscellaneous
        function r      = W(obj,sigma_E)
            % ----------------------------------------------------------
            % Apply Energy Resolution Function
            % ----------------------------------------------------------
            r = 0.5*(erf((obj.TeMax-obj.Te)/sqrt(2)/sigma_E./sqrt(obj.Te))...
                -erf((obj.Te-obj.TeMin)/sqrt(2)/sigma_E./sqrt(obj.Te)));
        end
        
    end %methods
    
    methods(Static)
        
        % Fit
        function chi2 = MyChi2(p,X)
            % Gaussian Xi2
            % X = [ "column vector of bin centers"  | "column vector of data" | "column vector of uncertainties" | "... extra info" ]
            % p = [ p1; p2; ... ] column vector of parameter to fit
            nn = X(:,2)~=0; % indexes of non-null bins
            x = X(nn,1);    % bin centers of non-null bins
            y = X(nn,2);    % bin content of non-null bins
            z = X(nn,3);    % uncertainties bin content of non-null bins
            m =  MyModel(obj,p);
            chi2 = sum(( y(1:numel(x)) - m(1:numel(x))).^2 ./ z(1:numel(x)).^2);
            %chi2 = -2*sum( y - m + y.*log(m./y) );                      % Poisson deviance chi square.
            
        end
        
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
