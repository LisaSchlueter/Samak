
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
% Th. Lasserre
% CEA/Saclay - TUM - IAS - MPP
%
% Last Updated: Jan. 2018
%
% ----------------------------------------------------------------------- %


classdef Kr < handle & FPD & WGTSMACE
    
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
        e_vel   = 9.9223e+07;% [m/s] velocity of the electron at 30474 keV...

    end
    
    properties (Access=public, Hidden = true)
        % 83mKr Parameters
        Z = 36;                     % Krypton
        
        % L3_32 - single pixel
        L3_32_E_i     = 30472.5; L3_32_E;     % Mean Energy, eV
        L3_32_W_i     = 1.11;    L3_32_W;     % Width, eV
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
        % K32 - all pixels
        K32_Phi0allPixels_i  = 100*ones(1,148);  % Intensity All Pixels
        K32_Phi0allPixels;                       % Intensity All Pixels
        
        tmpline;                              % all lines stored
        KrLine;                               % Kr Line
        KrLineTail;                           % Kr Line Offset Tail
        KrLineTailallPixels;                  % Kr Line Offset Tail all Pixels
        
        % Pixel-wise qU vector - HVL332 data
        qUpixel;
        Tepixel;
        TeAllPixelsInc; % Increment matrix for Te
        
        % Test
        TMPTF;
        TMPTFallPixels;
        
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
        %nPixels;
        
        % Doppler Effect
        DopplerEffectFlag;% turns on and off the Doppler Effect
        Eminus;         % Extra bins added to the end of the spectrum to account for convolution edge effect
        Eplus;          % Extra bins added to the beginning of the spectrum to account for conv. edge effect
        e_velParal;     % [m/s] electron vel. component parallel to emission direction
        cosO;           % cosine of the emission angle
        recomputeKernel;% switch to recompute kernel
        DE_sigma_e_vel; % [eV] standard deviation for the kernel in terms of the electron velocity
        DE_sigma;       % [eV] standard deviation for the kernel in terms of the electron energy (usually around 130 meV)
        DE_Bias;
        
        % Broadening kernel g parameters
        BulkVelKr;      % [m/s] bulk velocity of krypton
        GaussianKernelDE;% broadening kernel
        ProbKernDE;     % probability under the curve of the kernel (should = 1)
        
        % Binning
        nTeBinningFactor;
        
        % Fit: Bias - Single Pixel
        E_Bias;
        W_Bias;
        B_Bias;
        Phi0_Bias;     
    end
    
    methods
        function obj = Kr(varargin)
            fprintf(2,'Processing Kr Constructor ...\n');
            
            p = inputParser;

            % Krypton
            p.addParameter('L3_32_Phi0_i',100,@(x)isfloat(x) && x>0);
            p.addParameter('K32_Phi0_i',100,@(x)isfloat(x) && x>0);
            p.addParameter('L3_32_Phi0allPixels',100*ones(1,148));
            p.addParameter('K32_Phi0allPixels',100*ones(1,148));
            p.addParameter('BKG_RateSecallPixels',1*ones(1,148));
            
            % Doppler Effect Flag
            p.addParameter('DopplerEffectFlag','OFF',...
                @(x)ismember(x,{'OFF','matConv','Voigt','numConv'})); 
            p.addParameter('Eminus',5,@(x)isfloat(x) && x>0);
            p.addParameter('Eplus',5,@(x)isfloat(x) && x>0);
            p.addParameter('recomputeKernel','OFF',@(x)strcmp(x,'ON') || strcmp(x,'OFF'));
            
            % Broadening kernel parameters
            p.addParameter('E_range',1,@(x)isfloat(x) && x>0);
            p.addParameter('BulkVelKr',0,@(x)isfloat(x));
            p.addParameter('DE_sigma',0.058,@(x)isfloat(x) && x>0); % 0.134429846285963   
  
            % Binning
            p.addParameter('nTeBinningFactor',4,@(x)isfloat(x) && x>0);
            
            % Parse unmatched parameters to WGTSMACE.m
            p.KeepUnmatched=1;
            p.parse(varargin{:});
            if( isempty(fieldnames(p.Unmatched))); Unmatched={}; else
                Unmatched = reshape(...
                    [fieldnames(p.Unmatched),struct2cell(p.Unmatched)]',...
                    [1,length(fieldnames(p.Unmatched))*2]);
            end
            obj=obj@WGTSMACE(Unmatched{:}); %Parse to superclass WGTSMACE.m
            
            % KATRIN: GENERAL SETTINGS
            obj.nTeBinningFactor    = p.Results.nTeBinningFactor;

            % Doppler Effect Flag
            obj.DopplerEffectFlag   = p.Results.DopplerEffectFlag;
            obj.Eminus              = p.Results.Eminus;
            obj.Eplus               = p.Results.Eplus;
            obj.recomputeKernel     = p.Results.recomputeKernel;
            
            % Krypton
            obj.K32_Phi0_i          = p.Results.K32_Phi0_i;
            obj.L3_32_Phi0_i        = p.Results.L3_32_Phi0_i;
            
            % Krypton all Pixels
           % obj.nPixels             = p.Results.nPixels;
            obj.L3_32_Phi0allPixels = p.Results.L3_32_Phi0allPixels;
            obj.K32_Phi0allPixels   = p.Results.K32_Phi0allPixels;
            obj.BKG_RateSecallPixels= p.Results.BKG_RateSecallPixels;
            
            % Broadening kernel parameters
            obj.BulkVelKr           = p.Results.BulkVelKr;
            obj.DE_sigma            = p.Results.DE_sigma;
          
            % Binning - To get from the TD's
            obj.SetKrDSBinning; %Change Binning granularity with BinningFactor
           
            % Init
            SetFitBias(obj,0); %Set Fit Bias to 0
            
            % for HV data multi-fit only
            if obj.FPD_Pixel==0
                FillTableqUpixelHVdata(obj);
                SetKrDSBinningAllPixels(obj);
                obj.TeAllPixelsInc = -(obj.Tepixel - [obj.Tepixel(2:end,:);obj.Tepixel(end,:)])  ; 
            end
            
            % FPD:   Initialization
            InitializeBackground(obj);
            
            switch obj.FPD_Segmentation %Build and Store Transmission Functions in Matrices
                case 'OFF'
                    obj.KTF = obj.GetMeanTFAllFPD(); %function handle: Transmission Function (TF)
                    for i=1:1:obj.nqU 
                        tmp = obj.KTF(obj.Te,obj.qU(i)); %TF for one qU
                        tmp(isnan(tmp))=0;
                        obj.TMPTF(i,:) = tmp; %Matrix with all TF (--> 1 array for each qU)
                    end
                case 'RING'
                    
                case 'PIXEL'
                    SetPixel_MACE_BaEa(obj,'ReadConfigFile','OFF')
                    if obj.FPD_Pixel == 0
                        fprintf(2,'Processing All Pixels\n');
                        obj.KTFallPixels = obj.GetTFPixelWiseFPD_All();
                        for p=1:1:obj.nPixels
                            for i=1:1:obj.nqU
                                tfloc  = obj.KTFallPixels{p}; % function handle for Pixel p
                                switch obj.HVRipples
                                    case 'OFF'

                                tmp    = tfloc(obj.Tepixel(:,p),obj.qUpixel(i,p));
                                tmp(isnan(tmp))=0;
                                obj.TMPTFallPixels(p,i,:) = tmp;
                                    case 'ON'
                               x = -pi/2:pi/10:pi/2; %
                               tmp    = tfloc(obj.Tepixel(:,p),obj.qUpixel(i,p)+obj.HVRipplesP2PV/2.*sin(x));
                               tmp = mean(tmp');
                               tmp(isnan(tmp))=0;
                               obj.TMPTFallPixels(p,i,:) = tmp; % 3D Matrix: TF(Pixel,qU,:)
                                end                              
                            end
                        end  
                        
                        obj.L3_32_Phi0allPixels  = obj.L3_32_Phi0allPixels_i;
                        obj.K32_Phi0allPixels    = obj.K32_Phi0allPixels_i;
                        obj.BKG_RateSecallPixels = obj.BKG_RateSecallPixels_i;
                    else
                        obj.KTF = obj.GetTFPixelWiseFPD(obj.FPD_Pixel);
                        for i=1:1:obj.nqU
                            tmp = obj.KTF(obj.Te,obj.qU(i));
                            tmp(isnan(tmp))=0; obj.TMPTF(i,:) = tmp;
                        end
                        fprintf(2,'Processing Pixels %g\n',obj.FPD_Pixel);
                    end
            end     
            
            % Compute Normalization Factor
            ComputeNormFactorKrDS(obj);
            
            % for HV data multi-fit only
            %FillTableqUpixelHVdata(obj); 
           
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

        function SetKrDSBinning(obj)
            % Rebinning Te Vector with increased granularity
            % nTeBinningFactor=1 - Te = qU 
            % nTeBinningFactor=2 - Te = 2 x qU
            % nTeBinningFactor=3 - Te = 4 x qU
            % nTeBinningFactor=4 - Te = 8 x qU

            if obj.nTeBinningFactor>0
                obj.Te            = obj.qU;
                obj.nTe           = obj.nqU;
            end
            if obj.nTeBinningFactor>=2
                te2               = obj.Te + (obj.Te(2)-obj.Te(1))/2;
                obj.Te            = sort([obj.Te ; te2],1);
                obj.nTe           = 2*obj.nqU;
            end
            if obj.nTeBinningFactor>=3
                te3               = obj.Te + (obj.Te(2)-obj.Te(1))/2;
                obj.Te            = sort([obj.Te ; te3 ],1);
                obj.nTe           = 4*obj.nqU;
            end
            if obj.nTeBinningFactor>=4
                te4               = obj.Te + (obj.Te(2)-obj.Te(1))/2;
                obj.Te            = sort([obj.Te ; te4],1);
                obj.nTe           = 8*obj.nqU;
            end
            if obj.nTeBinningFactor>4
                obj.Te = interp1(linspace(1,obj.nqU,obj.nqU)',obj.qU,linspace(1,obj.nqU,obj.nqU*obj.nTeBinningFactor-(obj.nTeBinningFactor-1))');
            end

            % Truncating Te values > qUmax
            obj.TeMin             = min(obj.qU);
            obj.TeMax             = max(obj.qU);
            obj.Te                = obj.Te(obj.Te<=obj.TeMax);
            obj.nTe               = length(obj.Te);
            obj.TeStep            = obj.Te(2)-obj.Te(1);

        end
        
        function SetKrDSBinningAllPixels(obj)
            % Rebinning Te Vector with increased granularity
            % nTeBinningFactor=1 - Te = qU
            % nTeBinningFactor=2 - Te = 2 x qU
            % nTeBinningFactor=3 - Te = 4 x qU
            % nTeBinningFactor=4 - Te = 8 x qU
            
            if obj.nTeBinningFactor>0
                obj.Tepixel            = obj.qUpixel;
                obj.nTe                = obj.nqU;
            end
            if obj.nTeBinningFactor>=2
                te2               = obj.Tepixel + [diff(obj.Tepixel); zeros(1,obj.nPixels)]/2;
                obj.Tepixel       = sort([obj.Tepixel ; te2],1);
                obj.nTe           = 2*obj.nqU;
            end
            if obj.nTeBinningFactor>=3
                te3               = obj.Tepixel + [diff(obj.Tepixel); zeros(1,obj.nPixels)]/2;
                obj.Tepixel       = sort([obj.Tepixel ; te3 ],1);
                obj.nTe           = 4*obj.nqU;
            end
            if obj.nTeBinningFactor>=4
                te4               = obj.Tepixel + [diff(obj.Tepixel); zeros(1,obj.nPixels)]/2;
                obj.Tepixel       = sort([obj.Tepixel ; te4],1);
                obj.nTe           = 8*obj.nqU;
            end
            if obj.nTeBinningFactor>4
                % fprintf('ERROR: INVALID Te Binning Scale Factor! \n');
                 return;
                %obj.Tepixel = 
                %interp1(linspace(1,obj.nqU,obj.nqU)',obj.qU,linspace(1,obj.nqU,obj.nqU*obj.nTeBinningFactor-(obj.nTeBinningFactor-1))');
            end
            
            % Truncating Te values > qUmax
%             obj.TeMin             = min(obj.qU);
%             obj.TeMax             = max(obj.qU);
%             obj.Te                = obj.Te(obj.Te<=obj.TeMax);
%             obj.nTe               = length(obj.Te);
%             obj.TeStep            = obj.Te(2)-obj.Te(1);
        end
        
        
        %         function c = FindqUTeMatch(obj,tolerance)
        %             % Find match between Te and qU arrays
        %             % With a tolerance (1e-2 recommanded, 1 meV);
        %             % Return:
        %             %  - r: qU index array
        %             %  - c: qU matches in Te array
        %             out = abs(bsxfun(@minus,obj.qU(:),obj.Te(:).')) < tolerance;
        %             [r,c] = find(out([1:1:obj.nqU],:,1)==1);
        %         end
        
        % Compute Kr DS/IS Spectra
        function SetFitBias(obj,flag)
            if flag == 0
                obj.E_Bias     = 0;
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
                    obj.L3_32_Phi0  = obj.L3_32_Phi0_i + obj.Phi0_Bias; % L3_32
                    obj.K32_Phi0    = obj.K32_Phi0_i   + obj.Phi0_Bias; % K32
                end
                
                % L3_32
                obj.L3_32_E     = obj.L3_32_E_i    + obj.E_Bias;
                obj.L3_32_W     = obj.L3_32_W_i    + obj.W_Bias; 
                
                % K32
                obj.K32_E       = obj.K32_E_i      + obj.E_Bias;
                obj.K32_W       = obj.K32_W_i      + obj.W_Bias;      
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
                obj.Phi0_Bias  = zeros(obj.nPixels,1)'; % row vector
                obj.B_Bias     = zeros(obj.nPixels,1)';
            end
            
            if flag == 1
                obj.L3_32_E     = obj.L3_32_E_i    + obj.E_Bias;
                obj.L3_32_W     = obj.L3_32_W_i    + obj.W_Bias;
                obj.L3_32_Phi0allPixels(1:obj.nPixels)  = obj.L3_32_Phi0allPixels_i(1:obj.nPixels)  + obj.Phi0_Bias(1:obj.nPixels);
                obj.BKG_RateSecallPixels(1:obj.nPixels) = obj.BKG_RateSecallPixels_i(1:obj.nPixels) + obj.B_Bias(1:obj.nPixels);
                obj.K32_E       = obj.K32_E_i      + obj.E_Bias;
                obj.K32_W       = obj.K32_W_i      + obj.W_Bias;
            end
        end
        
        function SetLineInitParallPixels(obj)
            % Initialize Line Amplitude / Background
            % Based on Pixel-by-Pixel Fit
            switch obj.TD
                case {'HVL332'}
                    % Tail Good V1
                    %obj.L3_32_Phi0allPixels    = [44.6232      47.1541      46.3954      46.6741       45.761      46.6343      47.0666      45.9541      46.2233      45.9166      45.9511      45.5604      47.8608      45.5325      45.3468      46.2513      47.6844      48.7293      47.6418      46.6697      46.9427      50.1574      46.6698       44.659      46.4393       46.455      47.6318      47.1173       46.598      45.5871      46.1051       47.319      46.3033      46.8243      46.2841      47.1165      45.7393       46.505      44.4374      47.6845 zeros(1,108)];
                    % V2
                    %obj.L3_32_Phi0allPixels    = [45.1084      46.8574      46.3138      47.2566      45.7832      46.5892       49.269      44.6646      46.2282      45.8147      47.8512      47.2613      49.0241      46.7111      44.9991      46.0905      47.4269      47.5757      46.5872      44.9864      47.8579      50.1055      44.5695       44.799       44.742      45.6321      47.7861      47.6578      46.5218      44.9097      46.1765      46.2518      45.3604       46.354      48.7797      45.6433      44.2685      44.6916      45.1398       47.957 zeros(1,108)];
                    
                    % Tail Good V1
                    obj.BKG_RateSecallPixels_i   = [21.646      20.1148      20.1998      20.3904      21.3103      20.2357      19.2365      20.2849      20.4684      20.9596      20.6376       21.125      19.9332      20.5295      21.5831      20.7355      20.4928      19.1473      19.2369      19.6851       19.459      18.6784       20.078      20.5546      20.7442      20.4938       19.193      20.4008      19.9845      21.0239      19.9944      19.3057      20.1331      19.8468      20.8075      19.8936      20.3485      20.5098      21.5158      19.6338 zeros(1,108)];
                    %V2
                    %obj.BKG_RateSecallPixels   = [21.2143      20.1126      19.8434        19.85      21.4745      20.1822      17.6302      20.7845      20.5311      20.6216      19.3095      19.4966      19.0372      19.3769      21.5415      20.3688      20.4191      19.2642      19.8555      20.6339      18.5487      18.4958      21.3516      19.7295      21.3221      20.6244      18.9621      19.6201      19.6839      21.0746      19.4884       19.949       20.851      19.9666      19.4846      20.9079      21.0319      21.5459      20.7995      18.4634 zeros(1,108)];
                    
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
        
        function FillTableqUpixelHVdata(obj)
            % Build qU vector for all pixel
            % patch for HV data - C. Weinheimer cuisine ...
            switch obj.TD
                case {'HVL332'}
                    [~, hvdatamap] = read_hvdata('Display','OFF','MK35',1972.4531);
                    obj.qUpixel = squeeze(hvdatamap(1:obj.nPixels,3,1:end))';
            end
                        %                     obj.qUpixel(1,:) =    [30462.9055      30462.8976      30463.4316      30463.9247       30464.434       30464.949      30465.4961      30465.9465       30466.439      30466.9254      30467.4224      30467.9736       30468.439      30468.9663      30469.4845      30470.0545      30470.4993      30471.0092      30471.4719      30472.0163      30472.4638      30472.9619       30473.501      30473.9789      30474.4787      30474.9679      30475.4807      30476.0231       30476.533      30476.9967      30477.5425      30478.0252];
                        %                     obj.qUpixel(2,:) =    [30462.9012      30462.8933      30463.4272      30463.9204      30464.4294      30464.9447      30465.4916      30465.9421      30466.4344      30466.9211      30467.4181       30467.969      30468.4345      30468.9619      30469.4799        30470.05      30470.4947      30471.0046      30471.4676       30472.012      30472.4593      30472.9573      30473.4966      30473.9745      30474.4742      30474.9633      30475.4762      30476.0186      30476.5285      30476.9924      30477.5382      30478.0206];
                        %                     obj.qUpixel(3,:) =    [30462.8923      30462.8844      30463.4184      30463.9115      30464.4208      30464.9358      30465.4829      30465.9332      30466.4258      30466.9122      30467.4092      30467.9603      30468.4258      30468.9531      30469.4712      30470.0413      30470.4861      30470.9959      30471.4587      30472.0031      30472.4506      30472.9487      30473.4877      30473.9657      30474.4655      30474.9547      30475.4675      30476.0099      30476.5198      30476.9835      30477.5293       30478.012];
                        %                     obj.qUpixel(4,:) =    [30462.8955      30462.8876      30463.4215      30463.9146      30464.4237      30464.9389      30465.4859      30465.9364      30466.4287      30466.9153      30467.4124      30467.9633      30468.4288      30468.9562      30469.4742      30470.0442       30470.489      30470.9989      30471.4618      30472.0062      30472.4536      30472.9516      30473.4909      30473.9688      30474.4684      30474.9576      30475.4705      30476.0129      30476.5228      30476.9867      30477.5325      30478.0149];
                        %                     obj.qUpixel(5,:) =    [30462.9548      30462.947       30463.4809      30463.974      30464.4831      30464.9983      30465.5453      30465.9958      30466.4881      30466.9747      30467.4718      30468.0227      30468.4882      30469.0156      30469.5336      30470.1036      30470.5484      30471.0583      30471.5212      30472.0656       30472.513       30473.011      30473.5503      30474.0282      30474.5278       30475.017      30475.5298      30476.0722      30476.5821       30477.046      30477.5918      30478.0743];
                        %                     obj.qUpixel(6,:) =    [30462.9612      30462.9533      30463.487      30463.9801      30464.4894      30465.0046      30465.5516      30466.0019      30466.4944       30466.981      30467.4781       30468.029      30468.4945      30469.0217      30469.5399      30470.1099      30470.5547      30471.0646      30471.5275      30472.0719      30472.5193      30473.0173      30473.5564      30474.0345      30474.5341      30475.0233      30475.5361      30476.0786      30476.5884      30477.0522      30477.5979      30478.0806];
                        %                     obj.qUpixel(7,:) =    [30462.96        30462.9521      30463.4858      30463.9789      30464.4882      30465.0034      30465.5504      30466.0007      30466.4932      30466.9798      30467.4769      30468.0278      30468.4933      30469.0205      30469.5387      30470.1087      30470.5535      30471.0634      30471.5263      30472.0707      30472.5181      30473.0161      30473.5552      30474.0333      30474.5329      30475.0221       30475.535      30476.0774      30476.5873      30477.0512      30477.5968      30478.0794];
                        %                     obj.qUpixel(8,:) =    [30462.9525      30462.9446      30463.4783      30463.9714      30464.4807      30464.9959      30465.5429      30465.9932      30466.4857      30466.9721      30467.4692      30468.0203      30468.4858       30469.013      30469.5312      30470.1012       30470.546      30471.0559      30471.5186       30472.063      30472.5106      30473.0086      30473.5477      30474.0256      30474.5254      30475.0146      30475.5275      30476.0699      30476.5798      30477.0435      30477.5893      30478.0719];
                        %                     obj.qUpixel(9,:) =    [30462.9428      30462.9349      30463.4687      30463.9618      30464.4711      30464.9863      30465.5332      30465.9835      30466.4761      30466.9627      30467.4597      30468.0106      30468.4761      30469.0034      30469.5215      30470.0916      30470.5364      30471.0462      30471.5092      30472.0536      30472.5009       30472.999       30473.538      30474.0162      30474.5158       30475.005      30475.5178      30476.0602      30476.5701       30477.034      30477.5796      30478.0623];
                        %                     obj.qUpixel(10,:) =   [30462.9333      30462.9255      30463.4594      30463.9525      30464.4618      30464.9768       30465.524      30465.9743      30466.4668      30466.9532      30467.4503      30468.0014      30468.4669      30468.9941      30469.5123      30470.0823      30470.5271       30471.037      30471.4997      30472.0441      30472.4917      30472.9897      30473.5288      30474.0067      30474.5065      30474.9957      30475.5085      30476.0509      30476.5608      30477.0245      30477.5703       30478.053];
                        %                     obj.qUpixel(11,:) =   [30462.9282      30462.9203      30463.4541      30463.9472      30464.4565      30464.9717      30465.5186      30465.9689      30466.4615      30466.9481      30467.4451       30467.996      30468.4615      30468.9888      30469.5069       30470.077      30470.5218      30471.0316      30471.4946       30472.039      30472.4863      30472.9844      30473.5234      30474.0016      30474.5012      30474.9904      30475.5032      30476.0456      30476.5555      30477.0192       30477.565      30478.0477];
                        %                     obj.qUpixel(12,:) =   [30462.9264      30462.9185      30463.4525      30463.9456      30464.4547      30464.9699      30465.5169      30465.9674      30466.4597      30466.9463      30467.4434      30467.9943      30468.4598      30468.9872      30469.5052      30470.0752        30470.52      30471.0299      30471.4928      30472.0372      30472.4846      30472.9826      30473.5219      30473.9998      30474.4994      30474.9886      30475.5014      30476.0438      30476.5537      30477.0176      30477.5634      30478.0459];
                        %                     obj.qUpixel(13,:) =   [30462.9274      30462.9195      30463.4535      30463.9466      30464.4559      30464.9709      30465.5178      30465.9684      30466.4607      30466.9473      30467.4443      30467.9952      30468.4607      30468.9882      30469.5061      30470.0762       30470.521      30471.0309      30471.4938      30472.0382      30472.4855      30472.9836      30473.5229      30474.0008      30474.5004      30474.9896      30475.5024      30476.0448      30476.5547      30477.0186      30477.5644      30478.0469];
                        %                     obj.qUpixel(14,:) =   [30462.9304      30462.9225      30463.4564      30463.9495      30464.4586      30464.9738      30465.5208      30465.9713      30466.4636      30466.9502      30467.4473      30467.9982      30468.4637      30468.9911      30469.5091      30470.0791      30470.5239      30471.0338      30471.4967      30472.0411      30472.4885      30472.9865      30473.5256      30474.0037      30474.5034      30474.9925      30475.5054      30476.0478      30476.5577      30477.0216      30477.5672      30478.0498];
                        %                     obj.qUpixel(15,:) =   [30462.9373      30462.9294      30463.4631      30463.9563      30464.4655      30464.9807      30465.5277       30465.978      30466.4705      30466.9571      30467.4542      30468.0051      30468.4706      30468.9978       30469.516       30470.086      30470.5308      30471.0407      30471.5037       30472.048      30472.4954      30472.9934      30473.5325      30474.0106      30474.5103      30474.9994      30475.5123      30476.0547      30476.5646      30477.0285      30477.5741      30478.0567];
                        %                     obj.qUpixel(16,:) =   [30462.946       30462.9381      30463.4718      30463.9649      30464.4742      30464.9894      30465.5364      30465.9867      30466.4792      30466.9658      30467.4629      30468.0138      30468.4793      30469.0065      30469.5247      30470.0947      30470.5395      30471.0494      30471.5123      30472.0567      30472.5041      30473.0021      30473.5412      30474.0193      30474.5189      30475.0081      30475.5209      30476.0634      30476.5732      30477.0372      30477.5828      30478.0654];
                        %                     obj.qUpixel(17,:) =   [30463.027       30463.0191      30463.5531      30464.0462      30464.5553      30465.0705      30465.6175       30466.068      30466.5603      30467.0469       30467.544      30468.0949      30468.5604      30469.0876      30469.6058      30470.1758      30470.6206      30471.1305      30471.5934      30472.1378      30472.5851      30473.0832      30473.6223      30474.1004         30474.6      30475.0892       30475.602      30476.1444      30476.6543      30477.1182      30477.6638      30478.1465];
                        %                     obj.qUpixel(18,:) =   [30463.0314      30463.0235      30463.5574      30464.0505      30464.5598      30465.0748       30465.622      30466.0723      30466.5648      30467.0512      30467.5483      30468.0994      30468.5649      30469.0921      30469.6103      30470.1803      30470.6251       30471.135      30471.5977      30472.1421      30472.5897      30473.0877      30473.6268      30474.1047      30474.6045      30475.0937      30475.6066       30476.149      30476.6589      30477.1226      30477.6684       30478.151];
                        %                     obj.qUpixel(19,:) =   [30463.0241      30463.0162      30463.5499       30464.043      30464.5523      30465.0673      30465.6145      30466.0648      30466.5573      30467.0437      30467.5408      30468.0919      30468.5574      30469.0846      30469.6028      30470.1728      30470.6176      30471.1275      30471.5902      30472.1346      30472.5822      30473.0802      30473.6193      30474.0972       30474.597      30475.0862      30475.5991      30476.1415      30476.6514      30477.1151      30477.6609      30478.1435];
                        %                     obj.qUpixel(20,:) =   [30463.0099      30463.002      30463.5359       30464.029      30464.5383      30465.0533      30465.6005      30466.0508      30466.5433      30467.0297      30467.5268      30468.0779      30468.5434      30469.0706      30469.5888      30470.1588      30470.6036      30471.1135      30471.5762      30472.1206      30472.5682      30473.0662      30473.6053      30474.0832       30474.583      30475.0722      30475.5851      30476.1275      30476.6374      30477.1011      30477.6469      30478.1295];
                        %                     obj.qUpixel(21,:) =   [30462.9953      30462.9874      30463.5211      30464.0142      30464.5235      30465.0387      30465.5857       30466.036      30466.5285      30467.0151      30467.5122      30468.0631      30468.5286      30469.0558       30469.574       30470.144      30470.5888      30471.0987      30471.5616       30472.106      30472.5534      30473.0514      30473.5905      30474.0686      30474.5683      30475.0574      30475.5703      30476.1127      30476.6226      30477.0863      30477.6321      30478.1147];
                        %                     obj.qUpixel(22,:) =   [30462.9854      30462.9775      30463.5115      30464.0046      30464.5137      30465.0289      30465.5758      30466.0263      30466.5187      30467.0053      30467.5023      30468.0532      30468.5187      30469.0462      30469.5641      30470.1342       30470.579      30471.0888      30471.5518      30472.0962      30472.5435      30473.0416      30473.5806      30474.0588      30474.5584      30475.0476      30475.5604      30476.1028      30476.6127      30477.0766      30477.6222      30478.1049];
                        %                     obj.qUpixel(23,:) =   [30462.9827      30462.9748      30463.5087      30464.0018      30464.5109      30465.0261      30465.5731      30466.0236      30466.5159      30467.0025      30467.4996      30468.0505       30468.516      30469.0434      30469.5614      30470.1314      30470.5762      30471.0861       30471.549      30472.0934      30472.5408      30473.0388      30473.5781       30474.056      30474.5556      30475.0448      30475.5576      30476.1001      30476.6099      30477.0739      30477.6196      30478.1021];
                        %                     obj.qUpixel(24,:) =   [30462.9838      30462.9759      30463.5099       30464.003      30464.5121      30465.0273      30465.5743      30466.0248      30466.5171      30467.0037      30467.5008      30468.0517      30468.5172      30469.0446      30469.5626      30470.1326      30470.5774      30471.0873      30471.5502      30472.0946      30472.5419        30473.04      30473.5793      30474.0572      30474.5568       30475.046      30475.5588      30476.1012      30476.6111       30477.075      30477.6206      30478.1033];
                        %                     obj.qUpixel(25,:) =   [30462.9852      30462.9773      30463.5113      30464.0044      30464.5137      30465.0287      30465.5758      30466.0261      30466.5187      30467.0051      30467.5021      30468.0532      30468.5187       30469.046      30469.5641      30470.1342       30470.579      30471.0888      30471.5516       30472.096      30472.5435      30473.0416      30473.5806      30474.0586      30474.5584      30475.0476      30475.5604      30476.1028      30476.6127      30477.0764      30477.6222      30478.1049];
                        %                     obj.qUpixel(26,:) =   [30462.9909      30462.983      30463.5168      30464.0099      30464.5192      30465.0344      30465.5814      30466.0317      30466.5242      30467.0108      30467.5079      30468.0588      30468.5243      30469.0515      30469.5697      30470.1397      30470.5845      30471.0944      30471.5573      30472.1017      30472.5491      30473.0471      30473.5862      30474.0643      30474.5639      30475.0531      30475.5659      30476.1083      30476.6182      30477.0821      30477.6277      30478.1104];
                        %                     obj.qUpixel(27,:) =   [30463.0006      30462.9927      30463.5267      30464.0198      30464.5291      30465.0441      30465.5912      30466.0415      30466.5341      30467.0205      30467.5175      30468.0686      30468.5341      30469.0614      30469.5795      30470.1496      30470.5943      30471.1042       30471.567      30472.1114      30472.5589       30473.057       30473.596       30474.074      30474.5738      30475.0629      30475.5758      30476.1182      30476.6281      30477.0918      30477.6376      30478.1202];
                        %                     obj.qUpixel(28,:) =   [30463.0138      30463.0059      30463.5397      30464.0328      30464.5421      30465.0573      30465.6042      30466.0546      30466.5471      30467.0337      30467.5307      30468.0816      30468.5471      30469.0744      30469.5925      30470.1626      30470.6074      30471.1172      30471.5802      30472.1246      30472.5719        30473.07       30473.609      30474.0872      30474.5868       30475.076      30475.5888      30476.1312      30476.6411       30477.105      30477.6506      30478.1333];
                        %                     obj.qUpixel(29,:) =   [30463.0813      30463.0734      30463.6071      30464.1002      30464.6095      30465.1247      30465.6717       30466.122      30466.6145      30467.1011      30467.5982      30468.1491      30468.6146      30469.1418        30469.66        30470.23      30470.6748      30471.1847      30471.6476       30472.192      30472.6394      30473.1374      30473.6765      30474.1546      30474.6542      30475.1434      30475.6563      30476.1987      30476.7086      30477.1725      30477.7181      30478.2007];
                        %                     obj.qUpixel(30,:) =   [30463.0955      30463.0876      30463.6215      30464.1146      30464.6239      30465.1389      30465.6861      30466.1364      30466.6289      30467.1153      30467.6124      30468.1635       30468.629      30469.1562      30469.6744      30470.2444      30470.6892      30471.1991      30471.6618      30472.2062      30472.6538      30473.1518      30473.6909      30474.1688      30474.6686      30475.1578      30475.6707      30476.2131       30476.723      30477.1867      30477.7325      30478.2151];
                        %                     obj.qUpixel(31,:) =   [30463.0927      30463.0848      30463.6186      30464.1117       30464.621       30465.136      30465.6831      30466.1334       30466.626      30467.1124      30467.6094      30468.1605       30468.626      30469.1533      30469.6714      30470.2415      30470.6863      30471.1961      30471.6589      30472.2033      30472.6508      30473.1489      30473.6879      30474.1659      30474.6657      30475.1549      30475.6677      30476.2101        30476.72      30477.1837      30477.7295      30478.2122];
                        %                     obj.qUpixel(32,:) =   [30463.0771      30463.0692      30463.6032      30464.0963      30464.6054      30465.1206      30465.6676      30466.1181      30466.6104       30467.097      30467.5941       30468.145      30468.6105      30469.1379      30469.6559      30470.2259      30470.6707      30471.1806      30471.6435      30472.1879      30472.6352      30473.1333      30473.6726      30474.1505      30474.6501      30475.1393      30475.6521      30476.1945      30476.7044      30477.1683      30477.7139      30478.1966];
                        %                     obj.qUpixel(33,:) =   [30463.0596      30463.0517      30463.5856      30464.0787      30464.5878       30465.103        30465.65      30466.1005      30466.5928      30467.0794      30467.5765      30468.1274      30468.5929      30469.1203      30469.6383      30470.2083      30470.6531       30471.163      30471.6259      30472.1703      30472.6177      30473.1157       30473.655      30474.1329      30474.6326      30475.1217      30475.6346       30476.177      30476.6869      30477.1508      30477.6966       30478.179];
                        %                     obj.qUpixel(34,:) =   [30463.0436      30463.0357      30463.5697      30464.0628      30464.5721      30465.0871      30465.6342      30466.0845      30466.5771      30467.0635      30467.5605      30468.1116      30468.5771      30469.1044      30469.6225      30470.1926      30470.6373      30471.1472        30471.61      30472.1544      30472.6019         30473.1       30473.639       30474.117      30474.6168      30475.1059      30475.6188      30476.1612      30476.6711      30477.1348      30477.6806      30478.1632];
                        %                     obj.qUpixel(35,:) =   [30463.0385      30463.0306      30463.5645      30464.0576      30464.5669      30465.0819      30465.6291      30466.0794      30466.5719      30467.0583      30467.5554      30468.1065       30468.572      30469.0992      30469.6174      30470.1874      30470.6322      30471.1421      30471.6048      30472.1492      30472.5968      30473.0948      30473.6339      30474.1118      30474.6116      30475.1008      30475.6137      30476.1561       30476.666      30477.1297      30477.6755      30478.1581];
                        %                     obj.qUpixel(36,:) =   [30463.0395      30463.0316      30463.5653      30464.0584      30464.5677      30465.0829      30465.6299      30466.0802      30466.5727      30467.0593      30467.5564      30468.1073      30468.5728         30469.1      30469.6182      30470.1882       30470.633      30471.1429      30471.6058      30472.1502      30472.5976      30473.0956      30473.6347      30474.1128      30474.6124      30475.1016      30475.6144      30476.1569      30476.6667      30477.1307      30477.6762      30478.1589];
                        %                     obj.qUpixel(37,:) =   [30463.0428      30463.0349      30463.5687      30464.0618      30464.5711      30465.0863      30465.6332      30466.0835      30466.5761      30467.0627      30467.5597      30468.1106      30468.5761      30469.1034      30469.6215      30470.1916      30470.6364      30471.1462      30471.6092      30472.1536      30472.6009       30473.099       30473.638      30474.1162      30474.6158       30475.105      30475.6178      30476.1602      30476.6701       30477.134      30477.6796      30478.1623];
                        %                     obj.qUpixel(38,:) =   [30463.0434      30463.0355      30463.5695      30464.0626      30464.5719      30465.0869       30465.634      30466.0843      30466.5769      30467.0633      30467.5603      30468.1114      30468.5769      30469.1042      30469.6223      30470.1924      30470.6372       30471.147      30471.6098      30472.1542      30472.6017      30473.0998      30473.6388      30474.1168      30474.6166      30475.1057      30475.6186       30476.161      30476.6709      30477.1346      30477.6804       30478.163];
                        %                     obj.qUpixel(39,:) =   [30463.0507      30463.0428      30463.5768      30464.0699      30464.5792      30465.0942      30465.6413      30466.0916      30466.5842      30467.0706      30467.5676      30468.1187      30468.5842      30469.1115      30469.6296      30470.1997      30470.6444      30471.1543      30471.6171      30472.1615       30472.609      30473.1071      30473.6461      30474.1241      30474.6239       30475.113      30475.6259      30476.1683      30476.6782      30477.1419      30477.6877      30478.1703];
                        %                     obj.qUpixel(40,:) =   [30463.0649      30463.057       30463.591      30464.0841      30464.5932      30465.1084      30465.6553      30466.1058      30466.5982      30467.0848      30467.5818      30468.1327      30468.5982      30469.1257      30469.6436      30470.2137      30470.6585      30471.1683      30471.6313      30472.1757       30472.623      30473.1211      30473.6603      30474.1383      30474.6379       30475.127      30475.6399      30476.1823      30476.6922      30477.1561      30477.7019      30478.1843];
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
            
            p.parse(varargin{:});
            obj.E_Bias                 = p.Results.E_Bias;
            obj.W_Bias                 = p.Results.W_Bias;
            obj.B_Bias                 = p.Results.B_Bias;
            obj.Phi0_Bias              = p.Results.Phi0_Bias;
            
            % Initialization
            SetFitBias(obj,1);
            
            % Kinematics
            obj.KrDS = zeros(obj.nTe,1); obj.KrDSE = obj.KrDS;
      
            % Lorentzian Line Definition
            obj.KrLine = @(Phi0,GammaH,Eh,Ec) Phi0*(GammaH/(2*pi))./((Ec-Eh).^2+GammaH.^2/4);
            
            switch obj.TD
                case 'KrK32'
                   E_temp = obj.K32_E; W_temp = obj.K32_W;
                   Phi0_temp = obj.K32_Phi0; %B_temp = obj.K32_B;
                case {'HVL332','KrL3_32','KrL3_32_HS','MoSFitterJuly2017','DummyKrdc1'}
                   E_temp = obj.L3_32_E; W_temp = obj.L3_32_W;
                   Phi0_temp = obj.L3_32_Phi0; %B_temp = obj.L3_32_B;
            end
            
            switch obj.DopplerEffectFlag
                case 'OFF'
                    obj.KrDS   = (obj.KrLine(Phi0_temp,W_temp,E_temp,obj.Te));
                case {'matConv','numConv'}
                    obj.enlargeSpectrum();
                    if strcmp(obj.recomputeKernel,'ON') || strcmp(obj.DopplerEffectFlag,'numConv')
                        obj.computeKernel(E_temp); end
                    % Enlarge and restore spectrum to avoid convolution edge effects
                    obj.KrDS = obj.KrLine(Phi0_temp,W_temp,E_temp,obj.Te);
                    obj.KrDS = obj.KernelSpectrumConv(obj.KrDS,E_temp);
                    obj.KrDS = obj.restoreSpectrum(obj.KrDS);
                case 'Voigt'
                    if strcmp(obj.recomputeKernel,'ON'); obj.computeKernel(E_temp); end
                    z = (obj.Te-E_temp + 1i.*W_temp/2)/(obj.DE_sigma*1.4142135623731); %2^0.5
                    obj.KrDS = Phi0_temp*real(fadf(z))/(obj.DE_sigma*2.506628274631); %(pi)^0.5
            end
            
            obj.KrDSE  = obj.KrDS.^0.5;
            
        end
        
        function ComputeKrDSallPixels(obj,varargin)
            
            % Compute Differential Spectrum of a Krypton Line
            % Simultaneously for all pixels
            
            % Inputs
            p = inputParser;
            p.addParameter('E_Bias',0,@(x)isfloat(x));
            p.addParameter('W_Bias',0,@(x)isfloat(x));
            p.addParameter('Phi0_Bias',zeros(obj.nPixels,1)',@(x)isfloat(x)); 
            p.addParameter('B_Bias',zeros(obj.nPixels,1)',@(x)isfloat(x));       
            
            p.parse(varargin{:});
            obj.E_Bias                 = p.Results.E_Bias;
            obj.W_Bias                 = p.Results.W_Bias;
            obj.Phi0_Bias              = p.Results.Phi0_Bias; 
            obj.B_Bias                 = p.Results.B_Bias;

            % Initialization
            SetFitBiasallPixels(obj,1);
            %SetLineInitParallPixels(obj); %Sets Phi0 and/or BKG to single pixel fit value
                    
            % Line Definition
            obj.KrLine = @(Phi0,GammaH,Eh,Ec) Phi0*(GammaH/(2*pi))./((Ec-Eh).^2+GammaH.^2/4);
            
            % Loop on All Pixels
            switch obj.TD
                case 'KrK32'
                    E_temp = obj.K32_E; W_temp = obj.K32_W;
                    Phi0_temp = obj.K32_Phi0allPixels(1:obj.nPixels); %B_temp = obj.K32_B;
                case {'KrL3_32','KrL3_32_HS','MoSFitterJuly2017','DummyKrdc1','HVL332'}
                    E_temp = obj.L3_32_E; W_temp = obj.L3_32_W;
                    Phi0_temp = obj.L3_32_Phi0allPixels(1:obj.nPixels); %B_temp = obj.L3_32_B;
            end

            switch obj.DopplerEffectFlag
                case 'OFF'
                    obj.KrDSallPixels = arrayfun(@(x,y) obj.KrLine(x,W_temp,E_temp,y),repmat(Phi0_temp,obj.nTe,1),obj.Tepixel,'UniformOutput',false);
                case {'matConv','numConv'}
                    if strcmp(obj.recomputeKernel,'ON'); obj.computeKernel(); end
                    % Enlarge and restore spectrum to avoid convolution edge effects
                    obj.enlargeSpectrum();
                    obj.KrDSallPixels = arrayfun(@(x,y) obj.KrLine(x,W_temp,E_temp,y),repmat(Phi0_temp,obj.nTe,1),obj.Tepixel,'UniformOutput',false);
                    obj.KrDSallPixels = cellfun(@(x) obj.KernelSpectrumConv(x),obj.KrDSallPixels,'UniformOutput',false);
                    obj.KrDSallPixels = cellfun(@(x) obj.restoreSpectrum(x),obj.KrDSallPixels,'UniformOutput',false);
                    obj.restoreProperties();
                case 'Voigt'
                    if strcmp(obj.recomputeKernel,'ON'); obj.computeKernel(); end
                    z = (obj.Tepixel-E_temp + 1i.*W_temp/2)/(obj.DE_sigma*1.4142135623731); %2^0.5
                    %obj.KrDSallPixels = arrayfun(@(x) x*real(fadf(z))/(obj.DE_sigma*2.506628274631),repmat(Phi0_temp,obj.nTe,1),'UniformOutput',false); %(pi)^0.5
                    obj.KrDSallPixels = real(fadf(z))/(obj.DE_sigma*2.506628274631) .* repmat(Phi0_temp,obj.nTe,1);
            end
            
            %obj.KrDSallPixels = cell2mat(obj.KrDSallPixels); %Convert cell to matrix
            obj.KrDSEallPixels = obj.KrDSallPixels.^.5;
            
            %Truncate KrDS    
%             obj.KrDSallPixels = obj.KrDSallPixels(:,1:obj.nPixels);
%             obj.KrDSEallPixels = obj.KrDSEallPixels(:,1:obj.nPixels);
            end
        
        function ComputeKrIS(obj,varargin)
            p = inputParser;
            p.addParameter('IStype','TRAPEZ',@(x)ismember(x,{'TRAPEZ','ACCURATE'}));
            p.parse(varargin{:});
            IStype           = p.Results.IStype;
            
            obj.KrIS = zeros(obj.nqU,1); obj.KrISE = obj.KrIS;
            if strcmp(obj.HVRipples,'OFF')
                    switch IStype
                        case 'ACCURATE' % Real integral --> slow but accurate         
                            switch obj.TD
                                case 'KrK32'
                                    obj.KrLineTail = ComputeLorentzianTail(obj,obj.K32_Phi0,obj.K32_E,obj.K32_W,100e4,obj.KTF);
                                    
                                case {'HVL332','KrL3_32','KrL3_32_HS','MoSFitterJuly2017','DummyKrdc1'}
                                    obj.KrLineTail = ComputeLorentzianTail(obj,obj.L3_32_Phi0,obj.L3_32_E,obj.L3_32_W,100e4,obj.KTF);
                            end
            
                            for qUi = 1:1:obj.nqU
                                KrDSf = @(e) interp1(obj.Te,obj.KrDS.*obj.KTF(obj.Te,obj.qU(qUi)),e);
                                obj.KrIS(qUi) = ...
                                    (integral(KrDSf,obj.qU(qUi),obj.qUmax)+obj.BKG_RateSec+obj.KrLineTail(qUi))...
                                    *obj.TimeSec.*obj.qUfrac(qUi) ;
                                obj.KrISE(qUi) = obj.KrIS(qUi).^.5;
                            end
                        case 'TRAPEZ'  % Trapezoidal integral --> faster, less accurate
                            
                            switch obj.TD
                                case 'KrK32'
                                    obj.KrLineTail = ComputeLorentzianTail(obj,obj.K32_Phi0,obj.K32_E,obj.K32_W,1000e4,obj.KTF);
                                    
                                case {'HVL332','KrL3_32','KrL3_32_HS','MoSFitterJuly2017','DummyKrdc1'}
                                    obj.KrLineTail = ComputeLorentzianTail(obj,obj.L3_32_Phi0,obj.L3_32_E,obj.L3_32_W,1000e4,obj.KTF);
                            end
                            
%                             for qUi = 1:1:obj.nqU
                            %    qUieq = obj.qUTeMatch(qUi);
                                                             
                             %   if qUieq < obj.nTe
                            % obj.KrIS(qUi) = ...
                            %     (trapz(obj.Te(qUieq:end),obj.KrDS(qUieq:end).*obj.KTF(obj.Te(qUieq:end),obj.qU(qUi),obj.MACE_R_eV))+obj.BKG_RateSec+obj.KrLineTail(qUi)).*...
                            %     obj.TimeSec.*obj.qUfrac(qUi);
                             %obj.KrISE(qUi) = obj.KrIS(qUi).^.5;
%                              obj.KrIS(qUi) = ...
%                                  (trapz(obj.Te(qUi:end),obj.KrDS(qUi:end).*obj.KTF(obj.Te(qUi:end),obj.qU(qUi),obj.MACE_R_eV))+...
%                                  obj.BKG_RateSec+obj.KrLineTail(qUi)).*...
%                                  obj.TimeSec.*obj.qUfrac(qUi);
%                              obj.KrISE(qUi) = obj.KrIS(qUi).^.5;
%                              
                             obj.KrIS = ...
                                 (simpsons(obj.Te,obj.KrDS.*obj.TMPTF')+...
                                 obj.BKG_RateSec+obj.KrLineTail).*...
                                 obj.TimeSec.*obj.qUfrac;
                             obj.KrISE = obj.KrIS.^.5;                         
                             %  else
                                %    obj.KrIS(qUi)  = (obj.BKG_RateSec+obj.KrLineTail(qUi))*obj.TimeSec.*obj.qUfrac(qUi);
                                %    obj.KrISE(qUi) = obj.KrIS(qUi).^.5;
                               % end
                            %end
                    end
                    
            else % Trapez Integration Only
                
                switch obj.TD
                    case 'KrK32'
                        obj.KrLineTail = ComputeLorentzianTail(obj,obj.K32_Phi0,obj.K32_E,obj.K32_W,1000e4,obj.KTF);
                    case {'HVL332','KrL3_32','KrL3_32_HS','MoSFitterJuly2017','DummyKrdc1'}
                        obj.KrLineTail = ComputeLorentzianTail(obj,obj.L3_32_Phi0,obj.L3_32_E,obj.L3_32_W,1000e4,obj.KTF);
                end
                
                x = -pi/2:pi/10:pi/2;
                for qUi = 1:1:obj.nqU
                    %qUieq = obj.qUTeMatch(qUi);
                    %if qUieq <= obj.nTe
                        tmp = (simpsons(obj.Te,(obj.KrDS).*obj.KTF(obj.Te,obj.qU(qUi)+obj.HVRipplesP2PV/2.*sin(x)))+...
                            obj.BKG_RateSec+obj.KrLineTail(qUi)).*...
                            obj.TimeSec.*obj.qUfrac(qUi);
                        obj.KrIS(qUi)  = mean(tmp);% + obj.KrLineTail*obj.TimeSec.*obj.qUfrac(qUi);
                        obj.KrISE(qUi) = obj.KrIS(qUi).^.5;
                   % else
                      %  obj.KrIS(qUi)  = (obj.BKG_RateSec+obj.KrLineTail(qUi))*obj.TimeSec.*obj.qUfrac(qUi);
                      %  obj.KrISE(qUi) = obj.KrIS(qUi).^.5;
                    %end
                end
            end
            
            switch obj.CPS
                case 'ON'
                    obj.KrISE  = obj.KrISE  ./ (obj.TimeSec*obj.qUfrac ); % BEFORE KrIS!!!
                    obj.KrIS   = obj.KrIS   ./ (obj.TimeSec*obj.qUfrac );
            end
            
        end
        
        function f      = ComputeKrISf(obj,qu,eb,wb,pb,bb,nb)
            % Function for Matlab Fit
            ComputeKrDS(obj,...
                'E_bias',eb,...
                'W_bias',wb,...
                'Phi0_bias',pb,...
                'B_bias',bb);
            
            ComputeKrIS(obj,'IStype','TRAPEZ');
            
            qUieq = []; qUieq = (find(abs(obj.qU-qu)<1e-3));
            if qUieq>0
                f = (1+nb).*obj.KrIS(qUieq);
            else
                f = -1;
                return;
            end
        end
        
        function   ComputeKrISallPixels(obj,varargin)
            %
            % Compute Kr Integral Spectra
            % (simultaneously for all pixels)
            % Vectorization along qU and Pixels
            %
            % Th. Lasserre
            % Last Updated: March 8 2018
            %                         obj.qU = obj.qUpixel(i,:);
            
            switch obj.TD
                case 'KrK32'
                    obj.KrLineTailallPixels = ...
                        arrayfun(@(p) obj.ComputeLorentzianTailallPixels(...
                        obj.K32_Phi0allPixels(p),obj.K32_E,obj.K32_W,1e7,p),...
                        1:1:obj.nPixels,'UniformOutput',false);
                case {'HVL332','KrL3_32','KrL3_32_HS','MoSFitterJuly2017','DummyKrdc1'}
                    obj.KrLineTailallPixels = ...
                        arrayfun(@(p) obj.ComputeLorentzianTailallPixels(...
                        obj.L3_32_Phi0allPixels(p),obj.L3_32_E,obj.L3_32_W,1e7,p),...
                        1:1:obj.nPixels,'UniformOutput',false);
            end
            obj.KrLineTailallPixels = reshape(cell2mat(obj.KrLineTailallPixels),obj.nqU,obj.nPixels);
            
            % Old method
            % .*obj.TeAllPixelsInc
%                         obj.KrISallPixels = ((...
%                             squeeze(trapz(...
%                             ones(numel(obj.Tepixel(:,1)),1),permute((repmat(obj.TeAllPixelsInc,1,1,obj.nqU)),[1 3 2]).*reshape(repmat(obj.KrDSallPixels,obj.nqU,1,1),[obj.nTe,obj.nqU,obj.nPixels])...
%                             .*permute(obj.TMPTFallPixels,[3 2 1])))...
%                             +repmat(obj.BKG_RateSecallPixels,obj.nqU,1)+obj.KrLineTailallPixels).*obj.TimeSec.*obj.qUfrac);

%

% TRAPZ VECTORIZED - WORKS!
Y=(reshape(repmat(obj.KrDSallPixels,obj.nqU,1,1),[obj.nTe,obj.nqU,obj.nPixels]).*permute(obj.TMPTFallPixels,[3 2 1]));
obj.KrISallPixels = (squeeze(sum(abs(permute(repmat(diff(obj.Tepixel,1,1),1,1,obj.nqU),[1 3 2])).*0.5.*( Y(1:end-1,:,:) + Y(2:end,:,:))))+repmat(obj.BKG_RateSecallPixels,obj.nqU,1)+obj.KrLineTailallPixels).*obj.TimeSec.*obj.qUfrac;

            % CA MARCHE
%             for p = 1:1:obj.nPixels
%                 telocal = obj.Tepixel(:,p);
%                 tmp = ((...
%                     squeeze(trapz(telocal,reshape(repmat(obj.KrDSallPixels,obj.nqU,1,1),[obj.nTe,obj.nqU,obj.nPixels])...
%                     .*permute(obj.TMPTFallPixels,[3 2 1]),1))...
%                     +repmat(obj.BKG_RateSecallPixels,obj.nqU,1)+obj.KrLineTailallPixels).*obj.TimeSec.*obj.qUfrac);
%                 obj.KrISallPixels(:,p)=tmp(:,p);
%             end
            
        end
       
%        function ComputeKrISallPixels(obj,varargin)
%             
%             %Compute Krypton Integral Spectra
%             %Simultaneously for all pixels
%                         
%             %            Inputs
%             p = inputParser;
%             p.addParameter('IStype','TRAPEZRipples',@(x)ismember(x,{'TRAPEZ','ACCURATE','TRAPEZRipples'}));
%             p.addParameter('BKG_RateSecallPixels',1*ones(1,obj.nPixels));
%             
%             p.parse(varargin{:});
%             IStype  = p.Results.IStype;
%             %obj.BKG_RateSecallPixels   = p.Results.BKG_RateSecallPixels;
%             %SetLineInitParallPixels(obj);
%             
%             switch IStype
%                 case 'TRAPEZ'  % Trapezoidal integral --> faster, less accurate
%                     for i=1:1:obj.nPixels
%                         localKTF  = obj.KTFallPixels{i}; %function handle
%                         localKrDS = obj.KrDSallPixels{i};
%                         
%                         % Careful patch - read pixel wise qU
%                         obj.qU = obj.qUpixel(i);
%                         
%                         % Compute Lorentzian Tail Offset
%                         switch obj.TD
%                             case 'KrK32'
%                                 obj.KrLineTail = ComputeLorentzianTail(obj,obj.K32_Phi0allPixels(i),obj.K32_E,obj.K32_W,10000e4,localKTF);
%                                 
%                             case {'HVL332','KrL3_32','KrL3_32_HS','MoSFitterJuly2017','DummyKrdc1'}
%                                 obj.KrLineTail = ComputeLorentzianTail(obj,obj.L3_32_Phi0allPixels(i),obj.L3_32_E,obj.L3_32_W,10000e4,localKTF);
%                         end
%                         
%                         for qUi = 1:1:obj.nqU
%                             qUieq = obj.qUTeMatch(qUi);
%                             
%                             if qUieq < obj.nTe
%                                 obj.KrISallPixels{qUi,i} = ...
%                                     (trapz(obj.Te(qUieq:end),localKrDS(qUieq:end).*localKTF(obj.Te(qUieq:end),obj.qU(qUi),obj.MACE_R_eV))+obj.BKG_RateSecallPixels(i)+obj.KrLineTail(qUi)).*...
%                                     obj.TimeSec.*obj.qUfrac(qUi);
%                                 obj.KrISEallPixels{qUi,i} = obj.KrISallPixels{qUi,i}.^.5;
%                             else
%                                 obj.KrISallPixels{qUi,i}  = (obj.BKG_RateSecallPixels(i)+obj.KrLineTail(qUi))*obj.TimeSec.*obj.qUfrac(qUi);
%                                 obj.KrISEallPixels{qUi,i} = obj.KrISallPixels{qUi,i}.^.5;
%                             end
%                             
%                             switch obj.CPS
%                                 case 'ON'
%                                     obj.KrISEallPixels{qUi,i}  = cell2mat(obj.KrISEallPixels(qUi,i))  ./ ((obj.TimeSec*obj.qUfrac(qUi))); %BEFORE KrIS!!!
%                                     obj.KrISallPixels{qUi,i}   = cell2mat(obj.KrISallPixels(qUi,i))   ./ (obj.TimeSec*obj.qUfrac(qUi));
%                                     fprintf('CPS - %.0f - %.0f - %.3f \n',qUi,i,(obj.KrISallPixels{qUi,i}));
%                             end
%                         end
%                         
% %                         switch obj.CPS
% %                             case 'ON'
% %                                 obj.KrISEallPixels(:,i)  = (cellfun(@(x,y) x./y,{cell2mat(obj.KrISEallPixels(:,i))},{obj.TimeSec*obj.qUfrac},'UniformOutput',false));
% %                                 obj.KrISallPixels(:,i)   = (cellfun(@(x,y) x./y,{cell2mat(obj.KrISallPixels(:,i))},{obj.TimeSec*obj.qUfrac},'UniformOutput',false));
% %                         end
%                         
%                     end % Pixel Loop
%                 case 'TRAPEZRipples'  % Trapezoidal integral + HV Ripples
%                     for i=1:1:obj.nPixels
%                         
%                         localKTF = obj.KTFallPixels{i};
%                         localKDS = obj.KrDSallPixels(:,i);
%                         % Careful patch - read pixel wise qU
%                         obj.qU = obj.qUpixel(i,:);
%                         
%                         %Compute Lorentzian Tail Offset
%                         switch obj.TD
%                             case 'KrK32'
%                                 obj.KrLineTail = ComputeLorentzianTail(obj,obj.K32_Phi0allPixels(i),obj.K32_E,obj.K32_W,10000e4,localKTF);
%                                 
%                             case {'HVL332','KrL3_32','KrL3_32_HS','MoSFitterJuly2017','DummyKrdc1'}
%                                 obj.KrLineTail = ComputeLorentzianTail(obj,obj.L3_32_Phi0allPixels(i),obj.L3_32_E,obj.L3_32_W,10000e4,localKTF);
%                         end
%                         
%                         x = -pi/2:pi/10:pi/2; %
% 
%                         for qUi = 1:1:obj.nqU
%                            % qUieq = obj.qUTeMatch(qUi);
% 
%                            % if qUieq <= obj.nTe
%                                 tmp = (trapz(obj.Te,localKDS.*localKTF(obj.Te,obj.qU(qUi)+obj.HVRipplesP2PV/2.*sin(x),obj.MACE_R_eV))+obj.BKG_RateSecallPixels(i)+obj.KrLineTail(qUi)).*...
%                                     obj.TimeSec.*obj.qUfrac(qUi);
%                                 obj.KrISallPixels(qUi,i)  = mean(tmp);
%                                 obj.KrISEallPixels(qUi,i) = obj.KrISallPixels(qUi,i).^.5;
%                             %else
%                             %    obj.KrISallPixels{qUi,i}  = (obj.BKG_RateSec+obj.KrLineTail(qUi))*obj.TimeSec.*obj.qUfrac(qUi);
%                             %    obj.KrISEallPixels{qUi,i} = obj.KrISallPixels{qUi,i}.^.5;
%                             %end
%                             switch obj.CPS
%                                 case 'ON'
%                                     obj.KrISEallPixels{qUi,i}  = cell2mat(obj.KrISEallPixels(qUi,i))  ./ ((obj.TimeSec*obj.qUfrac(qUi))); %BEFORE KrIS!!!
%                                     obj.KrISallPixels{qUi,i}   = cell2mat(obj.KrISallPixels(qUi,i))   ./ (obj.TimeSec*obj.qUfrac(qUi));
%                             end
%                         end
%                     end
%             end
%         % obj.KrISallPixels = cell2mat(obj.KrISEallPixels);
%         % obj.KrISEallPixels = cell2mat(obj.KrISEallPixels);
%         end
       
        % Normalize Spectrum
        function    ComputeNormFactorKrDS(obj)
            % Normalize Spectrum according to KATRIN-settings
            %obj.TimeSec = obj.TimeYear*obj.Year2Sec;
        end
        
        % Compute Lorentzian Offset Tail
        function   output = ComputeLorentzianTail(obj,a,e,w,qUtail,ktf)
            % Compute the integral of the Lorentzian tail, 
            % from the qUmax to a large qU, qUtail
            % 
            % Account for Doppler Effect broadening, if activated
            % (Integration using the voigt profile approximation)
            % 
            % Input:
            %  - a      : line amplitude, cps
            %  - e      : line position, eV
            %  - w      : line width; eV
            %  - qUtail : upper bound of integral, V
            % Output:
            %  - integral line qUend --> qUtail
            %  Set obj.KrLineTail Vector
            %
            % Th. Lasserre, Feb. 5 2018
            %
            switch obj.DopplerEffectFlag
                case 'OFF'
                    tmpf = @(x) a*(w/(2*pi))./((e-x).^2+w.^2/4);
                otherwise
                    z = @(x) (x-e + 1i.*w/2)/(obj.DE_sigma*1.4142135623731); %2^0.5
                    tmpf = @(x) a*real(fadf(z(x)))/(obj.DE_sigma*2.506628274631); %(pi)^0.5 
            end
            
%             %output = zeros(1,obj.nqU);
%             tmptff = @(x) tmpf(x).*ktf(x,obj.qU(1));
%             tmp = integral(tmptff,obj.qUmax,qUtail);
%             
%             % Loop over retarding potentials - coud be vectorized!!!
%             for i=1:1:obj.nqU
%                 if i<(obj.nqU-obj.nqU/4)
%                     output(i)  = tmp;
%                 else
%                     tmptff = @(x) tmpf(x).*ktf(x,obj.qU(i));
%                     output(i) = integral(tmptff,obj.qUmax,qUtail);
%                 end
%             end
            
            gammaFacApprox = 1.018;
            ktf = @(te,qu) real(((te > qu-obj.Pixel_MACE_Ba_VCorr(obj.FPD_Pixel)) & (te < qu-obj.Pixel_MACE_Ba_VCorr(obj.FPD_Pixel) + qu*obj.MACE_Ba_T*obj.Pixel_MACE_Ba_TCorr(obj.FPD_Pixel)/obj.MACE_Bmax_T*gammaFacApprox) ).* ...
                (1 - sqrt(1 - (te-qu-obj.Pixel_MACE_Ba_VCorr(obj.FPD_Pixel))./te * obj.WGTS_B_T./obj.MACE_Ba_T*obj.Pixel_MACE_Ba_TCorr(obj.FPD_Pixel) * (1/gammaFacApprox) ))./...
                        (1 - sqrt(1- obj.WGTS_B_T/obj.MACE_Bmax_T))) +...
                        1.*(te >= qu-obj.Pixel_MACE_Ba_VCorr(obj.FPD_Pixel) + qu*obj.MACE_Ba_T*obj.Pixel_MACE_Ba_TCorr(obj.FPD_Pixel)/obj.MACE_Bmax_T*gammaFacApprox);
            tmp = integral(@(x) tmpf(x).*ktf(x,obj.qU(1)),obj.qUmax,qUtail);
            output=ones(1,obj.nqU)*tmp;
            output(floor(obj.nqU-floor(obj.nqU/6)):end)=...
                arrayfun(@(y) integral(@(x) tmpf(x).*ktf(x,y),obj.qUmax,qUtail),obj.qU(floor(obj.nqU-floor(obj.nqU/6))));
            
        end
        
        % Add Statistical Fluctuation
        function          AddStatFluctKrDS(obj)
            obj.KrDS   = obj.KrDS + obj.KrISE.*randn(obj.nTe,1);
            %obj.KrDSE  = obj.KrDS.^0.5;
        end
        
        function          AddStatFluctKrIS(obj)
            obj.KrIS   = obj.KrIS + obj.KrISE.*randn(obj.nqU,1);
            %obj.KrISE  = (obj.KrIS.^0.5);
        end
        
        function          AddStatFluctKrISallPixels(obj)
            obj.KrISallPixels   =  (cell2mat((obj.KrISallPixels)) + cell2mat(obj.KrISEallPixels).*randn(obj.nqU,obj.nPixels));
            obj.KrISallPixels   = mat2cell(obj.KrISallPixels,obj.nqU,obj.nPixels);
            %obj.KrISE  = (obj.KrIS.^0.5);
        end

        
        % Compute Lorentzian Offset Tail
        function   output = ComputeLorentzianTailallPixels(obj,a,e,w,qUtail,pixel)
            % Compute the integral of the Lorentzian tail, 
            % from the qUmax to a large qU, qUtail
            % 
            % Account for Doppler Effect broadening, if activated
            % (Integration using the voigt profile approximation)
            % 
            % Input:
            %  - a      : line amplitude, cps
            %  - e      : line position, eV
            %  - w      : line width; eV
            %  - qUtail : upper bound of integral, V
            % Output:
            %  - integral line qUend --> qUtail
            %  Set obj.KrLineTail Vector
            %
            % Th. Lasserre, Feb. 5 2018
            %
            
            % DISABLE TEMPORARILY
            switch obj.DopplerEffectFlag
                case 'OFF'
                    tmpf = @(x) a*(w/(2*pi))./((e-x).^2+w.^2/4);
                otherwise
                    z = @(x) (x-e + 1i.*w/2)./(obj.DE_sigma*1.4142135623731); %2^0.5
                    tmpf = @(x) a*real(fadf(z(x)))./(obj.DE_sigma*2.506628274631); %(pi)^0.5 
            end
            
            % TF
            gammaFacApprox = 1.018;
            ktf = @(te,qu) real(((te > qu-obj.Pixel_MACE_Ba_VCorr(pixel)) & (te < qu-obj.Pixel_MACE_Ba_VCorr(pixel) + qu*obj.MACE_Ba_T*obj.Pixel_MACE_Ba_TCorr(pixel)/obj.MACE_Bmax_T*gammaFacApprox) ).* ...
                (1 - sqrt(1 - (te-qu-obj.Pixel_MACE_Ba_VCorr(pixel))./te * obj.WGTS_B_T./obj.MACE_Ba_T*obj.Pixel_MACE_Ba_TCorr(pixel) * (1/gammaFacApprox) ))./...
                        (1 - sqrt(1- obj.WGTS_B_T/obj.MACE_Bmax_T))) +...
                        1.*(te >= qu-obj.Pixel_MACE_Ba_VCorr(pixel) + qu*obj.MACE_Ba_T*obj.Pixel_MACE_Ba_TCorr(pixel)/obj.MACE_Bmax_T*gammaFacApprox);
            tmp = integral(@(x) tmpf(x).*ktf(x,obj.qU(1)),obj.qUmax,qUtail);
            output=ones(1,obj.nqU)*tmp;
            output(floor(obj.nqU-floor(obj.nqU/6)):end)=...
                arrayfun(@(y) integral(@(x) tmpf(x).*ktf(x,y),obj.qUmax,qUtail),obj.qU(floor(obj.nqU-floor(obj.nqU/6))));
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
                        nCMTrials,obj.TD,obj.TimeSec);
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
                case {'KrL3_32','KrL3_32_HS','MoSFitterJuly2017','DummyKrdc1','DummyKrdc_tail'}
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
                case {'Conv','Voigt','realConv'}
                    fprintf('   - Doppler Effect: ON                     \n');
                    fprintf('    - Doppler Effect: T = %g K               \n',obj.WGTS_Temp);
                    %fprintf('    - Doppler Effect: Real Convolution  = %s \n',obj.realConv);
                    fprintf('    - Doppler Effect: Energy Broadening = %g eV\n',obj.DE_sigma);
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
            fprintf('   - Time (seconds / years) = %g / %g \n',obj.TimeSec,obj.TimeSec);
            fprintf(2,'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n');
            fprintf(2,'--------------------------------------------------------------\n');
            disp(' ');
            
        end
        
        function computeKernel(obj,E_cms)
            % The kernel is calculate according to [citation]
            
            % The standard deviations of the kernel funtion depending on
            % velocity (DE_sigma_e_vel) and energy (DE_sigma).
            obj.DE_sigma_e_vel = sqrt(obj.kb*obj.WGTS_Temp/obj.M);
            obj.DE_sigma = obj.DE_sigma_e_vel*sqrt(2*E_cms*obj.me/obj.c^2);
            
            % range for the kernel            
            obj.e_velParal = (obj.Te - E_cms)/((obj.me/obj.c^2)*obj.e_vel);
            obj.cosO = linspace(cos(obj.WGTS_CosMaxAAngle),1,obj.nTe);
            [e_velParal_M,cosO_M] = meshgrid(obj.e_velParal,obj.cosO);
            prefac = (1-obj.WGTS_CosMaxAAngle).^(-1);
            
            % integration of the emission angle
            g_vm = prefac./(sqrt(2*pi)*obj.DE_sigma_e_vel)*...
                simpsons(obj.cosO,exp(-0.5*((e_velParal_M-cosO_M*obj.BulkVelKr)/obj.DE_sigma_e_vel).^2));
            
            obj.GaussianKernelDE = 1/((obj.me/obj.c^2)*obj.e_vel)*g_vm;
            obj.ProbKernDE = trapz(obj.Te - E_cms,obj.GaussianKernelDE);
            g_fwhm = obj.GaussianKernelDE - 0.5*max(obj.GaussianKernelDE);
            g_fwhm_index = find(0 < g_fwhm);
            x1 = g_fwhm_index(1); xf = g_fwhm_index(end);
            x0_left = obj.Te(x1) - g_fwhm(x1)*(obj.Te(x1-1)-obj.Te(x1))/(g_fwhm(x1-1)-g_fwhm(x1));
            x0_right = obj.Te(xf) - g_fwhm(xf)*(obj.Te(xf+1)-obj.Te(xf))/(g_fwhm(xf+1)-g_fwhm(xf));
            FWHM = x0_right - x0_left;
            obj.DE_sigma = FWHM/(2*sqrt(2*log(2)));
        end % computeKernel
        
        function plotKernel(obj)
            figure(1);
            hold on
            plot(obj.Te,obj.GaussianKernelDE);
            leg_cell = {['u = ',num2str(obj.BulkVelKr),' m/s']};
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
            
        end % enlargeSpectrum
        
        function convolutedspectrum = KernelSpectrumConv(obj,spec,E_temp)
            switch obj.DopplerEffectFlag
                case 'numConv' % 
                    convolutedspectrum = obj.TeStep*conv(spec,obj.GaussianKernelDE');
                    nKernel = length(obj.GaussianKernelDE);
                    nconvspec = length(convolutedspectrum);
                    quarterconv = round(nconvspec/4);
                    middlekernel = round(nKernel/2);
                    indexkernel = (1:nKernel)';
                    indexEcms = indexkernel(obj.Te == E_temp);
                    differenceposition = abs(middlekernel - indexEcms);
                    begin = quarterconv + differenceposition + 1;
                    endspec = begin + nKernel - 1;
                    convolutedspectrum = convolutedspectrum(begin:endspec);
                case 'matConv' % convolution using matrix
                    obj.mate = obj.eres(obj.Te,obj.TeMin,obj.TeMax,obj.nTe,obj.DE_sigma);
                    %eres(E,Emin,Emax,Nbins,sigma)
                    %E = obj.Te; Emin = obj.TeMin; Emax = obj.TeMax; Nbins = obj.nTe; sigma = obj.DE_sigma;
                    convolutedspectrum = obj.mate*spec;
            end
        end % KernelSpectrumConv
        
        function restoredspectrum = restoreSpectrum(obj,spec)
            % Spectrum is restored to the original size
            restoredspectrum = spec(obj.Te>=obj.qUmin & obj.Te<=obj.qUmax);
            if obj.FPD_Pixel ~= 0 
            obj.Te = obj.Te(obj.Te>=obj.qUmin & obj.Te<=obj.qUmax);
            obj.nTe = length(obj.Te);
            end
        end
        
        function restoreProperties(obj)
            obj.Te = obj.Te(obj.Te>=obj.qUmin & obj.Te<=obj.qUmax);
            
            % properties are recalculated
            obj.nTe = length(obj.Te);

        end
        
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
    methods % some Methods for Fitting
            function y=FitCalculateKrIS(obj,varargin)
            obj.ComputeKrDS(varargin{:});
            obj.ComputeKrIS;
            y = obj.KrIS;
            end
            
            function y=FitCalculateKrISallPixels(obj,varargin)
            obj.ComputeKrDSallPixels(varargin{:});
            obj.ComputeKrISallPixels;
            y = reshape(obj.KrISallPixels,obj.nPixels*obj.nqU,1);
            end

            function onestring = BuildFitStrings(obj) 
                %For Fitting Instructions
                %Build one string with all Names for fMinuit
                parnames = {'E_Bias', 'W_Bias'};
                for k=1:1:obj.nPixels
                    tmp         = sprintf(' Phi0_Bias%g',k);
                    parnames = {parnames{:} tmp};
                end
                for k=1:1:obj.nPixels
                    tmp         = sprintf(' B_Bias%g',k);
                    parnames    = {parnames{:} tmp};
                end
                onestring = cell2mat(strcat(parnames,{' '}));
            end
            
            function [ParIni , fix , npar] = FitInit(obj,fixPar)                          
                % Init (Bias)
                i_E        = 0; i_W        = 0;
                i_Phi0     = zeros(obj.nPixels,1)';
                i_Bkg      = zeros(obj.nPixels,1)';
                ParIni = [i_E i_W  i_Phi0 i_Bkg];
                fix = sprintf('');
                switch fixPar
                    case 'none'
                        npar       = 2+obj.nPixels*2;    
                    case 'Phi0'
                        npar       = 2+obj.nPixels;
                        for i= 3:2+obj.nPixels  
                        fix_tmp = sprintf(' fix %u ;',i);
                        fix = strcat(fix, fix_tmp);          
                        end
                    case 'BKG'
                        npar       = 2+obj.nPixels*2;   
                        for i=3+obj.nPixels:3+2*obj.nPixels
                        fix_tmp = sprintf(' fix %u;',i);
                        fix = strcat(fix, fix_tmp); 
                        end
                end
            end    
            
            function tmp_z = zforint(obj,E,W) % function for fadz
            tmp_z = @(E,W) ((obj.Te-E + 1i.*W/2)/(obj.DE_sigma*1.4142135623731)); 
            end
    end
   
end % class
