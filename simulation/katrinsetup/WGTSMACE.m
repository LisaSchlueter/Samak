% ----------------------------------------------------------------------- %
% Class describing the KATRIN Tritium
% Source Apparatus (WGTS) + Magnetic Adiabatic Spectrometer (MACE)
% ----------------------------------------------------------------------- %
%
% - WGTS parameters
%  - Electron Scattering
%  - Covariance Matrix for Scattering Probabilities
%
% - MACE parameters
%  - Transmission Function
%  - KATRIN Response Function
%  - Covariance Matrix for WGTS+MACE
%
% Th. Lasserre, 2018
% CEA - Saclay
% TUM - IAS
%
% ----------------------------------------------------------------------- %

classdef WGTSMACE < FPD & handle %!dont change superclass without modifying parsing rules!
    
    properties (Constant = true, Hidden = true)
        % Physics Constants
        NA              = 6.02214129e23 ;    % Avogadro number
        kboltz          = 1.380e-23;         % Boltzmann, J/K
        bohrR           = 52.917721092e-12;  % Bohr Radius, m
        Eryd            = 13.60569252e-3;    % Rydberg energy, keV
    end
    
    properties (Access=public)
        
        % WGTS
        
        WGTS_Tp            ;       % T-T purity
        WGTS_DTHTr         ;       % D/H ratio
        WGTS_FTR_cm        ;       % Radius Flux Tube in cm
        WGTS_CD_MolPerCm2  ;       % Column Density in molecules/cm^2
        WGTS_CD_MolPerCm2allPixels;% Column Density in molecules/cm^2 per pixel, vector
        WGTS_CD_MolPerCm2_SubRun;  % Column Density in molecules/cm^2 per subrun
        WGTS_CDDist        ;       % Column Density Distribution per Pixel, vector
        WGTS_B_T           ;       % Source B field, Tesla
        WGTS_ErrB_T        ;       % Source B field, Tesla
        WGTS_Temp          ;       % [K] T2 temperature
        
        % Isotopologues
        WGTS_fTT           ;       % T-T fraction in WGTS
        WGTS_fHT           ;       % H-T fraction in WGTS
        WGTS_fDT           ;       % D-T fraction in WGTS
        WGTS_fHD           ;       % H-D fraction in WGTS
        
        % LARA Input
        WGTS_MolFrac_TT    ;    % T-T molecular fraction in WGTS
        WGTS_MolFrac_TT_SubRun;
        WGTS_MolFrac_DT    ;    % D-T molecular fraction in WGTS
        WGTS_MolFrac_DT_SubRun;
        WGTS_MolFrac_HT    ;    % H-T molecular fraction in WGTS
        WGTS_MolFrac_Tm;        % T^2  ion fraction in WGTS
        WGTS_MolFrac_Tm_i;      % init value for fit
        WGTS_MolFrac_HT_SubRun;
        WGTS_epsT          ;    % Tritium purity (nomenclature consistant with LARA)
         
        
        WGTS_MolFracRelErr_TT;  % Uncertainty on T-T molecular fraction in WGTS
        WGTS_MolFracRelErr_DT;  % Uncertainty on D-T molecular fraction in WGTS
        WGTS_MolFracRelErr_HT;  % Uncertainty on H-T molecular fraction in WGTS
        
        
        % FSD's - Decay to final states
        TTFSD              ;    % T-T Final States, 0(no), >0 (different computations)
        TTGSTh             ;    % T-T GS/ES Separation (bin)
        TTexP;TTexE        ;    % T-T FSD excitation energy and probability vectors
        DTFSD              ;    % D-T Final States, 0(no), >0 (different computations)
        DTGSTh             ;    % D-T GS/ES Separation (bin)
        DTexP;DTexE        ;    % D-T FSD excitation energy and probability vectors
        HTFSD              ;    % H-T Final States, 0(no), >0 (different computations)
        HTGSTh             ;    % H-T GS/ES Separation (bin)
        HTexP;HTexE        ;    % H-T FSD excitation energy and probability vectors
        TmFSD              ;    % T minus ion (OFF or SAENZ)
        TmGSTh             ;    % T- GS/ES Separation (bin)
        TmexP;TmexE        ;    % T- FSD excitation energy and probability vectors
        
        % WGTS Scattering Probabilities - Nominal
        NIS;                    %Number of Inelastic Scatterings
        ISXsection;             % inelastic scattering cross section
        is_P0  = 0.413339;
        is_P1  = 0.292658;
        is_P2  = 0.167332;
        is_P3  = 0.079127;
        is_P4  = 0.031777;
        is_P5  = 0.011086;
        is_P6  = 0.003423;
        is_P7  = 0.000949;
        is_P8  = 0.000239;
        is_P9  = 0.000054762381840;
        is_P10 = 0.000011677055841;
        is_Pv; % Vector of Inelastic Scattering Probabilites
        ISCS; % Inelastic Scattering Cross Section
        
        % Energy Loss Function Parametrization  Aseev/Abdurashitov
        %  (Inelastic Scattering - Aseev)
        is_A1; is_A1e;
        is_A2; is_A2e;
        is_w1; is_w1e;
        is_w2; is_w2e;
        is_eps1;is_eps1e;
        is_eps2;is_eps2e;
        is_epsc;
        is_EOffset;
        
        % Energy Loss Function Parametrization - Ch. Weinheimer 
        parCW_GLT;
        parCW_G2LT
        
        % Energy Loss Function Parametrization - KATRIN D2 - TOF
        parKatrinD2;
        errKatrinD2; 
         % Energy Loss Function Parametrization - KATRIN T2 - TOF+Integral
        parKatrinT2;
        errKatrinT2; 
        Ei; % Ionization Energy (eV)
        
        % Energy Loss Function 
        ELossFlag; % Aseev, Abdurashitov
        fscat    = [];
        
        % MACE
        MACE_Bmax_T        ;    % Pinch B field, Tesla
        ErrMACE_Bmax_T     ;    % Uncertainty on Pinch B field, Tesla
        MACE_Ba_T          ;    % Analysis Plane B field, Tesla
        ErrMACE_Ba_T       ;    % Uncertainty on Analysis Plane B field, Tesla
        MACE_R_eV          ;    % Mac-E Energy resolution, eV
        ErrMACE_R_eV       ;    % Uncertainty on Mac-E Energy resolution, eV
        MACE_Sigma         ;    % Broadening transmission function in eV            
        
        % HV Ripples
        HVRipples          ;    % Flag HV Ripples
        HVRipplesP2PV      ;    % HV Ripple Value Peak-2-Peak, Volt, Sin
        
        % E and B fields values in Analysis Plane - Pixel Wise
        MACE_Ba_Setting     ; % 2.7G by default, as first trial
        Pixel_MACE_Ba_TCorr ; % Analysis Plane B field, Tesla
        Pixel_MACE_Ba_VCorr = zeros(1,148) ; % Analysis Plane E field, Tesla
        MACE_Ba_TallPixels;
        SynchrotronFlag; % ON/OFF
        recomputeRF;           % Flag to recompute Response Function, even if there is one saved in cache already
        
    end
    
    properties (Dependent=true,Access=public)
    end
    
    properties (Access=public)
        
        % KATRIN Response Function
        KTFFlag             ;    % Flag defining Transmission Function
        KTF                 ;    % (KATRIN) Transmission Function
    end
    
    methods
        function obj = WGTSMACE(varargin)
            
            p = inputParser;
            % WGTS Parameters
            p.addParameter('WGTS_Tp',0.95,@(x)isfloat(x) && x>0); %0.95
            p.addParameter('WGTS_DTHTr',0,@(x)isfloat(x) && x>=0);
            p.addParameter('WGTS_FTR_cm',4.5,@(x)isfloat(x) && x>0);
            p.addParameter('WGTS_CD_MolPerCm2',5e17,@(x)isfloat(x) && x>0);
            p.addParameter('WGTS_CD_MolPerCm2_SubRun',[]);
            p.addParameter('WGTS_CDDist',0,@(x)isfloat(x) && sum(x)>0 && sum(x)<=148);
            p.addParameter('WGTS_B_T',3.6,@(x)isfloat(x) && x>0); % Source
            p.addParameter('WGTS_Temp',30,@(x)isfloat(x) && x>0);
            p.addParameter('NIS',7,@(x)isfloat(x) && x>0); %Numer inel. Scatterings
            p.addParameter('ISCS','Aseev',@(x)ismember(x,{'Aseev','Theory','Edep'})); %inel. scattering cross section
            p.addParameter('ISXsection',0,@(x)isfloat(x) || isa(x,'function_handle'));
            p.addParameter('is_EOffset',0,@(x)isfloat(x)); %energy loss offset
            
            % LARA Input
            % (Default values are first tritium values)
            p.addParameter('WGTS_MolFrac_TT',0,@(x)isfloat(x) && x>=0);
            p.addParameter('WGTS_MolFrac_TT_SubRun',[]);
            p.addParameter('WGTS_MolFrac_DT',1.4e-2,@(x)isfloat(x) && x>=0); % value from LARA-email 1.4e-2, 2e-2 is to get 1% tritium purity.
            p.addParameter('WGTS_MolFrac_DT_SubRun',[]);
            p.addParameter('WGTS_MolFrac_HT',0,@(x)isfloat(x) && x>=0);
            p.addParameter('WGTS_MolFrac_HT_SubRun',[]);
            p.addParameter('WGTS_MolFrac_Tm',0,@(x)isfloat(x) && x>=0);
            p.addParameter('WGTS_MolFracRelErr_TT',0,@(x)isfloat(x));
            p.addParameter('WGTS_MolFracRelErr_DT',0.06e-2,@(x)isfloat(x));
            p.addParameter('WGTS_MolFracRelErr_HT',0,@(x)isfloat(x));
            % WGTS: Flags for FSD: T-T / D-T / H-T
            p.addParameter('TTFSD','BlindingKNM2',@(x)ismember(x,{'OFF','DOSS','SAENZ','SAENZNOEE','ROLL','BlindingKNM1','WGTS100K','Sibille','Sibille0p5eV','SibilleFull','BlindingKNM2'}));
            p.addParameter('DTFSD','BlindingKNM2',@(x)ismember(x,{'OFF','DOSS','ROLL','HTFSD','TTFSD','BlindingKNM1','WGTS100K','Sibille','Sibille0p5eV','SibilleFull','BlindingKNM2'}));
            p.addParameter('HTFSD','BlindingKNM2',@(x)ismember(x,{'OFF','SAENZ','ROLL','BlindingKNM1','WGTS100K','Sibille','Sibille0p5eV','SibilleFull','BlindingKNM2'}));
            p.addParameter('TmFSD','SAENZ',@(x)ismember(x,{'OFF','SAENZ'}));
           
            % MACE Parameters
            p.addParameter('MACE_Ba_Setting','6.0G',@(x)ismember(x,{'2.7G','6.0G','Data'}));
            p.addParameter('MACE_Bmax_T',6,@(x)all(isfloat(x))); % Pinch
            p.addParameter('MACE_Ba_T',3e-4,@(x)all(isfloat(x))); % Analyzing plane 
            p.addParameter('MACE_Sigma',0,@(x)all(isfloat(x)));   
            p.addParameter('Pixel_MACE_Ba_TCorr',zeros(1,148),@(x)isfloat(x) && sum(x)>=0); % Analyzing plane all pixels
            p.addParameter('MACE_R_eV',2,@(x)isfloat(x) && x>0); %2
            p.addParameter('HVRipples','ON',@(x)ismember(x,{'OFF','ON'}));
            p.addParameter('HVRipplesP2PV',0.52,@(x)isfloat(x) && x>0); % V or eV
            p.addParameter('ELossFlag','Abdurashitov',@(x)ismember(x,{'KatrinT2','KatrinD2','Aseev','Abdurashitov','CW_GLT','CW_G2LT'}));
            % Both WGTS & MACE
            p.addParameter('KTFFlag','MACE',@(x)ismember(x,{'OFF','SSC',...
                'SSC2','SSCMW','MACE','MACER','Kle','SSCW_DCOfficial','SSC_BC','WGTSMACE'}));
            
            p.addParameter('recomputeRF','OFF',@(x)ismember(x,{'OFF','ON'}));
            p.addParameter('SynchrotronFlag','ON',@(x)ismember(x,{'OFF','ON'}));

            % Parse unmatched parameters to FPD.m
            p.KeepUnmatched=1;
            p.parse(varargin{:});
            if( isempty(fieldnames(p.Unmatched))); Unmatched={}; else
                Unmatched = reshape(...
                    [fieldnames(p.Unmatched),struct2cell(p.Unmatched)]',...
                    [1,length(fieldnames(p.Unmatched))*2]);
            end
            obj=obj@FPD(Unmatched{:}); %Parse to superclass FPD.m
            
            % WGTS Parameters
            obj.WGTS_Tp            = p.Results.WGTS_Tp;
            obj.WGTS_DTHTr         = p.Results.WGTS_DTHTr;
            obj.WGTS_FTR_cm        = p.Results.WGTS_FTR_cm;
            obj.WGTS_CD_MolPerCm2  = p.Results.WGTS_CD_MolPerCm2;
            obj.WGTS_CD_MolPerCm2_SubRun= p.Results.WGTS_CD_MolPerCm2_SubRun;
            obj.WGTS_CDDist        = p.Results.WGTS_CDDist;
            obj.WGTS_B_T           = p.Results.WGTS_B_T;
            obj.WGTS_Temp          = p.Results.WGTS_Temp;
            obj.is_Pv              = zeros(1,11); % 10 IS (first element is zero scattering)
            obj.NIS                = p.Results.NIS;
            obj.ISCS               = p.Results.ISCS;
            obj.ISXsection         = p.Results.ISXsection;
            
            % LARA Input
            obj.WGTS_MolFrac_TT    = p.Results.WGTS_MolFrac_TT;
            obj.WGTS_MolFrac_TT_SubRun    = p.Results.WGTS_MolFrac_TT_SubRun;
            obj.WGTS_MolFrac_DT    = p.Results.WGTS_MolFrac_DT;
            obj.WGTS_MolFrac_DT_SubRun    = p.Results.WGTS_MolFrac_DT_SubRun;
            obj.WGTS_MolFrac_HT    = p.Results.WGTS_MolFrac_HT;
            obj.WGTS_MolFrac_HT_SubRun    = p.Results.WGTS_MolFrac_HT_SubRun;
            obj.WGTS_MolFrac_Tm     = p.Results.WGTS_MolFrac_Tm;
            obj.WGTS_MolFracRelErr_TT = p.Results.WGTS_MolFracRelErr_TT;
            obj.WGTS_MolFracRelErr_DT = p.Results.WGTS_MolFracRelErr_DT;
            obj.WGTS_MolFracRelErr_HT = p.Results.WGTS_MolFracRelErr_HT;
            obj.is_EOffset            = p.Results.is_EOffset;
            
            % TBD: Flag FSD's
            obj.TTFSD              = p.Results.TTFSD;
            obj.DTFSD              = p.Results.DTFSD;
            obj.HTFSD              = p.Results.HTFSD;
            obj.TmFSD              = p.Results.TmFSD;
            
            % MACE Parameters
            obj.MACE_Ba_Setting    = p.Results.MACE_Ba_Setting;
            obj.MACE_Bmax_T        = p.Results.MACE_Bmax_T;
            obj.MACE_Ba_T          = p.Results.MACE_Ba_T;
            obj.MACE_R_eV          = p.Results.MACE_R_eV;
            obj.MACE_Sigma         = p.Results.MACE_Sigma;
            obj.HVRipples          = p.Results.HVRipples;
            obj.HVRipplesP2PV      = p.Results.HVRipplesP2PV;
            obj.ELossFlag          = p.Results.ELossFlag;
            % Both WGTS & MACE
            obj.KTFFlag            = p.Results.KTFFlag;
            obj.recomputeRF        = p.Results.recomputeRF;
            obj.SynchrotronFlag    = p.Results.SynchrotronFlag;
            
            obj.WGTS_MolFrac_Tm_i =  obj.WGTS_MolFrac_Tm;
            
            obj.GetISXsection; % inel. scattering cross section
            
            if( strcmp(obj.quiet, 'OFF') )
                fprintf('------------------- End WGTS Constructor----------------- \n');
            end
            
            
        end % constructor
        
    end % methods
    
    methods       
        function        InitializeColumnDensity(obj)
            % If there is information on the column density per pixel, it is
            % taken into account here. It should be given as a vector in
            % 'obj.WGTS_CDDist'
            switch obj.FPD_Segmentation
                case 'OFF'
                    % Do nothing
                case {'SINGLEPIXEL','MULTIPIXEL'}
                    if all(obj.WGTS_CDDist == 0)
                        obj.WGTS_CD_MolPerCm2allPixels = ...
                            obj.WGTS_CD_MolPerCm2.*ones(obj.nPixels,1)./148;
                    else
                        obj.WGTS_CD_MolPerCm2allPixels = ...
                            obj.WGTS_CD_MolPerCm2.*obj.WGTS_CDDist;
                    end
                case 'RING'
                    if all(obj.WGTS_CDDist == 0)
                        switch obj.FPD_Ring
                            case 1
                                obj.WGTS_CD_MolPerCm2allPixels = ...
                                    obj.WGTS_CD_MolPerCm2.*ones(obj.nRings,1)*(4./148);
                            otherwise
                                obj.WGTS_CD_MolPerCm2allPixels = ...
                                    obj.WGTS_CD_MolPerCm2.*ones(obj.nRings,1)*(12./148);
                        end
                    else
                        obj.WGTS_CD_MolPerCm2allPixels = ...
                            obj.WGTS_CD_MolPerCm2.*obj.WGTS_CDDist;
                    end
            end
        end    
        function        SetPixel_MACE_BaEa(obj)
            
            % Set Analysis Plane Pixel-wise E/B fields
            %
            % Simulation done for 2.7G --> in the analysis plane
            % E-field values, Volt --> Return E(2.7G,30 KV)/30 KV
            % B-field value, Tesla --> Return B(2.7G)/2.7G
            %
            % Simulation for 6.0G also available
            
            switch obj.MACE_Ba_Setting
                case 'Data'
                    %do nothing, taken care in the data
                    return
               % case '2.7G' %doesnt work anymore
                %    obj.Pixel_MACE_Ba_TCorr = importdata('MACE_Ba_TCorr_30kV_CWeinheimer.mat')';
                 %   obj.Pixel_MACE_Ba_VCorr = ones(1,148)*0;
                case '6.0G'
                    % Magnetic corrections for First Tritium
                    Ba_mpix = importdata([getenv('SamakPath'),'/inputs/WGTSMACE/MACE_Ba_TCorr6G_FT.mat']);
                    % Magnetic corrections for KNM1 TEMPORARY
                    %Ba_mpix = importdata([getenv('SamakPath'),'/inputs/WGTSMACE/MACE_Ba_TCorr6G_KNM1temporary.mat']);
                    obj.Pixel_MACE_Ba_TCorr = Ba_mpix./mean(Ba_mpix);
                    switch obj.FPD_Segmentation
                        case 'OFF'
                            obj.MACE_Ba_T    = mean(Ba_mpix(obj.FPD_PixList));
                            obj.ErrMACE_Ba_T = rms(obj.Pixel_MACE_Ba_TCorr(obj.FPD_PixList).*obj.MACE_Ba_T);
                        case {'MULTIPIXEL','SINGLEPIXEL'}
                            obj.MACE_Ba_T = Ba_mpix(obj.FPD_PixList);
                        case 'RING' 
                            obj.MACE_Ba_T    =  cell2mat(cellfun(@(x) mean(Ba_mpix(x)),obj.FPD_RingPixList,'UniformOutput',false)');
                            obj.ErrMACE_Ba_T = cell2mat(cellfun(@(x) rms(Ba_mpix(x)),obj.FPD_RingPixList,'UniformOutput',false)');          
                    end
            end
        end 
        function        ComputeIsotropologActivityWeight(obj)
             %% Decay to Molecular Excited States
            obj.WGTS_fTT = 2*obj.WGTS_Tp-1;
            obj.WGTS_fDT = 2.*(1-obj.WGTS_Tp).*obj.WGTS_DTHTr./(1+obj.WGTS_DTHTr);
            obj.WGTS_fHT = 2.*(1-obj.WGTS_Tp)./(1+obj.WGTS_DTHTr);
            obj.WGTS_fHD = 0;
            
            switch obj.TTFSD
                case 'OFF'
                    obj.WGTS_fTT = 0;
            end
            
            switch obj.DTFSD
                case 'OFF'
                    obj.WGTS_fDT = 0;
            end
            
            switch obj.HTFSD
                case 'OFF'
                    obj.WGTS_fHT = 0;
            end
        end    
        function        ComputeTritiumPurity(obj)
            obj.WGTS_epsT = obj.WGTS_MolFrac_TT + 0.5*obj.WGTS_MolFrac_DT + 0.5*obj.WGTS_MolFrac_HT + obj.WGTS_MolFrac_Tm;
        end
        
        %% Transmission Function (FPD-Segmentation dependent)
        function GetISXsection(obj)
            % define inel. scattering cross section according to ISCS Flag
            if strcmp(obj.ISCS,'Edep')
                % energy dependent inelastic scattering cross section
                % Unit: (m^2),  energy E in eV
                % Reference: High energy (Born approximation) electron inelastic total cross section
                %            (Bethe 1930, 1932; Inokuti71 eq. 4.41, 4.55, IKP67, Liu73, Liu87)
                MtotSq = 1.5356; % nuclear matrix element for T2
                dE = -0.0097; % relativistic and 1/E^^ correction near tritium end point
                cTot = 1.18;
                Ekin = @(E) 0.5*obj.me.*(1-obj.me.^2./(obj.me+E).^2);
                Eryd = obj.Eryd*1e3; % rydberg energy in eV
                obj.ISXsection = @(E) (4*pi*obj.bohrR^2)./(Ekin(E)./Eryd).*...
                    (MtotSq.*log(4*cTot.*Ekin(E)./Eryd)+dE);
            else
                if obj.ISXsection==0 || isempty(obj.ISXsection)
                    % inelastic scattering cross section
                    switch obj.ISCS
                        case 'Aseev'
                            obj.ISXsection   = 3.42e-22; % Troitzk measurement. Obsolete!
                        case 'Theory'
                            obj.ISXsection   = 3.637e-22;%3.642e-22;  % NEW FERENC T2 versus H2
                    end
                end
                
                if ~isa(obj.ISXsection,'function_handle')
                    obj.ISXsection = @(E) obj.ISXsection;
                end
            end
        end
        function out    = GetRF(obj,varargin)
            % Return TF of a specific pixel or all pixels
            switch obj.KTFFlag
                case 'OFF'
                    out = @(te,qu,p) 1;
                case 'WGTSMACE'
                    out = @obj.ComputeRF;
                case 'MACE'
                    out = @obj.ComputeMaceTF;
                otherwise
                    switch obj.KTFFlag
                        case 'SSC'
                            ssc  = importdata([getenv('SamakPath'),'/inputs/ResponseFunction/TFSSC2.txt']);
                        case 'SSCMW'
                            ssc  = importdata([getenv('SamakPath'),'/inputs/ResponseFunction/TFSSCMW.txt']);
                        case 'SSC_BC'
                            if isempty(varargin); qU = 18560;
                            else; qU = round(varargin{1}); end
                            ssc  = load(['RF_ssc_nScat6_DetailedOff_Xsection3p46_Bs3p6T_Ba5p83051G_BMax6T_ColumnDensity_5E21_at',num2str(qU),'.mat']);
                            ssc = ssc.ssc(:,2:3);
                        case 'SSCW_DCOfficial'
                            ssc  = importdata([getenv('SamakPath'),'/inputs/ResponseFunction/ResponseFunction_nScat6_nDetScat_off_qU18475_IS_X-section_3p46_bySSC.txt']);
                        case 'SSCW_DC'
                            % Response Functions calculated by SSC for the
                            % data challenge that took place in the 34th
                            % Collaboration Meeting.
                            %
                            % It can recieve as input the retarding
                            % potential value in V, if none is given, the
                            % default value is 18575.
                            %
                            % The Response Function includes 6 scatterings,
                            % detailed OFF, inelastic scattering cross-section of 3.46e-22
                            % Bs = 3.6T, Ba = 9G, Bmax 6T.
                            %
                            %
                            if isempty(varargin)
                                qU = 18575;
                            else
                                qU = round(varargin{1});
                                if qU < 18545; qU = 18545; end
                                if qU > 18580; qU = 18580; end
                            end
                            ssc  = load([getenv('SamakPath'),'/inputs/ResponseFunction/RFSSC/RF_ssc_nScat6_DetailedOff_Xsection3p46_Bs3p6T_Ba9G_BMax6T_qU',num2str(qU),'.mat']);
                            ssc = ssc.ssc(:,2:3);
                        case 'Kle' % Marco Response Function for Comparison
                            inRF = l([getenv('SamakPath'),'/inputs/ResponseFunction/RF_Marco_151202.txt']); % With ES
                            e_in = inRF(:,1);  r_in = inRF(:,2);  %aseev_-30eV_3Gauss
                            %                  r_in = in(:,3);    %aseev_-30eV_9Gauss
                            ssc = [e_in; r_in];
                        case 'Test'
                            ssctmp  = load([getenv('SamakPath'),'/inputs/ResponseFunction/samakRF_4.5e+17cm2_NIS11_Bm4.2T_Bs2.52T_Ba6.3018G.mat']);
                            ssc(:,1) = ssctmp.E; ssc(:,2) = ssctmp.RFip(:,1);
                    end
                    essc  = ssc(:,1); tfssc = ssc(:,2);
                    ssctftemp = @(e)(e>=0).*interp1(essc,tfssc,e,'spline',0);
                    out   = @(te,qu,p)(te>=qu).*(ssctftemp(te-qu));
            end
            
        end
        function out    = ComputeISProb(obj,varargin)
            % Compute:
            % - Inelastic Scattering Probabilities in WGTS
            % If no further arguments: Init KATRIN setting
            % Possible Inputs for variation of KATRIN setting:
            % ISXsection
            % MACE_Bmax_T, WGTS_B_T
            % WGTS_CD_MolPerCm2
            GetSamakPath;
            p = inputParser;
            p.addParameter('WGTS_DensityProfile','Flat',@(x)ismember(x,{'Table','File','Flat'}));
            p.addParameter('WGTS_ZCells',500,@(x)isfloat(x) && x>=1);
            p.addParameter('MACE_Bmax_T',obj.MACE_Bmax_T,@(x)isfloat(x));
            p.addParameter('WGTS_B_T',obj.WGTS_B_T,@(x)isfloat(x));
            p.addParameter('ISXsection',obj.ISXsection,@(x)isfloat(x) || isa(x,'function_handle'));
            p.addParameter('WGTS_CD_MolPerCm2',obj.WGTS_CD_MolPerCm2,@(x)isfloat(x) & x>0);
            p.addParameter('saveFile','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Method','Interp',@(x)ismember(x,{'Exact','Interp'}));
            p.addParameter('Energy',18575,@(x)all(isfloat(x)));
            
            p.parse(varargin{:});
            WGTS_DensityProfile        = p.Results.WGTS_DensityProfile;
            WGTS_ZCells                = p.Results.WGTS_ZCells;
            MACE_Bmax_T_local          = p.Results.MACE_Bmax_T;
            WGTS_B_T_local             = p.Results.WGTS_B_T;
            ISXsection_local           = p.Results.ISXsection;
            WGTS_CD_MolPerCm2_local    = p.Results.WGTS_CD_MolPerCm2;
            saveFile                   = p.Results.saveFile;
            Method                     = p.Results.Method;
            Energy                     = p.Results.Energy;
            
            % ---------------------------------------------------------
            % comment on New Method: (Lisa, Feb 2020)
            % -> no for loops and different integration style (less binning)
            % - when lambdaintegral already computed: Old method faster
            % - when lambdaintegral not calculated (e.g. in covariance matrix!), then Method 'New' is more than 2x faster
            % - both methods give the same result
            % ---------------------------------------------------------
            if strcmp(Method,'Interp')  
                try
                    %fprintf('Interpolation of inel. scattering probabilities ...')
                    out = ISProbInterp('WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2_local,...
                        'MACE_Bmax_T',MACE_Bmax_T_local,...
                        'WGTS_B_T',WGTS_B_T_local,...
                        'NIS',obj.NIS,...
                        'ISXsection', ISXsection_local(squeeze(Energy)),...
                        'SanityPlot','OFF');
                     %fprintf('succesful! \n')
                    return
                catch
                    fprintf(2,'failed - calculate exact inel. scattering probabilities \n')
                end  
            end
            
            switch WGTS_DensityProfile
                case 'File'
                    % KATRIN WGTS density profile - READ FILE
                    wgtsdp = importdata([getenv('SamakPath'),'/inputs/WGTSMACE/wgts_density_profile_for_thierry.dat']);
                case 'Table'
                    % KATRIN WGTS density profile - READ TABLE
                    wgtsdp(:,1) = (0:1:99)';
                    wgtsdp(:,2) = [5.16803e+19            8.76746e+19             1.2001e+20            1.49401e+20            1.76413e+20            2.01489e+20            2.24979e+20            2.47154e+20            2.68225e+20            2.88356e+20            3.07673e+20            3.26274e+20            3.44236e+20             3.6162e+20            3.78475e+20             3.9484e+20            4.10751e+20            4.26236e+20            4.41322e+20            4.56034e+20            4.70395e+20            4.84426e+20            4.98148e+20            5.11581e+20            5.24743e+20            5.37651e+20            5.50323e+20            5.62773e+20            5.75016e+20            5.87064e+20            5.98927e+20            6.10617e+20            6.22141e+20  6.335059999999999e+20            6.44719e+20            6.55785e+20            6.66708e+20            6.77492e+20             6.8814e+20  6.986540000000001e+20  7.090380000000001e+20            7.19295e+20            7.29427e+20             7.3944e+20  7.493379999999999e+20  7.591260000000001e+20            7.68812e+20            7.78401e+20            7.87904e+20  7.973299999999999e+20  7.973299999999999e+20            7.87904e+20            7.78401e+20            7.68812e+20  7.591260000000001e+20  7.493379999999999e+20             7.3944e+20            7.29427e+20            7.19295e+20  7.090380000000001e+20  6.986540000000001e+20             6.8814e+20            6.77492e+20            6.66708e+20            6.55785e+20            6.44719e+20  6.335059999999999e+20            6.22141e+20            6.10617e+20            5.98927e+20            5.87064e+20            5.75016e+20            5.62773e+20            5.50323e+20            5.37651e+20            5.24743e+20            5.11581e+20            4.98148e+20            4.84426e+20            4.70395e+20            4.56034e+20            4.41322e+20            4.26236e+20            4.10751e+20             3.9484e+20            3.78475e+20             3.6162e+20            3.44236e+20            3.26274e+20            3.07673e+20            2.88356e+20            2.68225e+20            2.47154e+20            2.24979e+20            2.01489e+20            1.76413e+20            1.49401e+20             1.2001e+20            8.76746e+19            5.16803e+19]';
                case 'Flat'
                    nBins = 1000;
                    wgtsdp(:,1) = linspace(0,99,nBins)';
                    wgtsdp(:,2) = ones(nBins,1)./wgtsdp(end,1)*5e22;
            end
            z_wgts = wgtsdp(:,1)/10;
            
            d = linspace(0,max(z_wgts),WGTS_ZCells); % Faster/Slower... 100 is good
            dir_lambda = [getenv('SamakPath'),'inputs/WGTSMACE/WGTS_CDProfileIntegral/'];
            if ~exist(dir_lambda,'dir')
                system(['mkdir ',dir_lambda]);
            end
            
            if strcmp(WGTS_DensityProfile,'Flat')
                rho_wgt = wgtsdp(1,2).*(WGTS_CD_MolPerCm2_local./5e17);
                lambdaInt = @(z) z.*rho_wgt;
            else
                rho_wgt = @(z) interp1(z_wgts,wgtsdp(:,2)/5e17*WGTS_CD_MolPerCm2_local,z); % uniform
                lambdaInt = @(z) arrayfun(@(z2) integral(rho_wgt,0,z2),z);
            end
            
            LambdaMax = lambdaInt(max(z_wgts)); % renormalize lambdaInt?
            lambdaInt = @(z) lambdaInt(z).*(1e4.*WGTS_CD_MolPerCm2_local./LambdaMax);
            
            % Effective Column Density
            lambda = @(z,theta) 1./cos(theta).*lambdaInt(z);
            thetamax    = asin(sqrt(WGTS_B_T_local/MACE_Bmax_T_local)); %maximum allowed starting angle
            
            % Mean scattering probability for electron starting at (z,theta)
            Pis_z_angle = @(i,z,theta) poisspdf(i,lambda(z,theta).*ISXsection_local(Energy));
            
            % Integration over angles
            ntmp = 720;
            f = @(i,z,theta) sin(theta).*Pis_z_angle(i,z,theta);
            Pis_z = @(i,z) 1./(1-cos(thetamax)).*simpsons(linspace(0,thetamax,ntmp)',f(i,z,linspace(0,thetamax,ntmp)')); %Integral Mean scattering probability over all angles
            
            % Integratio over z
            rho_wgt2 = @(z) arrayfun(@(z2) interp1(z_wgts,wgtsdp(:,2)/5e17*WGTS_CD_MolPerCm2_local,z2),z); % uniform
            Intfun = arrayfun(@(y) @(x) Pis_z(y-1,x).*rho_wgt2(x)./(WGTS_CD_MolPerCm2_local*1e4),1:obj.NIS+1,'UniformOutput',0);
            Pis_m = cellfun(@(y) 100*integral(y,0,max(z_wgts),'ArrayValued',1),Intfun,'UniformOutput',false)';
            Pis_m  = squeeze(cell2mat(Pis_m));
            obj.is_Pv(1:obj.NIS+1) = mean(Pis_m,2);
            
            % Save
            %Pis_mean = mean(Pis_m,2); % Averaged Probabilites
            ISProb_dir = [getenv('SamakPath'),sprintf('/inputs/WGTSMACE/WGTS_ISProb/')];
            MakeDir(ISProb_dir);
            
            if strcmp(obj.ISCS,'Edep')
                Estep = Energy(2)-Energy(1);
                IsXstr  = sprintf('Edep-Xsection-max%.0feV_Xstep%.1feV',max(Energy)-18575,Estep);
            else
                IsXstr = sprintf('%.5g-Xsection',ISXsection_local(E));
            end
            mystr = [ISProb_dir,sprintf('IS_%.5g-molPercm2_%s_%.0f-NIS_%.3g-Bmax_%.3g-Bs.mat',...
                WGTS_CD_MolPerCm2_local,IsXstr,obj.NIS+1,MACE_Bmax_T_local, WGTS_B_T_local)];
            if strcmp(saveFile,'ON')
                ISXsection_local = ISXsection_local(Energy);
                Energy = squeeze(Energy);
                save(mystr,'Pis_m','WGTS_CD_MolPerCm2_local','ISXsection_local','Energy','thetamax');
            end
            out = Pis_m;
            
        end
        function MaceTF    = ComputeMaceTF(obj,te,qu,varargin)
            % Compute MACE Filter Transmission Function - M. Slezak            
            % Switch depending on the segmentation 
            % (to include pixel-dependent magnetic field corrections)
           % if isscalar(varargin{1}); varargin = {'pixel',varargin{1}}; end
            p = inputParser;
            p.addParameter('MACE_Ba_T',obj.MACE_Ba_T, @(x)all(isfloat(x)));
            p.addParameter('MACE_Bmax_T',obj.MACE_Bmax_T,@(x)isfloat(x));
            p.addParameter('WGTS_B_T',obj.WGTS_B_T,@(x)isfloat(x));
            p.addParameter('pixel',1,@(x)isfloat(x) && x>0);
            p.addParameter('SynchrotronLoss',40e-03,@(x)isfloat(x) && x>=0);  % Renormalization, for Test
            p.addParameter('gammaFacApprox','',@(x)isfloat(x) && x>=0 || isempty(x));   % Approximation for Relativist Transmission (Gamma+1)/2=1.018

            p.parse(varargin{:});
            
            MACE_Ba_T_local     = p.Results.MACE_Ba_T;
            MACE_Bmax_T_local   = p.Results.MACE_Bmax_T;
            WGTS_B_T_local      = p.Results.WGTS_B_T;
            pixel               = p.Results.pixel; % pixel or a ring number
            SynchrotronLoss     = p.Results.SynchrotronLoss;
            gammaFacApprox      = p.Results.gammaFacApprox;
            
            if isempty(gammaFacApprox)
                gammaFacApprox = ((qu+obj.me)/obj.me+1)/2;
            end
            
            switch obj.FPD_Segmentation
                case 'OFF'
                    Ba = MACE_Ba_T_local;
                case {'RING','MULTIPIXEL','SINGLEPIXEL'}
                    Ba = MACE_Ba_T_local(pixel);
            end
            
            switch obj.SynchrotronFlag
                case 'ON'
                    % Samak Simulation Providing Synchrotron Loss in WGTS
                    PinchAngle =   [0.017453    0.034907     0.05236    0.069813    0.087266     0.10472     0.12217     0.13963     0.15708     0.17453     0.19199     0.20944     0.22689     0.24435      0.2618     0.27925     0.29671     0.31416     0.33161     0.34907     0.36652     0.38397     0.40143     0.41888     0.43633     0.45379     0.47124     0.48869     0.50615      0.5236     0.54105     0.55851     0.57596     0.59341     0.61087     0.62832     0.64577     0.66323     0.68068     0.69813     0.71558     0.73304     0.75049     0.76794      0.7854     0.80285      0.8203     0.83776     0.85521     0.87266];
                    SeLoss     =   [8.4048e-06   3.363e-05  7.5709e-05   0.0001347  0.00021067  0.00030374  0.00041402  0.00054166  0.00068685  0.00084979   0.0010307   0.0012299   0.0014476   0.0016842   0.0019401   0.0022157   0.0025114   0.0028278   0.0031655   0.0035251   0.0039073   0.0043129   0.0047428    0.005198   0.0056796   0.0061889   0.0067271   0.0072959   0.0078971   0.0085326   0.0092046   0.0099157    0.010669    0.011468    0.012316    0.013218    0.014179    0.015206    0.016306    0.017488    0.018763    0.020145     0.02165    0.023299    0.025119    0.027147    0.029428    0.032031    0.035047    0.038618];
                    DEsyn   = @(angleRad) interp1(PinchAngle,SeLoss,angleRad,'spline','extrap'); % data go till 50 deg, but computation needs 51 deg or so --> spline extrapolation
                    
                    % Define Edge Location of Transmission MACE Filter
                    qumax = qu*(1+Ba/MACE_Bmax_T_local*gammaFacApprox) + SynchrotronLoss;
                    
                    % Iterative Solving for Synch. Loss in Transmission Slope
                    ThetaTR_1 = real(((te >= qu) & (te <= qumax) ) .* (asin(sqrt((te-qu-DEsyn(0))./te*WGTS_B_T_local/Ba*gammaFacApprox))));
                    ThetaTR_2 = real(( (te >= qu) & (te <= qumax) ) .* asin(sqrt((te-qu-DEsyn(ThetaTR_1))./te*WGTS_B_T_local/Ba*gammaFacApprox)));
                    ThetaTR_3 = real(( (te >= qu) & (te <= qumax) ) .* asin(sqrt((te-qu-DEsyn(ThetaTR_2))./te*WGTS_B_T_local/Ba*gammaFacApprox)));
                    ThetaTR_4 = real(( (te >= qu) & (te <= qumax) ) .* asin(sqrt((te-qu-DEsyn(ThetaTR_3))./te*WGTS_B_T_local/Ba*gammaFacApprox)));
                    SynchrotronCorr = DEsyn(ThetaTR_4)./max(DEsyn(ThetaTR_4)).*SynchrotronLoss;
                    
                    %             tmpA=( (te >= qu) & (te <= qumax) ).*ThetaTR_4;tmpA=tmpA(tmpA>0);
                    %             tmpS=( (te >= qu) & (te <= qumax) ).*SynchrotronCorr;tmpS=tmpS(tmpS>0);
                    %             f = @(tdeg) interp1(tmpA,tmpS,tdeg*pi/180);
                case 'OFF'
                    qumax = qu*(1+Ba/MACE_Bmax_T_local*gammaFacApprox);
                    SynchrotronCorr = 0;
            end
            % Compute Transmission of MACE
            MaceTF = 0 + ((te >= qu) & (te <= qumax)) .* ...
                (1 - sqrt(1-(te-qu-SynchrotronCorr)./te * WGTS_B_T_local./Ba / gammaFacApprox )) ...
                ./ (1 - sqrt(1 - WGTS_B_T_local/MACE_Bmax_T_local)) + ...
                1 * (te > qumax);
            
            % if obj.MACE_Sigma>0
            %     out = TF_Convfun(MaceTF,qu',te,obj.MACE_Sigma);
            % else
            out = MaceTF;
            % end
            
            out(out>1)=1; 
        end
        %% Computating of empirical energy loss function and convolutions
        % Parameterization see SSC paper
        function [fscatn,fscatnE] = ComputeELossFunction(obj,varargin)
            % output: fscatn  -> energy loss function handle for all scatterings
            %         fscatnE -> energy loss function evaluated at energy E
            
            p = inputParser;
            p.addParameter('E','',@(x)isfloat(x));                                 % Energy (binning)
            p.addParameter('Debug','OFF',@(x)ismember(x,{'ON','OFF'}));
            % usefull for covariance matrix
            p.addParameter('is_A1',obj.is_A1,@(x)isfloat(x));                    % eloss parameters: Aseev or Abdurashitov
            p.addParameter('is_A2',obj.is_A2,@(x)isfloat(x));
            p.addParameter('is_w1',obj.is_w1,@(x)isfloat(x));
            p.addParameter('is_w2',obj.is_w2,@(x)isfloat(x));
            p.addParameter('is_epsc',obj.is_epsc,@(x)isfloat(x));
            p.addParameter('is_eps1',obj.is_eps1,@(x)isfloat(x));
            p.addParameter('is_eps2',obj.is_eps2,@(x)isfloat(x));
            p.addParameter('parKatrinD2',obj.parKatrinD2,@(x)all(isfloat(x)));     % eloss parameters: KATRIN D2
            p.addParameter('parKatrinT2',obj.parKatrinT2,@(x)all(isfloat(x)));     % eloss parameters: KATRIN T2
            p.addParameter('LoadOrSaveEloss','ON',@(x)ismember(x,{'ON','OFF'}));   % Only for covariance matrices has to be off
            p.addParameter('is_EOffset',obj.is_EOffset,@(x)isfloat(x));
            p.parse(varargin{:});
            
            is_A1_local       = p.Results.is_A1;
            is_A2_local       = p.Results.is_A2;
            is_w1_local       = p.Results.is_w1;
            is_w2_local       = p.Results.is_w2;
            is_eps1_local     = p.Results.is_eps1;
            is_eps2_local     = p.Results.is_eps2;
            is_epsc_local     = p.Results.is_epsc;
            parKatrinD2_local = p.Results.parKatrinD2;
            parKatrinT2_local = p.Results.parKatrinT2;
            E                 = p.Results.E;
            is_EOffset_local  = p.Results.is_EOffset;
            Debug             = p.Results.Debug;
            LoadOrSaveEloss   = p.Results.LoadOrSaveEloss;
            
            if isempty(E) % default binning if non specified
                maxE = 500;%9288;
                minE=-maxE; NbinE = (maxE-minE)/0.2;
                E = linspace(minE,maxE,NbinE);
            end
            Estep = E(2) - E(1);
            
            % labeling
            path_eloss = [getenv('SamakPath'),sprintf('inputs/WGTSMACE/WGTS_ELossConv/')];
            MakeDir(path_eloss);
            file_eloss = cell(numel((obj.NIS):10),1);
            file_eloss{1} = [path_eloss, sprintf('ELoss_%s_%geV-%geV-qU_%.3geV-Binning_%g-NIS.mat',...
                obj.ELossFlag,min(E),max(E),Estep,obj.NIS)];
            if is_EOffset_local~=0
              file_eloss{1} = strrep(file_eloss{1},obj.ELossFlag,...
                  [obj.ELossFlag,sprintf('_%.3fELOffset',is_EOffset_local)]); 
            end
            
            NIStmp = (obj.NIS+1):10;
            for i=1:numel(NIStmp)
                file_eloss{i+1} = strrep(file_eloss{1},sprintf('%.0f-NIS',obj.NIS),sprintf('%.0f-NIS',NIStmp(i)));
            end
            
            % load if already existing1.0932e+17
            file_logic = find(isfile(file_eloss)); % tells which file exists
            
            if ~isempty(file_logic)  %if one exists
                tmp = importdata(file_eloss{file_logic});
            end
            
            if ~isempty(file_logic)  && strcmp(obj.recomputeRF,'OFF') && strcmp(LoadOrSaveEloss,'ON') && numel(E)==numel(tmp.E) % if one exists and has correct binning
                fscatnE = tmp.fscatnE;
                fscatn  = tmp.fscatn;
                if strcmp(Debug,'ON')
                fprintf(2,'loading energy loss function from file \n')
                end
            else % compute energy loss otherwise
                
                fscatn = cell(obj.NIS,1);
                switch obj.ELossFlag
                    case {'KatrinD2'}
                        f1scat  = @(e) obj.eloss_katrin(e, parKatrinD2_local(1), parKatrinD2_local(2), parKatrinD2_local(3), parKatrinD2_local(4),...
                            parKatrinD2_local(5), parKatrinD2_local(6), parKatrinD2_local(7), parKatrinD2_local(8), parKatrinD2_local(9));
                        f1scatn = @(e) f1scat(e);%function is already properly normalized
                    case {'KatrinT2'}
                        f1scat  = @(e) obj.eloss_katrin(e, parKatrinT2_local(1), parKatrinT2_local(2), parKatrinT2_local(3), parKatrinT2_local(4),...
                            parKatrinT2_local(5), parKatrinT2_local(6), parKatrinT2_local(7), parKatrinT2_local(8), parKatrinT2_local(9));
                        f1scatn = @(e) f1scat(e);%function is already properly normalized
                    case {'Aseev','Abdurashitov'}
                        if max(E)<9288
                            tailFactor = GetELossTailFactor(obj,max(E));
                        else
                            tailFactor = 1;
                        end
                        f1scat  = @(e) (e<is_epsc_local).*is_A1_local.*exp(-2.*((e-is_eps1_local)/(is_w1_local)).^2)...
                            + (e>=is_epsc_local).*is_A2_local.*(is_w2_local.^2)./(is_w2_local.^2+4.*(e-is_eps2_local).^2);
                        f1scatn = @(e) f1scat(e)/simpsons(e,f1scat(e)).*tailFactor; %normalized
                    case 'CW_GLT' % eloss_gauss_single_line_with_tail
                        gauss  = @(x,pos,sigma) exp(-(x-pos).^2./(2.*sigma.^2));
                        tail   = @(x,square,cube) (square./(x+1e-99).^2+cube./(x+1e-99).^3);
                        y1     = obj.parCW_GLT(1).*gauss(obj.parCW_GLT(4), obj.parCW_GLT(2),obj.parCW_GLT(3));
                        y2     = tail(obj.parCW_GLT(5),obj.parCW_GLT(6),obj.parCW_GLT(7));
                        f1scat = @(x) ...
                            (x<=0).*(0) + ...
                            (x<obj.parCW_GLT(4)).*obj.parCW_GLT(1).*gauss(x,obj.parCW_GLT(2),obj.parCW_GLT(3)) + ...
                            (x>=obj.parCW_GLT(4)).*( (x<=obj.parCW_GLT(5)).*(y1.*(obj.parCW_GLT(5)-x) + y2.*(x-obj.parCW_GLT(4))./(obj.parCW_GLT(5)-obj.parCW_GLT(4))) + (x>obj.parCW_GLT(5)).*tail(x,obj.parCW_GLT(6),obj.parCW_GLT(7)));
                        f1scatn = @(e) f1scat(e)/simpsons(e,f1scat(e)); %normalized
                    case 'CW_G2LT' % eloss_gauss_double_line_with_tail
                        gauss  = @(x,pos,sigma) exp(-(x-pos).^2./(2.*sigma.^2));
                        tail   = @(x,square,cube) (square./(x+1e-99).^2+cube./(x+1e-99).^3);
                        
                        y1     = obj.parCW_G2LT(1).*gauss(obj.parCW_G2LT(4), obj.parCW_G2LT(2),obj.parCW_G2LT(3));
                        y2     = obj.parCW_G2LT(9);
                        y3     = tail(obj.parCW_G2LT(5),obj.parCW_G2LT(6),obj.parCW_G2LT(7));
                        
                        f1scat = @(x) ...
                            (x<=0).*(0) + ...
                            (x<obj.parCW_G2LT(4)).*(obj.parCW_G2LT(1).*gauss(x,obj.parCW_G2LT(2),obj.parCW_G2LT(3))) + ...
                            (x>=obj.parCW_G2LT(4)).*...
                            ( ...
                            (x<=obj.parCW_G2LT(8)).*...
                            (( y1*(obj.parCW_G2LT(8)-x) + y2.*(x-obj.parCW_G2LT(4)))/(obj.parCW_G2LT(8)-obj.parCW_G2LT(4))) + ...
                            (x>obj.parCW_G2LT(8)).*...
                            ((x<=obj.parCW_G2LT(5)).*...
                            ((y2*(obj.parCW_G2LT(5)-x) + y3.*(x-obj.parCW_G2LT(8)))./(obj.parCW_G2LT(5)-obj.parCW_G2LT(8))) + ...
                            (x>obj.parCW_G2LT(5)).*...
                            (tail(x,obj.parCW_G2LT(6),obj.parCW_G2LT(7))))...
                            );
                        f1scatn = @(e) f1scat(e)/simpsons(e,f1scat(e)); %normalized
                end
                
                %Convolutions for multiple scatterings
                %tmp = @(tau) f1scat(tau).*f1scat(E'-tau);
                %f2scat = @(e) integral(tmp,-inf,inf,'ArrayValued',1);
                %f2scatn = @(e) f2scat(e);%/simpsons(e,f2scat(e));
                
                 fscatnE      = zeros(7,numel(E)); % evaluate functions at E
                 fscatnE(1,:) = f1scatn(E+ is_EOffset_local);  
               
                if obj.NIS>1
                    f2scat = @(e) conv(f1scat(e),f1scat(e),'same');
                    f2scatn = @(e) f2scat(e);%/simpsons(e,f2scat(e));
                    fscatn(1:2) = {f1scatn; f2scatn};
                    fscatnE(2,:) = f2scatn(E+ is_EOffset_local);
                end
                
                if obj.NIS>2
                    f3scat = @(e) conv(f2scat(e),f1scat(e),'same');
                    f3scatn = @(e) f3scat(e);%/simpsons(e,f3scat(e));
                    fscatn(1:3) = {f1scatn; f2scatn; f3scatn;};
                    fscatnE(3,:) = f3scatn(E+ is_EOffset_local);
                end
                
                if obj.NIS>3
                    f4scat = @(e) conv(f3scat(e),f1scat(e),'same');
                    f4scatn = @(e) f4scat(e);%/simpsons(e,f4scat(e));
                    fscatn(1:4) = {f1scatn; f2scatn; f3scatn; f4scatn};
                    fscatnE(4,:) = f4scatn(E+ is_EOffset_local);
                end
                
                if obj.NIS>4
                    f5scat = @(e) conv(f4scat(e),f1scat(e),'same');
                    f5scatn = @(e) f5scat(e);%/simpsons(e,f5scat(e));
                    fscatn(1:5) = {f1scatn; f2scatn; f3scatn; f4scatn; f5scatn};
                    fscatnE(5,:) = f5scatn(E+ is_EOffset_local);
                end
                
                if obj.NIS>5
                    f6scat = @(e) conv(f5scat(e),f1scat(e),'same');
                    f6scatn = @(e) f6scat(e);%/simpsons(e,f6scat(e));
                    fscatn(1:6) = {f1scatn; f2scatn; f3scatn; f4scatn; f5scatn; f6scatn};
                    fscatnE(6,:) = f6scatn(E+ is_EOffset_local);
                end
                
                if obj.NIS>6
                    f7scat = @(e) conv(f6scat(e),f1scat(e),'same');
                    f7scatn = @(e) f7scat(e);%/simpsons(e,f7scat(e));
                    fscatn(1:7) = {f1scatn; f2scatn; f3scatn; f4scatn; f5scatn;f6scatn;f7scatn};
                    fscatnE(7,:) = f7scatn(E+ is_EOffset_local);
                end
                
                if obj.NIS>7
                    f8scat = @(e) conv(f7scat(e),f1scat(e),'same');
                    f8scatn = @(e) f8scat(e);%/simpsons(e,f8scat(e));
                    fscatn(1:8) = {f1scatn; f2scatn; f3scatn; f4scatn; f5scatn;f6scatn;f7scatn;f8scatn};
                end
                
                if obj.NIS>8
                    f9scat = @(e) conv(f8scat(e),f1scat(e),'same');
                    f9scatn = @(e) f9scat(e);%/simpsons(e,f9scat(e));
                    fscatn(1:9) = {f1scatn; f2scatn; f3scatn; f4scatn; f5scatn;f6scatn;f7scatn;f8scatn;f9scatn};
                end
                
                if obj.NIS>9
                    f10scat = @(e) conv(f9scat(e),f1scat(e),'same');
                    f10scatn = @(e) f10scat(e);%/simpsons(e,f10scat(e));
                    fscatn(1:10) = {f1scatn; f2scatn; f3scatn; f4scatn; f5scatn; f6scatn; f7scatn; f8scatn; f9scatn; f10scatn};
                end
                
                % output   
                % save
                if strcmp(LoadOrSaveEloss,'ON')
                    save(file_eloss{1},'fscatn','fscatnE','E');
                end
            end
        end
        
        function  InitializeELossFunction(obj)
            switch obj.ELossFlag
                case 'KatrinD2'
                    % Ionization energy for D2 [eV]
                    obj.Ei = 15.467;  
                    
                    % Global fit of parametrized energy loss function to measured integral and t.o.f. data
                    % Author: V. Hannen, Data preparation: C. Rodenbeck, R. Sack, L. Schimpf
                    % Last updated: 28/5/2019
                    amp1 	= 0.0296231;
                    pos1 	= 11.8467;
                    sig1 	= 0.182839;
                    amp2 	= 0.273946;
                    pos2 	= 12.785;
                    sig2 	= 0.486158;
                    amp3 	= 0.0749376;
                    pos3 	= 14.9051;
                    sig3 	= 1.0307;

                    amp1err 	= 0.00182837;
                    pos1err 	= 0.0177364;
                    sig1err 	= 0.0139258;
                    amp2err 	= 0.00127919;
                    pos2err 	= 0.00393521;
                    sig2err 	= 0.00393997;
                    amp3err 	= 0.000485682;
                    pos3err 	= 0.00812498;
                    sig3err 	= 0.0270375;    
                    
                    %amptof 	= 1;     % not used
                    mutof 	= 0.104426;  % not used
                    norm 	= 0.997937;  % not used
                    mu1 	= 0.0964532; % not used
                    mu2 	= 0.656861;  % not used
                    mu3 	= 1.64586;   % not used
                    %tail 	= 0.972224;  % not used
                    n1 	    = 1.0099;    % not used
                    
                    mutoferr 	= 0.000253665;  % not used
                    normerr 	= 0.000856045;  % not used
                    mu1err   	= 0.000144627; % not used
                    mu2err  	= 0.000195662;  % not used
                    mu3err  	= 0.000136446;   % not used
                    n1err 	    = 8.14*1e-05;    % not used
                    
                    obj.parKatrinD2 = [amp1, pos1, sig1, amp2, pos2, sig2, amp3, pos3, sig3, mutof, norm, mu1, mu2, mu3, n1];
                    obj.errKatrinD2 = [amp1err, pos1err, sig1err, amp2err, pos2err, sig2err, amp3err, pos3err, sig3err, mutoferr, normerr, mu1err, mu2err, mu3err, n1err];
                case 'KatrinT2'
                    % Ionization energy for T2 [eV]
                    % Global fit of parametrized energy loss function to measured integral and t.o.f. data
                    % Author: V. Hannen, Data preparation: C. Rodenbeck, R. Sack, L. Schimpf
                    % value from git repo: https://nuserv.uni-muenster.de:8443/katrin-git/KATRIN-eloss/blob/master/code/KNM1/
                    % Last updated: 18/02/2020
                    obj.Ei = 15.487;
                    [obj.parKatrinT2, obj.errKatrinT2] = GetElossPara;
                case 'Aseev'
                    obj.is_A1    = 0.204;   obj.is_A1e    = 0.001; % Normalization
                    obj.is_A2    = 0.0556;  obj.is_A2e    = 0.0003;% Normalization
                    obj.is_w1    = 1.85;    obj.is_w1e    = 0.02;  % Width of excitation peak
                    obj.is_w2    = 12.5;    obj.is_w2e    = 0.1;   % Width of ionization peak
                    obj.is_eps1  = 12.6;    obj.is_eps1e  = 0;     % Position of excitation peak
                    obj.is_eps2  = 14.3;    obj.is_eps2e  = 0.02;  % Position of ionization peak
                    obj.is_epsc  = obj.ComputeELossEpsilonC(obj.is_A1,obj.is_A2,obj.is_w1,obj.is_w2,obj.is_eps1,obj.is_eps2);
                case 'Abdurashitov' %arXiv: 1603.04243v1
                    obj.is_A1    = 0.070*3.60;  % notation in paper: Norm*A
                    obj.is_A1e   = 0.070*0.15;
                    obj.is_A2    = 0.070;  % == Norm
                    obj.is_A2e   = 0;      % no uncertainty given on second Norm
                    obj.is_w1    = 1.22;     obj.is_w1e    = 0.05; % Width of excitation peak
                    obj.is_w2    = 11.99;    obj.is_w2e    = 0.13; % Width of ionization peak
                    obj.is_eps1  = 12.695;   obj.is_eps1e  = 0.017;% Position of excitation peak
                    obj.is_eps2  = 13.29;    obj.is_eps2e  = 0.18; % Position of ionization peak
                    %obj.NIS = 8;
                    obj.is_epsc  = obj.ComputeELossEpsilonC(obj.is_A1,obj.is_A2,obj.is_w1,obj.is_w2,obj.is_eps1,obj.is_eps2);
                case 'CW_GLT' % eloss_gauss_single_line_with_tail
                    %        par[0] amplitude of gaussian excitation peak
                    par0   =    +1.84234012E+03;  % error 8.25624701E+00
                    %        par[1] position of gaussian excitation peak
                    par1   =    +1.27186562E+01;  % error +2.54038843E-03
                    %        par[2] sigma of gaussian excitation peak
                    par2   =    +5.61673258E-01;  % error +1.87857676E-03
                    %        par[3] transition point between gaussian peak and line
                    par3   =    +1.36732536E+01;  % error +6.64567708E-03
                    %        par[4] transition point between line and tail
                    par4   =    +1.48591947E+01;  % error +4.86608385E-02
                    %        par[5] amplitude of tail 1/(DeltaE)**2
                    par5   =    -3.14460383E+03;  % error +1.21982053E+03
                    %        par[6] amplitude of tail 1/(DeltaE)**3
                    par6   =    +1.69971770E+06;  % error +2.30037595E+04
                    %        par[7] amplitude of twofold scattering
                    par7   = 0;
                    obj.parCW_GLT = [par0 par1 par2 par3 par4 par5 par6 par7];
                case 'CW_G2LT' % eloss_gauss_double_line_with_tail
                    %         par[0] amplitude of gaussian excitation peak
                    par0   =    +1.84765404E+03; % +8.35634889E+00
                    %         par[1] position of gaussian excitation peak
                    par1   =    +1.27059057E+01;% +3.14140376E-03
                    %         par[2] sigma of gaussian excitation peak
                    par2   =    +5.53393287E-01;% +2.13304261E-03
                    %         par[3 ]transition point between gaussian peak and line1
                    par3   =    +1.34445690E+01; % +1.89733987E-02
                    %         par[7] transition point between line1 and line2
                    par7   =    +1.46988076E+01; % +3.48353214E-02
                    %         par[4] transition point between line2 and tail
                    par4   =    +1.38481042E+01; % +1.51379094E-02
                    %         par[5] amplitude of tail 1/(DeltaE)**2
                    par5   =    -2.39538346E+03; % +1.17824390E+03
                    %         par[6] amplitude of tail 1/(DeltaE)**3
                    par6   =    +1.68490678E+06; % +2.20745331E+04
                    %
                    par8   =    +4.05446329E+02; % +7.77796257E+00
                    obj.parCW_G2LT = [par0 par1 par2 par3 par4 par5 par6 par7 par8];
            end
        end
        
        function is_epsc = ComputeELossEpsilonC(obj,is_A1,is_A2,is_w1,is_w2,is_eps1,is_eps2)
            % Find is_epsc: continuity of Gaussian and Lorentian in Energy
            % Loss Function parameterization 
            eqn = @(x) is_A1.*exp(-2.*((x-is_eps1)/(is_w1)).^2)-is_A2.*(is_w2.^2)./(is_w2.^2+4.*(x-is_eps2).^2);
            is_epsc = fzero(eqn,14);  
             %syms x; % symbolic value for equation solver
            %eqn = is_A1.*exp(-2.*((x-is_eps1)/(is_w1)).^2)==is_A2.*(is_w2.^2)./(is_w2.^2+4.*(x-is_eps2).^2);
            % is_epsc = double(vpasolve(eqn,x,14)); % 14 is staring point (close to expected epsilon_c)   
            %  clear x; clear eqn;
        end
        
        %Compute Response Function
        function out  = ComputeRF(obj,te,qu,varargin)  
            % Compute Response Function
            p = inputParser;
            p.addParameter('Debug','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('ELossBinStep',0.04,@(x)isfloat(x)); % for ELoss Convolution
            p.addParameter('ELossRange',500,@(x)isfloat(x));   %  ELoss range: old default 9288
            p.addParameter('RFBinStep',0.04,@(x)isfloat(x));    % for Final RF Convolution 0.04
            p.addParameter('AdjustISProba','OFF',@(x)ismember(x,{'ON','OFF'}));    % for Final RF Convolution
            p.addParameter('ElossFunctions','',@(x)isfloat(x)); % for covariance matrix only
            p.addParameter('MACE_Bmax_T',obj.MACE_Bmax_T,@(x)isfloat(x)); 
            p.addParameter('WGTS_B_T',obj.WGTS_B_T,@(x)isfloat(x));
            p.addParameter('MACE_Ba_T',obj.MACE_Ba_T, @(x)all(isfloat(x)));
            p.addParameter('ISXsection',obj.ISXsection,@(x)isfloat(x) || isa(x,'function_handle'));
            p.addParameter('WGTS_CD_MolPerCm2',obj.WGTS_CD_MolPerCm2,@(x)isfloat(x) & x>0);
            p.addParameter('pixel',1,@(x)isfloat(x));
            
            p.parse(varargin{:});
            
            Debug               = p.Results.Debug;
            ELossBinStep        = p.Results.ELossBinStep;
            ELossRange          = p.Results.ELossRange;
            RFBinStep           = p.Results.RFBinStep;
            AdjustISProba       = p.Results.AdjustISProba;
            pixel               = p.Results.pixel;
            ElossFunctions      = p.Results.ElossFunctions;
            MACE_Bmax_T_local          = p.Results.MACE_Bmax_T;
            WGTS_B_T_local             = p.Results.WGTS_B_T;
            ISXsection_local           = p.Results.ISXsection;
            WGTS_CD_MolPerCm2_local    = p.Results.WGTS_CD_MolPerCm2;
            MACE_Ba_T_local            = p.Results.MACE_Ba_T;
            
            % binning for response function
            if strcmp(obj.TD,'RFcomparison')
                RFBinStep = 0.01;
                maxE_rf = 90;
                minE_rf = -maxE_rf;
                Estep_rf = RFBinStep;
                E_rf = minE_rf:Estep_rf:maxE_rf;
            else
                maxE_rf = (obj.qUmax-obj.qUmin)*1.5;
                minE_rf=-maxE_rf;
                NbinE_rf = (maxE_rf-minE_rf)/RFBinStep;
                E_rf = linspace(minE_rf,maxE_rf,NbinE_rf);
                Estep_rf = E_rf(2) - E_rf(1);
            end
            
            % binning for ELoss function -> very wide for convolution  
            maxE = ELossRange;
            minE=-maxE; NbinE = (maxE-minE)/ELossBinStep;
            E = linspace(minE,maxE,NbinE); Estep = E(2) - E(1);

            % binning for isProb E==te-qu
            IsProbBinStep = 2;
            maxEis = 500;
            Eiscs = 18575+(-maxEis:IsProbBinStep:maxEis);
            
            % IS Probabilities: labeling
            file_pis = cell(numel((obj.NIS+1):11),1);
            if strcmp(obj.ISCS,'Edep')
                file_pis{1} = [getenv('SamakPath'), sprintf('inputs/WGTSMACE/WGTS_ISProb/IS_%.5g-molPercm2_Edep-Xsection-max%.0feV_Xstep%.1feV_%.0f-NIS_%.3g-Bmax_%.3g-Bs.mat',...
                    WGTS_CD_MolPerCm2_local,maxEis,IsProbBinStep,obj.NIS+1,MACE_Bmax_T_local,WGTS_B_T_local)];
            else
                file_pis{1} = [getenv('SamakPath'), sprintf('inputs/WGTSMACE/WGTS_ISProb/IS_%.5g-molPercm2_%.5g-Xsection_%.0f-NIS_%.3g-Bmax_%.3g-Bs.mat',...
                               WGTS_CD_MolPerCm2_local,ISXsection_local(18575),obj.NIS+1,MACE_Bmax_T_local,WGTS_B_T_local)];
            end
            NIStmp = (obj.NIS+2):11;
            for i=1:numel(NIStmp)
                file_pis{i+1} = strrep(file_pis{1},sprintf('%.0f-NIS',obj.NIS+1),sprintf('%.0f-NIS',NIStmp(i)));
            end
            
            % IS Probabilities: load from lookup table or compute new
            if strcmp(obj.recomputeRF,'OFF')
                % file names
                file_logic = find(isfile(file_pis));
                if numel(file_logic)>1
                    file_logic = file_logic(1);
                end
                % import/calculation
                if ~isempty(file_logic) % if one of the files exist: load it
                    w = load(file_pis{file_logic});
                else % if not: compute scattering probabilities
                    w.Pis_m = obj.ComputeISProb('Energy',reshape(Eiscs,[1,1,numel(Eiscs)]),...
                        'MACE_Bmax_T',MACE_Bmax_T_local,...
                        'WGTS_B_T',WGTS_B_T_local,...
                        'ISXsection',ISXsection_local,...
                        'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2_local);
                end
            elseif strcmp(obj.recomputeRF,'ON')
                w.Pis_m = obj.ComputeISProb('Energy',reshape(Eiscs,[1,1,numel(Eiscs)]),...
                    'MACE_Bmax_T',MACE_Bmax_T_local,...
                    'WGTS_B_T',WGTS_B_T_local,...
                    'ISXsection',ISXsection_local,...
                    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2_local);
            end
            
            % IS Probabilities: Adjust Probabilities with SSC (1.1e17mol/cm2)
            ratioMasterSamak = [0.999906038165416   0.999952728590095   0.999967392247124   0.999983000052758   0.999998403728429   1.000013233420951   1.000027704218478   1.065797995861047  ones(1,numel(w.Pis_m)-8)];
            switch AdjustISProba
                case 'OFF'
                    obj.is_Pv = w.Pis_m/100;
                case 'ON'
                    obj.is_Pv = w.Pis_m /100 .* ratioMasterSamak;
            end
              
            %% Retrieve/Compute Energy Loss Functions
            obj.recomputeRF='OFF';
            if isempty(ElossFunctions)
                [~,ElossFunctions] = obj.ComputeELossFunction('E',E); % load if already exists, otherwise compute
            end
             obj.recomputeRF='ON';
            if numel(obj.is_Pv)<8 && obj.NIS<7
                obj.is_Pv(obj.NIS+2:end) = 0;
                obj.is_Pv =  [obj.is_Pv;zeros(8,1)]; %add some zeros
            end
            
            % cut eloss function to response function window
             EindexStart = find(E>=minE_rf,1)-5;
             EindexStop  = find(E>=maxE_rf,1)+5;
             Eel = E(EindexStart:EindexStop);
             
             if ~strcmp(obj.ISCS,'Edep')
                 obj.is_Pv = repmat(obj.is_Pv,[1,numel(Eiscs)]);
             end
             
             ISProb = interp1(Eiscs',obj.is_Pv',qu+Eel)';
             
            scatterings  = ISProb(2,:).*ElossFunctions(1,EindexStart:EindexStop) + ... % first scattering
                ISProb(3,:).*ElossFunctions(2,EindexStart:EindexStop).*(E(2)-E(1))^1 + ...
                ISProb(4,:).*ElossFunctions(3,EindexStart:EindexStop).*(E(2)-E(1))^2 + ...
                ISProb(5,:).*ElossFunctions(4,EindexStart:EindexStop).*(E(2)-E(1))^3 + ...
                ISProb(6,:).*ElossFunctions(5,EindexStart:EindexStop).*(E(2)-E(1))^4 + ...
                ISProb(7,:).*ElossFunctions(6,EindexStart:EindexStop).*(E(2)-E(1))^5 + ...
                ISProb(8,:).*ElossFunctions(7,EindexStart:EindexStop).*(E(2)-E(1))^6;
            
            obj.fscat       = @(e)interp1(Eel,scatterings,e);
            
            % Retreive Transmission Function
            TF  = @obj.ComputeMaceTF;
            ISProb0 = interp1(Eiscs',obj.is_Pv(1,:)',qu+E_rf);
            
            RF = TF(qu+E_rf,qu,'pixel',pixel).*ISProb0 + ...
                conv(TF(qu+E_rf,qu,'pixel',pixel,...
                'MACE_Bmax_T',MACE_Bmax_T_local,'WGTS_B_T',WGTS_B_T_local,'MACE_Ba_T',MACE_Ba_T_local),...
                obj.fscat(E_rf),'same').*Estep_rf;
            
            switch obj.FPD_Segmentation
                case {'SINGLEPIXEL','MULTIPIXEL'}
                    fprintf('WGTS:ComputeRF: RF, pixel = %d, qU = %d, rhoD=%.4g, Nis=%.0f, Bmax=%.2f T, Bwgts:%.2f T \n',pixel,qu,obj.WGTS_CD_MolPerCm2,obj.NIS+1,obj.MACE_Bmax_T,obj.WGTS_B_T);
                case 'RING'
                    fprintf('WGTS:ComputeRF: RF, for ring = %d, qU = %d, rhoD=%.4g, Nis=%.0f, Bmax=%.2f T, Bwgts:%.2f T \n',obj.FPD_RingList(pixel),qu,obj.WGTS_CD_MolPerCm2,obj.NIS+1,obj.MACE_Bmax_T,obj.WGTS_B_T);
                case 'OFF'
                    if numel(obj.FPD_PixList)==1
                        fprintf('WGTS:ComputeRF: RF, pixel %.0f, qU = %d,  rhoD=%.4g, Nis=%.0f, Bmax=%.2f T, Bwgts:%.2f T \n',qu,obj.FPD_PixList,obj.WGTS_CD_MolPerCm2,obj.NIS+1,obj.MACE_Bmax_T,obj.WGTS_B_T);
                    else
                        %h=msgbox(sprintf('WGTS:ComputeRF: RF, <pixel list>, qU = %d,  rhoD=%.4g, Nis=%.0f, Bmax=%.2f T, Bwgts:%.2f T \n',qu,obj.WGTS_CD_MolPerCm2,obj.NIS+1,obj.MACE_Bmax_T,obj.WGTS_B_T),'icon','help');
                        %delete(h);
                    end
            end
            
            % check is (te-qU) is already contained in E
            Etmp = round(E_rf,5);
            TeqUtmp = round(te-qu,5);
            Index  = ismember(Etmp,TeqUtmp);
            if sum(Index)==numel(te)
                out = RF(Index);
            else
                out = interp1(E_rf,RF,te-qu);
            end
        end
        
        function out  = ComputeRFTEST(obj,te,qu,varargin)  
            p = inputParser;
            p.addParameter('RFbinFactor',0.05,@(x)isfloat(x)); % SPEED UP or Less Precise...
            p.parse(varargin{:});
            RFbinFactor         = p.Results.RFbinFactor;
            
            % Binning
            maxE = (obj.qUmax-obj.qUmin)*10; minE=-maxE; NbinE = (maxE-minE)/0.25;
            E = linspace(minE,maxE,NbinE); Estep = E(2) - E(1);
            
            %% Inelastic Scattering Probabilities
            w = load( '../..//inputs/WGTSMACE/WGTS_ISProb/IS_1.1e+17-molPercm2_8-NIS_4.230-Bmax_2.520-Bs.mat');
            ratioMasterSamak = [0.999906038165416   0.999952728590095   0.999967392247124   0.999983000052758   0.999998403728429   1.000013233420951   1.000027704218478   1.065797995861047];
            obj.is_Pv = w.Pis_m /100 .* ratioMasterSamak;
            
            %% Energy Loss Function
            ElossFunctions = obj.ComputeELossFunction;
            f1scatn = ElossFunctions{1,1};
            f2scatn = ElossFunctions{2,1};
            f3scatn = ElossFunctions{3,1};
            f4scatn = ElossFunctions{4,1};
            f5scatn = ElossFunctions{5,1};
            f6scatn = ElossFunctions{6,1};
            f7scatn = ElossFunctions{7,1};
            tmpis1 = obj.is_Pv(2)*f1scatn(E);
            tmpis2 = obj.is_Pv(3)*f2scatn(E);
            tmpis3 = obj.is_Pv(4)*f3scatn(E);
            tmpis4 = obj.is_Pv(5)*f4scatn(E);
            tmpis5 = obj.is_Pv(6)*f5scatn(E);
            tmpis6 = obj.is_Pv(7)*f6scatn(E);
            tmpis7 = obj.is_Pv(8)*f7scatn(E);
            scatterings     = tmpis1 + tmpis2 + tmpis3 + tmpis4 + tmpis5 + tmpis6 + tmpis7;
            obj.fscat       = @(e)interp1(E,scatterings,e);
            
            % Retreive Transmission Function
            TF  = @obj.ComputeMaceTF;
            
            % Build Response Function - Change Binning
            maxE = (obj.qUmax-obj.qUmin)*10; minE=-maxE; NbinE = (maxE-minE)/0.05;
            E = linspace(minE,maxE,NbinE); Estep = E(2) - E(1);
            RF = TF(qu+E,qu,'pixel',1)*obj.is_Pv(1) + ...
                conv(TF(qu+E,qu,'pixel',1),obj.fscat(E),'same').*Estep;
            out = interp1(E,RF,te-qu);
            fprintf('WGTS:ComputeRF: RF qU = %d \n',qu);
        end
        
        function out = LoadOrSaveRF(obj,mode)
            % RF name structure: samak RF + column density + scaterrings (NIS) +
            % B max (Bm) + B source (Bs) + B analyzing plane (Ba) +
            % Binning Factor (Bin) + ELossModel
            % example for Segmentation OFF: '../../inputs/ResponseFunction/samakRF_5e+17cm2_NIS11_Bm4.2T_Bs3.6T_Ba6G_5Bin.mat'
            % example for Segmentation MULTIPIXEL: '../../inputs/ResponseFunction/samakRF_5e+17cm2_NIS11_Bm4.2T_Bs3.6T_Ba6G_5Bin_multipix.mat'
            
            % output: out=1 if loading was successful , out=0 if response function is not avaiable
            %% labeling
            if ~strcmp(obj.KTFFlag,'WGTSMACE')
                out = 0;
                return;
            end
            RFpath = [getenv('SamakPath'),'inputs/ResponseFunction/samakRF'];
            if strcmp(obj.ISCS,'Edep')
                IsXStr = 'Edep';
            else
                IsXStr = sprintf('%.5gm2',obj.ISXsection(18574));
            end
            
            switch obj.FPD_Segmentation
                case 'OFF'
                    RFfile = sprintf('_Uniform_%.5gcm2_IsX%s_NIS%.0d_Bm%0.3gT_Bs%0.3gT_Ba%0.4gG_Temin%.3f_Temax%.3f_Bin%.0d_%s.mat',...
                        obj.WGTS_CD_MolPerCm2,IsXStr,obj.NIS,obj.MACE_Bmax_T,obj.WGTS_B_T,...
                        obj.MACE_Ba_T*1e4,obj.TeMin,obj.TeMax,obj.nTeBinningFactor,obj.ELossFlag);
                case {'MULTIPIXEL','SINGLEPIXEL'}
                    RFfile = sprintf('_MultiPix_%0.5gcm2_IsX%s_NIS%.0d_Bm%0.3gT_Bs%0.3gT_Ba%0.4gG_Temin%3.f_Temax%3.f_Bin%.0d_%s.mat',...
                        obj.WGTS_CD_MolPerCm2,IsXStr,obj.NIS,obj.MACE_Bmax_T,obj.WGTS_B_T,...
                        mean(obj.MACE_Ba_T(obj.FPD_PixList))*1e4,obj.TeMin,obj.TeMax,obj.nTeBinningFactor,obj.ELossFlag);
                case 'RING'
                    RFfile = sprintf('_Ring_%s_%0.5gcm2_IsX%s_NIS%.0d_Bm%0.3gT_Bs%0.3gT_Ba%0.4gG_Temin%.3f_Temax%.3f_Bin%.0d_nring%.0d.mat',...
                        obj.FPD_RingMerge,obj.WGTS_CD_MolPerCm2,IsXStr,obj.NIS,obj.MACE_Bmax_T,obj.WGTS_B_T,...
                        mean(obj.MACE_Ba_T)*1e4,obj.TeMin,obj.TeMax,obj.nTeBinningFactor,obj.nRings);
            end
            if obj.is_EOffset ~=0
                RFfile = strrep(RFfile,obj.ELossFlag,[obj.ELossFlag,sprintf('_%.3fELOffset',obj.is_EOffset)]);
            end
            
            if  obj.MACE_Sigma>0
                if all(diff(obj.MACE_Sigma))==0 || numel(obj.MACE_Sigma)==1
                    RFfile = strrep(RFfile,'.mat',sprintf('_Sigma%.0fmeV.mat',obj.MACE_Sigma(1)*1e3));
                else
                    RFfile = strrep(RFfile,'.mat',sprintf('_MeanSigma%.0fmeV.mat',mean(obj.MACE_Sigma*1e3)));
                end
            end
            MakeDir(RFpath);
            RFName = [RFpath,RFfile];
            
            %% load
            switch mode
                case 'load'
                    if ~exist(RFName,'file')
                        fprintf(2,'response function %s not available in data bank --> compute \n',RFfile);
                        out = 0;
                        return
                    elseif exist(RFName,'file')==2
                        RFst = load(RFName);
                    end
                    
                    % sanity check 1: RF for correct Ba and PixList
                    firstLoadingCondition = any(round(RFst.MACE_Ba_TallPixels_save,5,'significant') ~= round(obj.MACE_Ba_T,5,'significant')) || ...
                        isfield(RFst,'pixlist') && any(~ismember(obj.FPD_PixList,RFst.pixlist));
                    if firstLoadingCondition
                        out = 0;
                        return;
                    end
                    
                    %sanity check 2: Correct energy binning
                    switch obj.FPD_Segmentation
                        case {'OFF','RING'}
                            secondLoadingCondition = all(size(obj.Te,1) == size(RFst.Te,1)) && ...
                                all(round(RFst.Te,7,'significant') == round(obj.Te,7,'significant'),'all') ...
                                && all(round(RFst.qU,7,'significant') == round(obj.qU,7,'significant'),'all');
%                         case 'RING'
%                             secondLoadingCondition = all(size(obj.Te,1) == size(RFst.Te,1)) && ...
%                                 all(all(RFst.Te(:,obj.FPD_RingList) == obj.Te));
                        case {'SINGLEPIXEL','MULTIPIXEL'}
                            secondLoadingCondition = all(size(obj.Te,1) == size(RFst.Te,1)) && ...
                                all(all(RFst.Te(:,obj.FPD_PixList) == obj.Te));
                    end
                   if ~secondLoadingCondition
                       out = 0;
                       return;
                   end

                    % load response function
                    switch obj.FPD_Segmentation
                        case {'OFF','RING'}
                            obj.RF = RFst.RF;
                        case {'SINGLEPIXEL','MULTIPIXEL'}
                            for p = 1:obj.nPixels
                                for ii = 1:obj.nqU
                                    obj.RF(:,ii,p) = ...
                                        RFst.RF(:,ii,obj.FPD_PixList(p));
                                end
                            end
                    end
                    out = 1;
                    fprintf('loading response function from file: %s\n',RFName);
                    
                    %% save
                case 'save'
                    switch obj.FPD_Segmentation
                        case {'MULTIPIXEL','SINGLEPIXEL'}
                            RF = nan(obj.nTe,obj.nqU,148); Te = nan(obj.nTe,148);
                            RF(:,:,obj.FPD_PixList) = obj.RF; Te(:,obj.FPD_PixList) = obj.Te;
                        case {'OFF','RING'}
                            RF = obj.RF; Te = obj.Te; 
                    end
                    qU=obj.qU;
                    MACE_Ba_TallPixels_save = obj.MACE_Ba_T;
                    pixlist = obj.FPD_PixList;
                    save(RFName,'Te','RF','MACE_Ba_TallPixels_save','qU','pixlist','-v7.3','-nocompression');
                    out = 1;
            end % switch mode
end %LoadorSaveRF
        
        % Plot / Display
        function DisplayWGTSMACEInfo(obj,varargin)
            p=inputParser;
            p.addParameter('Output','screen',@(x)ismember(x,{'file','screen'}));
            p.addParameter('filename','',@(x)ischar(x));
            p.parse(varargin{:});
            Output   = p.Results.Output;
            filename = p.Results.filename;
            
            switch Output
                case 'file'
                    try
                        fileID = fopen(filename,'w+');
                    catch
                        fprintf(fileID,'file doesnt exist. enter valid filename!\n');
                        return
                    end
                case 'screen'
                    fileID = 1;
            end
            
            fprintf(fileID,'--------------------------------------------------------------\n');
            fprintf(fileID,'KATRIN WGTS and MACE with following properties:\n');
            fprintf(fileID,'--------------------------------------------------------------\n');
            fprintf(fileID,'- - - - - WGTS \n');
            fprintf(fileID,'   - Radius Flux Tube in cm: %g \n',obj.WGTS_FTR_cm);
            fprintf(fileID,'   - Column Density in molecules/cm2: %g \n',obj.WGTS_CD_MolPerCm2);
            fprintf(fileID,'   - Temparature in K: %g \n',obj.WGTS_Temp);
            fprintf(fileID,'   - Bs in Tesla: %g +- %g \n',obj.WGTS_B_T,obj.WGTS_ErrB_T);
            fprintf(fileID,'   - Cos Maximum Acceptance Angle %g \n',cos(asin(sqrt(obj.WGTS_B_T./obj.MACE_Bmax_T))) );
            
            fprintf(fileID,'- - - - - Isotropologues \n');
            fprintf(fileID,'   - T-T purity: %g \n',obj.WGTS_epsT);
            fprintf(fileID,'   - T-T fraction in WGTS: %g \n',obj.WGTS_MolFrac_TT);
            fprintf(fileID,'   - D-T fraction in WGTS: %g \n',obj.WGTS_MolFrac_DT);
            fprintf(fileID,'   - H-T fraction in WGTS: %g \n',obj.WGTS_MolFrac_HT);
            fprintf(fileID,'- - - - - Scattering Probabilities (from no-scattering) \n');
            if obj.is_Pv(1) == 0  % check if IS probability already computed
                file_pis = [getenv('SamakPath'),sprintf('/inputs/WGTSMACE/WGTS_ISProb/IS_%g-molPercm2_%.0f-NIS_%.3f-Bmax_%.3f-Bs.mat',...
                    obj.WGTS_CD_MolPerCm2,obj.NIS+1,obj.MACE_Bmax_T,obj.WGTS_B_T)];
                if exist(file_pis,'file') == 2 % load from file
                    w = load(file_pis);
                else  % compute
                    w.Pis_m = 0;
                end
                % Scattering Probabilities
                obj.is_Pv = w.Pis_m /100;
            end
            fprintf(fileID,'%s\n',num2str(obj.is_Pv));
            fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n');
            fprintf(fileID,'- - - - - MACE \n');
            fprintf(fileID,'   - B pinch in Tesla: %g +- %g \n',obj.MACE_Bmax_T,obj.ErrMACE_Bmax_T);
            fprintf(fileID,'   - B analysis plane in Tesla: %g +- %g \n',obj.MACE_Ba_T,obj.ErrMACE_Ba_T);
            fprintf(fileID,'- - - - - Analysis Plane \n');
            fprintf(fileID,'   - MACE_Ba_Setting: %s \n',obj.MACE_Ba_Setting);
            fprintf(fileID,'   - Analysis Plane B-field Correction, Tesla: \n');
            if ~isempty(obj.Pixel_MACE_Ba_TCorr)
                fprintf(fileID,'%s\n',num2str(obj.Pixel_MACE_Ba_TCorr(1:13)));
            end
            fprintf(fileID,'   - Analysis Plane E-field Correction, Volt: \n');
            if ~isempty(obj.Pixel_MACE_Ba_VCorr)
                fprintf(fileID,'%s\n',num2str(obj.Pixel_MACE_Ba_VCorr(1:13)));
            end
            fprintf(fileID,'- - - - - Miscellaneous \n');
            fprintf(fileID,'   - Integral Spectrum - Time Distribution: %s \n',obj.TD);
            fprintf(fileID,'   - Transmission Function: %s \n',obj.KTFFlag);
            fprintf(fileID,'--------------------------------------------------------------\n');
            switch Output
                case 'file'
                    fclose(fileID);
            end
        end
        function PlotISXsection(obj,varargin)
            p = inputParser;
            p.addParameter('Emin',100,@(x)isfloat(x) && x>0);
            p.parse(varargin{:});
            Emin     = p.Results.Emin;
            Energy = 18575+(-Emin:0.1:50);
            f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
            plot(Energy-18574,1e22.*obj.ISXsection(18574).*ones(numel(Energy),1),':','LineWidth',2,'Color',rgb('SlateGray'));
            hold on;
            plot(Energy-18574 ,1e22.*obj.ISXsection(Energy),'-','LineWidth',2,'Color',rgb('DodgerBlue'));
            xlabel('Energy - 18574 (eV)');
            ylabel('Inel. cross section (10^{-22} m^2)');
            PrettyFigureFormat('FontSize',24)
            xlim([min(Energy) max(Energy)]-18574);
        end
        function PlotKTF(obj,varargin)
            % Plot Transmission Function
            p = inputParser;
            p.addParameter('fign',999,@(x)isfloat(x) && x>0);
            p.addParameter('pub',0,@(x)isfloat(x) && x>=0);
            p.addParameter('ring',1,@(x)isfloat(x) && x>=0);
            p.parse(varargin{:});
            fign      = p.Results.fign;
            pub       = p.Results.pub;
            ring      = p.Results.ring;
            
            if isequal(obj.KTFFlag,'OFF')
                fprintf(2,'WGTSMACE/PlotKTF: KTFFlag Flag = OFF - exit..\n');
                return;
            end
            
            fign = figure('Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
            maintitle=sprintf('qU-dependent Response Function - \\rho.d = %.4g mol/cm^2 - Bs = %.2f T - Ba = %.3g T',...
                obj.WGTS_CD_MolPerCm2,obj.WGTS_B_T,obj.MACE_Ba_T(ring));  
            if numel(obj.MACE_Ba_T)>1
                maintitle = [maintitle,sprintf(' (pseudo ring %.0f)',ring)];
            end
            a=annotation('textbox', [0 0.9 1 0.1], ...
                'String', maintitle, ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center');
            a.FontSize=20;a.FontWeight='bold';
            h = plot(obj.Te-obj.Q_i,squeeze(obj.RF(:,:,ring)),'LineWidth',2);            
            set(gca, 'YScale', 'lin');
            xlabel('T_e-E_0 (eV)','FontSize',18);
            ylabel('Transmission','FontSize',18);
            grid on
            PrettyFigureFormat
            set(gca,'FontSize',18);
            if pub>0
                MakeDir('./plots');
                export_fig(fign,'./plots/WGTSMACE_KTF.pdf','-painters')
            end
        end
        
        function PlotPIS(obj,varargin)
            % Plot KATRIN WGTS Inelastic Probabilities
            p = inputParser;
            p.addParameter('fign',999,@(x)isfloat(x) && x>0);
            p.addParameter('pub',0,@(x)isfloat(x) && x>=0);
            p.parse(varargin{:});
            fign     = p.Results.fign;
            pub      = p.Results.pub;
            is = importdata([getenv('SamakPath'),'/inputs/WGTSMACE/ProbMatrix_Trials-10000_NIS-11.mat']);
            %Pis_k = is.Pis_mean;
            Pis_k = obj.is_Pv'*100; % KATRIN Nominal
            
            %%
            figure(fign);
            subplot(3,2,1);
            nhist(is.Pis_m(1,:),'noerror');  hold on;
            line([Pis_k(1) Pis_k(1)],[0 1e3],'Color','Red','LineWidth',2); hold off;
            xlabel('%'); ylabel('P_{0}');
            
            subplot(3,2,2);
            nhist(is.Pis_m(2,:),'noerror');  hold on;
            line([Pis_k(2) Pis_k(2)],[0 1e3],'Color','Red','LineWidth',2); hold off;
            xlabel('%'); ylabel('P_{1}');
            
            subplot(3,2,3);
            nhist(is.Pis_m(3,:),'noerror'); hold on;
            line([Pis_k(3) Pis_k(3)],[0 1e3],'Color','Red','LineWidth',2); hold off;
            xlabel('%'); ylabel('P_{2}');
            
            subplot(3,2,4);
            nhist(is.Pis_m(4,:),'noerror'); hold on;
            line([Pis_k(4) Pis_k(4)],[0 1e3],'Color','Red','LineWidth',2); hold off;
            xlabel('%'); ylabel('P_{3}');
            
            subplot(3,2,5);
            nhist(is.Pis_m(5,:),'noerror'); hold on;
            line([Pis_k(5) Pis_k(5)],[0 1e3],'Color','Red','LineWidth',2); hold off;
            xlabel('%'); ylabel('P_{4}');
            
            subplot(3,2,6);
            nhist(is.Pis_m(6,:),'noerror'); hold on;
            line([Pis_k(6) Pis_k(6)],[0 1e3],'Color','Red','LineWidth',2); hold off;
            strtitle = sprintf('WGTS Inelastic Scattering Probabilites - %0.f Trials',numel(is.Pis_m(1,:)));
            set(gcf,'NextPlot','add'); axes; h = title(strtitle); set(gca,'Visible','off'); set(h,'Visible','on');
            xlabel('%'); ylabel('P_{5}');
            
            if pub>0
                publish_figure(fign,'./figures/WGTSMACE_PIS_0to5.eps')
            end
            
            %%
            figure(fign+1);
            subplot(3,2,1);
            nhist(is.Pis_m(7,:),'noerror'); hold on;
            line([Pis_k(7) Pis_k(7)],[0 1e3],'Color','Red','LineWidth',2); hold off;
            xlabel('%'); ylabel('P_{6}');
            
            subplot(3,2,2);
            nhist(is.Pis_m(8,:),'noerror'); hold on;
            line([Pis_k(8) Pis_k(8)],[0 1e3],'Color','Red','LineWidth',2); hold off;
            xlabel('%'); ylabel('P_{7}');
            
            subplot(3,2,3);
            nhist(is.Pis_m(9,:),'noerror'); hold on;
            line([Pis_k(9) Pis_k(9)],[0 1e3],'Color','Red','LineWidth',2); hold off;
            xlabel('%'); ylabel('P_{8}');
            
            subplot(3,2,4);
            nhist(is.Pis_m(10,:),'noerror'); hold on;
            line([Pis_k(10) Pis_k(10)],[0 1e3],'Color','Red','LineWidth',2); hold off;
            xlabel('%'); ylabel('P_{9}');
            
            subplot(3,2,5);
            nhist(is.Pis_m(11,:),'noerror'); hold on;
            line([Pis_k(11) Pis_k(11)],[0 1e3],'Color','Red','LineWidth',2); hold off;
            xlabel('%'); ylabel('P_{10}');
            
            strtitle = sprintf('WGTS Inelastic Scattering Probabilites - %0.f Trials',numel(is.Pis_m(1,:)));
            set(gcf,'NextPlot','add'); axes; h = title(strtitle); set(gca,'Visible','off'); set(h,'Visible','on');
            if pub>0
                publish_figure(fign+1,'./figures/WGTSMACE_PIS_6to10.eps')
            end
            
            %%
            figure(fign+2);
            j=1:1:11;
            subplot(2,1,1);
            h1 = errorbar(j-1,is.Pis_mean,sqrt(diag(is.PisCM)),'Marker','s','Color','Black'); hold on;
            %errorbar_tick(h1,200);
            h2 = plot(j-1,Pis_k(j),'LineStyle','--','LineWidth',1,'Color','Red'); hold off
            legend([h1 h2],'Covariance Matrix','KATRIN Nominal'); legend('boxoff','FontSize',12,'Location','northEast');
            set(gca,'Yscale','log'); grid on;
            ylabel('Scaterring Probability (%)');
            subplot(2,1,2)
            j=1:1:11;
            h3 = errorbar(j-1,Pis_k(j)./is.Pis_mean,sqrt(diag(is.PisCM))./is.Pis_mean,'Marker','s','Color','Red');
            %errorbar_tick(h1,200);
            legend(h3,'SSC Nominal / Thierry Covariance Computation'); legend('boxoff','FontSize',10,'Location','northEast');
            set(gca,'Yscale','lin'); grid on;
            xlabel('Number of Inelastic Scattering'); ylabel('Ratio');
            strtitle='Inelastic Probabilites in KATRIN WGTS';
            set(gcf,'NextPlot','add'); axes; h = title(strtitle); set(gca,'Visible','off'); set(h,'Visible','on');
            if pub>0
                publish_figure(fign+2,'./figures/WGTSMACE_PIS.eps')
            end
            
            figure(fign+3);
            imagesc((is.PisCM));
            colormap(hot);
            set(gca,'YDir','Normal');
            colorbar;
            title('KATRIN WGTS : Inelastic Scattering Covariance Matrix');
            if pub>0
                publish_figure(fign+3,'./figures/WGTSMACE_PIS_CovarianceMatrix.eps');
            end
        end   
        
        % --------------------------------------------------------------- %
        % Begin: KATRIN Energy Loss Function Parametrization
        function g = gaussian(~,x, pos, sigma)
            % KATRIN Energy Loss Parametrization: Gaussian Function
            g = exp( -0.5 .* ((x-pos)./sigma).^2 );
        end
        
        function f = dfdwH2(~,y)
            % KATRIN Energy Loss Parametrization: Differential dipole oscillator strength for H2; 
            % table I in  Kim, Rudd, Phys. Rev. A 50, 3954 (1994).
            c = 1.1262;
            d = 6.3982;
            e = -7.8055;
            f = 2.1440;
            f =  c .* y.^3 + d .* y.^4 + e .* y.^5 + f .* y.^6;            
        end
        
        function dsdW = Dsigma(obj,E, T, B)
            % KATRIN Energy Loss Parametrization: Unnormalized single differential ionization cross section.
            % E = electron energy loss [eV]
            % Ionization energy B for H2, D2 and T2 = 15.426 eV, 15.467 eV, 15.487 eV
            % Minimal energy loss: Elossmin = B
            % Maximal energy loss: Elossmax = (T - B)/2
            
            % Binary encounter dipole model, eq. 44 in Kim, Rudd, Phys. Rev. A 50, 3954 (1994).
            W  = E - B;   % secondary electron energy [eV]
            w  = W./B;
            t  = T./B;
            y  = 1./(w+1); % = B/E for expansion of oscillator strength
            N  = 2.;      % number of electrons in H2 molecule
            Ni = 1.173;   % eq. 39 and table I.
            
            dsdW = (Ni./N - 2)./(t+1) .* ( y + 1./(t-w) ) +...
                (2 - Ni./N) .* (y.^2 + 1./(t-w).^2) + ...
                log(t).*y./N .* obj.dfdwH2(y);
        end
        
        function f = en_loss_sc(obj,x, amp1, pos1, sig1, amp2, pos2, sig2, amp3, pos3, sig3)
            % KATRIN Energy Loss Parametrization: parametrized energy loss function, V. Hannen, March 2019, KATRIN CM36
            % Global fit of parametrized energy loss function to measured integral and t.o.f. data
            % Based on D2 scattering data

            T  = 18575;   % kinetic energy of incident electron [eV]
            
            f = (x <= 0)   .* 0  + ...
                (x>0 & x <  obj.Ei)  .*  (amp1 .* obj.gaussian(x, pos1, sig1)   + amp2 .* obj.gaussian(x, pos2, sig2)  + amp3 .* obj.gaussian(x, pos3, sig3)) + ...
                (x >= obj.Ei)  .* ((amp1 .* obj.gaussian(obj.Ei, pos1, sig1)  + amp2 .* obj.gaussian(obj.Ei, pos2, sig2) + amp3 .* obj.gaussian(obj.Ei, pos3, sig3)) .* obj.Dsigma(x, T, obj.Ei) ./ obj.Dsigma(obj.Ei, T, obj.Ei));
            f(isnan(f)) = 0;
        end
        
        function f = eloss_katrin(obj,e, amp1, pos1, sig1, amp2, pos2, sig2, amp3, pos3, sig3)
            % KATRIN Energy Loss Parametrization: parametrized energy loss function, V. Hannen, March 2019, KATRIN CM36
            % Global fit of parametrized energy loss function to measured integral and t.o.f. data
            % Based on D2 scattering data
            % Same as en_loss_sc but compactified
          
            
            w  = @(e) (e - obj.Ei)./obj.Ei;
            y  = @(e) 1./(w(e)+1);
            gaussian = @(x,pos,sigma) exp( -0.5 .* ((x-pos)./sigma).^2 );
            dfdwH2   = @(y) 1.1262 .* y.^3 + 6.3982 .* y.^4  -7.8055 .* y.^5 + 2.1440 .* y.^6;
            
            N = 2;      % numver of electron in T2 molecule
            Ni = 1.173; % eq. 39 and table I. (Volker's talk)
            E0 = 18575;

            Dsigma   = @(y,w) (Ni/N-2) ./ (E0/obj.Ei+1) .* ( y + 1./(E0/obj.Ei-w) ) + ...
                (2-Ni/N) .* (y.^2 + 1./(E0/obj.Ei-w).^2) + log(E0/obj.Ei)/2 .*y.* dfdwH2(y);
            
            %Dsigma   = @(y,w) (Ni/N-2) ./ (1.200943945173595e+03+1) .* ( y + 1./(1.200943945173595e+03-w) ) + ...
            %                  (2-Ni/N) .* (y.^2 + 1./(1.200943945173595e+03-w).^2) + 3.545431573764978.*y.* dfdwH2(y);
            
            f = (e <=  0) .* 0 + ...
                (e <  obj.Ei) .*   (amp1 .* gaussian(e , pos1, sig1)  + amp2 .* gaussian(e, pos2, sig2) + amp3 .* gaussian(e , pos3, sig3)) + ...
                (e >= obj.Ei) .*  ((amp1 .* gaussian(obj.Ei, pos1, sig1)  + ... 
                                    amp2 .* gaussian(obj.Ei, pos2, sig2) + ...
                                    amp3 .* gaussian(obj.Ei, pos3, sig3))...
                                    .* Dsigma(y(e),w(e)) ./ Dsigma(1,0));
            f(isnan(f)) = 0;
        end
        
        % End: KATRIN Energy Loss Function Parametrization
        % --------------------------------------------------------------- %
    end
  
end % class
