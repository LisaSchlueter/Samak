% ----------------------------------------------------------------------- %
% Fitting class for the KATRIN experiment
% ----------------------------------------------------------------------- %
% This class contains options to fit the data produced by the KATRIN
% experiment, as well as montecarlo data.
%
% Chi-squared functions included:
%   - chi2
%   - chi2 with covariance matrix
%   - chi2 Poissonian deviance
%
% Last update by Lisa: 04/07/2018 (maybe not accurate).
%
% T. Lasserre
% CEA - Saclay
%
% L. Schlueter
% TUM / MPP
%
% P. I. Morales Guzman
% TUM / MPP
%
%#ok<*ST2NM>
%
% ----------------------------------------------------------------------- %

classdef FITC < handle
    
    properties (Constant = true, Hidden = true)
        
    end
    
    properties (Access=protected)
        
    end
    
    properties (Dependent=true,Access=public)
        
    end
    
    properties (Access=public)
        
        %TODO: explain (comment) all properties
        
        % Main inputs (necessary)
        SO;             % Study Object
        SO2;            % 2nd Study Object for knm1+2 combi fit 
        SOCell;         % Study Objects in a cell (for use with the rhoD handle study)
        DATA;           % Data for fit, it could be a cell with three elemnts {qU,counts,errors} or a matrix 3*pixels columns [qU,counts,errors]
        
        % Main inputs (optional)
        COVMAT;         % Covariance Matrix for fit with the covariance matrix (optional)
        COVMATFrac;     % Fractional Covariance Matrix
        COVMATShape;    % Shape only part of Covariance Matrix
        COVMATFracShape;
        COVMATNorm;     % Normalization part of Covariance Matrix
        chi2name;       % chi2 function to be used for minimization
        fitter;         % fitter to use (matlab: single and multipixel; minuit: only single pixel)
        
        % options for fitter (optional)
        
        % options for analysis
        pixels;
        exclDataStart;  % remove data points from fit: exclDataStart=5, exclDataStop=15
        exclDataStop ;  %  -> used Data points [5:15]
        % options for
        
        pulls;
        pullFlag;      % 1 = nu mass pull, 2 = nu mass + FSD pull, 3 = FSD pull
        fixPar;        % fixed parameters
        nfixPar;       % number of fixed parameters
        fixPar_v;    % fixed parameters as vector and not string
        minuitOpt;     % minuit options
        parinitoriginal; % (not an input) temporal variable for the matlab minimizer
        
        % initialization of fitted paramaters
        i_mnu;
        i_Q;
        i_Q2; % second endpoint for knm1+2 combi fit
        i_B;
        i_N;
        i_BSlope;   %background qU slope
        i_BPTSlope; %background subrun time slope from interspec penning trap
        
        i_DTGS; % DT Ground  State Probability
        i_DTES; % DT Excited State Probability
        i_HTGS; % HT Ground  State Probability
        i_HTES; % HT Excited State Probability
        i_TTGS; % TT Ground  State Probability
        i_TTES; % TT Excited State Probability
        
        i_qUOffset; %retarding potential offset per ring or pixel
        i_mTSq;     % neutrino mass squared offset per ring or pixel
        i_FracTm;   % fraction T- ions

        i_mnu4Sq;
        i_sin2T4;

        parinit; % initial parameters for the fit
        
        % data
        qUdata; % qU bin values
        counts; % counts per qU
        c_error; % error on the counts (normally sqrt(counts))
        
        % option for autofit
        autofit; % autofit flag. Default behavior: 'ON' perform the fit when calling the class
        readdatadone = false;
        preparevarsdone = false;
        
        % Multipixel Fit Uncertainties
        UncMPixFlag; % Flag for calculating a better estimate for the uncertainties in the multipixel fit
        
        % output
        RESULTS;
        
        % CATS: Statistic Diagnostics
        DesignMatrix;                    % Design Matrix
        MaskFreeFitPar;                  % Masking the actual FREE fit parameters (no pulls)
        MaskFreeFitLabel;                % Labels of CATS free parameters       
        CATSstat;                        % CATS Diagnostic Sturcture
    end
    
    methods
        function obj = FITC(varargin)
            p = inputParser;
            % Main inputs (necessary)
            p.addParameter('SO',[],@(x) isa(x,'TBD') || isa(x,'RHOD')); %Study Object
            p.addParameter('SO2',[],@(x) isa(x,'TBD')); %2nd Study Object (optional! for knm1+2 test fit)
            p.addParameter('SOCell',{},@(x) iscell(x)); %Study Objects in cell
            p.addParameter('DATA',[],@(x) (iscell(x) || ismatrix(x)) && (size(x,2) <= 148*3));
            % Main inputs (optional)
            p.addParameter('COVMAT',[],@(x) ismatrix(x));
            p.addParameter('COVMATFrac',[],@(x) ismatrix(x));
            p.addParameter('COVMATShape',[],@(x) ismatrix(x));
            p.addParameter('COVMATFracShape',[],@(x) ismatrix(x));
            p.addParameter('COVMATNorm',[],@(x) ismatrix(x));
            p.addParameter('chi2name','chi2Stat', @(x)ismember(x,{'chi2Stat',...
                'chi2CM','chi2CMFrac','chi2P','chi2CMShape','chi2Nfix'}));
            p.addParameter('fitter', 'matlab' , @(x)ismember(x,{'matlab', 'minuit'}));
            p.addParameter('UncMPixFlag',true, @(x) islogical(x));
            
            
            % options for fitters (optional)
            % TO DO: add options for matlab and minuit one by one
            
            % options for analysis
            p.addParameter('pixels',[],@(x)isfloat(x) && x>0);
            p.addParameter('exclDataStart',1,@(x)isfloat(x)  && all(x>0));
            p.addParameter('exclDataStop',9999,@(x)isfloat(x) && x>0);
            % options for chi2
            p.addParameter('pulls',0,@(x)isfloat(x) && all(x>0));
            p.addParameter('pullFlag',1,@(x) isfloat(x));
            p.addParameter('fixPar','',@(x)ischar(x));
            p.addParameter('minuitOpt','min;migrad',@(x)ischar(x));
            % option for autofit
            p.addParameter('autofit','ON', @(x)ismember(x,{'ON','OFF'}));
            
            % initialization of fitted parameters
            p.addParameter('i_mnu',0,@(x)isfloat(x));
            p.addParameter('i_Q',0,@(x)isfloat(x));
             p.addParameter('i_Q2',0,@(x)isfloat(x));
            p.addParameter('i_DTGS',0,@(x)isfloat(x));
            p.addParameter('i_DTES',0,@(x)isfloat(x));
            p.addParameter('i_HTGS',0,@(x)isfloat(x));
            p.addParameter('i_HTES',0,@(x)isfloat(x));
            p.addParameter('i_TTGS',0,@(x)isfloat(x));
            p.addParameter('i_TTES',0,@(x)isfloat(x));
            p.addParameter('i_B',[],@(x)isfloat(x));
            p.addParameter('i_N',[],@(x)isfloat(x));
            p.addParameter('i_qUOffset',[],@(x)isfloat(x));
            p.addParameter('i_BSlope',0,@(x)isfloat(x));
            p.addParameter('i_BPTSlope',0,@(x)isfloat(x));
            p.addParameter('i_mTSq',[],@(x)isfloat(x));
            p.addParameter('i_FracTm',[],@(x)isfloat(x));
            p.addParameter('i_mnu4Sq',0,@(x)isfloat(x));
            p.addParameter('i_sin2T4',0,@(x)isfloat(x));
            
            p.parse(varargin{:});
            
            obj.SO              = p.Results.SO;
            obj.SO2             = p.Results.SO2;
            obj.DATA            = p.Results.DATA;
            obj.COVMAT          = p.Results.COVMAT;
            obj.COVMATFrac      = p.Results.COVMATFrac;
            obj.COVMATShape     = p.Results.COVMATShape;
            obj.COVMATFracShape = p.Results.COVMATFracShape;
            obj.COVMATNorm      = p.Results.COVMATNorm;
            obj.chi2name        = p.Results.chi2name;
            obj.fitter          = p.Results.fitter;
            
            obj.pixels          = p.Results.pixels;
            obj.exclDataStart   = p.Results.exclDataStart;
            obj.exclDataStop    = p.Results.exclDataStop;
            
            obj.pulls           = p.Results.pulls;
            obj.pullFlag        = p.Results.pullFlag;
            obj.fixPar          = p.Results.fixPar;
            obj.minuitOpt       = p.Results.minuitOpt;
            
            obj.autofit         = p.Results.autofit;
            
            obj.i_mnu           = p.Results.i_mnu;
            obj.i_Q             = p.Results.i_Q;
            obj.i_Q2             = p.Results.i_Q;
            obj.i_B             = p.Results.i_B;
            obj.i_N             = p.Results.i_N;
            obj.i_DTGS          = p.Results.i_DTGS;
            obj.i_DTES          = p.Results.i_DTES;
            obj.i_HTGS          = p.Results.i_HTGS;
            obj.i_HTES          = p.Results.i_HTES;
            obj.i_TTGS          = p.Results.i_TTGS;
            obj.i_TTES          = p.Results.i_TTES;
            obj.i_qUOffset      = p.Results.i_qUOffset;
            obj.i_mTSq          = p.Results.i_mTSq;
            obj.i_BSlope        = p.Results.i_BSlope;
            obj.i_BPTSlope      = p.Results.i_BPTSlope;
            obj.i_FracTm        = p.Results.i_FracTm;
            obj.UncMPixFlag     = p.Results.UncMPixFlag;
            obj.i_mnu4Sq        = p.Results.i_mnu4Sq;
            obj.i_sin2T4        = p.Results.i_sin2T4;
            
            if isempty(obj.SO) && isempty(obj.SOCell)
                error('You should provide a TBD object to use this class. Ex.: FIT(''SO'',TBDObject)');
            end
            
            if isempty(obj.DATA)
                error('You should provide qU and data with uncertainties to use this class (3 cells). Ex.: FIT(''DATA'',{qU(:),counts(:),errors(:)})');
            end
            
            if obj.exclDataStop == 9999 % if no boundary was set -> use all data points
                obj.exclDataStop = obj.SO.nqU;
            end
            
            obj.nfixPar = numel(str2num(strrep(strrep(obj.fixPar,'fix ',''),' ;',' ')));
            obj.fixPar_v = str2num(strrep(strrep(obj.fixPar,'fix ',''),' ;',' '));
            
            if strcmp(obj.autofit,'ON')
                readdata(obj);
                preparevars(obj);
                obj.RESULTS = startfit(obj);
                
                if isempty(obj.SO2)
                    displayresults(obj);
                end
            end
            
        end % constructor
    end % methods
    
    methods
        
        function readdata(obj)
            
            if iscell(obj.DATA)
                qU_tmp = obj.DATA{1};
                counts_tmp = obj.DATA{2};
                c_error_tmp = obj.DATA{3};
            else
                cut = size(obj.DATA,2)/3;
                qU_tmp = obj.DATA(:,1:cut);
                counts_tmp = obj.DATA(:,(cut+1):2*cut);
                c_error_tmp = obj.DATA(:,(2*cut+1):end);
            end
   
            if obj.SO.nPixels == size(counts_tmp,2)
                obj.qUdata = qU_tmp;
                obj.counts = counts_tmp;
                obj.c_error = c_error_tmp;
            else
               obj.qUdata = qU_tmp(:,obj.SO.FPD_Pixel);
               obj.counts = counts_tmp(:,obj.SO.FPD_Pixel);
               obj.c_error = c_error_tmp(:,obj.SO.FPD_Pixel);
%                 obj.qUdata = qU_tmp(:,obj.SO.nPixels);
%                 obj.counts = counts_tmp(:,obj.SO.nPixels);
%                 obj.c_error = c_error_tmp(:,obj.SO.nPixels);
            end
            
            obj.readdatadone = true;
        end
        
        function preparevars(obj)
            if obj.readdatadone
                
                % initiliazation of parameters for fit
                % i_mnu and i_Q (and i_DTES + i_DTGS) are already initialized to 0
                
                % background
                if any(ismember(obj.fixPar_v,3)) %background fixed
                    obj.i_B = 0;
                else
                    if isempty(obj.i_B)
                        obj.i_B = obj.counts(end,:)./...
                            (obj.SO.qUfrac(end,:).*obj.SO.TimeSec) - ...
                            obj.SO.BKG_RateSec_i;
                        bckzeroindex = ((obj.i_B + obj.SO.BKG_RateSec_i) == zeros(1,length(obj.counts(end,:))));
                        if any(bckzeroindex)
                            obj.i_B(bckzeroindex) = obj.i_B(bckzeroindex) + 1e-3;
                        end
                    end
                end
                %normalization
                if isempty(obj.i_N)
                    if contains(obj.fixPar,'4')
                        obj.i_N = zeros(1,length(obj.counts(end,:)));
                    else
                        obj.i_N = (sum(obj.counts) - obj.counts(end,:)*obj.SO.nqU)./...
                            (sum(obj.SO.TBDIS) - obj.SO.TBDIS(end,:)*obj.SO.nqU)- 1;
                    end
                end
                
                obj.parinit = [obj.i_mnu,obj.i_Q,obj.i_B,obj.i_N,...
                               obj.i_DTGS,obj.i_DTES,...
                               obj.i_HTGS,obj.i_HTES,...
                               obj.i_TTGS,obj.i_TTES,...
                               obj.i_qUOffset,...
                               obj.i_BSlope,...
                               obj.i_mTSq,...
                               obj.i_FracTm,...
                               obj.i_mnu4Sq,...
                               obj.i_sin2T4,...
                               obj.i_BPTSlope];
                           
                           if ~isempty(obj.SO2)
                               obj.i_N = [0,0];
                               %combi fit with knm1+2
                               obj.parinit = [obj.i_mnu,...
                                   obj.i_Q,...
                                   obj.i_Q2,...
                                   obj.i_B,...
                                   obj.i_N];
                           end
                % initilization for pulls (default all zero)
                if obj.pulls == 0
                    obj.pulls = inf(1,length(obj.parinit));
                end
                
                obj.preparevarsdone = true;
            else
                obj.readdata;
                obj.preparevars;
            end
        end
        
        function PullTerm = ComputePulls(obj,par)
            % PULLS TOLERANCE:
            FSDtol=1e-03;         % final state distribution PGS-PES
            Normtol=0.0017;        % normalization tolerance from ring to ring
            NormAbstol =  2;      % asbolute ringwise normalization tolerance (deviation from 0)
            qUOffsettol = 10;     % qUOffset absolute tolerance (eV)
            qUmeanOffsettol = 0.5;  % qUOffset tolerance from ring to ring (eV)
            BkgSlopetol = 4.74.*1e-06;%5*1e-6; % background slope constrain (cps/eV)
            BkgPTSlopetol24 = 3.*1e-06; % background time slope constrain (cps/s)
            BkgPTSlopetol25 = 2.5.*1e-06; % background time slope constrain (cps/s)
            mNuSqtol = 1;         % neutrino mass pull tolerance (eV^2)
            mTSqtol = 1;          % tachyonic neutrino mass (nu-mass offset) from ring to ring (eV^2)
            sin2T4tol = 0.5;
           
            
            PullTerm = 0;
            % 1: neutrino mass pull
            if any(ismember(obj.pullFlag,1))
                PullTerm = PullTerm + (par(1)-0)^2/mNuSqtol^2;
                %sum( ((par-zeros(1,length(par))).^2)./obj.pulls);
            end
            
            % 2: mTSq (neutrino mass offsets): small deviatin from mean
            if any(ismember(obj.pullFlag,2))
                ParmTSq = par(3*obj.SO.nPixels+11:4*obj.SO.nPixels+9);
                PullTerm =  PullTerm + ( (sum(ParmTSq)-numel(ParmTSq)*mean(ParmTSq)) /mTSqtol )^2;
            end
            
            % 3: FSD pulls: sum(GS+ES) should be close to init and GS shoud not vary by more than 1%
            if any(ismember(obj.pullFlag,3)) % FSD pull
                PullTerm = PullTerm + ...
                    (par(2*obj.SO.nPixels+3)+par(2*obj.SO.nPixels+4))^2/FSDtol^2 + ... % DT
                    (par(2*obj.SO.nPixels+5)+par(2*obj.SO.nPixels+6))^2/FSDtol^2 + ... % HT
                    (par(2*obj.SO.nPixels+7)+par(2*obj.SO.nPixels+8))^2/FSDtol^2 + ... % T2
                    par(2*obj.SO.nPixels+3)^2/(0.01^2)+ ...                % 1% on the GS DT probability
                    par(2*obj.SO.nPixels+5)^2/(0.01^2)+...                 % 1% on the GS HT probability
                    par(2*obj.SO.nPixels+7)^2/(0.01^2);                    % 1% on the GS T2 probability     
            end
            
            %4: Normalization: small deviation from mean normalization
            if any(ismember(obj.pullFlag,4))
                ParNorm = par(3+obj.SO.nPixels:3+2*obj.SO.nPixels-1); % -> ringwise normalizations
                meanNorm = repmat(mean(ParNorm),size(ParNorm));
                PullTerm = PullTerm + sum((ParNorm-meanNorm).^2./Normtol^2);
            end
            
            % 5: qU Offsets: small deviatin from mean qU Offset
            if any(ismember(obj.pullFlag,5))
                ParqUOffset = par(2*obj.SO.nPixels+9:3*obj.SO.nPixels+8);
                meanqU = repmat(mean(ParqUOffset),size(ParqUOffset));
                PullTerm =  PullTerm + sum((ParqUOffset-meanqU).^2./qUmeanOffsettol.^2);
                %sum((ParqUOffset-zeros(1,length(ParqUOffset))).^2  ./qUOffsettol^2);
            end
            
            % 6: qU Offsets: small absolute qU Offset
            if any(ismember(obj.pullFlag,6))
                ParqUOffset = par(2*obj.SO.nPixels+9:3*obj.SO.nPixels+8);
                PullTerm =  PullTerm + sum(ParqUOffset.^2./qUOffsettol.^2);
            end
            
            % 7: background slope constrain
            if any(ismember(obj.pullFlag,7))
                PullTerm = PullTerm + par(3*obj.SO.nPixels+9)^2/BkgSlopetol^2;
            end
            
            % 8: qU Offsets: follows U_RW .* (1 -
            % (1-offset)*(r./WGTS_Radius_cm).^2) dependence
            % WORKS ONLY FOR 4 RINGS
            if any(ismember(obj.pullFlag,8))
                % ParqUOffset = par(2*obj.SO.nPixels+9:3*obj.SO.nPixels+8);
                % o1index=par(2*obj.SO.nPixels+9);
                o2index=par(2*obj.SO.nPixels+10); %ring2
                o3index=par(2*obj.SO.nPixels+11); %ring3
                o4index=par(2*obj.SO.nPixels+12); %ring4
                PullTerm =  PullTerm + ...
                    (o4index/0.76-o3index/0.39)^2/Normtol.^2 + ...
                    (o4index/0.76-o2index/0.14)^2/Normtol.^2;
            end

            % 9: mixing angle pull
            if any(ismember(obj.pullFlag,9)) % sterile pull I 
            range =  ceil(abs(obj.SO.qU(min(obj.exclDataStart))-18574));
                    PullTerm = PullTerm + ...
                    exp(1e3*(par(4*obj.SO.nPixels+12)-(0.5+1e-02))) + ...  % sin2t4 upper bound 
                    exp(-1e3*(par(4*obj.SO.nPixels+12)+1e-02)) + ...       % sin2t4 lower bound 
                    exp(1e3*(par(4*obj.SO.nPixels+11)-range^2)) + ...      % mSq4 lower bound
                    exp(-1e3*(par(4*obj.SO.nPixels+11)+1e-02)) + ...       % mSq4 upper bound
                    (par(1)-0)^2/1^2;                                      % nu-mass
            end
            
            if any(ismember(obj.pullFlag,10))   %  background slope constrain with variable sigma
                PullTerm = PullTerm + par(3*obj.SO.nPixels+9)^2/obj.pulls.^2;
            end
            
            if any(ismember(obj.pullFlag,11)) % gaussian sterile pull II
                sin4Central = 0; sin4Tol = 10;
                m4Central = 0;   m4Tol   = 1e4;
                PullTerm = PullTerm + ...
                    (par(4*obj.SO.nPixels+12)-sin4Central).^2./sin4Tol^2+...  % sin2t4
                    (par(4*obj.SO.nPixels+11)-m4Central).^2./m4Tol^2;         % m4
            end
            
            
            if any(ismember(obj.pullFlag,12)) % nu-mass pull: mainz&troitzk
                mNuSqtol_MT =  1.94; % eV^2
                PullTerm = PullTerm + (par(1)-0)^2/mNuSqtol_MT^2;
            end
            
            if  any(ismember(obj.pullFlag,13)) % E0 pull: 1 eV
                E0tol1eV =  1; % eV
                PullTerm = PullTerm + (par(2)-0)^2/E0tol1eV^2;
            end
            
            if  any(ismember(obj.pullFlag,14)) % E0 pull: 2 eV
                E0tol2eV =  2; % eV
                PullTerm = PullTerm + (par(2)-0)^2/E0tol2eV^2;
            end
            
            if any(ismember(obj.pullFlag,15)) % nu-mass pull: 1 eV^2
                mNuSqtol =  1; % eV^2
                PullTerm = PullTerm + (par(1)-0)^2/mNuSqtol^2;
            end
            
            if any(ismember(obj.pullFlag,16)) % nu-mass pull: 2 eV^2
                mNuSqtol =  2; % eV^2
                PullTerm = PullTerm + (par(1)-0)^2/mNuSqtol^2;
            end
            
            if any(ismember(obj.pullFlag,17)) % nu-mass pull: 3 eV^2
                mNuSqtol =  3; % eV^2
                PullTerm = PullTerm + (par(1)-0)^2/mNuSqtol^2;
            end
            
            if any(ismember(obj.pullFlag,18)) % nu-mass pull: 0.5 eV^2
                mNuSqtol =  0.5; % eV^2
                PullTerm = PullTerm + (par(1)-0)^2/mNuSqtol^2;
            end
            
            if any(ismember(obj.pullFlag,19)) % non-sense sterile pull
                sin4Up = 0.125;
                m4Up   = 10;
                PullTerm = PullTerm + ...
                    exp(1e2*(par(4*obj.SO.nPixels+12)-(sin4Up+1e-02))) + ...  % sin2t4 upper bound
                    exp(-1e2*(par(4*obj.SO.nPixels+12)+1e-02)) + ...          % sin2t4 lower bound -> 0
                    exp(1e2*(par(4*obj.SO.nPixels+11)-m4Up-0.1)) + ...        % mSq4 lower bound
                    exp(-1e2*(par(4*obj.SO.nPixels+11)+0.1));                % mSq4 upper bound
            end
            
            if any(ismember(obj.pullFlag,20)) % non-sense sterile pull II
                sin4Central = 0; sin4Tol = 0.125;
                m4Central = 0;   m4Tol   = 10;
                PullTerm = PullTerm + ...
                    (par(4*obj.SO.nPixels+12)-sin4Central).^2./sin4Tol^2+...  % sin2t4
                    (par(4*obj.SO.nPixels+11)-m4Central).^2./m4Tol^2;         % m4
            end
            
            if any(ismember(obj.pullFlag,21)) % Neutrino 4 sterile pull
                sin4Central = 0.38; sin4Tol = 0.11;
                m4Central   = 7.26; m4Tol   = 0.7;
                PullTerm = PullTerm + ...
                    (par(4*obj.SO.nPixels+12)-sin4Central).^2./sin4Tol^2+...  % sin2t4
                    (par(4*obj.SO.nPixels+11)-m4Central).^2./m4Tol^2;         % m4
            end
            
            %22: Normalization: small abs deviation from 0
            if any(ismember(obj.pullFlag,22))
                ParNorm = par(3+obj.SO.nPixels:3+2*obj.SO.nPixels-1); % -> ringwise normalizations
                PullTerm = PullTerm + sum((ParNorm).^2./NormAbstol^2);
            end
            
            if any(ismember(obj.pullFlag,23)) % loose nu-mass pull
                mNuSqtol_Loose =  100; % eV^2
                PullTerm = PullTerm + (par(1)-0)^2/mNuSqtol_Loose^2;
            end
            
            if any(ismember(obj.pullFlag,24)) % background subrun time slope from penning trap pull
                PullTerm = PullTerm + (par(4*obj.SO.nPixels+13)-0)^2/BkgPTSlopetol24^2;
            end
            
            if any(ismember(obj.pullFlag,25)) % background subrun time slope from penning trap pull
                PullTerm = PullTerm + (par(4*obj.SO.nPixels+13)-0)^2/BkgPTSlopetol25^2;
            end
            
             if any(ismember(obj.pullFlag,26)) % nu-mass pull: mainz&troitzk
                mNuSqtol_KNM1 =  1.1; % eV^2 KATRIN KNM-1
                PullTerm = PullTerm + (par(1)-0)^2/mNuSqtol_KNM1^2;
             end
             
           % 27: exponential mixing angle + m4Sq pull
           if any(ismember(obj.pullFlag,[27,28,29,30])) % sterile pull KSN-2 I
               range =  ceil(abs(obj.SO.qU(min(obj.exclDataStart))-18574));
               if obj.pullFlag == 27
                   sinU = 0.5;     % sin2t4 upper bound
                   sinD = 0;       % sin2t4 lower bound
                   m4U   = range^2; % mSq4 upper bound
                   m4D   = 0;       % mSq4 lower bound
               elseif obj.pullFlag == 28
                   sinU = 1;     % sin2t4 upper bound
                   sinD = 0;       % sin2t4 lower bound
                   m4U   = range^2; % mSq4 upper bound
                   m4D   = 0;       % mSq4 lower bound
               elseif obj.pullFlag == 29
                   sinU = 1;     % sin2t4 upper bound
                   sinD = -1;       % sin2t4 lower bound
                   m4U   = range^2; % mSq4 upper bound
                   m4D   = -range^2;       % mSq4 lower bound
               elseif obj.pullFlag == 30
                   sinU = 0;     % sin2t4 upper bound
                   sinD = -1;       % sin2t4 lower bound
                   m4U   = range^2; % mSq4 upper bound
                   m4D   = 0;       % mSq4 lower bound
               end
               PullTerm = PullTerm + ...
                   exp(1e3*(par(4*obj.SO.nPixels+12)-sinU-0.01)) + ...    % sin2t4 upper bound
                   exp(-1e3*(par(4*obj.SO.nPixels+12)-sinD+0.01)) + ...     % sin2t4 lower bound
                   exp(-1e3*(par(4*obj.SO.nPixels+11)-m4D+0.01)) + ... % mSq4 lower bound
                   exp(1e3*(par(4*obj.SO.nPixels+11)-m4U-0.01));       % mSq4 upper bound
           end
           
        end
        function chi2 = chi2function(obj,par)
            
            % Fix parameters manually when the fitter is matlab
            switch obj.fitter
                case 'matlab'
                    fixpar = obj.fixPar_v;
                    partmp = obj.parinitoriginal;
                    partmp(~ismember(1:length(obj.parinitoriginal),fixpar)) = par;
                    par = partmp;
            end
            
            if ~isempty(obj.SO2)
                % Distribute the fitting parameters
                mnu_fit   = par(1);
                Q_fit     = par(2);
                Q2_fit    = par(3);
                backs_fit = par(4:5);
                norms_fit = par(6:7);
                obj.SO.nPixels = 1;
                obj.SO.ComputeTBDDS(...
                    'mSq_bias',mnu_fit,...
                    'E0_bias',Q_fit,...
                    'B_bias',backs_fit(1),...
                    'N_bias',norms_fit(1))
                obj.SO.ComputeTBDIS;
                
                obj.SO2.ComputeTBDDS(...
                    'mSq_bias',mnu_fit,...
                    'E0_bias',Q2_fit,...
                    'B_bias',backs_fit(2),...
                    'N_bias',norms_fit(2))
                obj.SO2.ComputeTBDIS;
            else
                % Distribute the fitting parameters
                mnu_fit   = par(1);
                Q_fit     = par(2);
                backs_fit = par(3:obj.SO.nPixels+2);
                norms_fit = par((obj.SO.nPixels+3):(2*obj.SO.nPixels+2));
                DTGS_fit  = par(2*obj.SO.nPixels+3);
                DTES_fit  = par(2*obj.SO.nPixels+4);
                HTGS_fit  = par(2*obj.SO.nPixels+5);%par(4*obj.SO.nPixels+3);
                HTES_fit  = par(2*obj.SO.nPixels+6);%par(4*obj.SO.nPixels+4);
                TTGS_fit  = par(2*obj.SO.nPixels+7);%par(6*obj.SO.nPixels+3);
                TTES_fit  = par(2*obj.SO.nPixels+8);%par(6*obj.SO.nPixels+4);
                qUOffset_fit = par((2*obj.SO.nPixels+9):(3*obj.SO.nPixels+8));
                BSlope_fit  = par(3*obj.SO.nPixels+9);
                mTSq_fit    =par(3*obj.SO.nPixels+10:4*obj.SO.nPixels+9);
                FracTm_fit = par(4*obj.SO.nPixels+10);
                mnu4Sq_fit   = par(4*obj.SO.nPixels+11);
                sin2T4_fit   = par(4*obj.SO.nPixels+12);
                BPTSlope_fit = par(4*obj.SO.nPixels+13);
                
                obj.SO.ComputeTBDDS(...
                    'mSq_bias',mnu_fit,...
                    'E0_bias',Q_fit,...
                    'B_bias',backs_fit,...
                    'N_bias',norms_fit,...
                    'DTGS_bias',DTGS_fit,...
                    'DTES_bias',DTES_fit,...
                    'HTGS_bias',HTGS_fit,...
                    'HTES_bias',HTES_fit,...
                    'TTGS_bias',TTGS_fit,...
                    'TTES_bias',TTES_fit,...
                    'qUOffset_bias',qUOffset_fit,...
                    'BSlope_bias',BSlope_fit,...
                    'BPTSlope_Bias',BPTSlope_fit,...
                    'mTSq_bias',mTSq_fit,...
                    'FracTm_bias',FracTm_fit,...
                    'mnu4Sq_bias',mnu4Sq_fit,...
                    'sin2T4_bias',sin2T4_fit);
                  obj.SO.ComputeTBDIS();
            end
            
 
            % exclude data points: data and model TBDIS
            if numel(obj.exclDataStart)==1
                IncludeIndices = obj.exclDataStart:obj.exclDataStop;
            else
                IncludeIndices = obj.exclDataStart;
            end
            
            m = obj.SO.TBDIS(IncludeIndices,:); %model
            y = obj.counts(IncludeIndices,:);   %data
            z = obj.c_error(IncludeIndices,:);  %error
            
            % exclude points in covariance matrix:
            % a bit complicated, because of multiring case
            nrings = numel(obj.SO.MACE_Ba_T); % always gives correct number of pseudo rings
            nqU_used = size(m,1);             % number of subruns, which are NOT excluded
            exclIndex = sort(reshape(repmat(IncludeIndices,[nrings,1])+[0:nrings-1]'.*obj.SO.nqU,[nqU_used*nrings,1]));
            CovMatFit = obj.COVMAT(exclIndex,exclIndex);
            CovMatFracFit = obj.COVMATFrac(exclIndex,exclIndex);
            CovMatShapeFit = obj.COVMATShape(exclIndex,exclIndex);
            
            %reshape for multiring
            if strcmp(obj.SO.FPD_Segmentation,'RING')
                 m = reshape(m,[nqU_used*obj.SO.nRings,1]);
                 y = reshape(y,[nqU_used*obj.SO.nRings,1]);
                 z = reshape(z,[nqU_used*obj.SO.nRings,1]);
            end
            
            if any(y==0)
                sprintf('Cannot fit data with zero counts \n');
                return
            end
            
            PullTerm = obj.ComputePulls(par);
 
            switch obj.chi2name
                case {'chi2Stat','chi2CM'}
                        chi2 = ((y - m)')*(CovMatFit \ (y - m))  + PullTerm;
                case 'chi2CMFrac'
                    % use model to renormalize covariance matrix in each iteration
                    mNoBKG =  obj.SO.TBDIS(IncludeIndices,:)-...
                        (obj.SO.BKG_RateSec.*obj.SO.qUfrac(IncludeIndices,:).*obj.SO.TimeSec); %model without background
                    mNoBKG(obj.SO.qU>obj.SO.Q,:) = 0; mNoBKG(mNoBKG<0) = 0;
                    CM_tmp = CovMatFracFit.*mNoBKG.*mNoBKG';
                    chi2 = ((y - m)')*(CM_tmp \ (y - m)) + PullTerm;
                case 'chi2CMShape'
                    chi2 = ((y - m)')* (CovMatShapeFit  \ (y - m)) + PullTerm;
                case 'chi2Nfix'
                    normFactor = sum(y(obj.qUdata<18575,:))/...
                        sum(m(obj.qUdata<18575,:))- 1;
                    chi2 = (y - normFactor*m)' * ...
                        (obj.COVMATShape(obj.exclDataStart:obj.exclDataStop,obj.exclDataStart:obj.exclDataStop)...
                        \ (y - normFactor*m))  + sum(((par-zeros(1,length(par))).^2)./obj.pulls);
                case 'chi2P'
                    chi2 = 0;
                    for ii = 1:length(y)
                        if y(ii) == 0
                            chi2temp = 2*m(ii);
                        else
                            chi2temp = -2*(y(ii) - m(ii) + y(ii)*log(m(ii)/y(ii)));
                        end
                        chi2 = chi2 + chi2temp;
                    end
                    chi2 = chi2 + PullTerm;
                case 'chi2meanSys'
                    chi2 = ((obj.DATA - par)')*obj.COVMAT\(obj.DATA - par);
            end
            
            if isa(obj.SO,'RHOD')
                obj.SO.FPD_Segmentation = obj.SO.TBDCell{1}.FPD_Segmentation;
                obj.SO.Q = obj.SO.TBDCell{1}.Q;
                bcki_temp = 0;
                for mm = 1:obj.SO.nModels
                    bcki_temp = bcki_temp + obj.SO.TBDCell{mm}.BKG_RateSec_i;
                end
                obj.SO.BKG_RateSec_i = bcki_temp;
            end
            
        end % chi2function
        
        function chi2 = chi2function2(obj,par)
            % chi2 function for combi knm1+2 fit
            % Distribute the fitting parameters
            mnu_fit   = par(1);
            Q_fit     = par(2);
            Q2_fit    = par(3);
            backs_fit = par(4:5);
            norms_fit = par(6:7);
            obj.SO.nPixels = 1;
            
            % model 1 (knm1)
            obj.SO.ComputeTBDDS(...
                'mSq_bias',mnu_fit,...
                'E0_bias',Q_fit,...
                'B_bias',backs_fit(1),...
                'N_bias',norms_fit(1))
            obj.SO.ComputeTBDIS;
            
            % model 2 (knm2)
            obj.SO2.ComputeTBDDS(...
                'mSq_bias',mnu_fit,...
                'E0_bias',Q2_fit,...
                'B_bias',backs_fit(2),...
                'N_bias',norms_fit(2))
            obj.SO2.ComputeTBDIS;

            % exclude data points: data and model TBDIS
            % only 40 eV range possible
            range = 40;   
            IncludeIndices = obj.qUdata-18574>=-range;
            m1tmp = obj.SO.TBDIS(2:end);
            m1 = m1tmp(IncludeIndices(:,1));%obj.SO.TBDIS(IncludeIndices(:,1)); %model
            y1 = obj.counts(IncludeIndices(:,1),1);   %data
            z1 = obj.c_error(IncludeIndices(:,1),1);  %error
            
            m2 = obj.SO2.TBDIS(IncludeIndices(:,2)); %model
            y2 = obj.counts(IncludeIndices(:,2),2);   %data
            z2 = obj.c_error(IncludeIndices(:,2),2);  %error
            
            m = [m1;m2];
            y = [y1;y2];
            z = [z1;z2];

            % exclude points in covariance matrix:
            IncludeIndicesAll = [IncludeIndices(:,1);IncludeIndices(:,2)];
            CovMatFit = obj.COVMAT(IncludeIndicesAll,IncludeIndicesAll);
            CovMatShapeFit = obj.COVMATShape(IncludeIndicesAll,IncludeIndicesAll);
             
            if any(y==0)
                sprintf('Cannot fit data with zero counts \n');
                return
            end
            
            PullTerm = obj.ComputePulls(par);
 
            switch obj.chi2name
                case {'chi2Stat','chi2CM'}
                        chi2 = ((y - m)')*(CovMatFit \ (y - m))  + PullTerm;
                case 'chi2CMShape'
                    chi2 = ((y - m)')* (CovMatShapeFit  \ (y - m)) + PullTerm;
            end
        end
        function results = startfit(obj)
            if obj.preparevarsdone
                switch obj.fitter
                    case 'matlab'
                        cprintf('blue','----------BEGIN FIT MATLAB-------------- \n');
                        
                        chi2matlab = @obj.chi2function;
                        nparinit = length(obj.parinit);
                        obj.parinitoriginal = obj.parinit;
                        partmp = obj.parinit;
                        errtmp = zeros(1,nparinit);
                        obj.parinit(obj.fixPar_v) = [];
                        options = optimoptions('fminunc','Algorithm','quasi-newton',...
                            'OptimalityTolerance',1e-6,'StepTolerance',5e-7,...
                            'FiniteDifferenceType','central',...
                            'MaxFunctionEvaluations',1e6,'UseParallel',true,...
                            'Display','notify','ObjectiveLimit',-1e20);
                        
                        if obj.UncMPixFlag
                            [par,chi2min,exitflag,~,~,Hessian] = fminunc(chi2matlab,obj.parinit,options);
                            errmat = 0.5*Hessian;
                            varcov = inv(errmat);
                            err = sqrt(diag(varcov))';
                            errtmp(~ismember(1:nparinit,obj.fixPar_v)) = err;
                            err = errtmp;
                        else
                            [par,chi2min] = fminunc(chi2matlab,obj.parinit,options);
                            err = zeros(1,nparinit);
                            errmat = zeros(nparinit);
                        end
                        
                        obj.parinit = obj.parinitoriginal;
                        partmp(~ismember(1:nparinit,obj.fixPar_v)) = par;
                        par = partmp;
                        
%                         if strcmp(obj.SO.FPD_Segmentation,'MULTIPIXEL') && obj.UncMPixFlag
%                             savefixPar = obj.fixPar;
%                             err(2) = obj.getCorrectUncertaintiesE0MultiPix(par,chi2min,err);
%                             obj.UncMPixFlag = true;
%                             obj.fixPar = savefixPar;
%                         end
                        
                        dof = obj.SO.nPixels*(obj.exclDataStop-obj.exclDataStart+1) - ...
                            nparinit+obj.nfixPar;
                        
                        results = {par, err, chi2min, errmat, dof, obj.COVMAT};

                    case 'minuit'
                        cprintf('blue','----------BEGIN FIT MINUIT-------------- \n');

                       % set lim 16 0 0.5 ; set lim 15 0 10000;
                        tmparg = sprintf(['%s set pri 1; set STRATEGY 2;'...
                            '%s'],obj.fixPar,obj.minuitOpt);
                        Args = {obj.parinit, obj, '-c',tmparg};
                        minuitOutputStr = [tempname(pwd),'.txt']; % random string 
                        diary(minuitOutputStr)
                        if isempty(obj.SO2)
                            [par, err, chi2min, errmat] = fminuit('chi2minuit',Args{:});
                        else
                            [par, err, chi2min, errmat] = fminuit('chi2minuitCombi',Args{:});
                        end
                        diary off
                        fixpar_local = strrep(obj.fixPar,';','');
                        fixpar_local = strrep(fixpar_local,'fix','');
                        
                        if numel(obj.exclDataStart)==1
                            nqU = obj.exclDataStop-obj.exclDataStart+1;
                        else
                            nqU = numel(obj.exclDataStart);
                        end
                        dof = obj.SO.nPixels*nqU - ...
                            length(obj.parinit)+numel(str2num(fixpar_local));
                        results = {par, err, chi2min, errmat, dof};
                        
                        if contains(obj.minuitOpt,'minos')
                            % read asymmetric uncertainties from text file
                            [errNeg, errPos] = GetAsymmErrMinos(minuitOutputStr);
                            results = {results{:},errNeg, errPos};
                        end
                        system(['rm -f ',minuitOutputStr]); % delete txt file 
                end
              
            else
                obj.preparevars;
                obj.startfit;
            end
        end
        
        function displayresults(obj)
            par     = obj.RESULTS{1}; 
            err     = obj.RESULTS{2};
            chi2min = obj.RESULTS{3};
            dof     = obj.RESULTS{5};
            
            % Normalization (pixelwise) (uniform 4)
            norms_fit = par(3+obj.SO.nPixels:3+2*obj.SO.nPixels-1);
            norms_fit_err = err(3+obj.SO.nPixels:3+2*obj.SO.nPixels-1);
            
            % Background (pixelwise) (uniform 3)
            bcks_fit     = par(3:2+obj.SO.nPixels);
            bcks_fit_err = err(3:2+obj.SO.nPixels);
            
            % Neutrino Mass squared  (pixelwise)
            % Begin WARNING: Changed 19/09/2019
            if isempty(strfind(obj.fixPar,'fix 1'))
                mnuSq_report = obj.parinit(1) + par(1) + obj.SO.mnuSq_i;
            else
                mnuSq_report = par(1)+ obj.SO.mnuSq_i;
            end
            % End Warning
           
            % Neutrino Mass  (pixelwise)
            if mnuSq_report>=0
            mnu_report   = sqrt(abs(mnuSq_report));
            err_mnu      = sqrt(abs((mnuSq_report + err(1)))) -  mnu_report;
            else
            mnu_report   = 0;
            err_mnu      = sqrt(abs(err(1))); 
            end

            % Endpoint
            e0_fit = par(2);
 
            %qUOffset 
            qUoffset_fit     = par(2*obj.SO.nPixels+9:3*obj.SO.nPixels+8);
            qUoffset_fit_err = err(2*obj.SO.nPixels+9:3*obj.SO.nPixels+8);
            
            % background slope
            Bslope_fit     = par(3*obj.SO.nPixels+9)*1e6;
            Bslope_fit_err = err(3*obj.SO.nPixels+9)*1e6;
            
            % tachynoc neutrino mass squared
            mTSq_fit      = par(3*obj.SO.nPixels+10:4*obj.SO.nPixels+9);
            mTSq_fit_err  = err(3*obj.SO.nPixels+10:4*obj.SO.nPixels+9);
            
            % T minus fraction
            FracTm_fit      = par(4*obj.SO.nPixels+10);
            FracTm_fit_err  = err(4*obj.SO.nPixels+10);
            
            % Sterile Neutrino Parameters
            mnu4Sq_fit      = par(4*obj.SO.nPixels+11);
            mnu4Sq_fit_err  = err(4*obj.SO.nPixels+11);
            sin2T4_fit      = par(4*obj.SO.nPixels+12);
            sin2T4_fit_err  = err(4*obj.SO.nPixels+12);
            
            % background PT slope
            BPTslope_fit     = (par(4*obj.SO.nPixels+13)+obj.SO.BKG_PtSlope_i).*1e6;
            BPTslope_fit_err = err(4*obj.SO.nPixels+13).*1e6;
            
            cprintf('blue','===============================================\n');
            cprintf('blue','  m^2       = %g +/- %g eV^2\n', mnuSq_report,err(1));
            cprintf('blue','  m         = %g +/- %g eV\n', mnu_report,err_mnu);
            cprintf('blue',' - - - - - - - - - - - - - - - - - - - - - - - \n');
            if (mnu4Sq_fit_err~=0) || (sin2T4_fit_err~=0)
                cprintf('blue','  mnu4Sq    = %g +/- %g eV\n', mnu4Sq_fit + obj.SO.mnu4Sq_i,mnu4Sq_fit_err);
                cprintf('blue','  sin2T4    = %g +/- %g \n', sin2T4_fit + obj.SO.sin2T4_i,sin2T4_fit_err);
                cprintf('blue',' - - - - - - - - - - - - - - - - - - - - - - - \n');
            end
            cprintf('blue','  (E0)eff   = %0.9g +/- %g eV\n', obj.SO.Q_i+e0_fit,err(2));
            cprintf('blue',' - - - - - - - - - - - - - - - - - - - - - - - \n');
            
            switch obj.SO.FPD_Segmentation
                case 'RING'
                    cprintf('blue','  B (peudo ring %.0f)       = %g +/- %g mcps\n', [1:numel(bcks_fit); (obj.SO.BKG_RateSec_i + bcks_fit)*1e3;bcks_fit_err*1e3]);
                    cprintf('blue',' - - - - - - - - - - - - - - - - - - - - - - - \n');
                    cprintf('blue','  N (pseudo ring %.0f)         = %g +/- %g\n', [1:numel(bcks_fit);norms_fit+1;norms_fit_err]);
              
                case 'OFF'
                    cprintf('blue','  B         = %g +/- %g mcps (%.0f pixels) \n', (obj.SO.BKG_RateSec_i + bcks_fit)*1e3,bcks_fit_err*1e3,numel(obj.SO.FPD_PixList));
                    cprintf('blue',' - - - - - - - - - - - - - - - - - - - - - - - \n');
                    cprintf('blue','  N         = %g +/- %g\n', norms_fit+1,norms_fit_err);
            end

            [~,path] = system('pwd');
            if contains(path,'lasserre')

                % LATEX FORMAT: START
                t = PrintTable(sprintf('Fit Results  - Fitter = %s',obj.fitter));
                t.addRow('Parameter','Value','Uncertainty','Unit');
                t.addRow('m^2_\beta',sprintf('%.2f',mnuSq_report),sprintf('%.2f',err(1)),'eV^2');
                t.addRow('E_{0,eff}',sprintf('%.2f',obj.SO.Q_i+e0_fit),sprintf('%.2f',err(2)),'eV');
                t.addRow('Background',sprintf('%.2f',(obj.SO.BKG_RateSec_i + bcks_fit)*1e3),sprintf('%.2f',bcks_fit_err*1e3),'mcps');
                t.addRow('Normalization',sprintf('%.2f',norms_fit+1),sprintf('%.2f',norms_fit_err),'');
                t.addRow('','','','')
                t.addRow('\chi^2','dof','','pvalue')
                t.addRow(sprintf('%.2f',chi2min),sprintf('%.2f',dof),'',sprintf('%.3f',chi2pvalue(chi2min,dof)));
                t.display;
                t.HasHeader = true;
                t.Format = 'tex';
                t.Caption = sprintf('Fit Results  - Fitter = %s - Range = [%.1f - %.1f] eV',...
                    obj.fitter,obj.SO.qU(obj.exclDataStart),obj.SO.qU(obj.exclDataStop));
                t.print;
                %LATEX FORMAT: END

            end
            %FSD
            if all(~ismember(char(strjoin(string([2*obj.SO.nPixels+3,2*obj.SO.nPixels+4]))),obj.fixPar)) %when DT FSD NOT fixed
                cprintf('blue',' - - - - - - - - - - - - - - - - - - - - - - - \n');
                cprintf('blue','DT-FSD Ground  State Prob = %.3f +/- %.3f \n',obj.SO.DTNormGS_i + par(2*obj.SO.nPixels+3),err(2*obj.SO.nPixels+3));
                cprintf('blue','DT-FSD Excited State Prob = %.3f +/- %.3f \n',obj.SO.DTNormES_i + par(2*obj.SO.nPixels+4),err(2*obj.SO.nPixels+4));
            end
            if all(~ismember(char(strjoin(string([2*obj.SO.nPixels+5,2*obj.SO.nPixels+6]))),obj.fixPar))
                cprintf('blue',' - - - - - - - - - - - - - - - - - - - - - - - \n');
                cprintf('blue','HT-FSD Ground  State Prob = %.3f +/- %.3f \n',obj.SO.HTNormGS_i + par(2*obj.SO.nPixels+5),err(2*obj.SO.nPixels+5));
                cprintf('blue','HT-FSD Excited State Prob = %.3f +/- %.3f \n',obj.SO.HTNormES_i + par(2*obj.SO.nPixels+6),err(2*obj.SO.nPixels+6));
            end
            if ~(contains(obj.fixPar,char(string(2*obj.SO.nPixels+7))) && contains(obj.fixPar,char(string(2*obj.SO.nPixels+8))))
                cprintf('blue',' - - - - - - - - - - - - - - - - - - - - - - - \n');
                cprintf('blue','TT-FSD Ground  State Prob = %.3f +/- %.3f \n',obj.SO.TTNormGS_i + par(2*obj.SO.nPixels+7),err(2*obj.SO.nPixels+7));
                cprintf('blue','TT-FSD Excited State Prob = %.3f +/- %.3f \n',obj.SO.TTNormES_i + par(2*obj.SO.nPixels+8),err(2*obj.SO.nPixels+8));
            end
            
            %qUOffset
            if  ~(contains(obj.fixPar,char(string(2*obj.SO.nPixels+10))))
                 cprintf('blue',' - - - - - - - - - - - - - - - - - - - - - - - \n');
                switch obj.SO.FPD_Segmentation
                    case 'RING'      
                        cprintf('blue','  qU offset (pseudo ring %.0f)         = %g +/- %g eV \n', [1:numel(bcks_fit);qUoffset_fit; qUoffset_fit_err]);
                    case 'OFF' 
                        cprintf('blue','  qU offset                     = %g +/- %g eV \n', qUoffset_fit, qUoffset_fit_err);        
                end
            end
            % background qU slope
            if  ~(contains(obj.fixPar,char(string(2*obj.SO.nPixels+10))))
                  cprintf('blue',' - - - - - - - - - - - - - - - - - - - - - - - \n');
                  cprintf('blue','  B slope  = %.3g +/- %.3g mcps/keV \n',  Bslope_fit, Bslope_fit_err);
            end
             % background time slope
            if  ~(contains(obj.fixPar,char(string(4*obj.SO.nPixels+13))))
                  cprintf('blue',' - - - - - - - - - - - - - - - - - - - - - - - \n');
                  cprintf('blue','  B PT slope  = %.3g +/- %.3g mucps/s \n',  BPTslope_fit, BPTslope_fit_err);
            end
            % neutrino mass offsets
            if ~(contains(obj.fixPar,char(string(4*obj.SO.nPixels+9))))
                 cprintf('blue',' - - - - - - - - - - - - - - - - - - - - - - - \n');
                 cprintf('blue',' mT^2 (pseudo ring %.0f)       = %g +/- %g eV^2 \n', [1:numel(mTSq_fit);mTSq_fit; mTSq_fit_err]);
            end
            % T-minus ions
            if ~(contains(obj.fixPar,char(string(4*obj.SO.nPixels+10))))
                cprintf('blue',' - - - - - - - - - - - - - - - - - - - - - - - \n');
                cprintf('blue',' Fraction T^- ion         = %g +/- %g   \n', FracTm_fit, FracTm_fit_err);
            end
            
            cprintf('blue',' - - - - - - - - - - - - - - - - - - - - - - - \n');
            cprintf('blue','  chi2/dof  = %g/%g\n', chi2min, dof);
            cprintf('blue','  p-value   = %g\n', chi2pvalue(chi2min,dof));
            cprintf('blue','===============================================\n');
  
        end
        
        function err = getCorrectUncertaintiesE0MultiPix(obj,par,chi2min,errHessian)
            % E0 Scan to get a better estimate of the uncertainties on E0
            % when doing a multipixel fit, as the uncertainties obtained by
            % the Hessian matrix might be innacurate due to the numerical
            % derivatives with many parameters
            
            
            %% E0 Scan
            obj.parinit = par;
            mnu_fit = par(1);
            E0_best = par(2);
            backs_fit = par(3:3+obj.SO.nPixels-1);
            norms_fit = par(3+obj.SO.nPixels:end-2);
            
            % WARNING
            DTGS_fit = par(end-1);
            DTES_fit = par(end);
            % End WARNING
            
            obj.fixPar = [obj.fixPar,' 2'];
            
            steps = 4;
 
            switch obj.chi2name
                case {'chi2Stat','chi2P'}  
                    uncRangeStat = 0.1; 
                    e0frac_list = [linspace(-uncRangeStat,-uncRangeStat/steps,steps),0,linspace(uncRangeStat/steps,uncRangeStat,steps)];
                    e0frac = e0frac_list  + E0_best;
                otherwise
                    uncRangeSys = 0.5;
                    e0frac_list = [linspace(-uncRangeSys,-uncRangeSys/steps,steps),0,linspace(uncRangeSys/steps,uncRangeSys,steps)];
                    e0frac = e0frac_list + E0_best;
            end

            chi2min_E0 = zeros(1,length(e0frac));
            
            for e0 = 1:length(e0frac)
                
                E0_fit = e0frac(e0);
                obj.parinit(2) = E0_fit;
                
                obj.SO.ComputeTBDDS(...
                    'mSq_bias',mnu_fit,...
                    'E0_bias',E0_fit,...
                    'B_bias',backs_fit,...
                    'N_bias',norms_fit,...
                    'DTGS_bias',DTGS_fit,...
                    'DTES_bias',DTES_fit);
                obj.SO.ComputeTBDIS();
                
                obj.UncMPixFlag = false;
                results = obj.startfit();
                
                chi2min_E0(:,e0) = results{3};
            end
            
            %% Find asymmetric errors
            chi2err = chi2min + 1;
            chi2min_E0LeftofMin = chi2min_E0(e0frac_list <= 0);
            E0LeftofMin = e0frac(e0frac_list <= 0);
            chi2min_E0RightofMin = chi2min_E0(e0frac_list >= 0);
            E0RightofMin = e0frac(e0frac_list >= 0);
            E0lowerUnc = interp1(chi2min_E0LeftofMin,E0LeftofMin,chi2err,'spline');
            E0upperUnc = interp1(chi2min_E0RightofMin,E0RightofMin,chi2err,'spline');
            
            err = (E0upperUnc-E0lowerUnc)/2;
            errH = errHessian(2);
            errDiff = err-errH;
            a = 1;
        end
        
        function f = model2designmatrix(obj,p)
            %
            % Function for Model used for Design Matrix
            % Input:
            %  p: fit parameters
            %
            % Output:
            %  function
            %
            % T. Lasserre, June 2018
            
            if numel(p)<5
                p(5)=0; p(6)=0;p(7)=0; p(8)=0;p(9)=0; p(10)=0;
            end
            %             obj.SO.ComputeTBDDS(...
            %                 'mSq_bias',p(1),...
            %                 'E0_bias',p(2),...
            %                 'B_bias',p(3),...
            %                 'N_bias',p(4),...
            %                 'DTGS_bias',p(5),...
            %                 'DTES_bias',p(6),...
            %                 'HTGS_bias',p(7),...
            %                 'HTES_bias',p(8),...
            %                 'TTGS_bias',p(9),...
            %                 'TTES_bias',p(10));
            %        end
            % %             obj.SO.ComputeTBDDS(...
            % %                 'mSq_bias',p(1),...
            % %                 'E0_bias',p(2),...
            % %                 'B_bias',p(3:6),...
            % %                 'N_bias',p(7:10),...
            % %                 'DTGS_bias',p(11),...
            % %                 'DTES_bias',p(12));
            
            % Distribute the fitting parameters
            mnu_fit   = p(1);
            Q_fit     = p(2);
            backs_fit = p(3:obj.SO.nPixels+2);
            norms_fit = p((obj.SO.nPixels+3):(2*obj.SO.nPixels+2));
            DTGS_fit  = p(2*obj.SO.nPixels+3);
            DTES_fit  = p(2*obj.SO.nPixels+4);
            HTGS_fit  = p(2*obj.SO.nPixels+5);%p(4*obj.SO.nPixels+3);
            HTES_fit  = p(2*obj.SO.nPixels+6);%p(4*obj.SO.nPixels+4);
            TTGS_fit  = p(2*obj.SO.nPixels+7);%p(6*obj.SO.nPixels+3);
            TTES_fit  = p(2*obj.SO.nPixels+8);%p(6*obj.SO.nPixels+4);
            qUOffset_fit = p((2*obj.SO.nPixels+9):(3*obj.SO.nPixels+8));
            BSlope_fit  = p(3*obj.SO.nPixels+9);
            mTSq_fit    =p(3*obj.SO.nPixels+10:4*obj.SO.nPixels+9);
            FracTm_fit = p(4*obj.SO.nPixels+10);
            
            if obj.SO.nPixels <2
                obj.SO.ComputeTBDDS(...
                    'mSq_bias',mnu_fit,...
                    'E0_bias',Q_fit,...
                    'B_bias',backs_fit,...
                    'N_bias',norms_fit,...
                    'DTGS_bias',DTGS_fit,...
                    'DTES_bias',DTES_fit,...
                    'HTGS_bias',HTGS_fit,...
                    'HTES_bias',HTES_fit,...
                    'TTGS_bias',TTGS_fit,...
                    'TTES_bias',TTES_fit);
            else
                obj.SO.ComputeTBDDS(...
                    'mSq_bias',mnu_fit,...
                    'E0_bias',Q_fit,...
                    'B_bias',backs_fit,...
                    'N_bias',norms_fit,...
                    'DTGS_bias',DTGS_fit,...
                    'DTES_bias',DTES_fit,...
                    'HTGS_bias',HTGS_fit,...
                    'HTES_bias',HTES_fit,...
                    'TTGS_bias',TTGS_fit,...
                    'TTES_bias',TTES_fit,...
                    'qUOffset_bias',qUOffset_fit,...
                    'BSlope_bias',BSlope_fit,...
                    'mTSq_bias',mTSq_fit,...
                    'FracTm_bias',FracTm_fit);
            end
            
            obj.SO.ComputeTBDIS();
            ftemp = obj.SO.TBDIS(obj.exclDataStart:obj.exclDataStop,:);
            f = ftemp(:);
        end
        
        function designmatrix(obj,varargin)
            %
            % Compute Design Matrix https://en.wikipedia.org/wiki/Design_matrix
            % Input:
            %  Epsilon for derivative computation - 1e-4 by default
            %  FitPar: fit parameters
            %  FitParList
            %
            % Output:
            %  DesignMatrix
            %
            % T. Lasserre, June 2018
            %
            
            % Parser
            p = inputParser;
            p.addParameter('epsilon',1e-4,@(x)isfloat(x));
            p.addParameter('FitParList',1:numel(obj.RESULTS{1}));
            p.parse(varargin{:});
            epsilon       = p.Results.epsilon;
            FitParList    = p.Results.FitParList;
            
            % Init
            nFitPar      = length(obj.RESULTS{1});
            f0           = feval(@obj.model2designmatrix,obj.RESULTS{1});
            nBins        = length(f0);
            obj.DesignMatrix = zeros(nBins,nFitPar);
            
            % Compute Design Matrix
            for i=FitParList
                Variation = [zeros(i-1,1); epsilon; zeros(nFitPar-i,1)];
                if i<11
                obj.DesignMatrix(:,i) = ( feval(@obj.model2designmatrix,obj.RESULTS{1}+Variation') ...
                    - feval(@obj.model2designmatrix,obj.RESULTS{1}-Variation')) / (2*epsilon);
                end
                
            end
            obj.DesignMatrix = obj.DesignMatrix(:,(obj.MaskFreeFitPar));
            feval(@obj.model2designmatrix,obj.RESULTS{1});
        end
        
    end 
    methods % cats + sanity checks
        function catss(obj)
            %
            % Compute CATS diagnostics
            %
            % T. Lasserre, June 2018
            %
            
            nFitPar       = length(obj.RESULTS{1});
            AllParameters = 1:nFitPar;
            if nFitPar    == 17  % Uniform - m2,E0,N,B
            AllLabels     = {'m^2 (eV^2)' ,'E_0 (eV)' ,'B (cps)' ,'N','Pgs_{DT}','Pes_{DT}','Pgs_{HT}','Pes_{HT}','Pgs_{TT}','Pes_{TT}','qUoffset','12','13','14','15'};
            elseif nFitPar    == 16  % Uniform - m2,E0,N,B
            AllLabels     = {'m^2 (eV^2)' ,'E_0 (eV)' ,'B (cps)' ,'N','Pgs_{DT}','Pes_{DT}','Pgs_{HT}','Pes_{HT}','Pgs_{TT}','Pes_{TT}','qUoffset','12','13','14','m_4^2 (eV^2)','sint4^2'};
            elseif nFitPar== 18  % MultiRing 4 - m2,E0,N,B
            AllLabels     = {'m^2 (eV^2)' ,'E_0 (eV)', 'B1 (cps)', 'B2 (cps)' , 'N1', 'N2', 'Pgs_{DT}','Pes_{DT}','Pgs_{HT}','Pes_{HT}','Pgs_{TT}','Pes_{TT}','qUoffset1','qUoffset2','Bslope','mTSq1','mTSq2','T-'};
            elseif nFitPar== 26  % MultiRing 4 - m2,E0,N,B
            AllLabels     = {'m^2 (eV^2)' ,'E_0 (eV)', 'B1 (cps)', 'B2 (cps)' ,'B3 (cps)' ,'B4 (cps)' , 'N1', 'N2', 'N3', 'N4','Pgs_{DT}','Pes_{DT}','Pgs_{HT}','Pes_{HT}','Pgs_{TT}','Pes_{TT}','qUoffset1','qUoffset2','qUoffset3','qUoffset4','Bslope','mTSq1','mTSq2','mTSq3','mTSq4','T-'};
            end
                
            % obj.nfixPar = numel(str2num(strrep(strrep(obj.fixPar,'fix ',''),' ;',' ')));
            % obj.fixPar_v = str2num(strrep(strrep(obj.fixPar,'fix ',''),' ;',' '));
            [obj.MaskFreeFitPar , ilabels] = setdiff(AllParameters,obj.fixPar_v);
            obj.MaskFreeFitLabel           = AllLabels(ilabels);
            
%             % Define Number of Free Parameters
%             switch  obj.fixPar
%                 case '1 5 6 7 8 9 10'         % nu mass, FSD Fixed
%                     obj.MaskFreeFitPar     = [2 3 4];
%                     
%                 case ''                       % all free
%                     obj.MaskFreeFitPar     = 1:4;
%                 case '1'                      % nu mass fixed
%                     obj.MaskFreeFitPar     = 2:4;
%                 case '1 3 4'                  % nu mass, background, norm fixed, FSD Fixed
%                     obj.MaskFreeFitPar     = 2;
%                 case '5 6 7 8 9 10'           %  FSD Fixed
%                     obj.MaskFreeFitPar     = 1:4;
%                 
%                 case '1 11 12'                % ?
%                     obj.MaskFreeFitPar     = [2:10];
%             end
            
            % Data - before 27/3/20
            % D=obj.DATA(obj.exclDataStart:obj.exclDataStop,2);
            % % D=obj.DATA(obj.exclDataStart:obj.exclDataStop,obj.SO.nPixels+[1:4]);
            % D = D(:);

            % Data - 27/3/20
            m         = obj.SO.TBDIS(obj.exclDataStart:obj.exclDataStop,:); %model
            nrings    = numel(obj.SO.MACE_Ba_T); % always gives correct number of pseudo rings
            nqU_used  = size(m,1);             % number of subruns, which are NOT excluded
            exclIndex = sort(reshape(repmat(obj.exclDataStart:obj.exclDataStop,[nrings,1])+[0:nrings-1]'.*obj.SO.nqU,[nqU_used*nrings,1]));            
            D         = obj.counts(exclIndex);
                            
            % Model at best fit
            Mbf=obj.model2designmatrix(obj.RESULTS{1});
            
            % Design Matrix
            obj.designmatrix('epsilon',1e-4);
            
            % Covariance Matrix
            switch obj.chi2name
                case 'chi2Stat'
                    CovMat = diag(D);
                case 'chi2CM'
                    %CovMat = obj.COVMAT(obj.exclDataStart:obj.exclDataStop,obj.exclDataStart:obj.exclDataStop);
                    CovMat = obj.COVMAT(exclIndex,exclIndex);
                case 'chi2CMFrac'
                    %CovMat = obj.COVMATFrac(obj.exclDataStart:obj.exclDataStop,obj.exclDataStart:obj.exclDataStop).*m.*m';
                    CovMat = obj.COVMATFrac(exclIndex,exclIndex).*m.*m';
                case 'chi2CMShape'
                    %CovMat = obj.COVMATShape(obj.exclDataStart:obj.exclDataStop,obj.exclDataStart:obj.exclDataStop);
                    CovMat = obj.COVMATShape(exclIndex,exclIndex);
            end
            
            % CATS call (in tools)
            obj.CATSstat = cats( D - Mbf + obj.DesignMatrix*obj.RESULTS{1}(obj.MaskFreeFitPar)', ...
                obj.DesignMatrix , ...
                'Vinv' , inv(CovMat),...
                'RenormalizeErrors',false);
            feval(@obj.model2designmatrix,obj.RESULTS{1});
            
             % Label
            %obj.CATSstat.names = [round(obj.DATA(obj.exclDataStart:end,1)-obj.SO.Q_i,1) ; obj.CATSstat.nsys];
            obj.CATSstat.names = [round(obj.qUdata(exclIndex)-obj.SO.Q_i,0) ; obj.CATSstat.nsys];

        end
        
        function CATSplot(obj)
            %
            % Plot CATS diagnostics
            % -Residuals
            % -Leverages
            % -Cook's Distance
            % -Residual/Leverage/Cook's
            %
            % T. Lasserre, June 2018
            %
            
            % Residual Plot
            figure(10001)
            resplot(obj.CATSstat);
            PrettyFigureFormat;%XXX
            
            % Leverage Plot
            figure(10002)
            levplot(obj.CATSstat);
            PrettyFigureFormat;
            
            % Cook Plot
            figure(10003)
            cookdplot(obj.CATSstat);
            PrettyFigureFormat;
            ylabel('qU (measurement points)')
            % summary plot
            f1 = figure(10004)
            set(f1, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.8, 0.8]);
            stdrlevplot(obj.CATSstat)
            PrettyFigureFormat;
            set(gca,'FontSize',16);
            %title(sprintf('%s Fit to Stacked Runs (40531-40693) - %.0f eV below E0',obj.chi2name,obj.SO.Q_i-obj.SO.qU(obj.exclDataStart)));
            figure(10005)
            standresplot(obj.CATSstat);
            PrettyFigureFormat
            figure(10006)
            covratioplot(obj.CATSstat);
            PrettyFigureFormat
            figure(10007)
            dmseplot(obj.CATSstat);
            %            PrettyFigureFormat
            %             figure(10008)
            %             dfxsplot(obj.CATSstat);
            %             figure(100081)
            %             dfxiplot(obj.CATSstat,2);
            PrettyFigureFormat
            figure(10009)
            plotcorh(obj.CATSstat);
            PrettyFigureFormat
            figure(10010)
            subplot(1,2,1)
            dffitplot(obj.CATSstat);
            PrettyFigureFormat
            subplot(1,2,2)
            dffitsplot(obj.CATSstat);
            PrettyFigureFormat
            figure(100011)
            deleteparplot(obj.CATSstat);
            figure(100012)
            deletevarianceplot(obj.CATSstat);
            figure(100013)
            sclchinparplot(obj.CATSstat);
        end
        
        function Samakcats_StandResidual_Leverage(obj,varargin)
            % Standardized residual plot of the fit performed by CATS
            % Plot the standardized residuals of the fit performed by CATS
            % Plot Leverage  of the fit performed by CATS
            % Copyright 2005-2011 Guillaume MENTION, CEA Saclay
            % $Revision: 0.99$  $Date: 2011/11/15$
            
            p=inputParser;
            p.addParameter('savePlot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            savePlot   = p.Results.savePlot;
            
            nobs = obj.CATSstat.nobs;
            nsys = obj.CATSstat.nsys;
            npar = obj.CATSstat.npar;
            
            % Figure definition
            myMainTitle=sprintf('Samak: Standardized Residuals & Leverages (CATS)');
            savefile=sprintf('plots/cats/CATS_StandResidual_Leverage.pdf');
            fig = figure('Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 1000 1000]);
            a=annotation('textbox', [0 0.88 1 0.1], ...
                'String', myMainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');
            a.FontSize=20;a.FontWeight='bold';
            
            subplot(1,2,1)
            barh(obj.CATSstat.standres(1:nobs),...
                'FaceColor',rgb('Amethyst'),...
                'EdgeColor',rgb('SteelBlue'),...
                'LineWidth',2);
            colormap(lines(nsys));
            cmap = [.5*ones(1,3); circshift(lines(nsys),-1)];
            hold on;
            for i=1:nsys
                barh(nobs+i,obj.CATSstat.standres(nobs+i),'FaceColor',cmap(i,:));
            end
            hold off;
            box on;
            axis([min(0,min(obj.CATSstat.standres)*1.2) max(obj.CATSstat.standres)*1.2 0 length(obj.CATSstat.standres)+1 ]);
            grid off;
            xlabel('Standardized Residuals');
           ylabel(sprintf('Retarding energy - 18574 (eV)'));
            % set(gca,'Xtick',round(linspace(1,length(obj.CATSstat.leverage),min(10,nobs+nsys))));
            if isempty(obj.CATSstat.names)
                set(gca,'YDir','Reverse','Ytick',1:(nobs+nsys),'YtickLabel',1:(nobs+nsys));
                grid on;
                set(gca,'YGrid','off');
                PrettyFigureFormat('FontSize',24);
            else
                set(gca,'YDir','Reverse','Ytick',1:(nobs+nsys),'YtickLabel',obj.CATSstat.names(1:(nobs+nsys)));
                PrettyFigureFormat('FontSize',24);
                grid on;
                set(gca,'YGrid','off');
            end
           
            if nobs>40
                 a = get(gca,'YTickLabel');
                set(gca,'YTickLabel',a,'fontsize',8);
            end
            
            subplot(1,2,2)
            barh(obj.CATSstat.leverage(1:nobs),...
                'FaceColor',rgb('Amethyst'),...
                'EdgeColor',rgb('SteelBlue'),...
                'LineWidth',2);
            if nsys==0
                colormap('default')
            else
                colormap(lines(nsys));
            end
            cmap = [.5*ones(1,3); circshift(lines(nsys),-1)];
            hold on;
            for i=1:nsys
                barh(nobs+i,obj.CATSstat.leverage(nobs+i),'FaceColor',cmap(i,:));
            end
            hold off;
            box on;
            axis([ 0 max(obj.CATSstat.leverage)*1.2 0 length(obj.CATSstat.leverage)+1]);
            grid off;
            ylabel(sprintf('Retarding energy - 18574 (eV)'));
            xlabel('Leverages');
            % set(gca,'Xtick',round(linspace(1,length(obj.CATSstat.leverage),min(10,nobs+nsys))));
            if isempty(obj.CATSstat.names)
                set(gca,'YDir','Reverse','Ytick',1:(nobs+nsys),'YtickLabel',1:(nobs+nsys));
                grid on;
                set(gca,'YGrid','off');
                PrettyFigureFormat('FontSize',24);
            else
                set(gca,'YDir','Reverse','Ytick',1:(nobs+nsys),'YtickLabel',obj.CATSstat.names(1:(nobs+nsys)));
                grid on;
                set(gca,'YGrid','off');
                PrettyFigureFormat('FontSize',24);
            end
             ylabel(sprintf('Retarding energy - 18574 (eV)'));
            xlabel('Leverages');
            if nobs>40
                a = get(gca,'YTickLabel');
                set(gca,'YTickLabel',a,'fontsize',8);
            end
            line(2*(nsys+npar)/(nobs+nsys)*[1 1],[1-1 length(obj.CATSstat.leverage)+1],...
                'LineWidth',2,'LineStyle','--','Color',.5*ones(3,1));
            line((nsys+npar)/(nobs+nsys)*[1 1],[1-1 length(obj.CATSstat.leverage)+1],...
                'LineWidth',2,'LineStyle','-.','Color',.7*ones(3,1));
            
            if strcmp(savePlot,'ON')
                if exist('./plots','dir')~=7
                    mkdir plots
                end
                if exist('./plots/cats','dir')~=7
                    mkdir plots/cats
                end
                publish_figurePDF(fig,savefile);
            end
        end
        
        function Samakcats_StandResidualLeverage2D(obj,varargin)
            % Residual versus leverage plot of the fit performed by CATS to identify
            % indluential points.
            %
            % stdrvslplot(s,RainbowFlag,NameFlag) plots the studentized residuals versus
            % the leverages of the fit performed by CATS where s is the output result of CATS
            % ( s = cats(...) ). The Cook's distance is computed for each point and and
            % Cook's distance level cruve are plotted (D = 0.5 and D = 1). By default, RainbowFlag is
            % set to true and colors the points from blue to red to identify begining to end of the
            % spectrum in this plot as well as systematic bins. By default, NameFlag is
            % set to true and indicates names of bins where the Cook's distance is
            % above 1.
            %
            % Copyright 2005-2012 Guillaume MENTION, CEA Saclay
            % $Revision: 1.1$  $Date: 2012/02/08$
            
            p=inputParser;
            p.addParameter('savePlot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            savePlot      = p.Results.savePlot;
            
            % Figure definition
            %            myMainTitle=sprintf('Samak: Standardized Residuals Versus Leverages (CATS)');
            savefile=sprintf('plots/cats/CATS_StandResidualLeverage2D.pdf');
            fig = figure('Units','normalized','Position',[0.1,0.1,0.7,0.45]);%'Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 800 800]);
            %             a=annotation('textbox', [0 0.88 1 0.1], ...
            %                 'String', myMainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');
            %             a.FontSize=20;a.FontWeight='bold';
            % Build uniform / multi-ring indexes
            m         = obj.SO.TBDIS(obj.exclDataStart:obj.exclDataStop,:); %model
            nrings    = numel(obj.SO.MACE_Ba_T); % always gives correct number of pseudo rings
            nqU_used  = size(m,1);             % number of subruns, which are NOT excluded
            exclIndex = sort(reshape(repmat(obj.exclDataStart:obj.exclDataStop,[nrings,1])+[0:nrings-1]'.*obj.SO.nqU,[nqU_used*nrings,1]));
            obj.CATSstat.names = [round(obj.qUdata(exclIndex)-18574,1) ; obj.CATSstat.nsys];
            
            stdrlevplot(obj.CATSstat);
            %            obj.CATSstat.names = [round(obj.DATA(obj.exclDataStart:end,1)-obj.SO.Q_i,1) ; obj.CATSstat.nsys];
            
            
            PrettyFigureFormat('FontSize',18);

            if strcmp(savePlot,'ON')
                if exist('./plots','dir')~=7
                    mkdir plots
                end
                if exist('./plots/cats','dir')~=7
                    mkdir plots/cats
                end
                publish_figurePDF(fig,savefile);
            end
            
        end
        
        function Samakcats_Dffits(obj,varargin)
            % Scaled change in fitted values. DFFITS (with studentized residuals)
            % The deletion influence of ith observation on the predicted or fitted value 
            % can be investigated by using diagnostic by Belsley, Kuh and Welsch
            % DFFITS_i is the number of standard deviations that the fitted value 
            % y_i changes of i-th observation is removed.
            
            p=inputParser;
            p.addParameter('savePlot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            savePlot      = p.Results.savePlot;
            
            % Figure definition
            myMainTitle=sprintf('Samak: DFFITS estimator (scaled change in fitted values, CATS)');
            savefile=sprintf('plots/cats/CATS_Dffits.pdf');
            fig = figure('Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
            a=annotation('textbox', [0 0.9 1 0.1], ...
                'String', myMainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');
            a.FontSize=20;a.FontWeight='bold';
            
            bar(obj.CATSstat.dffit(1:obj.CATSstat.nobs),...
                'FaceColor',rgb('Amethyst'),...
                'EdgeColor',rgb('SteelBlue'),...
                'LineWidth',2);
            colormap(lines(obj.CATSstat.nsys+obj.CATSstat.npar));
            cmap = [.5*ones(1,3); circshift(lines(obj.CATSstat.nsys+obj.CATSstat.npar),-1)];
            hold on;
            for i=1:obj.CATSstat.nsys
                bar(obj.CATSstat.nobs+i,obj.CATSstat.dffit(obj.CATSstat.nobs+i),'FaceColor',cmap(i,:));
            end
            hold off;
            box on;
            axis([0 length(obj.CATSstat.dffit)+1 1.2*max(abs(obj.CATSstat.dffit))*[-1 1]]);
            grid off;
            set(gca,'Xtick',...
                1:length(obj.CATSstat.dffit),...
                'XtickLabel',obj.CATSstat.names(1:(obj.CATSstat.nobs+obj.CATSstat.nsys)));
            xlabel(sprintf('retarding potential - %.1f (eV)',obj.SO.Q_i));
            set(gca,'FontSize',10);
            xtickangle(45);
            ylabel('dffits (standard deviations)');
            grid on;
            PrettyFigureFormat;

            if strcmp(savePlot,'ON')
                if exist('./plots','dir')~=7
                    mkdir plots
                end
                if exist('./plots/cats','dir')~=7
                    mkdir plots/cats
                end
                publish_figurePDF(fig,savefile);
            end
            
        end
        
        function Samakcats_Dfbetas_Rel(obj,varargin)
            % Delete-1 parameters - CATS function deleteparplot
            % Measures the relative influence of ith observation
            % if it is removed from the sample, applied to each
            % relevant fitted parameters
            % DFBETAS which indicates that how much the regression
            % coefficient changes if the ith observation were deleted.
            % Large (in magnitude) value of DFBETASj,i , indicates that ith
            % observation has considerable influence on the
            % jth regression coefficient.
            % Rel --> Such change is measured in terms of standard deviation units.
            
            p=inputParser;
            p.addParameter('savePlot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            savePlot      = p.Results.savePlot;
            
            % Figure definition
            myMainTitle=sprintf('Samak: DFBETAS estimator (scaled change in fitted parameters, CATS)');
            savefile=sprintf('plots/cats/CATS_Dfbetas_Rel.pdf');
            fig = figure('Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
            a=annotation('textbox', [0 0.9 1 0.1], ...
                'String', myMainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');
            a.FontSize=20;a.FontWeight='bold';
            
            XD=1:size(obj.CATSstat.x_i,2);
            YD=1:size(obj.CATSstat.x_i,1);
            diffxi = obj.CATSstat.x_i-repmat(obj.CATSstat.xhat,[1 size(obj.CATSstat.x_i,2)]);
            imagesc(XD,YD,diffxi./repmat(max(abs(diffxi),[],2),[1 size(obj.CATSstat.x_i,2)]));
            set(gca,'Xtick',round(linspace(1,length(XD),max(8,obj.CATSstat.nobs+obj.CATSstat.nsys))));
%           set(gca,'YTick',round(linspace(1,length(YD),max(8,obj.CATSstat.nsys+obj.CATSstat.npar))));
            set(gca,'YTick',round(linspace(1,length(YD),obj.CATSstat.nsys+obj.CATSstat.npar)));
            caxis([-1 1]);
            cm=0:1/9:1;
            my_cm_blue=[1-cm;1-cm;ones(size(cm))]';
            my_cm_red=[ones(size(cm));1-cm;1-cm]';
            colormap([flipud(my_cm_blue);[1 1 1];my_cm_red]);
            colorbar;
            set(gcf,'Color',[1 1 1]);
            colorbar;
            set(gca,'XAxisLocation','bottom');
            grid on;
            set(gca,'YTickLabel',obj.MaskFreeFitLabel)
            %PrettyFigureFormat;
            xtickangle(90);
            if obj.CATSstat.nobs>40
            a = get(gca,'XTickLabel');set(gca,'XTickLabel',a,'fontsize',8);
            end
            
            if strcmp(savePlot,'ON')
                if exist('./plots','dir')~=7
                    mkdir plots
                end
                if exist('./plots/cats','dir')~=7
                    mkdir plots/cats
                end
                publish_figurePDF(fig,savefile);
            end
            
        end
        
        function Samakcats_Dfbetas_Abs(obj,varargin)
            % Delete-1 parameters - CATS function deleteparplot
            % Measures the relative influence of ith observation
            % if it is removed from the sample, applied to each
            % relevant fitted parameters
            % DFBETAS which indicates that how much the regression
            % coefficient changes if the ith observation were deleted.
            % Large (in magnitude) value of DFBETASj,i , indicates that ith
            % observation has considerable influence on the
            % jth regression coefficient.
            % Abs --> Such change is measured in terms of absolute change
            % of the jth regression coefficient.
            
            p=inputParser;
            p.addParameter('savePlot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            savePlot      = p.Results.savePlot;
            
            % Figure definition
            myMainTitle=sprintf('Samak: DFBETAS estimator (change in fitted parameters, CATS)');
            savefile=sprintf('plots/cats/CATS_Dfbetas_Abs.pdf');
            fig = figure('Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
            a=annotation('textbox', [0 0.9 1 0.1], ...
                'String', myMainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');
            a.FontSize=20;a.FontWeight='bold';
            
            diffxi = obj.CATSstat.x_i-repmat(obj.CATSstat.xhat,[1 size(obj.CATSstat.x_i,2)]);
            labelYaxis = obj.MaskFreeFitLabel;
            BestFitCoeff = [obj.SO.mnuSq  obj.SO.Q ...
                obj.SO.BKG_RateSec   (1+obj.SO.normFit) ...
                obj.SO.DTNormGS(1)*100   obj.SO.DTNormES(1)*100 ...
                obj.SO.HTNormGS(1)*100   obj.SO.HTNormES(1)*100 ...
                obj.SO.TTNormGS(1)*100   obj.SO.TTNormES(1)*100];
            
            m         = obj.SO.TBDIS(obj.exclDataStart:obj.exclDataStop,:); % model
            nrings    = numel(obj.SO.MACE_Ba_T); % always gives correct number of pseudo rings
            nqU_used  = size(m,1);               % number of subruns, which are NOT excluded
            exclIndex = sort(reshape(repmat(obj.exclDataStart:obj.exclDataStop,[nrings,1])+[0:nrings-1]'.*obj.SO.nqU,[nqU_used*nrings,1]));
            
            sv=[];counter=0;
            for k=obj.MaskFreeFitPar
                counter=counter+1;
                s(counter)=subplot(numel(obj.MaskFreeFitPar),1,counter);
                %                p=plot([round(obj.DATA(obj.exclDataStart:end,1)-obj.SO.Q_i,1)],diffxi(counter,:)+BestFitCoeff(k),...
                %                    's:','MarkerSize',10,'LineWidth',2,'MarkerFaceColor',rgb('Amethyst'));
                p=plot(round(obj.qUdata(exclIndex)-obj.SO.Q_i,1),diffxi(counter,:)+BestFitCoeff(k),'.:','MarkerSize',20,'LineWidth',2);
                hold on
%                 l=line([ min(round(obj.DATA(obj.exclDataStart:end,1)-obj.SO.Q_i,1)) max(round(obj.DATA(obj.exclDataStart:end,1)-obj.SO.Q_i,1))],...
%                     [BestFitCoeff(k),BestFitCoeff(k)],...
%                     'LineStyle','--','LineWidth',2,'Color',rgb('CadetBlue'));
                l=line([ min(round(obj.qUdata(exclIndex)-obj.SO.Q_i,1)) max(round(obj.qUdata(exclIndex)-obj.SO.Q_i,1))],...
                    [BestFitCoeff(k),BestFitCoeff(k)],...
                    'LineStyle','--','LineWidth',2,'Color',rgb('CadetBlue'));
                hold off
                ylabel(labelYaxis{counter});
                if k==3
                    legend([p l],'delete 1 data - best fit','all data - best fit','Location','southwest') ;
                    legend('boxoff')
                end
                PrettyFigureFormat('FontSize',20); %set(gca,'FontSize',18);
                sv = [sv s(counter)];
            end
            xlabel(sprintf('Retarding energy - %.1f (eV)',obj.SO.Q_i));
            for k=1:4
                linkaxes(sv,'x');
            end
            
            if strcmp(savePlot,'ON')
                if exist('./plots','dir')~=7
                    mkdir plots
                end
                if exist('./plots/cats','dir')~=7
                    mkdir plots/cats
                end
                export_fig(fig,savefile);
            end
            
            %% plot only mnu
            fnew = figure('Units','normalized','Position',[0.1,0.1,0.6,0.4]);
            counter = 1; k=1;
            l=line([min(round(obj.qUdata(exclIndex)-18574,1))-5 5+max(round(obj.qUdata(exclIndex)-18574,1))],...
                [BestFitCoeff(k),BestFitCoeff(k)],...
                'LineStyle','-','LineWidth',1.5,'Color',rgb('Black'));
            hold on;
            p=plot(round(obj.qUdata(exclIndex)-18574,1),diffxi(counter,:)+BestFitCoeff(k),...
                '.:','MarkerSize',26,'LineWidth',2,'Color',rgb('DodgerBlue'));
            hold off 
            ylabel(sprintf('{\\itm}_\\nu^{ 2} (eV^{ 2})'));
            xlabel('Retarding energy - 18574 (eV)');
            xlim([min(round(obj.qUdata(exclIndex)-18574,1))-4 4+max(round(obj.qUdata(exclIndex)-18574,1)) ])
            ylim([min(diffxi(counter,:)+BestFitCoeff(k))-0.05,max(diffxi(counter,:)+BestFitCoeff(k))+0.05]);
            PrettyFigureFormat('FontSize',19);
           % set(gca,'FontSize',18);
            grid off;
            leg = legend([l,p],'Best-fit with all scan steps in analysis range',sprintf('Fit without 1 scan step'));
            legend boxoff
            leg.FontSize = get(gca,'FontSize')+4;
            leg.Location = 'southeast';
            if strcmp(savePlot,'ON')
                export_fig(fnew,strrep(savefile,'.pdf','_mnu.pdf'));
            end
            %% plot only E0 
             fnew2 = figure('Units','normalized','Position',[0.1,0.1,0.6,0.4]);
            counter = 2; k=2;
            l=line([min(round(obj.qUdata(exclIndex)-18574,1))-5 5+max(round(obj.qUdata(exclIndex)-18574,1))],...
                [0,0],...
                'LineStyle','-','LineWidth',1.5,'Color',rgb('Black'));
            hold on;
            p=plot(round(obj.qUdata(exclIndex)-18574,1),diffxi(counter,:),...
                '.:','MarkerSize',24,'LineWidth',2,'Color',rgb('DodgerBlue'));
            hold off 
            ylabel(sprintf('{\\itE}_0^{fit} - {\\itE}_0^{bf} (eV)'));
            xlabel('Retarding energy - 18574 (eV)');
            xlim([min(round(obj.qUdata(exclIndex)-18574,1))-4 4+max(round(obj.qUdata(exclIndex)-18574,1)) ])
           ylim([min(diffxi(counter,:))-0.01,max(diffxi(counter,:))+0.01]);
          
            PrettyFigureFormat('FontSize',19);
          
           leg=   legend([l,p],'Best-fit with all scan steps in analysis range',sprintf('Fit without 1 scan step'));
            legend boxoff
              leg.FontSize = get(gca,'FontSize')+4;
                 leg.Location = 'southeast';
            if strcmp(savePlot,'ON')
                export_fig(fnew2,strrep(savefile,'.pdf','_E0.pdf'));
            end
        end
        
        function Samakcats_DMSE(obj,varargin)
            % Change in Mean Squared Errors (chi2min/ndof).
            % dmseplot(s) plots the mse change in the fit performed
            % when observation are removed in turn
            % Copyright 2005-2011 Guillaume MENTION, CEA Saclay
            % $Revision: 0.99$  $Date: 2011/11/15$
            p=inputParser;
            p.addParameter('savePlot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            savePlot      = p.Results.savePlot;
            
            % Figure definition
            myMainTitle=sprintf('Samak: DMSE estimator (CATS)');
            savefile=sprintf('plots/cats/CATS_DMSE.pdf');
            fig = figure('Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 600 1000]);
            a=annotation('textbox', [0 0.9 1 0.1], ...
                'String', myMainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');
            a.FontSize=20;a.FontWeight='bold';
            
            nobs = obj.CATSstat.nobs;
            nsys = obj.CATSstat.nsys;
            npar = obj.CATSstat.npar;
            
            h1 = barh(obj.CATSstat.s2_i(1:nobs),...
                'FaceColor',rgb('Amethyst'),...
                'EdgeColor',rgb('SteelBlue'),...
                'LineWidth',2);
            set(h1,'BaseValue',obj.CATSstat.mse);
            colormap(lines(nsys));
            cmap = [.5*ones(1,3); circshift(lines(nsys),-1)];
            hold on;
            for i=1:nsys
                barh(nobs+i,obj.CATSstat.s2_i(nobs+i),'FaceColor',cmap(i,:));
            end
            hold off;
            box on;
            axis([min(1,min(obj.CATSstat.s2_i)/1.1) max(obj.CATSstat.s2_i)*1.1 0 length(obj.CATSstat.s2_i)+1 ]);
            grid off;
            % set(gca,'Xtick',round(linspace(1,length(obj.CATSstat.s2_i),min(10,nobs+nsys))));
            if isempty(obj.CATSstat.names)
                set(gca,'YDir','Reverse','Ytick',1:(nobs+nsys),'YtickLabel',1:(nobs+nsys));
                grid on;
                set(gca,'YGrid','off');
                PrettyFigureFormat;
            else
                set(gca,'YDir','Reverse','Ytick',1:(nobs+nsys),'YtickLabel',obj.CATSstat.names(1:(nobs+nsys)));
                grid on;
                set(gca,'YGrid','off');
                PrettyFigureFormat;
            end
            
            xlabel('change in mse (  \chi^2_{min} / ndof )');
            
            line([1 1]*obj.CATSstat.mse,[0 length(obj.CATSstat.s2_i)+1],'LineWidth',2,'LineStyle','-',...
                'Color',.5*ones(3,1));
            line([1 1]*obj.CATSstat.mse*obj.CATSstat.ndof/(obj.CATSstat.ndof-1),...
                [0 length(obj.CATSstat.s2_i)+1],'LineWidth',2,'LineStyle','--',...
                'Color',.5*ones(3,1));
            line([1 1]*(obj.CATSstat.mse*obj.CATSstat.ndof-1)/(obj.CATSstat.ndof-1),...
                [0 length(obj.CATSstat.s2_i)+1],'LineWidth',2,'LineStyle',':',...
                'Color',.5*ones(3,1));
            
            if (max(obj.CATSstat.s2_i>5))
                set(gca,'XLim',[0 5]);
            end
            
            ylabel(sprintf('retarding potential - %.1f (eV)',obj.SO.Q_i));
            PrettyFigureFormat;
            
            if strcmp(savePlot,'ON')
                if exist('./plots','dir')~=7
                    mkdir plots
                end
                if exist('./plots/cats','dir')~=7
                    mkdir plots/cats
                end
                publish_figurePDF(fig,savefile);
            end
            
        end
        
        function Samakcats_Delete1Variance(obj,varargin)
            % Delete-1 variance
            % The delete-1 variance shows how the mean squared error
            % changes when an observation is removed from the data set.
            p=inputParser;
            p.addParameter('savePlot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            savePlot      = p.Results.savePlot;
            
            % Figure definition
            myMainTitle=sprintf('Samak: Delete-1 Variance (CATS)');
            savefile=sprintf('plots/cats/CATS_Delete1Variance.pdf');
            fig = figure('Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 1400 400]);
            a=annotation('textbox', [0 0.9 1 0.1], ...
                'String', myMainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');
            a.FontSize=20;a.FontWeight='bold';
            
            plot(obj.CATSstat.s2_i,'s:','MarkerSize',10,'MarkerFaceColor',rgb('Amethyst'),'LineWidth',2);
            axis([0 length(obj.CATSstat.s2_i)+1 obj.CATSstat.mse-6/sqrt(obj.CATSstat.nobs-obj.CATSstat.npar-1) obj.CATSstat.mse+6/sqrt(obj.CATSstat.nobs-obj.CATSstat.npar-1)]);
            line([1-1 obj.CATSstat.nobs+obj.CATSstat.nsys+1],[obj.CATSstat.mse obj.CATSstat.mse],'LineWidth',2,'Color',.5*ones(3,1));
            grid off;
            set(gca,'Xtick',round(linspace(1,length(obj.CATSstat.s2_i),min(8,obj.CATSstat.nobs+obj.CATSstat.nsys))));
            ylabel('\sigma_{i}');
            if obj.CATSstat.nobs+obj.CATSstat.nsys-obj.CATSstat.npar>5
                line([1-1 length(obj.CATSstat.s2_i)+1],(obj.CATSstat.mse+2/sqrt(obj.CATSstat.nobs-obj.CATSstat.npar-1))*[1 1],...
                    'LineWidth',2,'LineStyle','--','Color',.5*ones(3,1));
                line([1-1 length(obj.CATSstat.s2_i)+1],(obj.CATSstat.mse-2/sqrt(obj.CATSstat.nobs-obj.CATSstat.npar-1))*[1 1],...
                    'LineWidth',2,'LineStyle','--','Color',.5*ones(3,1));
            end
            if isempty(obj.CATSstat.names)
                set(gca,'Xtick',1:(obj.CATSstat.nobs+obj.CATSstat.nsys),'XtickLabel',1:(obj.CATSstat.nobs+obj.CATSstat.nsys));
                grid on;
                xtickangle(gca,45)
                set(gca,'XGrid','off');
                xlabel(sprintf('retarding potential - %.1f (eV)',obj.SO.Q_i));
                PrettyFigureFormat;
            else
                set(gca,'Xtick',1:(obj.CATSstat.nobs+obj.CATSstat.nsys),'XtickLabel',obj.CATSstat.names(1:(obj.CATSstat.nobs+obj.CATSstat.nsys)));
                grid on;
                set(gca,'XGrid','off');
                xtickangle(gca,45)
                xlabel(sprintf('retarding potential - %.1f (eV)',obj.SO.Q_i));
                PrettyFigureFormat;
            end
            
            if strcmp(savePlot,'ON')
                if exist('./plots','dir')~=7
                    mkdir plots
                end
                if exist('./plots/cats','dir')~=7
                    mkdir plots/cats
                end
                publish_figurePDF(fig,savefile);
            end
            
          end
         
    end % methods
end % class
