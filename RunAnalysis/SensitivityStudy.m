classdef SensitivityStudy < handle
    % ----------------------------------------------------------
    % Class for Neutrino Mass Sensitivity Studies
    % L.Schlueter (November 2018)
    %-----------------S-------------------------------------------
   
    properties  (Access=public)
        % Sensitivity properties
        mNu90;                % neutrino mass squared sensitivity 90% C.L.
        mNumin;               % neutrino mass for which chi2 was smallest (should be zero)
        FitObj;               % Model Object for FIT 
        SimObj;               % Model Object for simulated sample spectra
        MultimNu90;           % sensitivity of single effects: stat + TC, stat + FSD,....
        MultimNu90Stack;      % sensitivities N+1: stat, stat+tc, stat + TC + FSD, stat +.....
        nSys = 4;             % number of stacked sys
        nSingleSys = 9;       % numbe of single sys
        
        % properties for TBD object
        TimeSec;               % measurement time
        range;                 % energy range below endpoint
        FPD_MeanEff;           % Focal plane detector efficiency
        MACE_Ba_T;             % magnetic field analyzing plane
        WGTS_B_T;              % source magnetic field
        MACE_Bmax_T;           % magnetic field of pinch magnet
        Q_i;                   % endpoint
        WGTS_CD_MolPerCm2;     % column density
        ELossFlag;
        recomputeRF;
        
        %TD Options
        TD;                    % measurement time distribution (name)
        TD_NuMassFactor;       % only needed when MTCcreator is used
        TD_BkgFraction;
        TD_BqU;     
        
        % Background Options
        BKG_RateSec;           % background for whole detector
        AnchorBkg6G;           % if BKG_RateSec and FPD_ROIlow not specified: BKG is scaled to Anchor6G level
        BkgCM;                 % Optional: Background covariance matrix (works also for chi2 = 'chi2Stat')
        
        % FPD
        FPD_SegOffPixList;     % List of Pixels for Uniform Analysis
        FPD_ROIlow;            % if BKG_RateSec not specified: BKG is scaled to First Tritium level

        % properties for Asimov Data Scan Method
        mNuStart;              % start value for scan (neutrino mass)
        mNuStop;               % stop value for scan (neutrino mass)
        ScanPrcsn;             % tolerated delta chi2
        Scan;                  % ON = Do fit with neutrino mass vector - OFF = only 1 fit, fast, but sometimes fails
        nFitMax=20;            % maximal number of fits
        ScanResults;           % struct with scan Results
        RecomputeFlag; 
        
        % Systematics Settings
        chi2;                  % 'chi2Stat' or 'chi2CM' or 'chi2CMShape'
        SysEffect;             % systematic effect used in covariance matrix
        SysBudget;             % different predefined set of systematic budgets (in ref_CovarianceMatrix_NominalKATRIN)
        MultiCM;               % Covariance Matrix
        MultiCMFrac;           % Fractional Covariance Matrix
        MultiCMShape;          % Shape Only Covariance Matrix
        % single sys: used in NuMassScan_SystematicsLoop
        SysEffectsAll = {'TC','FSD','RF','all','RF_EL','RF_BF','RF_RX','RF_BFRX'};  
        % N plus 1: used in NuMassScan_SensitivityNominal_Nplus1Systematics_Loop
        SysEffectsStack = {'TC',...%TC only
                             struct('TC','ON','FSD','ON'),...%TC + FSD
                             'all'}; %TC + RF + FSD = everything   % N+1: stat, stat+tc, stat + TC + FSD, stat +.....
                  
        % Display Options
        PlotFit;               % plot neutrino mass scan
        PlotCM;                % display covariance matrix used for fit
    end
    methods
        function obj = SensitivityStudy(varargin) % constructor
            fprintf('-----------------Start SensitivityStudy contructor ----------------------\n');
            p = inputParser;
            p.addParameter('TimeSec',3*365*24*60*60,@(x)isfloat(x) && x>0);
            p.addParameter('range',30,@(x)isfloat(x));
            p.addParameter('TD','DR30',@(x)ischar(x)); %only for non SensitivityBa....
            p.addParameter('TD_NuMassFactor',0.60,@(x)isfloat(x)); 
            p.addParameter('TD_BkgFraction',0.15,@(x)isfloat(x)); 
            p.addParameter('TD_BqU',[5 10 20],@(x)all(isfloat(x))); 
        
            % Background
            p.addParameter('BKG_RateSec','',@(x)isfloat(x) || isempty(x)); %if empty, filled according to Ba
            p.addParameter('AnchorBkg6G','');

            % FPD
            p.addParameter('FPD_SegOffPixList',1:119,@(x)all(isfloat(x)) && all(x)>0);
            p.addParameter('FPD_MeanEff',0.95,@(x)isfloat(x) && x>=0);
            p.addParameter('FPD_ROIlow',14,@(x)isfloat(x) && x>=0);

            % WGTSMACE Model Parameter
            p.addParameter('WGTS_CD_MolPerCm2',5e17,@(x)isfloat(x) && x>=0);
            p.addParameter('MACE_Ba_T',3e-04,@(x)isfloat(x) && x>=0);
            p.addParameter('WGTS_B_T',3.6,@(x)isfloat(x) && x>=0);
            p.addParameter('MACE_Bmax_T',6,@(x)isfloat(x) && x>=0);
            p.addParameter('ELossFlag','Abdurashitov',@(x)ismember(x,{'Aseev','Abdurashitov','KatrinD2'}));
            p.addParameter('recomputeRF','OFF',@(x)ismember(x,{'ON','OFF'}));

            %TBD Model Paramater
            p.addParameter('mNuStart',0.2,@(x)all(isfloat(x))); % Start mNu for neutrino mass scan (eV)
            p.addParameter('mNuStop',0.6,@(x)all(isfloat(x))); % Stop mNu for neutrino mass scan (eV)
            p.addParameter('ScanPrcsn',0.02,@(x)all(isfloat(x))); % Scan Precision (Acceptable Delta chi2)
            p.addParameter('Q_i',18575,@(x)all(isfloat(x))); % To give Bias of Work Function (in sample spectra)
           
            %Fit Settings
            p.addParameter('chi2','chi2CM',@(x)ismember(x,{'chi2Stat','chi2CM','chi2CMShape'}));
            p.addParameter('SysEffect','FSD'); %1 sys effect e.g. FSD or struct for more than 1
            p.addParameter('PlotFit','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            
            %Systematics settings
            p.addParameter('PlotCM','OFF',@(x)ismember(x,{'ON','OFF'})); %produces a set of plots and a text file with information
            p.addParameter('SysBudget','03',@(x)ischar(x)); % defines the covmat reference file
            p.addParameter('Scan','ON',@(x)ismember(x,{'ON','OFF'}));   % Off = no scan, 1 asimov fit with free nu-mass
            p.addParameter('BkgCM','OFF',@(x)ismember(x,{'ON','OFF'})); % Background Covariance matrix
            p.parse(varargin{:});
            
            obj.TD              = p.Results.TD;
            obj.TimeSec         = p.Results.TimeSec;
            obj.range           = p.Results.range;
            obj.TD_NuMassFactor = p.Results.TD_NuMassFactor;
            obj.TD_BkgFraction  = p.Results.TD_BkgFraction;
            obj.TD_BqU          = p.Results.TD_BqU;
            
            % WGTS MACE          
            obj.WGTS_CD_MolPerCm2 = p.Results.WGTS_CD_MolPerCm2;
            obj.MACE_Ba_T         = p.Results.MACE_Ba_T;
            obj.WGTS_B_T          = p.Results.WGTS_B_T;
            obj.MACE_Bmax_T       = p.Results.MACE_Bmax_T;
            obj.ELossFlag         = p.Results.ELossFlag;
            obj.recomputeRF       = p.Results.recomputeRF;
            obj.Q_i               = p.Results.Q_i;
            % Systematic settings
            obj.chi2         = p.Results.chi2;
            obj.SysEffect    = p.Results.SysEffect;
            obj.SysBudget    = p.Results.SysBudget;
            obj.BkgCM        = p.Results.BkgCM;
            obj.PlotCM       = p.Results.PlotCM;
            % Scan Options
            obj.mNuStart     = p.Results.mNuStart;
            obj.mNuStop      = p.Results.mNuStop;
            obj.ScanPrcsn    = p.Results.ScanPrcsn;
            obj.Scan         = p.Results.Scan;
            obj.PlotFit      = p.Results.PlotFit;
            obj.RecomputeFlag= p.Results.RecomputeFlag;
            % Background / FPD
            obj.FPD_SegOffPixList = p.Results.FPD_SegOffPixList;
            obj.FPD_MeanEff  = p.Results.FPD_MeanEff;
            obj.FPD_ROIlow   = p.Results.FPD_ROIlow;
            obj.AnchorBkg6G  = p.Results.AnchorBkg6G;
            obj.BKG_RateSec  = p.Results.BKG_RateSec;
            
            if strcmp(obj.chi2,'chi2Stat')
                obj.SysEffect = '';
            end
            % Compute Background
            if isempty(obj.BKG_RateSec)
                obj.BKG_RateSec = GetBackground('MACE_Ba_T',obj.MACE_Ba_T,...
                    'WGTS_B_T',obj.WGTS_B_T,'FPD_ROIlow',obj.FPD_ROIlow,...
                    'AnchorBkg6G',obj.AnchorBkg6G, ...
                    'FPD_SegOffPixList',obj.FPD_SegOffPixList);
            end
            
            % Init TD if neccessary
            obj.SetTD;
            % Compute Models
            obj.InitializeModels;
            % Compute Covariance Matrix
             obj.ComputeCM;
           
            fprintf('-----------------End SensitivityStudy contructor ----------------------\n');
        end
    end
    methods % Initialization methods
        function InitializeModels(obj)
            %Simulation Model
            ModelArg = {...
                'FPD_SegOffPixList',obj.FPD_SegOffPixList,...
                'WGTS_CD_MolPerCm2',obj.WGTS_CD_MolPerCm2,...
                'BKG_RateAllFPDSec',obj.BKG_RateSec,...
                'MACE_Ba_T',obj.MACE_Ba_T,...
                'WGTS_B_T',obj.WGTS_B_T,...
                'MACE_Bmax_T',obj.MACE_Bmax_T,...
                'TimeSec',obj.TimeSec,...
                'TD',obj.TD,...
                'mnuSq_i',0,...
                'Q_i',obj.Q_i,...
                'FPD_MeanEff',obj.FPD_MeanEff,...
                'ELossFlag',obj.ELossFlag,...
                'recomputeRF',obj.recomputeRF,...
                'FPD_ROIlow',obj.FPD_ROIlow};
            obj.SimObj = ref_TBD_NominalKATRIN(ModelArg{:});
            obj.SimObj.ComputeTBDDS;
            obj.SimObj.ComputeTBDIS;
            %Fit Model
            obj.FitObj = ref_TBD_NominalKATRIN(ModelArg{:});
            obj.FitObj.ComputeTBDDS; obj.FitObj.ComputeTBDIS;
        end
        function ComputeCM(obj)
            % Compute/ Load Covariance Matrices
            if ~strcmp(obj.chi2,'chi2Stat')
                [~, obj.MultiCM, obj.MultiCMFrac, obj.MultiCMShape, ~] = ...
                    ref_CovarianceMatrix_NominalKATRIN('RecomputeFlag','OFF','ModelObj',obj.FitObj,...
                    'SysEffect',obj.SysEffect,'PlotCM',obj.PlotCM,'SysBudget',obj.SysBudget);
                
                %Init Fit Model again
                obj.InitializeModels;
                
                % Optional: Background Covariance Matrix
                if strcmp(obj.BkgCM,'ON')
                    [BkgCM, ~, BkgCMShape, ~] = ComputeCM_Background('StudyObject',obj.FitObj,...
                        'plotFit',obj.PlotCM,'nTrials',1000);
                    obj.MultiCM      = obj.MultiCM + BkgCM;
                    obj.MultiCMShape = obj.MultiCMShape + BkgCMShape;
                end
                
            elseif strcmp(obj.chi2,'chi2Stat') % statistics only
                if strcmp(obj.BkgCM,'OFF')
                    obj.MultiCM      = zeros(obj.FitObj.nqU);
                    obj.MultiCMFrac  = zeros(obj.FitObj.nqU);
                    obj.MultiCMShape = zeros(obj.FitObj.nqU);
                elseif strcmp(obj.BkgCM,'ON')
                    [BkgCM, ~, BkgCMShape, ~] = ComputeCM_Background('StudyObject',obj.FitObj,...
                        'plotFit',obj.PlotCM,'nTrials',1000);
                    obj.MultiCM      = BkgCM;
                    obj.MultiCMFrac  = zeros(obj.FitObj.nqU);
                    obj.MultiCMShape = BkgCMShape;
                end
            end
        end
    end
    methods % Sensitivity Computation
        function [mnuSq_i_Fit, par, err, chi2min, dof,mNu90,mNumin] = NuMassScan(obj)
            obj.ScanResults = struct('mnuSq_i_Fit',zeros(obj.nFitMax,1),...% neutrino mass squared scan vector
                'par',zeros(10,obj.nFitMax),... % Fit Parameter
                'err',zeros(10,obj.nFitMax),... % Error on Fit Parameter
                'chi2min',zeros(obj.nFitMax,1),...
                'dof',zeros(obj.nFitMax,1));
            
            obj.ScanResults.mnuSq_i_Fit(2)           = obj.mNuStart^2;    % start value for scan
            obj.ScanResults.mnuSq_i_Fit(obj.nFitMax) = obj.mNuStop^2;   % stop value for scan
            F                                        = cell(length(obj.nFitMax),1); % Fit Class Objects
            
            
            Data = [obj.SimObj.qU, obj.SimObj.TBDIS, sqrt(obj.SimObj.TBDIS)];
            if strcmp(obj.Scan,'ON')
                % Test of mNu=0 and Stop Value
                % mNu=0 should yield chi2=0
                % Stopvalue should yield chi2>2.7
                TestVal = [1,obj.nFitMax];
                for i=1:numel(TestVal)
                    obj.FitObj.mnuSq_i = obj.ScanResults.mnuSq_i_Fit(TestVal(i)); %zero
                    F{TestVal(i)} = FITC('SO',obj.FitObj,'DATA',Data,'fitter','minuit',...
                        'chi2name',obj.chi2,...
                        'COVMAT', obj.MultiCM+diag(obj.SimObj.TBDIS),...
                        'COVMATFrac', obj.MultiCMFrac+diag(1./obj.SimObj.TBDIS),...
                        'COVMATShape', obj.MultiCMShape,...
                        'fixPar','1 5 6 7 8 9 10',...
                        'exclDataStart',1); % at the moment: take all qU
                    obj.ScanResults.par(:,TestVal(i))    = F{TestVal(i)}.RESULTS{1};
                    obj.ScanResults.par(1,TestVal(i))    = obj.ScanResults.par(1,TestVal(i))+obj.FitObj.mnuSq_i;
                    obj.ScanResults.err(:,TestVal(i))   = F{TestVal(i)}.RESULTS{2};
                    obj.ScanResults.chi2min(TestVal(i)) = F{TestVal(i)}.RESULTS{3};
                    obj.ScanResults.dof(TestVal(i))     = F{TestVal(i)}.RESULTS{5};
                    
                    if  obj.ScanResults.chi2min(1)>1e-05 %if chi2 isnt 0, problem!
                        fprintf(2,'Error: \chi^2 is not 0 for neutrino mass 0! \n');
                        return
                    elseif TestVal(i)==obj.nFitMax && obj.ScanResults.chi2min(obj.nFitMax)<2.7
                        fprintf(2,'Error: Neutrino Mass Scan Range is too small. Enlarge mNuStop! \n');
                        return
                    end
                end
                
                % Do Scan
                for i=2:obj.nFitMax
                    obj.FitObj.mnuSq_i = obj.ScanResults.mnuSq_i_Fit(i);
                    F{i} = FITC('SO',obj.FitObj,'DATA',Data,'fitter','minuit',...
                        'chi2name',obj.chi2,...
                        'COVMAT', obj.MultiCM+diag(obj.SimObj.TBDIS),'COVMATFrac', obj.MultiCMFrac+diag(1./obj.SimObj.TBDIS),...
                        'COVMATShape', obj.MultiCMShape+diag(obj.SimObj.TBDIS),...
                        'fixPar','1 5 6 7 8 9 10',...
                        'exclDataStart',1);
                    obj.ScanResults.par(:,i)   = F{i}.RESULTS{1};
                    obj.ScanResults.par(1,i)   = obj.ScanResults.par(1,i)+obj.FitObj.mnuSq_i;
                    obj.ScanResults.err(:,i)   = F{i}.RESULTS{2};
                    obj.ScanResults.chi2min(i) = F{i}.RESULTS{3};
                    obj.ScanResults.dof(i)     = F{i}.RESULTS{5};
                    if  abs(obj.ScanResults.chi2min(i)-2.7)<obj.ScanPrcsn % exit condition
                        fprintf('90 C.L. limit reached \n');
                        obj.ScanResults.mnuSq_i_Fit(i+1:end)=NaN; % set all not fitted to NaN
                        obj.ScanResults.chi2min(i+1:end) = NaN;
                        obj.ScanResults.dof(i+1:end)     = NaN;
                        obj.ScanResults.par(:,i+1:end)   = NaN;
                        obj.ScanResults.err(:,i+1:end)   = NaN;
                        fprintf('Exist Scan, 90%% C.L. reached \n');
                        break
                    elseif obj.ScanResults.chi2min(i)<2.7
                        mNuSqLarger= min(obj.ScanResults.mnuSq_i_Fit(obj.ScanResults.mnuSq_i_Fit>obj.ScanResults.mnuSq_i_Fit(i))); % find next higher mNuSq
                        obj.ScanResults.mnuSq_i_Fit(i+1) =  0.5*(obj.ScanResults.mnuSq_i_Fit(i)+mNuSqLarger); % find middle between current and next higher one
                    elseif obj.ScanResults.chi2min(i)>2.7
                        mNuSqSmaller= max(obj.ScanResults.mnuSq_i_Fit(obj.ScanResults.mnuSq_i_Fit<obj.ScanResults.mnuSq_i_Fit(i))); % find next smaller mNuSq
                        obj.ScanResults.mnuSq_i_Fit(i+1) =  0.5*(obj.ScanResults.mnuSq_i_Fit(i)+mNuSqSmaller);  % find middle between current and next smaller one
                    end
                end
                
                % Compute Sensitivity
                obj.mNu90 =interp1(obj.ScanResults.chi2min(~isnan(obj.ScanResults.chi2min)),obj.ScanResults.mnuSq_i_Fit(~isnan(obj.ScanResults.chi2min)),min(obj.ScanResults.chi2min(~isnan(obj.ScanResults.chi2min)))+2.7,'spline'); % sensitivity on mnu^2 (90% C.L.)
                obj.mNumin = interp1(obj.ScanResults.chi2min(~isnan(obj.ScanResults.chi2min)),obj.ScanResults.mnuSq_i_Fit(~isnan(obj.ScanResults.chi2min)),min(obj.ScanResults.chi2min(~isnan(obj.ScanResults.chi2min))),'spline');   % mass with minimal chi2
                if strcmp(obj.PlotFit,'ON')
                    obj.PlotScan;
                end
            elseif strcmp(obj.Scan,'OFF')
                F = FITC('SO',obj.FitObj,'DATA',Data,'fitter','minuit',...
                    'chi2name',obj.chi2,...
                    'COVMAT', obj.MultiCM+diag(obj.SimObj.TBDIS),'COVMATFrac', obj.MultiCMFrac+diag(1./obj.SimObj.TBDIS),...
                    'COVMATShape', obj.MultiCMShape+diag(obj.SimObj.TBDIS),...
                    'fixPar','5 6 7 8 9 10',...
                    'exclDataStart',1);
                obj.ScanResults.par         = squeeze(repmat(F.RESULTS{1},1,1,obj.nFitMax));
                obj.ScanResults.par(1,:)    = obj.ScanResults.par(1,:)+obj.FitObj.mnuSq_i;
                obj.ScanResults.err         = squeeze(repmat(F.RESULTS{2},1,1,obj.nFitMax));
                obj.ScanResults.chi2min     = squeeze(repmat(F.RESULTS{3},1,1,obj.nFitMax));
                obj.ScanResults.dof         = squeeze(repmat(F.RESULTS{5},1,1,obj.nFitMax));
                obj.ScanResults.mnuSq_i_Fit = zeros(obj.nFitMax,1);
                obj.mNu90 = obj.ScanResults.err(1)*1.64; % sensitivity on mnu^2 (90% C.L.) central, two-sided interval
                obj.mNumin = obj.ScanResults.par(1);
            end
            fprintf(2,'Neutrino Mass Sensitivity =  %.1f meV (90%% C.L.) \n', sqrt(obj.mNu90)*1e3);
            % output
            mnuSq_i_Fit = obj.ScanResults.mnuSq_i_Fit;
            par =  obj.ScanResults.par;
            err=  obj.ScanResults.err;
            chi2min= obj.ScanResults.chi2min;
            dof = obj.ScanResults.dof;
            mNu90 = obj.mNu90;
            mNumin = obj.mNumin;
        end
        function NuMassScan_SystematicsLoop(obj)
            % --------------------------------------------------------------------------------
            % Compute Sensitivity for stat + 1sys
            % loop over all systematic effects:
            % Theoretical Corrections, Final States Distribution, Response
            % Function, all three, Response Function (energy loss only), Response Function (magnetic fields),
            %Response Function (rhod + cross section), Response Function
            %(rhod + cross section + magnetic fields)
            % save results for later plotting
            % -----------------------------------------------------------------------
            % labeling
            TDlabel = strrep(obj.TD,sprintf('_Ba%.0fG',obj.MACE_Ba_T*1e4),'');
            TDlabel = strrep(TDlabel,'Sensitivity_','');
            if strcmp(obj.SysBudget,'08') || strcmp(obj.SysBudget,'09')
                SysBudget_local = [obj.SysBudget,'_NoELoss'];
            else
                SysBudget_local = obj.SysBudget;
            end
            save_name = sprintf('./results/%s_%0.0fd_SensitivityNominal_ResultsNuMassScan_MTD-%s_Ba%.0fG_Bs%.2fT_FPDeff%.2f_BKG-%.0fmcps_Systematics_%s.mat',...
                SysBudget_local,round(obj.TimeSec/(86400)),TDlabel,obj.MACE_Ba_T*1e4,obj.WGTS_B_T,obj.FPD_MeanEff*obj.FitObj.FPD_Coverage,obj.BKG_RateSec*1e3,obj.ELossFlag);
            % if doesnt exist but RecomputeFlag OFF -> try without ELoss label
            if strcmp(obj.RecomputeFlag,'OFF')
                if~exist(save_name,'file') % if file doesnt exist, try without ELossFlag
                    save_name = strrep(save_name,sprintf('_%s',obj.ELossFlag),'');
                    if ~strcmp(obj.ELossFlag,'Aseev')
                        fprintf(2,'Warning: Loading Results with Aseev ELoss \n');
                    end
                end
            end
            if exist(save_name,'file') && strcmp(obj.RecomputeFlag,'OFF')
                fprintf('Loading Results from File: %s \n', save_name);
                d = load(save_name);
                obj.MultimNu90.Stat    = d.mNu90(1);
                obj.MultimNu90.TC      = d.mNu90(2);
                obj.MultimNu90.FSD     = d.mNu90(3);
                obj.MultimNu90.RF      = d.mNu90(4);
                obj.MultimNu90.all     = d.mNu90(5);
                obj.MultimNu90.RF_EL   = d.mNu90(6);
                obj.MultimNu90.RF_BF   = d.mNu90(7);
                obj.MultimNu90.RF_RX   = d.mNu90(8);
                obj.MultimNu90.RF_BFRX = d.mNu90(9);
            elseif  strcmp(obj.RecomputeFlag,'OFF')
                fprintf('File: %s \ndoes not exist.Do you want to recompute? \n', save_name);
                return;
            elseif  strcmp(obj.RecomputeFlag,'ON')
                %Init: gather results
                chi2minScan  = zeros(numel(obj.SysEffectsAll)+1,obj.nFitMax); % chi2min distribution
                mNu_Fit      = zeros(numel(obj.SysEffectsAll)+1,obj.nFitMax);
                mNu90        = zeros(numel(obj.SysEffectsAll)+1,1);                  % sensitivity on mnu^2 (90% C.L.)
                mNumin       = zeros(numel(obj.SysEffectsAll)+1,1);                  % mass with minimal chi2
                parScan      = zeros(numel(obj.SysEffectsAll)+1,10,obj.nFitMax);
                errScan      = zeros(numel(obj.SysEffectsAll)+1,10,obj.nFitMax);
                
                %% Do Scan
                chi2prev      = obj.chi2;
                SysEffectprev = obj.SysEffect;
                %stat
                obj.chi2 = 'chi2Stat';
                obj.ComputeCM; % Reset CM to 0
                obj.NuMassScan;
                [mNu_Fit(1,:), parScan(1,:,:), errScan(1,:,:), chi2minScan(1,:,:),~,mNu90(1),mNumin(1)] = ...
                    obj.NuMassScan;
                
                % Systematics
                obj.chi2 = 'chi2CM';
                for i=1:numel(obj.SysEffectsAll)
                    obj.SysEffect = obj.SysEffectsAll{i};
                    obj.ComputeCM;
                    [mNu_Fit(i+1,:), parScan(i+1,:,:), errScan(i+1,:,:), chi2minScan(i+1,:,:), dof ,mNu90(i+1),mNumin(i+1)] = ...
                        obj.NuMassScan;
                end
                           
                % save
                TD = obj.TD; ScanPrcsn = obj.ScanPrcsn;TimeSec = obj.TimeSec;
                WGTS_B_T = obj.WGTS_B_T; MACE_Ba_T = obj.MACE_Ba_T; BKG_RateSec = obj.BKG_RateSec;
                mySysEffects = obj.SysEffectsAll;
                save(save_name,'parScan','errScan','chi2minScan','mNu_Fit','ScanPrcsn','mNu90', 'mNumin','TD','mySysEffects','dof','WGTS_B_T','MACE_Ba_T','BKG_RateSec','TimeSec');
                obj.chi2      = chi2prev;
                obj.SysEffect = SysEffectprev;
                obj.MultimNu90.Stat    = mNu90(1);
                obj.MultimNu90.TC      = mNu90(2);
                obj.MultimNu90.FSD     = mNu90(3);
                obj.MultimNu90.RF      = mNu90(4);
                obj.MultimNu90.all     = mNu90(5);
                obj.MultimNu90.RF_EL   = mNu90(6);
                obj.MultimNu90.RF_BF   = mNu90(7);
                obj.MultimNu90.RF_RX   = mNu90(8);
                obj.MultimNu90.RF_BFRX = mNu90(9);
            end
        end
        function NuMassScan_SensitivityNominal_Nplus1Systematics_Loop(obj)
            TDlabel = strrep(obj.TD,sprintf('_Ba%.0fG',obj.MACE_Ba_T*1e4),''); %
            TDlabel = strrep(TDlabel,'Sensitivity_','');
            if strcmp(obj.SysBudget,'08') || strcmp(obj.SysBudget,'09')
                SysBudget_local = [obj.SysBudget,'_NoELoss'];
            else
                SysBudget_local = obj.SysBudget;
            end
            save_name = sprintf('./results/%s_%0.0fd_SensitivityNominal_ResultsNuMassScan_MTD-%s_Ba%.0fG_Bs%.2fT_FPDeff%.2f_BKG-%.0fmcps_SystematicsNplus1_%s.mat',...
                SysBudget_local,round(obj.TimeSec/(86400)),TDlabel,obj.MACE_Ba_T*1e4,obj.WGTS_B_T,obj.FPD_MeanEff*obj.FitObj.FPD_Coverage,obj.BKG_RateSec*1e3,obj.ELossFlag);
            if strcmp(obj.RecomputeFlag,'OFF')
                if~exist(save_name,'file') % if file doesnt exist, try without ELossFlag
                    save_name = strrep(save_name,sprintf('_%s',obj.ELossFlag),'');
                    if ~strcmp(obj.ELossFlag,'Aseev')
                        fprintf(2,'Warning: Loading Results with Aseev ELoss \n');
                    end
                end
            end
            
            if exist(save_name,'file') && strcmp(obj.RecomputeFlag,'OFF')
                fprintf('Loading Results from File: %s \n', save_name);
                d = load(save_name);
                obj.MultimNu90Stack.Stat = d.mNu90(1);
                obj.MultimNu90Stack.TC   = d.mNu90(2);
                obj.MultimNu90Stack.FSD  = d.mNu90(3);
                obj.MultimNu90Stack.RF   = d.mNu90(4);
            elseif  strcmp(obj.RecomputeFlag,'OFF')
                fprintf('File: %s \ndoes not exist. Do you want to recompute? \n', save_name);
                return;
            elseif  strcmp(obj.RecomputeFlag,'ON')             
                %Init: gather results
                chi2minScan  = zeros(numel(obj.SysEffectsStack)+1,obj.nFitMax);     % chi2min distribution
                mNu90        = zeros(numel(obj.SysEffectsStack)+1,1);               % sensitivity on mnu^2 (90% C.L.)
                mNumin       = zeros(numel(obj.SysEffectsStack)+1,1);               % mass with minimal chi2
                parScan      = zeros(numel(obj.SysEffectsStack)+1,10,obj.nFitMax);
                errScan      = zeros(numel(obj.SysEffectsStack)+1,10,obj.nFitMax);
                mNu_Fit      = zeros(numel(obj.SysEffectsStack)+1,obj.nFitMax);
                
                chi2prev = obj.chi2;
                SysEffectprev = obj.SysEffect;
                
                % Stat
                obj.chi2      = 'chi2Stat';
                obj.SysEffect = '';
                obj.ComputeCM; %reset CM to 0
                [mNu_Fit(1,:), parScan(1,:,:), errScan(1,:,:), chi2minScan(1,:,:),~,mNu90(1),mNumin(1)] = ...
                    obj.NuMassScan;
                
                % Systematics
                obj.chi2 = 'chi2CM';
                for i=1:numel(obj.SysEffectsStack)
                    obj.SysEffect = obj.SysEffectsStack{i};
                    obj.ComputeCM;
                    [mNu_Fit(i+1,:), parScan(i+1,:,:), errScan(i+1,:,:), chi2minScan(i+1,:,:), dof ,mNu90(i+1),mNumin(i+1)] = ...
                        obj.NuMassScan;
                end
                
                % save
                TD = obj.TD; ScanPrcsn = obj.ScanPrcsn;TimeSec = obj.TimeSec;
                WGTS_B_T = obj.WGTS_B_T; MACE_Ba_T = obj.MACE_Ba_T; BKG_RateSec = obj.BKG_RateSec;
                mySysEffects = obj.SysEffectsAll;
                save(save_name,'parScan','errScan','chi2minScan','mNu_Fit','ScanPrcsn','mNu90', 'mNumin','TD','mySysEffects','dof','WGTS_B_T','MACE_Ba_T','BKG_RateSec','TimeSec');
                obj.chi2      = chi2prev;
                obj.SysEffect = SysEffectprev;
                obj.MultimNu90Stack.Stat = mNu90(1);
                obj.MultimNu90Stack.TC   = mNu90(2);
                obj.MultimNu90Stack.FSD  = mNu90(3);
                obj.MultimNu90Stack.RF   = mNu90(4);
            end
        end
    end
    
    methods % Plot Options
        function PlotScan(obj)
            if isempty(obj.ScanResults)%all(obj.ScanResults.chi2min(2:end))
                fprintf(2,'No Fit Results - Run NuMassScan first! \n');
                return
            end
            % Plot chi2 curve
            mNuPlot = 0:0.01:obj.mNuStop^2;
            chi2func = interp1(obj.ScanResults.mnuSq_i_Fit(~isnan(obj.ScanResults.chi2min)),obj.ScanResults.chi2min(~isnan(obj.ScanResults.chi2min)),mNuPlot,'spline');
            
            % Plot chi2 curve
            f11 = figure('Name','NuMassScanChi2Curve','Renderer','opengl');
            set(f11, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
            x = [0,obj.mNu90];
            y = [min(obj.ScanResults.chi2min),min(obj.ScanResults.chi2min)+2.7];
            pInterp = plot(mNuPlot,chi2func,'-','LineWidth',3,'Color',rgb('CadetBlue'));
            hold on;
            pLine = plot(x,(min(obj.ScanResults.chi2min)+2.7).*ones(numel(x),1),'--','LineWidth',3,'Color',rgb('IndianRed'));
            pScan = plot(obj.ScanResults.mnuSq_i_Fit,obj.ScanResults.chi2min,'o','MarkerSize',10,'MarkerFaceColor',rgb('Navy'),'MarkerEdgeColor',rgb('Navy'));
            plot(obj.mNu90.*ones(numel(y),1),y,'--','LineWidth',3,'Color',rgb('IndianRed'));
            xlabel('m_{\nu}^2 (eV^2)')
            ylabel(['\chi^2 (',num2str(obj.ScanResults.dof(1)),' dof)']);
            scan_leg = ['m_{\nu}^2  = ',sprintf('%.3f \\pm %.3f eV^2 (90%% C.L.)\n',obj.mNumin,obj.mNu90),...
                'm_{\nu}  = ',sprintf('%.3f \\pm %.3f eV   (90%% C.L.)',sqrt(obj.mNumin),sqrt(obj.mNu90))];
            line_leg ='\chi^2_{min}+2.7';
            leg = legend([pInterp, pLine],scan_leg,line_leg,'Location','northwest'); legend boxoff;
            if ~strcmp(obj.chi2,'chi2Stat')
                leg.Title.String = ['stat + ',sprintf(strrep(obj.SysEffect,'_','-'))];
            else
                leg.Title.String = 'stat';
            end
            PrettyFigureFormat;
            title(sprintf('Samak Fit to Asimov \n - KATRIN %.0f years - MTD: %s - BKG %.0f mcps - B_{T},B_{max} = %.0f %%',...
                obj.TimeSec/(60*60*24*365),strrep(obj.TD,'_',' '),obj.BKG_RateSec*1e3, obj.WGTS_B_T/3.6*100));
            set(gca,'FontSize',18);
            xlim([0 obj.mNu90*1.5]);
            ylim([0,interp1(obj.ScanResults.mnuSq_i_Fit(~isnan(obj.ScanResults.chi2min)),obj.ScanResults.chi2min(~isnan(obj.ScanResults.chi2min)),obj.mNu90*1.5,'spline')]);
            save_name = sprintf('Sensitivity_NuMassScan_%s%s_%s_%.2fT-BT_%.2fyears_SysBudget%s',...
                obj.chi2,strrep(obj.SysEffect,'_','-'),obj.TD,obj.WGTS_B_T,obj.TimeSec/(60*60*24*365),obj.SysBudget);
            if ~exist('../sensitivity_nominalKATRIN/plots/png/NuMassScanChi2Curve/','dir')
                mkdir ../sensitivity_nominalKATRIN/plots/png/NuMassScanChi2Curve/;
                mkdir ../sensitivity_nominalKATRIN/plots/pdf/NuMassScanChi2Curve/;
                mkdir ../sensitivity_nominalKATRIN/plots/fig/NuMassScanChi2Curve/;
            end
            export_fig(f11,['../sensitivity_nominalKATRIN/plots/png/NuMassScanChi2Curve/',save_name,'.png']);
            savefig(f11,['../sensitivity_nominalKATRIN/plots/fig/NuMassScanChi2Curve/',save_name,'.fig'],'compact');
            publish_figurePDF(f11,['../sensitivity_nominalKATRIN/plots/pdf/NuMassScanChi2Curve/',save_name,'.pdf']);
        end
        function PlotSysBreakdownBars(obj,varargin)
            % Multibar plot for systematic breakdown
            % mandatory input: Ranges
            p = inputParser;
            p.addParameter('Ranges','', @(x)all(isfloat(x)) && all(x)>0);
            p.addParameter('SingleSys','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('TDRSys','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SysInfoBox','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SavePlot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('DispTitle','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('xLimits','',@(x)all(isfloat(x)) || @(x)isempty(x)); % when empty, y-limit automatically set
            p.addParameter('RFBreakdown','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            
            Ranges     = p.Results.Ranges;
            SingleSys  = p.Results.SingleSys;
            TDRSys     = p.Results.TDRSys;
            SysInfoBox = p.Results.SysInfoBox;
            SavePlot   = p.Results.SavePlot;
            DispTitle  = p.Results.DispTitle;
            xLimits    = p.Results.xLimits;
            RFBreakdown = p.Results.RFBreakdown;
            
            Nranges = numel(Ranges);
            
            % Load Sensitivities (Stacked and Single Contributions)   
            % always use the same systematic budget for all ranges here (usefull in Plot_runtime)
            SysBudgets  = {obj.SysBudget};
            for i=2:Nranges
                SysBudgets = {SysBudgets{:}, obj.SysBudget};
            end
            [PlotmNu90Stack, PlotmNu90, StackBarX, RangeLabel] = obj.LoadSensitivitiesBarPlot(...
                'Ranges',Ranges,'TimeSecs',repmat(obj.TimeSec,1,Nranges),'SysBudgets',SysBudgets,...
                'ImportOpt','Stack+Single');
            
            % Stacked Sensitivity: Prepare Variables for barh plot
            StackBarY = zeros(Nranges,obj.nSys);
            for i=1:Nranges % Neutrino Mass sensitivity (meV)
                StackBarY(i,:) = Convert2Stack(sqrt(PlotmNu90Stack(i,:))*1e3,sqrt(PlotmNu90Stack(i,1))*1e3);
                if StackBarY(i,3)<0 % then TC+FSD is smaller than TC
                    StackBarY(i,3)=0;
                end
            end
            SingleBarStat = sqrt(PlotmNu90Stack(:,1))'*1e3;
            
            % Single Contribution Sensitivity: Prepare Variables for barh plot
            SingleBarY = zeros(Nranges,obj.nSingleSys);
            for i=1:obj.nSingleSys
                SingleBarY(:,i) = (sqrt(PlotmNu90(:,i))'*1e3-SingleBarStat)';
            end
            [SingleBarX, TDRBarX] = obj.GetRangeSingleSys('PlotRanges',StackBarX'); % Single Bar range values
            
            % In case only 1 range: add one invisible row (workaround barh)
            if Nranges==1
                StackBarY = [StackBarY;NaN.*zeros(1,obj.nSys)];
                StackBarX =  [StackBarX,StackBarX+1];               
                SingleBarX = [SingleBarX;SingleBarX+1];
                SingleBarY = [SingleBarY;NaN.*zeros(1,obj.nSingleSys)];
                SingleBarStat = [SingleBarStat,NaN];
                if strcmp(obj.TD,'DR30')
                    RangeLabel = {RangeLabel{:};'TDR'};
                else
                    RangeLabel = {RangeLabel{:};' '};
                end
            end
           
            % MultiBar Plot (Stacked)
            f55 = figure('Name','MultiBarPlot','Renderer','opengl');
            set(f55, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
            b = barh(StackBarX,StackBarY,'stacked');
            b(1).LineStyle = '--';
            b(2).LineStyle = 'none';b(3).LineStyle = 'none'; b(4).LineStyle = 'none';
            b(1).FaceColor = rgb('White');
            b(2).FaceColor = rgb('FireBrick');
            b(3).FaceColor = rgb('GoldenRod');
            b(4).FaceColor = rgb('CadetBlue');
            leg_str = {' Statistical';'+ Theoretical Corrections';'+ Final State Distribution';'+ Response Function'};
         
            % Plot Single Contributions
            switch SingleSys
                case 'ON'
                    if Nranges~=1
                        b(1).BarWidth = b(1).BarWidth/2;
                    end
                    hold on;
                    SingleSysBarWidth = b(1).BarWidth/(2*3);
                    pnone = plot(NaN*[1 1],NaN*[1 1],'Color','w');
                    bTC  = barh(SingleBarX(:,1)',[SingleBarStat;SingleBarY(:,2)']','stacked','BarWidth',SingleSysBarWidth);
                    bFSD = barh(SingleBarX(:,2)', [SingleBarStat;SingleBarY(:,3)']','stacked','BarWidth',SingleSysBarWidth);
                    bTC(2).FaceColor =rgb('FireBrick'); bTC(2).LineStyle = 'none'; bTC(2).FaceAlpha = 0.5;
                    bFSD(2).FaceColor =rgb('GoldenRod'); bFSD(2).LineStyle = 'none'; bFSD(2).FaceAlpha = 0.5;
                    bTC(1).FaceColor = 'w';bFSD(1).FaceColor = 'w'; 
                    bTC(1).LineStyle = 'none'; bFSD(1).LineStyle = 'none';  bTC(1).FaceAlpha = 0;
                    bRF  = barh(SingleBarX(:,3)', [SingleBarStat;SingleBarY(:,4)']','stacked','BarWidth',SingleSysBarWidth);
                    bRF(2).FaceColor =rgb('CadetBlue'); bRF(2).LineStyle = 'none'; bRF(2).FaceAlpha = 0.5;
                    bRF(1).LineStyle = 'none';bRF(1).FaceColor = 'w';
                     leg_str = {leg_str{:},'Single Contributions','Stat + Theoretical Corrections', 'Stat + Final State Distribution','Stat + Response Function'};
                    if strcmp(RFBreakdown,'OFF')
                    leg = legend([b,pnone,bTC(2),bFSD(2),bRF(2)],leg_str{:});
                    leg.NumColumns = 2;
                    elseif strcmp(RFBreakdown,'ON')
                        SingleBarX2 = [11.43 10.1 8.77; 31.43 30.1 28.77];
                        bRX  = barh(SingleBarX2(:,2)', [SingleBarStat;SingleBarY(:,8)']','stacked','BarWidth',SingleSysBarWidth);
                        bBF  = barh(SingleBarX2(:,1)', [SingleBarStat;SingleBarY(:,7)']','stacked','BarWidth',SingleSysBarWidth);
                        bRX(1).FaceColor = 'w'; bBF(1).FaceColor = 'w'; 
                        bRX(1).LineStyle = 'none'; bBF(1).LineStyle = 'none';  bRX(1).FaceAlpha = 0;
                        bRX(2).FaceColor =rgb('CornFlowerBlue'); bRX(2).LineStyle = 'none'; bRX(2).FaceAlpha = 1;
                        bBF(2).FaceColor =rgb('Navy'); bBF(2).LineStyle = 'none'; bBF(2).FaceAlpha = 1;
                        
                        if ~strcmp(obj.SysBudget,'09')
                            bEL  = barh(SingleBarX2(:,3)', [SingleBarStat;SingleBarY(:,6)']','stacked','BarWidth',SingleSysBarWidth);
                            bEL(2).FaceColor =rgb('PowderBlue'); bEL(2).LineStyle = 'none'; bEL(2).FaceAlpha = 1;
                            bEL(1).LineStyle = 'none';bEL(1).FaceColor = 'w';
                            leg_str = {leg_str{:},'Stat + Magnetic Fields',...
                                'Stat + Column Density + Inel. Cross Section','Stat + Energy Loss Function'};
                             leg = legend([b,pnone,bTC(2),bFSD(2),bRF(2),bBF(2),bRX(2),bEL(2)],leg_str{:});
                        else
                            leg_str = {leg_str{:},'Stat + Magnetic Fields',...
                                'Stat + Column Density + Inel. Cross Section'};
                             leg = legend([b,pnone,bTC(2),bFSD(2),bRF(2),bBF(2),bRX(2)],leg_str{:});
                        end
                        leg.NumColumns = 3;
                       % ylim([2 60])
                    end
                    
                    if strcmp(TDRSys,'ON') && ~strcmp(obj.TD,'DR30')
                        hold on;
                        %TDR Systemaitcs for higher range: TDR + Marco Kleesiek's value due to rhod sigma and Ba
                        if Nranges==3 % 30, 45 ,60
                            TDRSys90 = sqrt(sqrt(PlotmNu90Stack(:,1).^2+(1.64.*[0.0213;0.0220;0.0226]).^2))*1e3-SingleBarStat';
                        elseif Nranges==2 % 30 and 60
                            TDRSys90 = sqrt(sqrt(PlotmNu90Stack(:,1).^2+(1.64.*[0.0213;0.0226]).^2))*1e3-SingleBarStat';
                        end
                        bTDR =   barh(TDRBarX,[SingleBarStat',TDRSys90, sqrt(PlotmNu90Stack(:,4))*1e3-TDRSys90],'stacked','FaceColor',rgb('SlateGray'),'FaceAlpha',1);
                        bTDR(1).BarWidth = bTC.BarWidth; bTDR(2).BarWidth = b(2).BarWidth/2;
                        bTDR(1).FaceColor = 'w'; bTDR(1).LineStyle='none';bTDR(2).LineStyle='none';bTDR(1).FaceAlpha = 0;
                        bTDR(3).FaceColor = 'w'; bTDR(3).LineStyle = 'none';
                        leg_str = {leg_str{1:4},'TDR + Neutrino16',leg_str{5:end}};
                        leg = legend([b,bTDR(2),pnone,bTC(2),bFSD(2),bRF(2)],leg_str{:});
                        leg.NumColumns = 2;
                    end
                    
                case 'OFF'
                    leg = legend(leg_str{:});
            end
            
            % special treatment for TDR
            if strcmp(obj.TD,'DR30')
                hold on;
                mNu90DRStat = sqrt(0.018*1.64)*1e3;
                mNu90DRSys  = sqrt(sqrt(0.018^2+0.017^2)*1.64)*1e3-mNu90DRStat;
                bTDR = barh([21;NaN],[mNu90DRStat,mNu90DRSys;NaN,NaN],'Stacked');
                bTDR(1).FaceColor = 'w'; bTDR(1).LineStyle = '--';
                bTDR(2).FaceColor = rgb('SlateGray'); bTDR(2).LineStyle = 'none';
                leg_str = {leg_str{:},'TDR04 Statistical    (\sigma = 0.018eV^2)','TDR04 Systematics (\sigma = 0.017eV^2)'}; %_{m²_\nu}
                if strcmp(SingleSys,'ON')
                    leg = legend([b,pnone,bTC(2),bFSD(2),bRF(2),bTDR(1),bTDR(2)],leg_str{:});
                else
                    leg = legend([b,bTDR(1),bTDR(2)],leg_str{:});
                end
                ylim([min(StackBarX)-1.5; max(StackBarX)+1.3]); 
                xlim([160 300]);
                xticks((160:5:300));
                leg.NumColumns = 3;
                yticklabels({'30 eV','30 eV TDR'});
            elseif Nranges==1
                ylim([min(StackBarX)-1.5; max(StackBarX)]);
            end
            
            % legend and labels, display options 
            leg.FontSize = 12;
            legend boxoff
            yticklabels(RangeLabel);
            ytickangle(90)
            xlabel('neutrino mass sensitivity 90% C.L. (meV)');
            ylabel('MTD range');
            
            if ~strcmp(obj.TD,'DR30')
                if ~isempty(xLimits)
                    xmin = xLimits(1);
                    xmax = xLimits(2);
                else
                    xmin =min(floor(min(SingleBarStat)/10)*10-10);
                    xmax = sqrt(max(PlotmNu90Stack(:,end)))*1e3+10;
                end
                xlim([200 xmax*2]);
                xticks((200:20:xmax*2));
                if contains(obj.TD,'Optimize')
                    xlim([xmin, xmax]);
                    xticks((xmin:10:xmax));
                elseif contains(obj.TD,'MTDcreator')
                    xlim([xmin, xmax]);
                    if xmax-xmin>200
                        xticks((xmin:40:xmax));
                    elseif xmax-xmin>100
                        xticks((xmin:20:xmax));
                    else
                         xticks((xmin:10:xmax)); 
                    end
                end
            end
            

            grid on;
            PrettyFigureFormat;
            set(gca,'FontSize',18);
            
            if contains(obj.TD,'MTDcreator') || contains(obj.TD,'Optimize')
                labelTime = obj.TimeSec/(24*60*60*124/148); % Correct number of pixels
            else
                labelTime = obj.TimeSec/(24*60*60);
            end
            %TDlabel1 = strrep(obj.TD,sprintf('_Ba%.0fG',obj.MACE_Ba_T*1e4),''); %
            %TDlabel2 = strrep(strrep(TDlabel1,sprintf('_%.0feV',Ranges(end)),'Opt'),'_',' ');
            TDlabel = strrep(obj.TD,sprintf('_E018575.0_%.0feV_B35_Ba7.0_RedF0.7_NuMF0.60_BkgF0.15_B364',Ranges(1)),'');
            TDlabel = strrep(TDlabel,sprintf('_E018575.0_%.0feV_B35_Ba7.0_RedF0.7_NuMF0.60_BkgF0.15_B364',Ranges(end)),'');
            if strcmp(DispTitle,'ON')
%                 title(sprintf('Systematic Uncertainty Breakdown \n Nominal KATRIN %.0f days, B_a = %.0fG, B_{T/max} = %.0f %%, Background = %.0f mcps',...
%                     labelTime, obj.MACE_Ba_T*1e4,obj.WGTS_B_T/3.6*100,obj.BKG_RateSec*1e3));
                title(sprintf('KATRIN - Systematic Uncertainty Breakdown \n Column density = %.2g mol/cm^2, %.0f days, B_a = %.0fG, B_{T/max} = %.0f %%, Background = %.0f mcps',...
                    obj.SimObj.WGTS_CD_MolPerCm2, labelTime, obj.MACE_Ba_T*1e4,obj.WGTS_B_T/3.6*100,obj.BKG_RateSec*1e3));
            end
            % Annotation Box with Systematic Info
            if strcmp(SysInfoBox,'ON')
                SysUncertainties = obj.GetSysInfo('SysFlag','SysAll');
                a=annotation('textbox', [0.14 0.1 1 0.12], ...
                    'String', SysUncertainties, ...
                    'EdgeColor', 'none', ...
                    'HorizontalAlignment', 'left');
                a.FontSize=11;a.FontWeight='bold';
            end
            
            % Save Plot
            if strcmp(SavePlot,'ON')
                if ~exist('../sensitivity_nominalKATRIN/plots/pdf/SysBreakdown/','dir')
                    mkdir ../sensitivity_nominalKATRIN/plots/pdf/SysBreakdown/
                    mkdir ../sensitivity_nominalKATRIN/plots/png/SysBreakdown/
                    mkdir ../sensitivity_nominalKATRIN/plots/fig/SysBreakdown/
                end
                if strcmp(obj.TD,'DR30')
                    save_name = sprintf('%s_%0.0fd_SysBreakdownBar_%s_Ba%.0fG_B%.0fSingleSys%s_%s_%.2fFPDeff',...
                        obj.SysBudget,round(obj.TimeSec/(86400)),strrep(TDlabel,sprintf('_%.0feV',Ranges(end)),'Opt'),...
                        obj.MACE_Ba_T*1e4,obj.WGTS_B_T/3.6*100,SingleSys,obj.ELossFlag,obj.FPD_MeanEff*obj.FitObj.FPD_Coverage);
                else
                    save_name = sprintf('%s_%0.0fd_SysBreakdownBar_%s_Ba%.0fG_B%.0fSingleSys%s_TDRSys%s_%s_%.2fFPDeff',...
                        obj.SysBudget,round(obj.TimeSec/(86400)),strrep(TDlabel,sprintf('_%.0feV',Ranges(end)),'Opt'),...
                        obj.MACE_Ba_T*1e4,obj.WGTS_B_T/3.6*100,SingleSys,TDRSys,obj.ELossFlag,obj.FPD_MeanEff*obj.FitObj.FPD_Coverage);
                end
                if strcmp(RFBreakdown,'ON')
                    save_name = [save_name,'_RF'];
                end
                print(f55,['../sensitivity_nominalKATRIN/plots/png/SysBreakdown/',save_name,'.png'],'-dpng');
                savefig(f55,['../sensitivity_nominalKATRIN/plots/fig/SysBreakdown/',save_name,'.fig'],'compact');
                publish_figurePDF(f55,['../sensitivity_nominalKATRIN/plots/pdf/SysBreakdown/',save_name,'.pdf']);
                
            end
        end
        function PlotRFBreakdownBars(obj,varargin)
            % Multibar plot for systematic breakdown
            % mandatory input: Ranges
            p = inputParser;
            p.addParameter('Ranges','', @(x)all(isfloat(x)) && all(x)>0);
            p.addParameter('SingleSys','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SysInfoBox','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SavePlot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('DispTitle','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('xLimits','',@(x)all(isfloat(x)) || @(x)isempty(x)); % when empty, y-limit automatically set 
          
            p.parse(varargin{:});
            
            Ranges     = p.Results.Ranges;
            SingleSys  = p.Results.SingleSys;
            SysInfoBox = p.Results.SysInfoBox;
            SavePlot   = p.Results.SavePlot;
            DispTitle  = p.Results.DispTitle;
            xLimits    = p.Results.xLimits;
            
            Nranges = numel(Ranges);
            
            % always use the same systematic budget for all ranges here
            % (usefull in Plot_runtime)
            SysBudgets  = {obj.SysBudget};
            for i=2:Nranges
                SysBudgets = {SysBudgets{:}, obj.SysBudget};
            end
            
            % Load Sensitivities (Stacked and Single Contributions)
            [~, PlotmNu90, StackBarX, RangeLabel] = obj.LoadSensitivitiesBarPlot(...
                'Ranges',Ranges,'TimeSecs',repmat(obj.TimeSec,1,Nranges),'SysBudgets',SysBudgets,...
                'ImportOpt','Single'); 
            
            % Prepare Stacked Parameter
            StackBarY = zeros(Nranges,obj.nSys);
            RFStack_tmp = zeros(Nranges,obj.nSys);
            RFStack_tmp(:,1) = 1e3*sqrt(PlotmNu90(:,1)); % stat
            RFStack_tmp(:,2) = 1e3*sqrt(PlotmNu90(:,8)); % RF - rho d
            RFStack_tmp(:,3) = 1e3*sqrt(PlotmNu90(:,9)); % RF - rho d and b-fields
            RFStack_tmp(:,4) = 1e3*sqrt(PlotmNu90(:,4)); % RF - all
            for i=1:Nranges % Neutrino Mass sensitivity (meV)
                StackBarY(i,:) = Convert2Stack(RFStack_tmp(i,:),RFStack_tmp(i,1));
            end
            
            %Prepare Single Contribution Parameters
            SingleBarY(:,1) = 1e3*sqrt(PlotmNu90(:,8))-1e3*sqrt(PlotmNu90(:,1)); % rhod only
            SingleBarY(:,2) = 1e3*sqrt(PlotmNu90(:,7))-1e3*sqrt(PlotmNu90(:,1)); % bfields only
            SingleBarY(:,3) = 1e3*sqrt(PlotmNu90(:,6))-1e3*sqrt(PlotmNu90(:,1)); % eloss only
            
            % Get x values for single contributions
            [SingleBarX, ~] = obj.GetRangeSingleSys('PlotRanges',StackBarX'); % Single Bar range values
            
            % In case only 1 range: add one invisible row (workaround barh)
            if Nranges==1
                StackBarY = [StackBarY;NaN.*zeros(1,obj.nSys)];
                StackBarX =  [StackBarX,StackBarX+1];
                RangeLabel = {RangeLabel{:};' '};
                SingleBarX = [SingleBarX;SingleBarX+1];
                SingleBarY = [SingleBarY;NaN.*zeros(1,3)];   
            end
            
            % MultiBar plot
            f55 = figure('Name','MultiBar_RFBreakdown','Renderer','opengl');
            set(f55, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
            if strcmp(obj.SysBudget,'08') || strcmp(obj.SysBudget,'09')
                b = barh(StackBarX,StackBarY(:,1:3),'stacked','BaseValue',0);
            else
                b = barh(StackBarX,StackBarY,'stacked','BaseValue',0);
            end
            b(1).FaceColor = rgb('White');
            b(3).FaceColor = rgb('DarkCyan');
            b(2).FaceColor = rgb('PowderBlue');
            b(1).LineStyle = '--'; b(2).LineStyle = 'none';b(3).LineStyle = 'none';
            % legends
            if strcmp(obj.SysBudget,'08') || strcmp(obj.SysBudget,'09')
                legStack = {'Statistical','+ Column Density and Inel. Cross Section','+ Magnetic Fields'};
            else
                b(4).FaceColor = rgb('Orange');b(4).LineStyle = 'none';
                legStack = {'Statistical','+ Column Density and Inel. Cross Section','+ Magnetic Fields','+ Energy Loss Functions'};
            end
            leg = legend(legStack{:},'Location','northwest');
            if strcmp(SingleSys,'ON')
                if Nranges~=1
                    b(1).BarWidth = b(1).BarWidth/2;
                end
                hold on;
                pnone = plot(NaN*[1 1],NaN*[1 1],'Color','w');
                bRX = barh(SingleBarX(:,1),[StackBarY(:,1),SingleBarY(:,1)],'stacked','BarWidth',b(1).BarWidth/6);
                bBF = barh(SingleBarX(:,2),[StackBarY(:,1),SingleBarY(:,2)],'stacked','BarWidth',b(1).BarWidth/6);
                bRX(1).FaceAlpha = 0; bRX(1).LineStyle = 'none';
                bRX(2).LineStyle = 'none'; bRX(2).FaceAlpha =0.5; bRX(2).FaceColor = rgb('PowderBlue');
                bBF(1).FaceAlpha = 0;  bBF(1).LineStyle = 'none';
                bBF(2).LineStyle = 'none';  bBF(2).FaceAlpha =0.5; bBF(2).FaceColor = rgb('DarkCyan');
                if ~strcmp(obj.SysBudget,'08') && ~strcmp(obj.SysBudget,'09')
                    bEL = barh(SingleBarX(:,3),[StackBarY(:,1),SingleBarY(:,3)],'stacked','BarWidth',b(1).BarWidth/6);
                    bEL(1).FaceAlpha = 0; bEL(1).LineStyle = 'none';
                    bEL(2).LineStyle = 'none'; bEL(2).FaceAlpha =0.5; bEL(2).FaceColor = rgb('Orange');
                    legSingle = {'Single Contributions','Stat + Column Density and Inel. Cross Section',...
                        'Stat + Magnetic Fields','Stat + Energy Loss Functions'};
                    leg = legend([b,pnone,bRX(2),bBF(2),bEL(2)], legStack{:},legSingle{:});
                else
                    legSingle = {'Single Contributions','Stat + Column Density and Inel. Cross Section',...
                        'Stat + Magnetic Fields'};
                    leg = legend([b,pnone,bRX(2),bBF(2)], legStack{:},legSingle{:});
                end
                leg.NumColumns = 2;
            end
            leg.FontSize = 12;
            legend boxoff
            yticks(StackBarX);
            yticklabels(RangeLabel);
            ytickangle(90)
            xlabel('neutrino mass sensitivity 90% C.L. (meV)');
            ylabel('MTD range');
            xmax = ceil(max(max(sqrt(PlotmNu90(:,4))*1e3))*1.01);
            if Nranges==3
                ylim([min(min((SingleBarX)))-20 max(StackBarX)+30]);
                xmin = 200;%ceil(min(min([mNu90Stack60 mNu90Stack30]))-10);
            elseif Nranges==2
                ylim([min(min((SingleBarX)))-8 max(StackBarX)+15]);
                xmin = floor(min(min(sqrt(PlotmNu90)*1e3))/10)*10-10;
            elseif Nranges==1
                ylim([min(min((SingleBarX)))-0.5 max(StackBarX)]);
                xmin = floor(min(min(sqrt(PlotmNu90)*1e3))/10)*10;
            end
            %xmax = 900;xmin = 540;
            if ~isempty(xLimits)
                xmin = xLimits(1);
                xmax = xLimits(2);
            end
            xlim([xmin xmax]);
            if xmax-xmin>200
                xticks((xmin:40:xmax));
            elseif xmax-xmin>100
                xticks((xmin:20:xmax));
            else
                xticks((xmin:10:xmax));
            end
            if contains(obj.TD,'MTDcreator') || contains(obj.TD,'Optimize')
                labelTime = obj.TimeSec/(24*60*60*124/148);
            elseif strcmp(obj.TD,'DR30')
                labelTime = obj.TimeSec/(24*60*60);
                xticks(xmin:5:xmax);
            else
                labelTime = obj.TimeSec/(24*60*60);
            end
            if strcmp(DispTitle,'ON')
                title(sprintf('Response Function Uncertainty Breakdown \nNominal KATRIN %.0f days, B_a = %.0fG, B_{T/max} = %.0f %%, Background = %.0f mcps',...
                    labelTime, obj.MACE_Ba_T*1e4,obj.WGTS_B_T/3.6*100,1e3*obj.BKG_RateSec));
            end
            PrettyFigureFormat;
            PlotFontSize = 18;
            set(gca,'FontSize',PlotFontSize);
            grid on;
            
            % Annotation Box with Systematic Info
            if strcmp(SysInfoBox,'ON')
                SysUncertainties = obj.GetSysInfo('SysFlag','SysRF');
                a=annotation('textbox', [0.14 0.1 1 0.12], ...
                    'String', SysUncertainties, ...
                    'EdgeColor', 'none', ...
                    'HorizontalAlignment', 'left');
                a.FontSize=11;a.FontWeight='bold';
            end
           
            if strcmp(SavePlot,'ON')
               % TDlabel1 = strrep(obj.TD,sprintf('_Ba%.0fG',obj.MACE_Ba_T*1e4),''); %
                %TDlabel_tmp = strrep(TDlabel1,['_',num2str(Ranges(end)),'eV'],'');
                TDlabel = strrep(obj.TD,sprintf('_E018575.0_%.0feV_B35_Ba7.0_RedF0.7_NuMF0.60_BkgF0.15_B364',Ranges(1)),'');
                TDlabel = strrep(TDlabel,sprintf('_E018575.0_%.0feV_B35_Ba7.0_RedF0.7_NuMF0.60_BkgF0.15_B364',Ranges(end)),'');
                 save_name = sprintf('%s_%0.0fd_%s_RFBreakdownBar_Bfields%.0fpercent_Ba%.0fG_SingleSys%s_%s_%.2fFPDeff_BKG-%.0fmcps',...
                    obj.SysBudget,obj.TimeSec/86400,TDlabel,obj.WGTS_B_T/3.6*100,obj.MACE_Ba_T*1e4,SingleSys,obj.ELossFlag,obj.FPD_MeanEff*obj.FitObj.FPD_Coverage,obj.BKG_RateSec*1e3);
                if ~exist('../sensitivity_nominalKATRIN/plots/png/RFBreakdown/','dir')
                    mkdir ../sensitivity_nominalKATRIN/plots/png/RFBreakdown/
                    mkdir ../sensitivity_nominalKATRIN/plots/pdf/RFBreakdown/
                    mkdir ../sensitivity_nominalKATRIN/plots/fig/RFBreakdown/
                end
                print(f55,['../sensitivity_nominalKATRIN/plots/png/RFBreakdown/',save_name,'.png'],'-dpng');
                savefig(f55,['../sensitivity_nominalKATRIN/plots/fig/RFBreakdown/',save_name,'.fig'],'compact');
                publish_figurePDF(f55,['../sensitivity_nominalKATRIN/plots/pdf/RFBreakdown/',save_name,'.pdf']);
            end
        end
        function Plot_RunTimeEvolution(obj,varargin)
            % Multibar plot for systematic breakdown
            % mandatory input: SysBudgets, TimeSecs
            p = inputParser;
            p.addParameter('SysBudgets',{'09_NoELoss','07','06'},@(x)iscell(x));
            p.addParameter('TimeSecs',(124/148*24*60*60).*[900, 300, 42],@(x)all(isfloat(x)));
            p.addParameter('SingleSys','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SavePlot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SysInfoBox','OFF',@(x)ismember(x,{'ON','OFF'}));
            
            p.parse(varargin{:});
            
            SysBudgets = p.Results.SysBudgets;
            TimeSecs   = p.Results.TimeSecs;
            SingleSys  = p.Results.SingleSys;
            SavePlot   = p.Results.SavePlot;    
            SysInfoBox = p.Results.SysInfoBox;
            Ntimes = numel(TimeSecs);
           
            % Load Sensitivities (Stacked and Single Contributions)
            [PlotmNu90Stack, PlotmNu90, StackBarX, TimeLabel] = obj.LoadSensitivitiesBarPlot(...
                'TimeSecs',TimeSecs,'SysBudgets',SysBudgets,'Ranges',repmat(obj.range,Ntimes,1),...
                'ImportOpt','Stack+Single');
            
            % Stacked Sensitivity: Prepare Variables for barh plot
            StackBarY = zeros(Ntimes,obj.nSys);
            for i=1:Ntimes % Neutrino Mass sensitivity (meV)
                StackBarY(i,:) = Convert2Stack(sqrt(PlotmNu90Stack(i,:))*1e3,sqrt(PlotmNu90Stack(i,1))*1e3);
                if StackBarY(i,3)<0 % then TC+FSD is smaller than TC
                    StackBarY(i,3)=0;
                end
            end
            SingleBarStat = sqrt(PlotmNu90Stack(:,1))'*1e3;
            
            % Single Contribution Sensitivity: Prepare Variables for barh plot
            SingleBarY = zeros(Ntimes,obj.nSingleSys);
            for i=1:obj.nSingleSys
                SingleBarY(:,i) = (sqrt(PlotmNu90(:,i))'*1e3-SingleBarStat)';
            end
            [SingleBarX, ~] = obj.GetRangeSingleSys('PlotRanges',StackBarX'); % Single Bar range values
            
            % In case only 1 range: add one invisible row (workaround barh)
            if Ntimes==1
                StackBarY = [StackBarY;NaN.*zeros(1,obj.nSys)];
                StackBarX =  [StackBarX,StackBarX+1];
                RangeLabel = {RangeLabel{:};' '};
                SingleBarX = [SingleBarX;SingleBarX+1];
                SingleBarY = [SingleBarY;NaN.*zeros(1,obj.nSingleSys)];
                SingleBarStat = [SingleBarStat,NaN];
            end
            
            % Multibar Plot (Stacked)
            f51 = figure('Name','MultiBar_RunTime','Renderer','opengl');
            set(f51, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
            b = barh(StackBarX,StackBarY,'stacked','BaseValue',0);
            b(1).FaceColor = rgb('White');
            b(2).FaceColor = rgb('FireBrick');
            b(3).FaceColor = rgb('GoldenRod');
            b(4).FaceColor = rgb('CadetBlue');
            b(1).LineStyle = '--'; b(2).LineStyle = 'none';b(3).LineStyle = 'none';b(4).LineStyle = 'none';
            leg_str = {' Statistical';'+ Theoretical Corrections';'+ Final State Distribution';'+ Response Function'};
     
            % Plot Single Contributions
            switch SingleSys
                case 'ON'
                    if Ntimes~=1
                        b(1).BarWidth = b(1).BarWidth/2;
                    end
                    hold on;
                    SingleSysBarWidth = b(1).BarWidth/(2*3);
                    pnone = plot(NaN*[1 1],NaN*[1 1],'Color','w');
                    bTC  = barh(SingleBarX(:,1)',[SingleBarStat;SingleBarY(:,2)']','stacked','BarWidth',SingleSysBarWidth);
                    bFSD = barh(SingleBarX(:,2)', [SingleBarStat;SingleBarY(:,3)']','stacked','BarWidth',SingleSysBarWidth);
                    bRF  = barh(SingleBarX(:,3)', [SingleBarStat;SingleBarY(:,4)']','stacked','BarWidth',SingleSysBarWidth);
                    bTC(2).FaceColor =rgb('FireBrick'); bTC(2).LineStyle = 'none'; bTC(2).FaceAlpha = 0.5;
                    bFSD(2).FaceColor =rgb('GoldenRod'); bFSD(2).LineStyle = 'none'; bFSD(2).FaceAlpha = 0.5;
                    bRF(2).FaceColor =rgb('CadetBlue'); bRF(2).LineStyle = 'none'; bRF(2).FaceAlpha = 0.5;
                    bTC(1).FaceColor = 'w';bFSD(1).FaceColor = 'w'; bRF(1).FaceColor = 'w';
                    bTC(1).LineStyle = 'none'; bFSD(1).LineStyle = 'none'; bRF(1).LineStyle = 'none'; bTC(1).FaceAlpha = 0;
                    leg_str = {leg_str{:},'Single Contributions','Stat + Theoretical Corrections', 'Stat + Final State Distribution','Stat + Response Function'};
                    leg = legend([b,pnone,bTC(2),bFSD(2),bRF(2)],leg_str{:});                      
                case 'OFF'
                    leg = legend(leg_str{:});
            end
            
            % legend and labels, display options
            leg.NumColumns = 2;
            leg.FontSize = 12;
            legend boxoff
            yticklabels(TimeLabel);
            ytickangle(90)
            xlabel('neutrino mass sensitivity 90% C.L. (meV)');
            ylabel('measurement time (days)');
            
            if ~strcmp(obj.TD,'DR30')
                xmax = 10*ceil(0.1*(sqrt(max(PlotmNu90Stack(:,end)))*1e3+10));
                xmin = 10*floor(0.1*(sqrt(min(PlotmNu90Stack(:,1)))*1e3-10));
                %xmin = 540; xmax = 900;
                xlim([xmin xmax]);
                xticks((xmin:50:xmax));
                if Ntimes==1
                    ylim([min(StackBarX)-1.5; max(StackBarX)]);
                elseif Ntimes==2
                    ylim([min(min((SingleBarX)))-8 max(StackBarX)+15])
                elseif Ntimes==3
                    ylim([min(min((SingleBarX)))-15 max(StackBarX)+25]);
                elseif Ntimes==5
                    ylim([min(min((SingleBarX)))-15 max(StackBarX)+35]); 
                end
            end
            grid on;
            PrettyFigureFormat;
            set(gca,'FontSize',14);
            if contains(obj.TD,'Optimize') || contains(obj.TD,'MTDcreator')
                titlestr = sprintf('Samak Optimized MTD - %.0feV scan range -  Background %.0f mcps',obj.range,obj.BKG_RateSec*1e3);
                title(titlestr);
            else
                title(sprintf('MTD %s - Background %.0f mcps ',strrep(obj.TD,'_',' '),obj.BKG_RateSec*1e3));
            end
            % Annotation Box with Systematic Info
            if strcmp(SysInfoBox,'ON')
                SysUncertainties = obj.GetSysInfo('SysFlag','SysAll');
                a=annotation('textbox', [0.14 0.1 1 0.12], ...
                    'String', SysUncertainties, ...
                    'EdgeColor', 'none', ...
                    'HorizontalAlignment', 'left');
                a.FontSize=11;a.FontWeight='bold';
            end
            
            if strcmp(SavePlot,'ON')
                if ~exist('../sensitivity_nominalKATRIN/plots/png/RunTime/','dir')
                    mkdir ../sensitivity_nominalKATRIN/plots/png/RunTime/
                    mkdir ../sensitivity_nominalKATRIN/plots/pdf/RunTime/
                    mkdir ../sensitivity_nominalKATRIN/plots/fig/RunTime/
                end
                TDlabel = strrep(obj.TD,sprintf('_E018575.0_%.0feV_B35_Ba7.0_RedF0.7_NuMF0.60_BkgF0.15_B364',obj.range),'');
                save_name = sprintf('RunTime_%.0feV_%s_BKG-%.0f_FPDeff%.2f_%.0fBars_%s',...
                    obj.range,TDlabel,obj.BKG_RateSec*1e3,obj.FPD_MeanEff*obj.FitObj.FPD_Coverage,numel(TimeSecs),obj.ELossFlag);
                print(f51,['../sensitivity_nominalKATRIN/plots/png/RunTime/',save_name,'.png'],'-dpng');
                savefig(f51,['../sensitivity_nominalKATRIN/plots/fig/RunTime/',save_name,'.fig'],'compact');
                publish_figurePDF(f51,['../sensitivity_nominalKATRIN/plots/pdf/RunTime/',save_name,'.pdf']);
                
            end
        end
    end
    methods % Auxillary methods
        function [PlotmNu90Stack, PlotmNu90, PlotRange, XLabel] = LoadSensitivitiesBarPlot(obj,varargin)
            % Load sensitivites for MultiBarPlot for all 'Ranges' if they exist
            p = inputParser;
            p.addParameter('Ranges',[30 45 60], @(x)all(isfloat(x)) && all(x)>0);
            p.addParameter('SysBudgets',{'09_NoELoss','07','06'}',@(x)iscell(x)); %
            p.addParameter('TimeSecs',(124/148*24*60*60).*[900, 300, 42],@(x)all(isfloat(x)));
            p.addParameter('ImportOpt','Stack+Single',@(x)ismember(x,{'Stack','Single','Stack+Single'}));
            
            p.parse(varargin{:});
            
            Ranges     = p.Results.Ranges;
            ImportOpt  = p.Results.ImportOpt;
            SysBudgets = p.Results.SysBudgets;
            TimeSecs   = p.Results.TimeSecs;
            
            if ~isequal(numel(Ranges),numel(TimeSecs),numel(SysBudgets))
                fprintf(2,'Ranges, SysBudgets and TimeSecs need to have the same size! \n');
                return
            end
            
            % save previous setting
            range_prev         = obj.range;
            RecomputeFlag_prev = obj.RecomputeFlag;
            SysBudget_prev     = obj.SysBudget;
            TimeSec_prev       = obj.TimeSec;
            
            % Init
            Nranges = numel(Ranges);
            PlotmNu90Stack  = zeros(Nranges,obj.nSys);       % stat, tc, tc+fsd, tc+fsd+rf
            PlotmNu90       = zeros(Nranges,obj.nSingleSys); % stat, tc, fsd, rf, all, rf_el, rf_bf, rf_rx, rf_bfrx
            PlotRange = zeros(Nranges,1);
            RangeLabel = cell(Nranges,1);
            TimeLabel = cell(Nranges,1);
            % Import
            for i=1:Nranges
                obj.range = Ranges(i);
                obj.SysBudget = SysBudgets{i};
                obj.TimeSec = TimeSecs(i);
                obj.SetTD;
                obj.RecomputeFlag = 'OFF';
                if strcmp(ImportOpt,'Stack') || strcmp(ImportOpt,'Stack+Single')
                    obj.NuMassScan_SensitivityNominal_Nplus1Systematics_Loop;
                    PlotmNu90Stack(i,:) = structfun(@(x) x,obj.MultimNu90Stack); % convert to sensitivity on neutrino mass (meV)
                end
                if strcmp(ImportOpt,'Single') || strcmp(ImportOpt,'Stack+Single')
                    obj.NuMassScan_SystematicsLoop;
                    PlotmNu90(i,:) =  structfun(@(x) x,obj.MultimNu90);
                end
                RangeLabel{i} = sprintf('%.0f eV',Ranges(i));
                TimeLabel{i} = sprintf('%.0f',TimeSecs(i)/(60*60*24*124/148));
                PlotRange(i) = i*20;
            end
            
            if ~strcmp(obj.TD,'DR30')
                if isequal(Ranges(1),Ranges(2))
                    XLabel = TimeLabel;
                else
                    XLabel = RangeLabel;
                end
            elseif strcmp(obj.TD,'DR30')
                XLabel = RangeLabel;
            end
            % reset to previous setting
            obj.range         = range_prev;
            obj.RecomputeFlag = RecomputeFlag_prev;
            obj.TimeSec       = TimeSec_prev;
            obj.SysBudget     = SysBudget_prev;
        end
        function y = Convert2Stack(x,startValue)
            %x is array to be converted
            y = zeros(1,numel(x));
            y(1)= startValue;
            for i=2:numel(y)
                y(i)= x(i)-sum(y(1:i-1));
            end
        end
        function SetTD(obj,varargin)
            p=inputParser;
            p.addParameter('MTDcreatorFlag','ON',@(x)ismember(x,{'ON','OFF'}));  
            p.parse(varargin{:});
            MTDcreatorFlag = p.Results.MTDcreatorFlag;
            
            % Make sure the correct TD is used
            % e.g. when range is changed --> TD has to be adapted
            if strcmp(obj.TD,'DR30') && obj.range~=30
                fprintf(2,'WARNING: TD = Design Report 30! Change range to 30!\n')
                return
            elseif contains(obj.TD,'Sensitivity')
                obj.TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',obj.range,obj.MACE_Ba_T*1e4);
            elseif contains(obj.TD,'Optimize')
                obj.TD  = sprintf('OptimizeMTD%.0f_NuMassFactor_BkgFraction_all_03',obj.range);
            elseif contains(obj.TD,'MTDcreator')
                
               % NuMassFactor = 0.6;
               % BkgFraction  = 0.15;
                % obj.TD  = sprintf('MTDcreator_E018575.0_%.0feV_B35_Ba7.0_RedF0.7_NuMF0.59_BkgF0.20',obj.range);
               % obj.TD  = sprintf('MTDcreator_E0%.1f_%.0feV_B35_Ba%.1f_RedF%.1f_NuMF%.2f_BkgF%.2f_B%0.f',...
                %    obj.Q_i,obj.range,obj.MACE_Ba_T*1e4,obj.WGTS_B_T/3.6,NuMassFactor,BkgFraction,obj.BKG_RateSec*1e3);
                obj.TD = sprintf('MTDcreator_E0%.1f_%.0feV_B%.0f_Ba%.1f_RedF%0.1f_NuMF%0.2f_BkgF%0.2f_B%.0f',...
                    obj.Q_i,obj.range,sum(obj.TD_BqU),obj.MACE_Ba_T/1e-4,obj.WGTS_B_T/3.6,...
                    obj.TD_NuMassFactor,obj.TD_BkgFraction,round(obj.BKG_RateSec*1e3));
                % Create MTD if it does not exist yet
                save_name = ['../../simulation/katrinsetup/TD_DataBank/' obj.TD '.mat'];
                if ~exist(save_name,'file') && strcmp(MTDcreatorFlag,'ON')
                    [qU, qUfrac, TD] = MTDcreator(...
                        'E0',obj.Q_i,...
                        'Range',obj.range,...
                        'BqU',obj.TD_BqU,...
                        'MACE_Ba_T',obj.MACE_Ba_T,...
                        'BsBmRed',obj.WGTS_B_T/3.6,...
                        'RunTime',obj.TimeSec,...
                        'NuMassFactor',obj.TD_NuMassFactor,...
                        'BkgFraction',obj.TD_BkgFraction,...
                        'Save','ON',...
                        'MTD_Plot','OFF',...
                        'NuMassSignal_Plot','OFF',...
                        'AnchorBkg6G',obj.AnchorBkg6G,...
                        'FPD_ROIlow',obj.FPD_ROIlow,...
                        'BKG_RateSec',obj.BKG_RateSec);
                end
            end
        end
        
        function [SingleBarX, TDRBarX] = GetRangeSingleSys(obj,varargin)
            % Get Range Vector for Single Contributions (MultiBarPlot)
            p = inputParser;
            p.addParameter('PlotRanges',[30, 45, 60],@(x)all(isfloat(x)));
            p.parse(varargin{:});
            PlotRanges = p.Results.PlotRanges;
            % Init Var
            nSingleSys = 3; %for plot only 3: tc, fsd, rf or el, bf, rx
            Nranges = numel(PlotRanges);
            SingleBarX = zeros(Nranges,nSingleSys);
            TDRBarX    = zeros(Nranges,1);
            % Set ranges (manually), found to be good
            for i=1:Nranges
                if Nranges==1
                    SingleBarX(i,:) = PlotRanges(i)-((0:0.135:(nSingleSys-1)/6)+0.47);
                    TDRBarX(i) = PlotRanges(i)+1.5;
                elseif Nranges==2
                    SingleBarX(i,:) = PlotRanges(i)-([0,1.32,2.64]+4.6);
                    TDRBarX(i) = PlotRanges(i)+2;
                elseif Nranges==3
                    SingleBarX(i,:) = PlotRanges(i)-((0:1.32:nSingleSys)+4.6);
                    TDRBarX(i) = PlotRanges(i)+2;
                elseif Nranges==4
                     SingleBarX(i,:) = PlotRanges(i)-((0:1.32:nSingleSys)+4.6);
                    TDRBarX(i) = PlotRanges(i)+2;
                elseif Nranges==5
                     SingleBarX(i,:) = PlotRanges(i)-((0:1.32:nSingleSys)+4.6);
                    TDRBarX(i) = PlotRanges(i)+2;
                end
            end
        end
        function SysUncertainties = GetSysInfo(obj,varargin)
            p = inputParser;
            p.addParameter('SysFlag','SysAll',@(x)ismember(x,{'SysAll','SysRF'}));
            p.parse(varargin{:});
            SysFlag = p.Results.SysFlag;
            
            s = GetSysBudget('SysBudget',obj.SysBudget,'ELossFlag',obj.ELossFlag);
            % Systematic Budget Annotation
            switch SysFlag
                case 'SysAll'
                    SysUncertainties = [sprintf('Response Function: '),...
                        sprintf('\\Delta \\rhod\\sigma = %.1f%%, ',100*sqrt(s.WGTS_CD_MolPerCm2_RelErr^2+s.ISXsection_RelErr^2)),...
                        sprintf('\\DeltaB_a = %s, \\DeltaB_{T/max} = %.1f%%, ',s.MACE_Ba_T_Err_str,100*s.MACE_Bmax_T_RelErr),...
                        s.ELoss,...
                        sprintf('\nFinal State Distribution:'),...
                        sprintf(' Normalization %.0f %%, ',s.FSDNorm_RelErr*100),...
                        sprintf('Bin-to-Bin uncorrelated %.0f %% (GS), %.0f %% (ES)',100*s.FSDShapeGS_RelErr,100*s.FSDShapeES_RelErr)];
                case 'SysRF'
                    SysUncertainties = [sprintf('Response Function: '),...
                        sprintf('\\Delta \\rhod\\sigma = %.1f%%, ',100*sqrt(s.WGTS_CD_MolPerCm2_RelErr^2+s.ISXsection_RelErr^2)),...
                        sprintf('\\DeltaB_a = %s, \\DeltaB_{T/max} = %.1f%%',s.MACE_Ba_T_Err_str,100*s.MACE_Bmax_T_RelErr),...
                        s.ELoss];
            end
        end
    end
end