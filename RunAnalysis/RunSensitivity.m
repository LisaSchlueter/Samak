classdef RunSensitivity < handle
    % ----------------------------------------------------------
    % Class for Sensitivity Studies (Neutrino Mass, Endpoint,...)
    % new class - coupled to RunAnalysis class
    % old class still avaiable under "SensitivityStudy.m"
    %
    % Lisa Schlueter (June 2019)
    %------------------------------------------------------------
   
    properties  (Access=public)
        % Sensitivity properties
        ConfLevel;         % confidence level: e.g. 0.9
        Lpar;              % current sensitivity (all parameters)
        LparMinos;         % current mNuSq sensitivity with MINOS errors (averaged)
        LparNeg
        LparPos
        MultiLpar;         % sensitivity of single effects: stat + TC, stat + FSD,....
        MultiLparMinos;
        MultiLparNeg;
        MultiLparPos;
        MultiLparStack;    % sensitivities N+1: stat, stat+SysEffectAll{1}, stat + +SysEffectAll{1} ++SysEffectAll{2}, .....
        MultiFitResults;   % fit results
        NPcomponent;       % non poissonian component
        AsymErr;           % use asymmetric uncertainties (only minos) or symmetrized 
        
        % method
        Mode % Asimov, Scan, MC, ..
        LimitFlag; % Central or UP %prelimiary! only used for asimov conversion at the moment
        
        % systematics / covariance matrices
        SysEffectsAll;     % cell with systematic effects (used in loop)
        SysEffect;         % struct with effects (used when only 1 sensitivity is calculated)
        SysEffectLeg;      % legend (has to match SysEffectsAll)
        nSys;              % number of systematic effects
        CovMatFrac;        % struct with covariance matrices from MultiLpar
        CovMatFracStack;   % struct with covariance matrices from MultiLparStack
        CovMatShape;
        CovMatShapeStack;
        chi2sys;           %CM or CMShape
        
        RunAnaObj;         % RunAnalysis Object
        RecomputeFlag;     % compute or load sensitivites
        
        % plot
        PlotColor;
        
        % neutrino mass scan
        ScanResults;
        ScanSide; %Pos or Neg
        
        %upper limit: Neyman  confidence interval construction
        NFitResults %Neyman fit results (normal fit with free neutrino mass to samples)
        %N_LimitUp;
        %N_LimitLow;
        N_LimitUp;
        N_mNuSq;
       
        FC_mNuSqFit  % fit neutrino masses
        FC_mNuSqTrue
        FC_DeltaChi2
        FC_DeltaChi2C; % critical delta chi2
        FC_Chi2True;
        FC_x1;
        FC_x2;
        FC_savename;
        FC_PDF;
        FC_CumProb;
    end
    methods
        function obj = RunSensitivity(varargin) % constructor
            fprintf('-----------------Start RunSensitivity contructor ----------------------\n');
        
            p = inputParser;
            p.addParameter('AsymErr','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('ConfLevel',0.9, @(x)isfloat(x));
            p.addParameter('RunAnaObj','', @(x) isa(x,'RunAnalysis') || isa(x,'MultiRunAnalysis'));
            p.addParameter('SysEffectsAll','',@(x)iscell(x));
            p.addParameter('SysEffectLeg','',@(x)iscell(x));
            p.addParameter('SysEffect',struct('FPDeff','ON'),@(x)isstruct(x));
            p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Mode','Asimov',@(x)ismember(x,{'Asimov','Scan'}));
            p.addParameter('chi2sys','chi2CMShape',@(x)ismember(x,{'CM','CMShape'}));
            p.addParameter('ScanSide','Pos',@(x)ismember(x,{'Pos','Neg'}));
            p.addParameter('LimitFlag','Central',@(x)ismember(x,{'Central','Up'}));
            p.parse(varargin{:});
            
            obj.ConfLevel      = p.Results.ConfLevel;
            obj.RunAnaObj      = p.Results.RunAnaObj;
            obj.SysEffectsAll  = p.Results.SysEffectsAll;
            obj.SysEffectLeg   = p.Results.SysEffectLeg;
            obj.SysEffect      = p.Results.SysEffect;
            obj.RecomputeFlag  = p.Results.RecomputeFlag;
            obj.Mode           = p.Results.Mode;
            obj.chi2sys        = p.Results.chi2sys;
            obj.ScanSide       = p.Results.ScanSide;
           
            obj.LimitFlag      = p.Results.LimitFlag;
            obj.AsymErr        = p.Results.AsymErr;
            
            if isempty(obj.RunAnaObj)
                fprintf(2,'Error: RunAnalysis Object is necessary input! \n');
                return
            end
            
            %define colors for plot (preliminary)
            obj.PlotColor = {rgb('White'),rgb('Navy'),rgb('GoldenRod'),rgb('PowderBlue'),...
                rgb('CadetBlue'),rgb('DarkOrange'),rgb('FireBrick'),rgb('DarkSlateGray'),...
                rgb('YellowGreen'),rgb('Magenta'),...
                rgb('SeaGreen'),rgb('DodgerBlue'),rgb('LightGreen')};
            
            if strcmp(obj.RunAnaObj.fitter,'matlab')
                % do 1 init fit
                obj.RunAnaObj.Fit;
            end
            
            if isempty(obj.SysEffectsAll)
                obj.GetDefSysEffect;
            end
            obj.nSys           = numel(obj.SysEffectsAll);
            fprintf('-----------------End RunSensitivity contructor ----------------------\n');
        end
        function GetData(obj)
            % call this function, when changing obj.RunAnaObj.DataType
%             if strcmp(obj.RunAnaObj.DataSet,'FirstTritium.katrin')
%                  obj.RunAnaObj.StackRuns;
%                 obj.RunAnaObj.SimulateStackRuns;
%                 obj.RunAnaObj.ModelObj.ComputeTBDDS; obj.RunAnaObj.ModelObj.ComputeTBDIS;
%                 obj.RunAnaObj.RunData.TBDIS = obj.RunAnaObj.ModelObj.TBDIS;
%                 obj.RunAnaObj.RunData.TBDISE = obj.RunAnaObj.ModelObj.TBDISE;
%             else
                if ismember(obj.RunAnaObj.DataType,'Twin')
                    fprintf('RunSensitivity: Retrieve/Compute Twins \n');
                    if isa(obj.RunAnaObj,'MultiRunAnalysis')
                        obj.RunAnaObj.StackRuns;
                        obj.RunAnaObj.SimulateStackRuns;
                    else
                        obj.RunAnaObj.ReadData;
                        obj.RunAnaObj.SimulateRun;
                    end
                elseif ismember(obj.RunAnaObj.DataType,'Real')
                    fprintf('RunSensitivity: Retrieve Real Data \n');
                    if isa(obj.RunAnaObj,'MultiRunAnalysis')
                        obj.RunAnaObj.StackRuns;
                        obj.RunAnaObj.SimulateStackRuns;
                    else
                        obj.RunAnaObj.ReadData;
                        obj.RunAnaObj.SimulateRun;
                    end
                end
            end
        %end
    end
    methods % sensitivity computation
        function ComputeSensitivity_Asimov(obj,varargin)
            p=inputParser;
            p.addParameter('cl',obj.ConfLevel,@(x)isfloat(x));
            p.addParameter('GetCM','ON',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            cl    = p.Results.cl;
            GetCM = p.Results.GetCM;
            
            fixPar_v = str2num(strrep(strrep(obj.RunAnaObj.fixPar,'fix ',''),' ;',' '));
            if any(ismember(fixPar_v,1))
                fprintf(2,'WARNING: neutrino mass is fixed! \n');
            end
            
            factor = obj.ConvertCl2Std('cl',cl);
            
            nTrials = obj.GetnTrials(obj.SysEffect);
            if ~strcmp(obj.RunAnaObj.chi2,'chi2Stat') && strcmp(GetCM,'ON')
                if strcmp(obj.SysEffect,'Bkg')
                    if strcmp(obj.RunAnaObj.DataSet,'Knm1')
                        BkgMode = 'SlopeFit';
                    else
                        BkgMode = 'Gauss';
                    end
                    obj.RunAnaObj.ComputeCM('SysEffect',struct('FSD','OFF'),...
                        'BkgCM','ON','BkgPtCM','OFF','nTrials',nTrials,'Mode',BkgMode,...
                        'BkgScalingOpt',2,'BkgRingCorrCoeff',0);
                    obj.RunAnaObj.NonPoissonScaleFactor = 1;
                    obj.RunAnaObj.ComputeCM('SysEffect',struct('FSD','OFF'),...
                        'BkgCM','ON','BkgPtCM','OFF','nTrials',nTrials,'Mode',BkgMode,...
                        'BkgScalingOpt',2,'BkgRingCorrCoeff',0);
                elseif strcmp(obj.SysEffect,'BkgPT')
                     [SysErr,~] = GetSysErr(obj.RunAnaObj.SysBudget);
                      obj.RunAnaObj.ComputeCM('SysEffect',struct('FSD','OFF'),...
                        'BkgCM','OFF','BkgPtCM','ON','nTrials',nTrials,...
                        'BKG_PtSlopeErr',SysErr.BKG_PtSlopeErr);
                else
                    [~,CmArg] = GetSysErr(obj.RunAnaObj.SysBudget);
                    obj.RunAnaObj.ComputeCM('SysEffects',obj.SysEffect,CmArg{:},...
                        'BkgCM','OFF','BkgPtCM','OFF''nTrials',nTrials);
                end
            end
            
            % compute sensitvity by fitting asimov dataset
            obj.RunAnaObj.Fit;
            obj.Lpar      = obj.RunAnaObj.FitResult.err.*factor;  %sensitivity all parameter at given confidence level
            if isfield(obj.RunAnaObj.FitResult,'errNeg')
                MeanErr =  0.5*(abs(obj.RunAnaObj.FitResult.errNeg(1))+obj.RunAnaObj.FitResult.errPos(1));
                obj.LparMinos = zeros(size(obj.Lpar));
                obj.LparNeg   = zeros(size(obj.Lpar));
                obj.LparPos   = zeros(size(obj.Lpar));
                %save asymmetric MINOS errors (only for mNu)
                obj.LparMinos(1) = MeanErr*factor; 
                obj.LparNeg(1) = obj.RunAnaObj.FitResult.errNeg(1);
                obj.LparPos(1) = obj.RunAnaObj.FitResult.errPos(1);
            else
                obj.LparMinos = zeros(size(obj.Lpar));
            end
        end
        function ComputeSensitivity_Scan(obj,varargin)
            fprintf('Neutrino Mass Scan - %s Side \n',obj.ScanSide)
            p=inputParser;
            p.addParameter('cl',obj.ConfLevel,@(x)isfloat(x));
            p.addParameter('ScanPrcsn',0.01,@(x)isfloat(x));
            p.addParameter('ScanStart',0.7,@(x)isfloat(x));
            p.addParameter('ScanStop',2,@(x)isfloat(x));
            p.addParameter('nFitMax',20,@(x)isfloat(x));
            p.addParameter('GetCM','ON',@(x)ismember(x,{'ON','OFF'}));
            
            p.parse(varargin{:});
            
            cl            = p.Results.cl;
            ScanPrcsn     = p.Results.ScanPrcsn; %scan precision (defines when loop is exited)
            ScanStart     = p.Results.ScanStart;
            ScanStop      = p.Results.ScanStop;
            nFitMax       = p.Results.nFitMax;
            GetCM         = p.Results.GetCM;
            
            DeltaChi2 = obj.ConvertCl2Std('cl',cl)^2;     % defines confidence level
            
            %init scan par
            ScanmNuSq          = zeros(nFitMax,1);
            ScanmNuSq(2)       = ScanStart;
            
            switch obj.ScanSide
                case 'Pos'
                    ScanmNuSq(nFitMax) = ScanStop;
                case 'Neg'
                    ScanmNuSq(nFitMax) = -ScanStop;
            end
            
            ScanChi2           = zeros(nFitMax,1);
            ScanPar          = zeros(nFitMax,obj.RunAnaObj.nPar);
            ScanErr          = zeros(nFitMax,obj.RunAnaObj.nPar);
            
            fixPar_i = obj.RunAnaObj.fixPar;
            obj.RunAnaObj.fixPar = '1 5 6 7 8 9 10 11'; %parameter of interest has to be fixed in scan
            obj.Lpar = zeros(obj.RunAnaObj.nPar,1);
            if strcmp(obj.RunAnaObj.DataSet,'Knm1')
                BkgMode = 'SlopeFit';
            else
                BkgMode = 'Gauss';
            end
            nTrials = obj.GetnTrials(obj.SysEffect);
            if ~strcmp(obj.RunAnaObj.chi2,'chi2Stat') && strcmp(GetCM,'ON')
                if strcmp(obj.SysEffect,'Bkg')
                    obj.RunAnaObj.ComputeCM('SysEffect',struct('FSD','OFF'),'BkgCM','ON','BkgMode',BkgMode,...
                        'nTrials',nTrials);
                else
                    [~,CmArg] = GetSysErr(obj.RunAnaObj.SysBudget);
                    obj.RunAnaObj.ComputeCM('SysEffects',obj.SysEffect,CmArg{:},...
                        'BkgCM','OFF','nTrials',nTrials,'BkgMode',BkgMode);
                end
            end
            
            %% Test if nu-mass=0 (expected MC value) gives small chi2
            obj.RunAnaObj.ModelObj.mnuSq_i = ScanmNuSq(1);
            obj.RunAnaObj.Fit;
            ScanChi2(1)  = obj.RunAnaObj.FitResult.chi2min;
            ScanPar(1,:) = obj.RunAnaObj.FitResult.par;
            ScanErr(1,:) = obj.RunAnaObj.FitResult.err;
            if ScanChi2(1)>1
                fprintf('ScanStart (%.2f) gives large chi2 \\chi2 = %.1f! \n',ScanStart,ScanChi2(1));
                return
            end
            
            %% Test if stop value for scan Stop is large enough
            obj.RunAnaObj.ModelObj.mnuSq_i = ScanmNuSq(end);
            obj.RunAnaObj.Fit;
            ScanChi2(end) = obj.RunAnaObj.FitResult.chi2min;
            ScanPar(end,:) = obj.RunAnaObj.FitResult.par;
            ScanErr(end,:) = obj.RunAnaObj.FitResult.err;
            if  ScanChi2(end)<(DeltaChi2+ScanChi2(1))
                fprintf('ScanStop (%.2f) is not large enough \\chi2 = %.1f! \n',ScanStop,ScanChi2(end));
                return
            end
            
            %% Start Scan
            for i=2:nFitMax
                if strcmp(obj.ScanSide,'Neg')
                    ScanmNuSq(i) = -ScanmNuSq(i);
                end
                obj.RunAnaObj.ModelObj.mnuSq_i = ScanmNuSq(i);
                obj.RunAnaObj.Fit;
                ScanChi2(i) = obj.RunAnaObj.FitResult.chi2min;
                ScanPar(i,:) = obj.RunAnaObj.FitResult.par;
                ScanErr(i,:) = obj.RunAnaObj.FitResult.err;
                if abs(ScanChi2(i)-(DeltaChi2+ScanChi2(1)))<ScanPrcsn %exit condition
                    fprintf('Exist Scan, %.0f%% C.L. reached \n',obj.ConfLevel*100);
                    ScanChi2(i+1:end-1) = NaN;
                    ScanmNuSq(i+1:end-1) = NaN;
                    ScanPar(i+1:end-1,:) = NaN;
                    ScanErr(i+1:end-1,:) = NaN;
                    break
                elseif ScanChi2(i)<(DeltaChi2+ScanChi2(1))
                    mNuSqLarger= min(abs(ScanmNuSq(abs(ScanmNuSq)>abs(ScanmNuSq(i))))); % find next higher mNuSq
                    ScanmNuSq(i+1) =  0.5*(abs(ScanmNuSq(i))+mNuSqLarger); % find middle between current and next higher one
                elseif ScanChi2(i)>(DeltaChi2+ScanChi2(1))
                    mNuSqSmaller= max(abs(ScanmNuSq(abs(ScanmNuSq)<abs(ScanmNuSq(i))))); % find next smaller mNuSq
                    ScanmNuSq(i+1) =  0.5*(abs(ScanmNuSq(i))+mNuSqSmaller);  % find middle between current and next smaller one
                elseif i==nFitMax
                    fprintf('Sensitivity not found after nFitMax fits - enlarge nFitMax or give different start value \n');
                    break
                end
            end
            
            % Compute Sensitivity
            obj.Lpar(1) =interp1(ScanChi2(~isnan(ScanChi2)),ScanmNuSq(~isnan(ScanChi2)),DeltaChi2+ScanChi2(1),'spline');  %sensitivity
            obj.Lpar = obj.Lpar';
            mNuSqmin = interp1(ScanChi2(~isnan(ScanChi2)),ScanmNuSq(~isnan(ScanChi2)),min(ScanChi2(~isnan(ScanChi2))),'spline'); % mass^2 with minimal chi2
            
            % fill structure
            obj.ScanResults = struct('mnuSq_i_Fit',ScanmNuSq,...% neutrino mass squared scan vector
                'par',ScanPar,... % Fit Parameter
                'err',ScanErr,... % Error on Fit Parameter
                'chi2min',ScanChi2,...
                'dof',repmat(obj.RunAnaObj.FitResult.dof,1,numel(ScanChi2)),...
                'mNuSqmin',mNuSqmin,...
                'DeltaChi2',DeltaChi2+ScanChi2(1),...
                'ConfLevel',cl,...
                'chi2',obj.RunAnaObj.chi2,...
                'SysEffect','',...
                'ScanPrcsn',ScanPrcsn,...
                'ScanSide',obj.ScanSide);
            
            obj.RunAnaObj.fixPar =  fixPar_i;
        end
    end
    methods % upper limit computation
        function [UpperLimit,mNuSq_Quantil90,mNuSq_t] = ComputeUpperLimit(obj,varargin)
            % calculate upper limit with Neyman confidence interval
            % construction
            p=inputParser;
            p.addParameter('mNuSq_t',0:0.2:2,@(x)isfloat(x));
            p.addParameter('nSamples',500,@(x)isfloat(x));
            p.addParameter('ReFit','OFF',@(x)isfmember(x,{'ON','OFF'}));
            p.addParameter('SanityPlot','OFF',@(x)isfmember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            mNuSq_t     = p.Results.mNuSq_t; % true neutrino masses in simulation
            nSamples    = p.Results.nSamples;
            ReFit       = p.Results.ReFit;
            SanityPlot  = p.Results.SanityPlot;
            
            % sanity check: neutrino mass should be free
            fixPar_v = str2num(strrep(strrep(obj.RunAnaObj.fixPar,'fix ',''),' ;',' '));
            if any(ismember(fixPar_v,1))
                fprintf(2,'WARNING: neutrino mass is fixed! \n');
            end
            
            obj.RunAnaObj.Fit;
            mNuSqFit = obj.RunAnaObj.FitResult.par(1);
            
            obj.NFitResults = obj.ComputeFitNuMassTwins('mNuSq',mNuSq_t,'nSamples',nSamples,'ReFit',ReFit,'FC','OFF');
            
            mNuSq_Quantil90 = zeros(numel(mNuSq_t),1);
            mNuSq = squeeze(obj.NFitResults.par(:,1,:))';
            DeltamNuSq = abs(mNuSq-repmat(mNuSq_t,1000,1));
            mNuSq_Corr = cell(numel(mNuSq_t),1);
            for i=1:numel(mNuSq_t)
                mNuSq_Corr{i} = mNuSq(DeltamNuSq(:,i)<50,i);
                mNuSq_Quantil90(i) = prctile(mNuSq_Corr{i},100*(1-obj.ConfLevel));
            end
            
            %  mNuSq_QuantilLow =  prctile(mNuSq,100*(1-obj.ConfLevel)./2);
            %  mNuSq_QuantilUp =  prctile(mNuSq,100*(1+obj.ConfLevel)./2);
            
            % nu-mass quantil: when experiment is repeated many times: 90% of the
            % measured neutrino masses are smaller than this value
            % mNuSq_Quantil90 = cell2mat(cellfun(@(x) x(obj.ConfLevel),obj.NFitResults.Percentiles,'UniformOutput',0));
            
            %             LimitUp = zeros(numel(mNuSq_t),1);
            %             LimitLow = zeros(numel(mNuSq_t),1);
            %             for i=1:numel(mNuSq_t)
            %                 if mNuSq_QuantilLow(i)>0
            %                     LimitLow(i) = mNuSq_QuantilLow(i);
            %                     LimitUp(i)  = mNuSq_QuantilUp(i);
            %                 else
            %                      LimitLow(i) = NaN;
            %                      LimitUp(i) = mNuSq_Quantil90(i);
            %                 end
            %             end
            %             obj.N_LimitLow = LimitLow;
            %             obj.N_LimitUp  = LimitUp;
            %             obj.N_mNuSq    = mNuSq_t;
            
            UpperLimit = interp1(mNuSq_Quantil90,mNuSq_t,mNuSqFit);
            
            obj.N_LimitUp = UpperLimit;
            %obj.N_CentralUp = interp1(mNuSq_QuantilUp,mNuSq_t,mNuSqFit);
            %obj.NCentralLow =  interp1(mNuSq_QuantilLow,mNuSq_t,mNuSqFit);
            if strcmp(SanityPlot,'ON')
                obj.PlotsUpperLimit
            end
            
        end
        function ComputeFC(obj,varargin)
            % test other FC implementation
            p=inputParser;
            p.addParameter('mNuSq_t',0.2,@(x)isfloat(x));
            p.addParameter('nSamples',3,@(x)isfloat(x));
            p.addParameter('ReFit','OFF',@(x)isfmember(x,{'ON','OFF'}));
            p.addParameter('SanityPlot','ON',@(x)isfmember(x,{'ON','OFF'}));
            p.addParameter('DummyGaussPdf','OFF',@(x)isfmember(x,{'ON','OFF'})); % use gaussian distribution instead of MC fits
            
            p.parse(varargin{:});
            mNuSq_t     = p.Results.mNuSq_t; % true neutrino masses in simulation
            nSamples    = p.Results.nSamples;
            ReFit       = p.Results.ReFit;
            SanityPlot  = p.Results.SanityPlot;
            DummyGaussPdf = p.Results.DummyGaussPdf;
            
            % 1. Get neutrino-mass pdfs for different model nu-masses
            switch DummyGaussPdf
                case 'OFF'
                    MClabel = 'TwinMC';
                    obj.ComputeFitNuMassTwins('mNuSq',mNuSq_t,'nSamples',nSamples,'ReFit',ReFit);
                    mNuSqFit = permute(obj.NFitResults.par(:,1,:),[2,3,1]);
                case 'ON'
                    MClabel = 'TwinDummyMC';
                    mNuSqFit = obj.ComputeDummyMC('mNuSq',mNuSq_t,'nSamples',nSamples);
            end
            
            % 2. calculate Delta Chi2 for each true mNu and each sample
            % label
            savedir = [getenv('SamakPath'),sprintf('tritium-data/fit/%s%s/',MClabel,obj.RunAnaObj.DataSet)];
            obj.CreateDir(savedir);
            
            % save previous setting
            tmp_fixPar = obj.RunAnaObj.fixPar;
            
            % init
            obj.RunAnaObj.fixPar = '1 5 6 7 8 9 10 11'; % fix neutrino mass
            Chi2TrueAll  = zeros(numel(mNuSq_t),nSamples);
            Chi2BestAll  = zeros(numel(mNuSq_t),nSamples);
            DeltaChi2All = zeros(numel(mNuSq_t),nSamples);
            
            LoadLogic = zeros(numel(mNuSq_t),1); % incidates which files were loaded and which new computed
            
            for i=1:numel(mNuSq_t)
                thisfile =  [savedir,obj.FC_savename{i}];
                
                % load if possible
                if exist(thisfile,'file') && strcmp(ReFit,'OFF')
                    d = importdata(thisfile);
                else
                    d='';
                end
                
                if isfield(d,'Chi2True')
                    Chi2TrueAll(i,:) = d.Chi2True;
                    Chi2BestAll(i,:) = d.Chi2Best;
                    DeltaChi2All(i,:) = d.DeltaChi2;
                    LoadLogic(i) = 1;
                else
                    
                    for s=1:nSamples
                        progressbar(s/nSamples)
                        % 2. compute Chi2(x_measured,mu_true):
                        obj.RunAnaObj.SimulateStackRuns('mNuSq_i',mNuSq_t(i))
                        obj.RunAnaObj.ModelObj.ComputeTBDDS; obj.RunAnaObj.ModelObj.ComputeTBDIS;
                        obj.RunAnaObj.RunData.TBDIS = obj.RunAnaObj.ModelObj.TBDIS;
                        obj.RunAnaObj.RunData.TBDISE = sqrt(obj.RunAnaObj.ModelObj.TBDIS);
                        obj.RunAnaObj.ModelObj.mnuSq_i = mNuSqFit(1,s,i);
                        obj.RunAnaObj.Fit;
                        Chi2TrueAll(i,s) = obj.RunAnaObj.FitResult.chi2min;
                        obj.RunAnaObj.ModelObj.SetFitBias(0);
                        
                        % compute chi2(x_measured, x_measured OR x=0 if negative)
                        if  obj.NFitResults.par(i,1,s)>0
                            Chi2BestAll(i,s) = 0; % very close to zero anyway
                        elseif obj.NFitResults.par(i,1,s)<0
                            obj.RunAnaObj.SimulateStackRuns('mNuSq_i',0)
                        end
                        obj.RunAnaObj.ModelObj.ComputeTBDDS; obj.RunAnaObj.ModelObj.ComputeTBDIS;
                        obj.RunAnaObj.RunData.TBDIS = obj.RunAnaObj.ModelObj.TBDIS;
                        obj.RunAnaObj.RunData.TBDISE = sqrt(obj.RunAnaObj.ModelObj.TBDIS);
                        obj.RunAnaObj.ModelObj.mnuSq_i = mNuSqFit(1,s,i);
                        obj.RunAnaObj.Fit;
                        Chi2BestAll(i,s) = obj.RunAnaObj.FitResult.chi2min;
                        obj.RunAnaObj.ModelObj.SetFitBias(0);
                    end
                    
                    DeltaChi2All(i,:) = Chi2TrueAll(i,:)-Chi2BestAll(i,:);
                    
                    Chi2True  = Chi2TrueAll(i,:);
                    Chi2Best  = Chi2BestAll(i,:);
                    DeltaChi2 = DeltaChi2All(i,:);
                    save(thisfile,'Chi2True','Chi2Best','DeltaChi2','-append');
                end
            end
            
            
            obj.RunAnaObj.fixPar = tmp_fixPar;
            %   plot(mNuSqFit,Chi2BestAll-Chi2TrueAll,'x');
            
            obj.FC_mNuSqFit     = squeeze(mNuSqFit)'; %squeeze(obj.NFitResults.par(:,1,:)); % fit neutrino masses
            obj.FC_mNuSqTrue    = mNuSq_t'; % true neutrino masses
            obj.FC_DeltaChi2    = squeeze(DeltaChi2All);
            
            % obj.FC_ComputeX1X2;
        end
        
        function mNuSqFitAll = ComputeDummyMC(obj,varargin)
            p = inputParser;
            p.addParameter('mNuSq_t',0,@(x)isfloat(x));
            p.addParameter('nSamples',1000,@(x)isfloat(x));
            
            p.parse(varargin{:});
            
            mNuSq_t    = p.Results.mNuSq_t;
            nSamples = p.Results.nSamples;
            
            tmpDataType = obj.RunAnaObj.DataType;
            tmpfixPar = obj.RunAnaObj.fixPar;
            tmpmNuSqBias = obj.RunAnaObj.TwinBias_mnuSq;
            
            obj.RunAnaObj.DataType = 'Twin';
            obj.RunAnaObj.fixPar = '5 6 7 8 9 10 11';
            
            %label
            MClabel ='TwinDummyMC';
            savedir = [getenv('SamakPath'),sprintf('tritium-data/fit/%s%s/',MClabel,obj.RunAnaObj.DataSet)];
            obj.CreateDir(savedir);
            savename = sprintf('%s%s_%s_%.0fbE0_fixPar%s_%.0fsamples.mat',MClabel,obj.RunAnaObj.RunData.RunName,obj.RunAnaObj.chi2,...
                obj.GetRange,strrep(strrep(strrep(obj.RunAnaObj.fixPar,'fix ',''),' ;',' '),' ',''),nSamples);
            mNuSqFitAll = zeros(1,nSamples,numel(mNuSq_t));
            
            for i=1:numel(mNuSq_t)
                thisfile =  [savedir,sprintf('%.3geV2mNuSq_',mNuSq_t(i)),savename];
                
                if exist(thisfile,'file')
                    mNuSqFitAll(1,:,i) = importdata(thisfile);
                else
                    % get 1 sigma error
                    obj.RunAnaObj.TwinBias_mnuSq = mNuSq_t(i);
                    obj.GetData; %retrieve correct twins
                    obj.RunAnaObj.Fit;
                    mNuSqFitAll(1,:,i) = obj.RunAnaObj.FitResult.err(1)*randn(nSamples,1)+mNuSq_t(i);
                    
                    mNuSqFit =  mNuSqFitAll(1,:,i);
                    save(thisfile,'mNuSqFit');
                end
            end
            
            %reset to init
            obj.RunAnaObj.DataType        = tmpDataType;
            obj.RunAnaObj.TwinBias_mnuSq  = tmpmNuSqBias;
            obj.RunAnaObj.fixPar          = tmpfixPar;
            
        end
        function out = ComputeFitNuMassTwins(obj,varargin)
            % generate sets of twins with neutrino mass
            % fit these with fluctuations
            % if not previously computed: take a lot of time -> run on server
            p=inputParser;
            p.addParameter('mNuSq_t',0:0.2:2,@(x)isfloat(x));
            p.addParameter('nSamples',500,@(x)isfloat(x));
            p.addParameter('ReFit','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            mNuSq_t    = p.Results.mNuSq_t;
            nSamples = p.Results.nSamples;
            ReFit    = p.Results.ReFit;
            
            % save previous settings and reset at the end
            tmpDataType = obj.RunAnaObj.DataType;
            tmpTwinBias = obj.RunAnaObj.TwinBias_mnuSq;
            
            % init
            parAll      = zeros(numel(mNuSq_t),obj.RunAnaObj.nPar,nSamples);
            errAll      = zeros(numel(mNuSq_t),obj.RunAnaObj.nPar,nSamples);
            chi2minAll  = zeros(numel(mNuSq_t),nSamples);
            dofAll      = zeros(numel(mNuSq_t),1);
            
            % label
            savedir = [getenv('SamakPath'),sprintf('tritium-data/fit/TwinMC%s/',obj.RunAnaObj.DataSet)];
            obj.CreateDir(savedir);
            fix_str = ConvertFixPar('freePar',obj.RunAnaObj.fixPar,'Mode','Reverse');%strrep(strrep(strrep(obj.RunAnaObj.fixPar,'fix ',''),' ;',' '),' ','');
           
            savename = sprintf('%s_%s_%.0fbE0_freePar%s_%.0fsamples.mat',obj.RunAnaObj.RunData.RunName,obj.RunAnaObj.chi2,...
                obj.GetRange,fix_str,nSamples);
            LoadLogic = zeros(numel(mNuSq_t),1); % incidates which files were loaded and which new computed
            
            obj.FC_savename = cell(numel(mNuSq_t),1);
            obj.RunAnaObj.DataType = 'Twin';
            
            for i=1:numel(mNuSq_t)
                obj.FC_savename{i} = [sprintf('TwinMC%.3geV2mNuSq_',mNuSq_t(i)),savename];
                thisfile =  [savedir,obj.FC_savename{i}];
                
                if exist(thisfile,'file') && strcmp(ReFit,'OFF')
                    load(thisfile,'par','err','chi2min','dof');
                    LoadLogic(i) = 1;
                    fprintf('Loading %s from file \n',thisfile);
                else
                    obj.RunAnaObj.TwinBias_mnuSq = mNuSq_t(i);
                    obj.GetData; %retrieve correct twins
                    [par, err, chi2min, dof] = obj.RunAnaObj.FitTwin('nSamples',nSamples);
                    
                    mNuSq = mNuSq_t(i);
                    save(thisfile,'par','err','chi2min','dof','mNuSq');
                end
                parAll(i,1:size(par,1),:)      = par;
                errAll(i,1:size(par,1),:)      = err;
                chi2minAll(i,:)    = chi2min;
                dofAll(i)          = dof;
            end
            
            % save in structure
            obj.NFitResults             = struct('par',parAll);
            obj.NFitResults.err         = errAll;
            obj.NFitResults.chi2min     = chi2minAll;
            obj.NFitResults.dof         = dofAll;
            obj.NFitResults.mNuSq_t     = mNuSq_t;         % true neutrino masses in simulations
            out = obj.NFitResults;
            
            if sum(LoadLogic)~=numel(mNuSq_t)
                %  if something was new computed: set back to init
                obj.RunAnaObj.TwinBias_mnuSq = tmpTwinBias;
                obj.RunAnaObj.DataType=tmpDataType;
                obj.GetData;
            end
        end
        function FC_ComputeX1X2(obj,varargin)
            % feldman cousins:
            % compute [x1,x2]
            x1 = zeros(numel(obj.FC_mNuSqTrue),1);
            x2 = zeros(numel(obj.FC_mNuSqTrue),1);
            % 1. find critical delta chi2
            DeltaChi2C = prctile(obj.FC_DeltaChi2',obj.ConfLevel*100);
            
            for i=1:numel(obj.FC_mNuSqTrue)
                % find lower acceptance region x1
                IndexLow = obj.FC_mNuSqFit(i,:)<=obj.FC_mNuSqTrue(i); % lower neutrino masses
                DeltaChi2Low = obj.FC_DeltaChi2(i,IndexLow);
                mNuSqLow = obj.FC_mNuSqFit(i,IndexLow);
                
                [DeltaChi2Low, index_new, ~] = unique(DeltaChi2Low,'stable');
                
                if numel(index_new)==1 || all(DeltaChi2Low<=DeltaChi2C(i)) % if all lower Delta chi2 are smaller than critical value
                    x1(i) = NaN;
                else
                    x1(i) = interp1(DeltaChi2Low ,mNuSqLow,DeltaChi2C(i),'spline');
                end
                
                % find higher acceptance region x2
                IndexUp = obj.FC_mNuSqFit(i,:)>obj.FC_mNuSqTrue(i); % upper neutrino masses
                DeltaChi2Up = obj.FC_DeltaChi2(i,IndexUp);
                mNuSqUp = obj.FC_mNuSqFit(i,IndexUp);
                x2(i) = interp1(DeltaChi2Up,mNuSqUp,DeltaChi2C(i),'spline');
            end
            obj.FC_x1 = x1;
            obj.FC_x2 = x2;
        end 
        function DeltaChi2 = FC_GetDeltaChi2(obj,varargin)
            
            savedir = [getenv('SamakPath'),'tritium-data/FC/DeltaChi2LookupTable/'];
            MakeDir(savedir)
            %   savename = sprintf('DeltaChi2LookupTables_TwinMC%.3geV2mNuSq_%s_%s_%.0fbE0_fixPar%s',mNuSq_t,obj.RunAnaObj.RunData.RunName,obj.RunAnaObj.chi2,...
            %      obj.GetRange,strrep(strrep(strrep(obj.RunAnaObj.fixPar,'fix ',''),' ;',' '),' ',''));
            DeltaChi2 = zeros(size(obj.FC_mNuSqFit));
            
            for i=1:numel(obj.FC_mNuSqTrue)
                savename = ['DeltaChi2LookupTables_',obj.FC_savename{i}];
                % find files
                mydir = arrayfun(@(x) x.name,dir(savedir),'UniformOutput',0);
                Index = cell2mat(cellfun(@(x) contains(x,savename),mydir,'UniformOutput',0));
                myfiles = mydir(Index);
                % check if file has enough samples
                nSamples = str2double(extractBetween(myfiles,[strrep(strrep(strrep(obj.RunAnaObj.fixPar,'fix ',''),' ;',' '),' ',''),'_'],'samples'));
                
                if max(nSamples)>=100
                    thisfile = [savedir,myfiles{nSamples==max(nSamples)}];
                    d = importdata(thisfile);
                    DeltaChi2(i,:)=interp1(d.mNuSqFit,d.DeltaChi2,obj.FC_mNuSqFit(i,:),'spline');
                else
                    %return
                    %calculate delta chi2 and save
                end
                
            end
        end
        function [mNuSq,DeltaChi2,Chi2True,Chi2Best] = FC_ComputeDeltaChi2LookupTables(obj,varargin)
            p=inputParser;
            p.addParameter('mNuSq_t',0,@(x)isfloat(x)); % only 1
            p.addParameter('nSamples',100,@(x)isfloat(x));
            p.addParameter('SanityPlot','ON',@(x)isfmember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            mNuSq_t    = p.Results.mNuSq_t;
            nSamples   = p.Results.nSamples;
            SanityPlot = p.Results.SanityPlot;
            
            % label
            fix_str = ConvertFixPar('freePar',obj.RunAnaObj.fixPar,'Mode','Reverse');%strrep(strrep(strrep(obj.RunAnaObj.fixPar,'fix ',''),' ;',' '),' ','');
            savedir = [getenv('SamakPath'),'tritium-data/FC/DeltaChi2LookupTable/'];
            if strcmp(obj.RunAnaObj.chi2,'chi2Stat')
                chiStr = sprintf('%s',obj.RunAnaObj.chi2);
            else
                 chiStr = sprintf('%s_SysBudget%.0f',obj.RunAnaObj.chi2,obj.RunAnaObj.SysBudget);
            end
            
            if strcmp(obj.RunAnaObj.DataSet,'Knm1')
                RunName = 'KNM1';
                
                if strcmp(obj.RunAnaObj.DataSet,'Knm1') && obj.RunAnaObj.TwinBias_BKG_PtSlope~=0
                     save_str = sprintf('AsimovDeltaChi2_mNuSq%.3geV2_%s_%s_%.0fbE0_BkgPT%.3gmuCpsS_freePar%s_%.0fsamples.mat',...
                        mNuSq_t,RunName,chiStr,obj.GetRange,1e6.*obj.RunAnaObj.TwinBias_BKG_PtSlope,fix_str,nSamples);
                %    save_str = strrep(save_str,'.mat',sprintf('_BkgPT%.3gmuCpsS.mat',1e6.*obj.RunAnaObj.TwinBias_BKG_PtSlope));
                else
                    save_str = sprintf('AsimovDeltaChi2_mNuSq%.3geV2_%s_%s_%.0fbE0_freePar%s_%.0fsamples.mat',...
                        mNuSq_t,RunName,chiStr,obj.GetRange,fix_str,nSamples);
                end
            else
                RunName = obj.RunAnaObj.RunData.RunName;
                save_str = sprintf('AsimovDeltaChi2_mNuSq%.3geV2_%s_%s_%.0fbE0_freePar%s_%.0fsamples.mat',...
                    mNuSq_t,RunName,chiStr,obj.GetRange,fix_str,nSamples);
            end
            
            
            if strcmp(obj.RunAnaObj.AnaFlag,'Ring')
                save_str = strrep(save_str,'.mat',sprintf('_Ring%s.mat',obj.RunAnaObj.RingMerge));
            end
            
           
            
            savefile = [savedir,save_str];
            
            % look for file with same or more samples than required
            tmp_dir = arrayfun(@(x) x.name,dir(savedir),'UniformOutput',0); % content of directoy
            Index = cell2mat(cellfun(@(x) contains(x,extractBefore(save_str,[fix_str,'_'])),tmp_dir,'UniformOutput',0)); % files that are similar to save_str
            myfiles = tmp_dir(Index);
            
            % check if file has enough samples
            FileSamples = str2double(extractBetween(myfiles,[fix_str,'_'],'samples'));
            
            if max(FileSamples)>=nSamples
                thisfile = [savedir,myfiles{FileSamples==max(FileSamples)}];
                d = importdata(thisfile);
                DeltaChi2 = d.DeltaChi2;
                Chi2True  = d.Chi2True;
                Chi2Best  = d.Chi2Best;
                mNuSq     = d.mNuSq;
                fprintf('Loading Delta Chi2 from file %s \n',thisfile)
            else
                % save previous setting
                tmp_fixPar = obj.RunAnaObj.fixPar;
                
                % Get Correct Normalization and Background
                obj.RunAnaObj.DataType = 'Real';
                obj.RunAnaObj.StackRuns;
                obj.RunAnaObj.InitModelObj_Norm_BKG;
                FitBkg =obj.RunAnaObj.ModelObj.BKG_RateSec;
                FitN  = obj.RunAnaObj.ModelObj.normFit;
                InitBkg = {'BKG_RateAllFPDSec',FitBkg};
                InitN = {'N_bias',FitN};
                
                %                 % get 1 sigma interval
                %                 obj.RunAnaObj.fixPar = '5 6 7 8 9 10 11';   % free neutrino mass
                %                 obj.RunAnaObj.SimulateStackRuns('mNuSq_i',mNuSq_t,InitBkg{:})
                %                 obj.RunAnaObj.ModelObj.ComputeTBDDS(InitN{:});
                %                 obj.RunAnaObj.ModelObj.ComputeTBDIS;
                %                 obj.RunAnaObj.RunData.TBDIS = obj.RunAnaObj.ModelObj.TBDIS;
                %                 obj.RunAnaObj.RunData.TBDISE = sqrt(obj.RunAnaObj.ModelObj.TBDIS);
                %                 if ~strcmp(obj.RunAnaObj.chi2,'chi2Stat')
                %                     obj.RunAnaObj.ComputeCM;
                %                 end
                %                 obj.RunAnaObj.Fit;
                obj.RunAnaObj.fixPar = 'E0 Norm Bkg'; % fix neutrino mass for all other fits
                obj.RunAnaObj.InitFitPar;
                % set up samples mNuSq: true nu-mass +- 15 sigma
                % stepsize increasing the closer to mNuSq_t
                %                mNuSqSigma = obj.RunAnaObj.FitResult.err(1);
                %                 mNuSq_4 = linspace(mNuSq_t-15*mNuSqSigma,mNuSq_t+15*mNuSqSigma,nSamples/4);
                %                 mNuSq_3 = linspace(mNuSq_t-7*mNuSqSigma,mNuSq_t+10*mNuSqSigma,nSamples/4);
                %                 mNuSq_2 = linspace(mNuSq_t-3*mNuSqSigma,mNuSq_t+3*mNuSqSigma,nSamples/4);
                %                 mNuSq_1 = linspace(mNuSq_t-1*mNuSqSigma,mNuSq_t+1*mNuSqSigma,nSamples/4);
                
                mNuSq = linspace(mNuSq_t-10,mNuSq_t+10,nSamples);
                %mNuSq_2 = linspace(mNuSq_t-5,mNuSq_t+5,nSamples/4);
                %mNuSq_1 = linspace(mNuSq_t-2,mNuSq_t+2,nSamples/4);
                %mNuSq = sort([mNuSq_1,mNuSq_2,mNuSq_3]);
                % true neutrino mass shall be inside test mNuSq
                TestmNuSqTrue = sum(mNuSq == mNuSq_t);
                if TestmNuSqTrue==0
                    tmpIndex = find(abs(mNuSq-mNuSq_t)==min(abs(mNuSq-mNuSq_t)),1);
                    mNuSq(tmpIndex) = mNuSq_t; %replace closest value with true neutrino mass
                end
                % get rid of redundant points
                [mNuSq, ~,~] = unique(mNuSq);
                % fill up again to nSamples with some outliers
                nDiff = nSamples-numel(mNuSq);
                mNuSqDiff_1 = linspace(mNuSq_t-15,mNuSq_t-19.01,nDiff/2);
                mNuSqDiff_2 = linspace(mNuSq_t+10.01,mNuSq_t+15,nDiff/2);
                mNuSq = sort([mNuSq,mNuSqDiff_1,mNuSqDiff_2]);
                
                Chi2True = zeros(nSamples,1);
                Chi2Best = zeros(nSamples,1);
                
                progressbar('Compute Delta Chi2 LookupTables (Asimov)')
                for s=1:nSamples
                    progressbar(s/nSamples)
                    % 2. compute Chi2(x_measured,mu_true):
                    obj.RunAnaObj.SimulateStackRuns('mNuSq_i',mNuSq_t,InitBkg{:})
                    obj.RunAnaObj.ModelObj.ComputeTBDDS(InitN{:});
                    obj.RunAnaObj.ModelObj.ComputeTBDIS;
                    obj.RunAnaObj.RunData.TBDIS = obj.RunAnaObj.ModelObj.TBDIS;
                    obj.RunAnaObj.RunData.TBDISE = sqrt(obj.RunAnaObj.ModelObj.TBDIS);
                    obj.RunAnaObj.ModelObj.mnuSq_i = mNuSq(s);
                    obj.RunAnaObj.Fit;
                    Chi2True(s) = obj.RunAnaObj.FitResult.chi2min;
                    obj.RunAnaObj.ModelObj.SetFitBias(0);
                    
                    % compute chi2(x_measured, x_measured OR x=0 if negative)
                    if  mNuSq(s)>0
                        Chi2Best(s) = 0; % very, very close to zero anyway
                    elseif mNuSq(s)<0
                        obj.RunAnaObj.SimulateStackRuns('mNuSq_i',0,InitBkg{:})
                        obj.RunAnaObj.ModelObj.ComputeTBDDS(InitN{:});
                        obj.RunAnaObj.ModelObj.ComputeTBDIS;
                        obj.RunAnaObj.RunData.TBDIS = obj.RunAnaObj.ModelObj.TBDIS;
                        obj.RunAnaObj.RunData.TBDISE = sqrt(obj.RunAnaObj.ModelObj.TBDIS);
                        obj.RunAnaObj.ModelObj.mnuSq_i = mNuSq(s);
                        obj.RunAnaObj.Fit;
                        Chi2Best(s) = obj.RunAnaObj.FitResult.chi2min;
                        obj.RunAnaObj.ModelObj.SetFitBias(0);
                    end
                end
                
                DeltaChi2 = Chi2True-Chi2Best;
                save(savefile,'Chi2True','Chi2Best','DeltaChi2','mNuSq','mNuSq_t');
                
                %reset to previous setting
                obj.RunAnaObj.fixPar = tmp_fixPar;
            end
        end
        function ComputeFC_Asimov(obj,varargin)
            p=inputParser;
            p.addParameter('mNuSq_t',0.2,@(x)isfloat(x));
            p.addParameter('nSamples',10000,@(x)isfloat(x));
            p.addParameter('nSamplesAsimov',300,@(x)isfloat(x));
            p.addParameter('ReFit','OFF',@(x)isfmember(x,{'ON','OFF'}));
            
            p.parse(varargin{:});
            mNuSq_t        = p.Results.mNuSq_t; % true neutrino masses in simulation
            nSamples       = p.Results.nSamples;
            nSamplesAsimov = p.Results.nSamplesAsimov;
            
            spline = 'ON';
            switch spline
                case 'ON'
                    interpStyle = 'spline';
                case 'OFF'
                    interpStyle = 'lin';
            end
            
            %init
            obj.FC_DeltaChi2  = zeros(numel(mNuSq_t),nSamples);
            obj.FC_DeltaChi2C = zeros(numel(mNuSq_t),1);
            obj.FC_mNuSqFit   = zeros(numel(mNuSq_t),nSamples);
            obj.FC_x1         = zeros(numel(mNuSq_t),1);
            obj.FC_x2         = zeros(numel(mNuSq_t),1);
            obj.FC_PDF        = zeros(numel(mNuSq_t),nSamples);
            obj.FC_CumProb    = zeros(numel(mNuSq_t),nSamples);
            
            for i=1:numel(mNuSq_t)
                % Get Delta Chi2 Curve
                [mNuSq,DeltaChi2,Chi2True,Chi2Best] = obj.FC_ComputeDeltaChi2LookupTables('mNuSq_t',mNuSq_t(i),'nSamples',nSamplesAsimov);
                
                % Convert Chi2True into Probability
                Prob_tmp = exp(-Chi2True./2);
%                 % tmp start
%                 NormProb1 = 0.5.*simpsons(mNuSq(mNuSq<=mNuSq_t(i)),Prob_tmp(mNuSq<=mNuSq_t(i)));
%                 NormProb2 = 0.5.*simpsons(mNuSq(mNuSq>mNuSq_t(i)),Prob_tmp(mNuSq>mNuSq_t(i)));
%                 Prob_tmp(mNuSq<=mNuSq_t(i)) = Prob_tmp(mNuSq<=mNuSq_t(i))./NormProb1;
%                 Prob_tmp(mNuSq>mNuSq_t(i)) = Prob_tmp(mNuSq>mNuSq_t(i))./NormProb2;
%                 Prob = Prob_tmp;
%                 %end
                Prob = Prob_tmp./simpsons(mNuSq,Prob_tmp); % normalization
                CumProb = GetCDF(mNuSq,Prob);%CumProb = cumsum(Prob); CumProb = CumProb./max(CumProb);
                
                % solve (numerically) after x1 and x2
                options = optimoptions(@fsolve,'Display','off');
                if mNuSq_t(i)~=0
                    a =  @(x) (interp1(mNuSq,CumProb,x(2),interpStyle)-interp1(mNuSq,CumProb,x(1),interpStyle))-obj.ConfLevel;
                    b =  @(x)  (interp1(mNuSq,DeltaChi2,x(1),interpStyle)-interp1(mNuSq,DeltaChi2,x(2),interpStyle));
                    fun = @(x) [a(x),b(x)];
                    x0 = [-1,1];
                    x1x2 = fsolve(fun,x0,options);
                    obj.FC_x1(i) = x1x2(1);
                    obj.FC_x2(i) = x1x2(2);
                    obj.FC_DeltaChi2C(i) = interp1(mNuSq,DeltaChi2,obj.FC_x1(i),interpStyle);
                elseif mNuSq_t(i)==0 % Delta Chi2 is always zero for negative masses
                    fun =  @(x) interp1(mNuSq,CumProb,x,interpStyle)-obj.ConfLevel;
                    x1x2 = fsolve(fun,1,options);
                    obj.FC_x1(i) = NaN;
                    obj.FC_x2(i) = x1x2;
                    obj.FC_DeltaChi2C(i) = interp1(mNuSq,DeltaChi2,obj.FC_x2(i),interpStyle);
                end
                % save delta chi2 for later
                obj.FC_mNuSqFit(i,:) = linspace(min(mNuSq),max(mNuSq),nSamples);
                TestmNuSqTrue = sum(obj.FC_mNuSqFit(i,:) == mNuSq_t(i));
                if TestmNuSqTrue==0
                    tmpIndex = find(abs(obj.FC_mNuSqFit(i,:)-mNuSq_t(i))==min(abs(obj.FC_mNuSqFit(i,:)-mNuSq_t(i))),1);
                    obj.FC_mNuSqFit(i,tmpIndex) = mNuSq_t(i);
                end
                
                obj.FC_DeltaChi2(i,:) = interp1(mNuSq,DeltaChi2,obj.FC_mNuSqFit(i,:),interpStyle);
                obj.FC_PDF(i,:)       = interp1(mNuSq,Prob,obj.FC_mNuSqFit(i,:),interpStyle);
                obj.FC_CumProb(i,:)   = interp1(mNuSq,CumProb,obj.FC_mNuSqFit(i,:),interpStyle);
                obj.FC_CumProb(i,:)   = obj.FC_CumProb(i,:)./max(obj.FC_CumProb(i,:));
                if mNuSq_t(i)==0
                    obj.FC_Chi2True = interp1(mNuSq,Chi2True,obj.FC_mNuSqFit(i,:),interpStyle);
                end
            end
            
            %save
            obj.FC_mNuSqTrue = mNuSq_t';
        end
        function ComputeLokhov_Asimov(obj,varargin)
            
            p=inputParser;
            p.addParameter('mNuSq_t',0.2,@(x)isfloat(x));
            p.addParameter('nSamples',10000,@(x)isfloat(x));
            p.addParameter('nSamplesAsimov',300,@(x)isfloat(x));
            p.addParameter('ReFit','OFF',@(x)isfmember(x,{'ON','OFF'}));
            
            p.parse(varargin{:});
            mNuSq_t     = p.Results.mNuSq_t; % true neutrino masses in simulation
            nSamples    = p.Results.nSamples;
            nSamplesAsimov = p.Results.nSamplesAsimov;
            
            %init
            obj.FC_DeltaChi2  = zeros(numel(mNuSq_t),nSamples);
            obj.FC_DeltaChi2C = zeros(numel(mNuSq_t),1);
            obj.FC_mNuSqFit   = zeros(numel(mNuSq_t),nSamples);
            obj.FC_x1         = zeros(numel(mNuSq_t),1);
            obj.FC_x2         = zeros(numel(mNuSq_t),1);
            obj.FC_PDF        = zeros(numel(mNuSq_t),nSamples);
            obj.FC_CumProb    = zeros(numel(mNuSq_t),nSamples);
            
            for i=1:numel(mNuSq_t)
                % Get Delta Chi2 Curve
                [mNuSq,DeltaChi2,Chi2True,~] = obj.FC_ComputeDeltaChi2LookupTables('mNuSq_t',mNuSq_t(i),'nSamples',nSamplesAsimov);
                
                % Convert Chi2True into Probability
                Prob_tmp = exp(-Chi2True./2);
                
                %                 % start new
                %                 NormProb1 = 0.5.*simpsons(mNuSq(mNuSq<=mNuSq_t(i)),Prob_tmp(mNuSq<=mNuSq_t(i)));
                %                 NormProb2 = 0.5.*simpsons(mNuSq(mNuSq>mNuSq_t(i)),Prob_tmp(mNuSq>mNuSq_t(i)));
                %                 Prob_tmp(mNuSq<=mNuSq_t(i)) = Prob_tmp(mNuSq<=mNuSq_t(i))./NormProb1;
                %                 Prob_tmp(mNuSq>mNuSq_t(i)) = Prob_tmp(mNuSq>mNuSq_t(i))./NormProb2;
                %                 Prob = Prob_tmp;
                %                 %end new
                Prob = Prob_tmp;%./simpsons(mNuSq,Prob_tmp); % normalization
                
                CumProb_tmp = GetCDF(mNuSq,Prob);
                % find acceptance region
                [CumProb,ia,~] = unique(CumProb_tmp);
                obj.FC_x1(i) = interp1(CumProb,mNuSq(ia),(1-obj.ConfLevel)/2,'spline');
                if i~=1
                    if obj.FC_x1(i)-obj.FC_x1(i-1)<0
                        % should always be more positive! try piecewise
                        % cubic interpolation instead
                        obj.FC_x1(i) = interp1(CumProb,mNuSq(ia),(1-obj.ConfLevel)/2,'pchip');
                    end
                end
                
                if obj.FC_x1(i)>0  % --> twosided
                    obj.FC_x2(i)= interp1(CumProb,mNuSq(ia),(1+obj.ConfLevel)/2,'spline');
                elseif obj.FC_x1(i)<0 % --> one-sided
                    obj.FC_x2(i) = interp1(CumProb,mNuSq(ia),obj.ConfLevel,'spline');
                end
                
                spline = 'ON';
                switch spline
                    case 'ON'
                        interpStyle = 'spline';
                    case 'OFF'
                        interpStyle = 'lin';
                end
                % save delta chi2 for later
                obj.FC_mNuSqFit(i,:) = linspace(min(mNuSq),max(mNuSq),nSamples);
                TestmNuSqTrue = sum(obj.FC_mNuSqFit(i,:) == mNuSq_t(i));
                if TestmNuSqTrue==0
                    tmpIndex = find(abs(obj.FC_mNuSqFit(i,:)-mNuSq_t(i))==min(abs(obj.FC_mNuSqFit(i,:)-mNuSq_t(i))),1);
                    obj.FC_mNuSqFit(i,tmpIndex) = mNuSq_t(i);
                end
                
                obj.FC_DeltaChi2(i,:) = interp1(mNuSq,DeltaChi2,obj.FC_mNuSqFit(i,:),interpStyle);
                obj.FC_PDF(i,:)       = interp1(mNuSq,Prob,obj.FC_mNuSqFit(i,:),interpStyle);
                obj.FC_CumProb(i,:)   = interp1(mNuSq(ia),CumProb,obj.FC_mNuSqFit(i,:),interpStyle);
                obj.FC_CumProb(i,:)   = obj.FC_CumProb(i,:)./max(obj.FC_CumProb(i,:));
                if mNuSq_t(i)==0
                    obj.FC_Chi2True = interp1(mNuSq,Chi2True,obj.FC_mNuSqFit(i,:),interpStyle);
                end
                
            end
            %save
            obj.FC_mNuSqTrue = mNuSq_t';
        end
    end
    methods % systematics
        function SensitivitySys(obj)
            % Compute Sensitivity for stat + 1sys
            % save result in MultiLpar
            NP_prev= obj.RunAnaObj.NonPoissonScaleFactor;
            obj.RunAnaObj.NonPoissonScaleFactor = 1;
             
            loadSuccess = obj.LoadSaveSensitivity('Mode','load','SysEffect','Stat');
            if loadSuccess==0
                obj.RunAnaObj.chi2                  = 'chi2Stat';
                obj.RunAnaObj.ComputeCM_StatPNP;
                switch obj.Mode
                    case 'Asimov'
                        obj.ComputeSensitivity_Asimov('cl',0,'GetCM','OFF');  % factor applied later
                    case 'Scan'
                        obj.ComputeSensitivity_Scan('cl',obj.ConfLevel,'GetCM','OFF'); %factor has to be applied now
                end
                
                obj.LoadSaveSensitivity('Mode','save','SysEffect','Stat');
            end
            
            % stat + 1 sys sensitivity
            [~,CmArg] = GetSysErr(obj.RunAnaObj.SysBudget);
            
            for i=1:obj.nSys
                % get covariance matrix
                nTrials = obj.GetnTrials(obj.SysEffectsAll{i});
                
                loadSuccess = obj.LoadSaveSensitivity('Mode','load','SysEffect',obj.SysEffectsAll{i});
                
                if loadSuccess==0  % calculate sensitivity
                    obj.RunAnaObj.chi2     = obj.chi2sys;
                    obj.RunAnaObj.NonPoissonScaleFactor=1;
                    obj.RunAnaObj.SetNPfactor;
                     
                    % get covariance matrix
                    if ~ismember(obj.SysEffectsAll{i},{'Bkg','NP','BkgPT'})
                        obj.RunAnaObj.NonPoissonScaleFactor=1;
                        obj.RunAnaObj.ComputeCM('SysEffect',struct(obj.SysEffectsAll{i},'ON'),...
                            'nTrials',nTrials,CmArg{:},'BkgCM','OFF','BkgPtCM','OFF');
                    elseif strcmp(obj.SysEffectsAll{i},'Bkg')
                        if strcmp(obj.RunAnaObj.DataSet,'Knm1')
                            BkgMode = 'SlopeFit';
                        else
                            BkgMode = 'Gauss';
                        end
                        % exception for background shape cov mat
                        obj.RunAnaObj.NonPoissonScaleFactor=NP_prev;
                        obj.RunAnaObj.SetNPfactor;
                        obj.RunAnaObj.ComputeCM('SysEffect',struct('FSD','OFF'),...
                            'BkgCM','ON','nTrials',nTrials,'BkgMode',BkgMode);
                        obj.RunAnaObj.NonPoissonScaleFactor=1;
                        obj.RunAnaObj.SetNPfactor;
                        obj.RunAnaObj.ComputeCM('SysEffect',struct('FSD','OFF'),...
                            'BkgCM','ON','BkgPtCM','OFF','nTrials',nTrials,'BkgMode',BkgMode);
                    elseif strcmp(obj.SysEffectsAll{i},'BkgPT')
                        obj.RunAnaObj.NonPoissonScaleFactor=1;
                        obj.RunAnaObj.ComputeCM('SysEffect',struct('FSD','OFF'),...
                            'nTrials',nTrials,CmArg{:},'BkgCM','OFF','BkgPtCM','ON');
                    elseif strcmp(obj.SysEffectsAll{i},'NP')
                        % exception for background rate (Non-Poiss)
                        obj.RunAnaObj.chi2          = 'chi2Stat';
                        obj.RunAnaObj.NonPoissonScaleFactor = NP_prev;
                        obj.RunAnaObj.SetNPfactor;  
                    end
 
                    %compute sensitivity
                    switch obj.Mode
                        case 'Asimov'
                            obj.ComputeSensitivity_Asimov('cl',0,'GetCM','OFF');  % factor applied later
                        case 'Scan'
                            obj.ComputeSensitivity_Scan('cl',obj.ConfLevel,'GetCM','OFF'); %factor has to applied now
                    end
                    
                    obj.RunAnaObj.NonPoissonScaleFactor= NP_prev;
                  
                    %save
                    obj.LoadSaveSensitivity('Mode','save','SysEffect',obj.SysEffectsAll{i});
                end
            end
            
            obj.RunAnaObj.NonPoissonScaleFactor=NP_prev;
            obj.RunAnaObj.chi2     = obj.chi2sys;
               
                switch obj.Mode
                    case 'Asimov'
                        factor          = obj.ConvertCl2Std;            %convert to confidence level
                        obj.MultiLpar   = structfun(@(x) x.*factor,obj.MultiLpar,'UniformOutput',0);
                end
            %end
        end
        function SensitivitySysStack(obj)
            % Compute Sensitivity for stat + sys1 + sys2 +...
            % save result in MultiLparStack
            
            %label
            savedir = [getenv('SamakPath'),'tritium-data/sensitivity/',obj.RunAnaObj.DataSet,'/'];
            obj.CreateDir(extractBefore(savedir,obj.RunAnaObj.DataSet));
            obj.CreateDir(savedir);
            switch obj.Mode
                case 'Scan'
                    clLabel = sprintf('_%.0fcl_%sSide',obj.ConfLevel*100,obj.ScanSide);
                case 'Asimov'
                    clLabel = '';
            end
            TDlabel = extractBefore(strrep(obj.RunAnaObj.RunData.RunName,'_',''),'E0');
            if isempty(TDlabel)
                TDlabel = obj.RunAnaObj.RunData.RunName;
            end
            savefile = [savedir,sprintf('SensitivitySysStack_%s_%s_%.0feV_%s_%s_budget%.0f%s%s.mat',obj.Mode,TDlabel,...
                obj.GetRange,[obj.SysEffectsAll{:}],obj.chi2sys,obj.RunAnaObj.SysBudget,clLabel,obj.GetFitLabel)];
            
            if exist(savefile,'file') && strcmp(obj.RecomputeFlag,'OFF')
                d = importdata(savefile);
                switch obj.Mode
                    case 'Asimov'
                        factor          = obj.ConvertCl2Std;         % convert to confidence level
                        obj.MultiLparStack   = structfun(@(x) x.*factor,d.MultiLparStack,'UniformOutput',0);
                    case 'Scan'
                        obj.MultiLparStack = d.MultiLparStack;
                end
                obj.CovMatFracStack = d.CovMatFracStack;
                try %not in all files yet
                    obj.CovMatShapeStack = d.CovMatShapeStack;
                catch
                end
                fprintf('Loading sensitivites from file %s \n',savefile)
                
                if ~isfield(d,'NPcomponent')
                    NPfactor = obj.RunAnaObj.NonPoissonScaleFactor;
                    chi2tmp = obj.RunAnaObj.chi2;
                    
                    obj.RunAnaObj.NonPoissonScaleFactor=1;
                    obj.RunAnaObj.chi2 = 'chi2Stat';
                    obj.RunAnaObj.Fit;
                    obj.NPcomponent = obj.RunAnaObj.FitResult.err;
                    NPcomponent = obj.NPcomponent;
                    save(savefile,'NPcomponent','-append');
                    
                    obj.RunAnaObj.chi2 = chi2tmp;
                    obj.RunAnaObj.NonPoissonScaleFactor=NPfactor;
                end
            else
                
                % statistical sensitivity
                obj.RunAnaObj.chi2       = 'chi2Stat';
                obj.RunAnaObj.ComputeCM_StatPNP;
                %compute sensitivity
                switch obj.Mode
                    case 'Asimov'
                        obj.ComputeSensitivity_Asimov('cl',0,'GetCM','OFF');  % factor applied later
                    case 'Scan'
                        obj.ComputeSensitivity_Scan('cl',obj.ConfLevel,'GetCM','OFF'); %factor has to applied now
                end
                
                obj.MultiLparStack.Stat       = obj.Lpar;
                obj.CovMatFracStack.Stat      = obj.RunAnaObj.FitCMFrac;
                obj.CovMatShapeStack.Stat     = obj.RunAnaObj.FitCMShape;
                
                % stat + 1 + N sys sensitivity
                [~,CmArg] = GetSysErr(obj.RunAnaObj.SysBudget);
                obj.RunAnaObj.chi2     = obj.chi2sys;
                SysEffect_local  =  struct(obj.SysEffectsAll{1},'ON');
                for i=1:obj.nSys
                    % get covariance matrix
                    nTrials = obj.GetnTrials(SysEffect_local);
                    if strcmp(obj.SysEffectsAll{i},'Bkg')
                        obj.RunAnaObj.ComputeCM('SysEffect',SysEffect_local,CmArg{:},...
                            'BkgCM','ON','BkgPtCM','OFF','nTrials',nTrials);
                    elseif strcmp(obj.SysEffectsAll{i},'BkgPT')
                        obj.RunAnaObj.ComputeCM('SysEffect',SysEffect_local,CmArg{:},...
                            'BkgCM','OFF','BkgPtCM','ON','nTrials',nTrials);
                    else
                        obj.RunAnaObj.ComputeCM('SysEffect',SysEffect_local,CmArg{:},...
                            'BkgCM','OFF','BkgPtCM','OFF','nTrials',nTrials);
                    end
                    
                    %compute sensitivity
                    switch obj.Mode
                        case 'Asimov'
                            obj.ComputeSensitivity_Asimov('cl',0,'GetCM','OFF');  % factor applied later
                        case 'Scan'
                            obj.ComputeSensitivity_Scan('cl',obj.ConfLevel,'GetCM','OFF'); %factor has to applied now
                    end
                    
                    %save
                    obj.MultiLparStack.(obj.SysEffectsAll{i})   = obj.Lpar;
                    obj.CovMatFracStack.(obj.SysEffectsAll{i})  = obj.RunAnaObj.FitCMFrac;
                    obj.CovMatShapeStack.(obj.SysEffectsAll{i}) = obj.RunAnaObj.FitCMShape;
                    
                    if i<obj.nSys
                        %set next sys effect ON (in addition to old ones!)
                        SysEffect_local.(obj.SysEffectsAll{i+1}) = 'ON';
                    end
                end
                %save to file
                MultiLparStack     = obj.MultiLparStack;
                SysEffectsAll      = obj.SysEffectsAll;
                CovMatFracStack    = obj.CovMatFracStack;
                CovMatShapeStack   = obj.CovMatShapeStack;
                save(savefile,'MultiLparStack','SysEffectsAll','CovMatFracStack','CovMatShapeStack');
                
                %convert to confidence level
                switch obj.Mode
                    case 'Asimov'
                        factor          = obj.ConvertCl2Std;            %convert to confidence level
                        obj.MultiLparStack   = structfun(@(x) x.*factor,obj.MultiLparStack,'UniformOutput',0);
                end
            end
        end
    end
    
    methods % plots
        function PlotSysBreakdownBars(obj,varargin)
            % Multibar plot for systematic breakdown
            % mandatory input: Ranges
            p = inputParser;
            p.addParameter('Ranges',12, @(x)all(isfloat(x)) && all(x)>0);      % give exclDataStarts
            p.addParameter('Parameter',1,@(x)isfloat(x) & x>0); %1==mNuSq, 2==E0
            p.addParameter('SingleSys','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SysInfoBox','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','pdf'}));
            p.addParameter('DispTitle','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('xLimits','',@(x)all(isfloat(x)) || @(x)isempty(x)); % when empty, x-limit automatically set
            p.addParameter('HoldOn','OFF',@(x)ismember(x,{'ON','OFF'}));    %only used for PlotCompareSensitivityTwinData()
            p.addParameter('NP','OFF',@(x)ismember(x,{'ON','OFF'}));    % plot NP-background uncertainty
            p.addParameter('StackSys','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            
            Ranges     = p.Results.Ranges;
            Parameter  = p.Results.Parameter;
            SingleSys  = p.Results.SingleSys;
            SysInfoBox = p.Results.SysInfoBox;
            SavePlot   = p.Results.SavePlot;
            DispTitle  = p.Results.DispTitle;
            xLimits    = p.Results.xLimits;
            HoldOn     = p.Results.HoldOn;
            NP         = p.Results.NP;
            StackSys   = p.Results.StackSys;
            
            LocalFontSize = 28;
            
            if Parameter==2 && strcmp(obj.Mode,'Scan')
                fprintf('Scan for endpoint not avaibale! Set Mode to Asimov \n');
                return
            end
            
            Nranges = numel(Ranges);
            RangeLabel = cell(Nranges,1);
            
            % load sensitivities: prepare variables for barh plot
            StackBarY  = zeros(Nranges,obj.nSys+1);
            SingleBarY = zeros(Nranges,obj.nSys+1);
            StackBarX  = zeros(Nranges,1);
            SingleBarStat = zeros(Nranges,1);
            
            for i = 1:Nranges
                obj.RunAnaObj.exclDataStart = Ranges(i);
                RangeLabel{i} = sprintf('%.0f eV',obj.GetRange);
                
                if strcmp(obj.RunAnaObj.DataType,'Twin')
                    StartX = 20;
                elseif strcmp(obj.RunAnaObj.DataType,'Real')
                    StartX = 22;
                end
                
                StackBarX(i) = i*StartX;
                %stacked sys
                obj.SensitivitySysStack;
               
                PlotVarStackTmp = struct2array(structfun(@(x)x(Parameter),obj.MultiLparStack,'UniformOutput',0));
                for s=1:obj.nSys %stacked sensitivity has to increase
                    if (PlotVarStackTmp(s+1)-PlotVarStackTmp(s))<0
                        PlotVarStackTmp(s+1) = PlotVarStackTmp(s);
                    end
                end
                
                StackBarY(i,:) = Convert2Stack(PlotVarStackTmp,obj.MultiLparStack.Stat(Parameter));
                
                %single sys
                obj.SensitivitySys;
                PlotVarTmp = struct2array(structfun(@(x)x(Parameter),obj.MultiLpar,'UniformOutput',0));
                PlotVarTmp(PlotVarTmp<PlotVarStackTmp(1))=PlotVarStackTmp(1); % no sensitivity can be smaller as stat only
                SingleBarStat(i) = obj.MultiLpar.Stat(Parameter)';
                SingleBarY(i,:)  = PlotVarTmp-SingleBarStat(i);
            end
            
            SingleBarX = obj.GetRangeSingleSys('PlotRanges',StackBarX'); % Single Bar range values
            
            % In case only 1 range: add one invisible row (workaround barh)
            if Nranges==1
                StackBarY = [StackBarY;NaN.*zeros(1,obj.nSys+1)];
                StackBarX =  [StackBarX,StackBarX+1];
                SingleBarX = [SingleBarX;SingleBarX+1];
                SingleBarY = [SingleBarY;NaN.*zeros(1,obj.nSys+1)];
                SingleBarStat = [SingleBarStat,NaN]';
                RangeLabel = {RangeLabel{:};' '};
            end
            
            % MultiBar Plot (Stacked)
            if strcmp(HoldOn,'OFF')
                f55 = figure('Name','MultiBarPlot','Renderer','painters');
                set(f55, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
            end
            
            b = barh(StackBarX,StackBarY,'stacked');
            b(1).LineStyle = ':'; b(1).FaceColor = obj.PlotColor{1};
            b(1).EdgeColor = rgb('DimGray'); b(1).LineWidth = 2;
            for i=2:numel(b)
                b(i).LineStyle = 'none';
                b(i).FaceColor = obj.PlotColor{i};
            end
            
             leg_str = {'Statistical';'Theoretical corrections';'Molecular final states';...
                'Tritium activity fluctuations';'Energy loss function';'Magnetic fields';...
                sprintf('Gas density and inel. scattering cross section');'HV fluctuations';'Detector efficiency';'Background slope'};
            
            
            if strcmp(NP,'ON')
                hold on;
                bNP = barh(StackBarX,[obj.NPcomponent(Parameter),NaN],'stacked');
                bNP.LineStyle = ':'; bNP.FaceColor = obj.PlotColor{1};
                b(1).FaceColor = rgb('Silver'); bNP.LineWidth=2;
                b(1).LineStyle = 'none';    
                leg_str = {leg_str{1},'Non Poissonian Background',leg_str{2:end}};
            end
            
           
            % Plot Single Contributions
            switch SingleSys
                case 'ON'
                    if Nranges~=1
                        b(1).BarWidth = b(1).BarWidth/1.6;
                    end
                    hold on;
                    SingleSysBarWidth = b(1).BarWidth/(1.6*obj.nSys);
                    pnone = plot(NaN*[1 1],NaN*[1 1],'Color','w');
                    bsingle = cell(obj.nSys,1);
                    
                    for i=1:obj.nSys
                        bsingle{i}  = barh(SingleBarX(:,i)',[SingleBarStat';SingleBarY(:,i+1)']','stacked','BarWidth',SingleSysBarWidth);
                        btmp= bsingle{i};
                        btmp(1).FaceColor = 'w';  btmp(1).LineStyle = 'none';  btmp(1).FaceAlpha = 0;
                        btmp(2).FaceColor = obj.PlotColor{i+1}; btmp(2).LineStyle ='none';
                        btmp(2).FaceAlpha=0.5;
                        
                    end
                    b2 =cellfun(@(x) x(2),bsingle,'UniformOutput',0);
                    
                    %  leg_str = {leg_str{:},'Single Contributions','Stat + Theoretical Corrections', ...
                    %     'Stat + Final State Distribution','Stat + Activity Fluctuations','Stat + Energy Loss Function',...
                    %     'Stat + Magnetic Fields',sprintf('Stat + \\rhod\\sigma'),'Stat + Stacking',...
                    %     'Stat + FPD Efficiency','Stat + Background'};
                    if strcmp(NP,'ON')
                        leg = legend([bNP,b],leg_str{:});
                        leg.NumColumns = 2; 
                    else
                        leg = legend([b],leg_str{:});
                        leg.NumColumns = 2;
                    end
                case 'OFF'
                    leg = legend(leg_str{:});
                    leg.NumColumns = 2;
            end
            
            % legend and labels, display options
            legend boxoff
            yticklabels(RangeLabel);
            ytickangle(90)
            if obj.ConfLevel==0
                cl = 0.683;
            else
                cl = obj.ConfLevel;
            end
            if Parameter==1
                if strcmp(obj.LimitFlag,'Up')
                    xlabel(sprintf('Sensitivity upper limit {\\itm}^2_\\nu at %.0f%% C.L. (eV^{ 2})',cl*100));
                else
                    if strcmp(obj.RunAnaObj.DataType,'Real')
                         xlabel(sprintf('Uncertainty {\\itm}^2_\\nu at %.0f%% C.L. (eV^{ 2})',cl*100));
                    else
                        xlabel(sprintf('Sensitivity {\\itm}^2_\\nu at %.0f%% C.L. (eV^{ 2})',cl*100));
                    end
                end
            elseif Parameter==2
                xlabel(sprintf('Sensitivity E_{0_{fit}} at %.0f%% C.L. (eV)',cl*100));
            end
            ylabel('MTD range');
            
            if Nranges==1
                if strcmp(SingleSys,'ON')
                    ymin = min(StackBarX)-1.2;
                    ymax = max(StackBarX)+0.1;
                else
                    ymin = min(StackBarX)-0.7;
                    ymax = max(StackBarX)+0.1;
                end
                RangeLabel{1} = ['Range  ',RangeLabel{1}];
                yticklabels(RangeLabel);
                ylabel('');
            elseif Nranges==2
                ymin = min(StackBarX)-20;
                ymax = max(StackBarX)+23;
            elseif Nranges==3
                ymin = min(StackBarX)-25;
                ymax = max(StackBarX)+38;
            elseif Nranges==4
                ymin = min(StackBarX)-35;
                ymax = max(StackBarX)+49;
            elseif Nranges==5
                ymin =min(StackBarX)-45;
                ymax = max(StackBarX)+60;
            elseif Nranges==6
                ymin = min(StackBarX)-43;
                ymax = max(StackBarX)+73;
            end
            
            if strcmp(SysInfoBox,'OFF')
                ylim([ymin*1.01, ymax]);
            else
                ylim([ymin, ymax]);
            end
            
            if ~isempty(xLimits)
                xmin = xLimits(1);
                xmax = xLimits(2);
            else
                if strcmp(NP,'ON')
                    xmin =min(obj.NPcomponent(Parameter))*0.999;
                else
                    xmin =min(SingleBarStat)*0.998;
                end
                xmax = max(PlotVarStackTmp)*1.001;
            end
            xlim([xmin, xmax]);
            
            switch Parameter
                case 1
                    xticks((round(xmin,2):0.005:round(xmax,2)));
                    %xticks((round(xmin,2):0.005:round(xmax,2)));
                case 2
                    xticks((round(xmin,2):0.02:round(xmax,2)));
            end
            
            PrettyFigureFormat;
            set(gca,'FontSize',LocalFontSize);
            leg.FontSize = LocalFontSize-4;
            set(gca,'YMinorTick','off');
            h = gca;
            h.YRuler.TickLength = [0,0];
            if numel(Ranges)==1
                myYticks = yticks;
                yticks(myYticks(1));
            end
            
            labelTime = obj.RunAnaObj.ModelObj.TimeSec/(24*60*60);
            labelTD = extractBefore(obj.RunAnaObj.RunData.RunName,'E0');
            if isempty(labelTD)
                labelTD = obj.RunAnaObj.RunData.RunName;
            end
            if strcmp(DispTitle,'ON')
                title(sprintf('KATRIN %s (%.2f days) - Systematic Uncertainty Breakdown - ',...
                    labelTD,labelTime));
            end
            
            % Annotation Box with Systematic Info
            if strcmp(SysInfoBox,'ON')
                SysUncertainties = obj.GetSysInfo;
                a=annotation('textbox', [0.14 0.1 1 0.12], ...
                    'String', SysUncertainties, ...
                    'EdgeColor', 'none', ...
                    'HorizontalAlignment', 'left');
                a.FontSize=LocalFontSize-5;a.FontWeight='normal';
            end
            
            %save plot
            if strcmp(SavePlot,'ON') || strcmp(SavePlot,'pdf')
                savedir = [getenv('SamakPath'),'tritium-data/sensitivity/',obj.RunAnaObj.DataSet,'/plots/'];
                obj.CreateDir(savedir);
                ranges_str = strrep(num2str(Ranges),' ','');
                savefile = [savedir,sprintf('SensitivitySys_%s_%.0feV_%s_Par%.0f_%s%s_ranges%s.png',obj.RunAnaObj.RunData.RunName,...
                    obj.GetRange,[obj.SysEffectsAll{:}],Parameter,obj.chi2sys,obj.GetFitLabel,ranges_str)];
                if strcmp(SavePlot,'ON')
                    obj.PlotWhiteSpace;
                    print(f55,savefile,'-dpng','-r100');
                else
                    export_fig(f55,strrep(savefile,'.png','.pdf'));
                end
                 fprintf('Save plot to %s \n',savefile);
            end
        end
        function PlotSysBreakdownBars2(obj,varargin)
            % Multibar plot for systematic breakdown
            % mandatory input: Ranges
            p = inputParser;
            p.addParameter('Ranges',12, @(x)all(isfloat(x)) && all(x)>0);      % give exclDataStarts
            p.addParameter('Parameter',1,@(x)isfloat(x) & x>0); %1==mNuSq, 2==E0
            p.addParameter('SingleSys','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SysInfoBox','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','pdf'}));
            p.addParameter('DispTitle','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('xLimits','',@(x)all(isfloat(x)) || @(x)isempty(x)); % when empty, x-limit automatically set
            p.addParameter('HoldOn','OFF',@(x)ismember(x,{'ON','OFF'}));    %only used for PlotCompareSensitivityTwinData()
            p.addParameter('NP','OFF',@(x)ismember(x,{'ON','OFF'}));    % plot NP-background uncertainty
            p.addParameter('StackSys','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            
            Ranges     = p.Results.Ranges;
            Parameter  = p.Results.Parameter;
            %SingleSys  = p.Results.SingleSys;
            %SysInfoBox = p.Results.SysInfoBox;
            SavePlot   = p.Results.SavePlot;
            %DispTitle  = p.Results.DispTitle;
            xLimits    = p.Results.xLimits;
            HoldOn     = p.Results.HoldOn;
            %NP         = p.Results.NP;
            %StackSys   = p.Results.StackSys;
            
            LocalFontSize = 28;
            
            if Parameter==2 && strcmp(obj.Mode,'Scan')
                fprintf('Scan for endpoint not avaibale! Set Mode to Asimov \n');
                return
            end
            
            Nranges = numel(Ranges);
            %RangeLabel = cell(Nranges,1);
            
            % load sensitivities: prepare variables for barh plot
            SingleBarY = zeros(Nranges,obj.nSys+1);
            %StackBarX  = 20+linspace(0,10,Nranges);
            SingleBarStat = zeros(Nranges,1);
            
            for i = 1:Nranges
                obj.RunAnaObj.exclDataStart = Ranges(i);
                obj.SensitivitySys;  % get single sys
                SingleBarStat(i) = obj.MultiLpar.Stat(Parameter)';
                PlotVarTmp = struct2array(structfun(@(x)x(Parameter),obj.MultiLpar,'UniformOutput',0));
                PlotVarTmp(PlotVarTmp<PlotVarTmp(1))=PlotVarTmp(1); % no sensitivity can be smaller as stat only
                SingleBarY(i,1)      = PlotVarTmp(1);
                SingleBarY(i,2:end)  = sqrt(PlotVarTmp(2:end).^2-SingleBarStat(i).^2);
                
                if strcmp(obj.AsymErr,'ON') && strcmp(obj.RunAnaObj.fitter,'minuit')
                    PlotVarTmpNeg = struct2array(structfun(@(x)abs(x(Parameter)),obj.MultiLparNeg,'UniformOutput',0));
                    %PlotVarTmpNeg(PlotVarTmpNeg<PlotVarTmpNeg(1))=PlotVarTmpNeg(1); % no sensitivity can be smaller as stat only
                    PlotVarTmpPos = struct2array(structfun(@(x)abs(x(Parameter)),obj.MultiLparPos,'UniformOutput',0));
                    
                    %PlotVarTmpPos(PlotVarTmpPos<PlotVarTmpPos(1))=PlotVarTmpPos(1); % no sensitivity can be smaller as stat only
                    SingleBarY(i,1)      =  SingleBarStat(i);
                    %SingleBarY(i,2:end)  = 0.5.*(sqrt(PlotVarTmpNeg(2:end).^2-PlotVarTmpNeg(1).^2)+...
                    %    sqrt(PlotVarTmpPos(2:end).^2-PlotVarTmpPos(1).^2));
                    
                    MeanErr = 0.5*(PlotVarTmpPos+PlotVarTmpNeg);
                    MeanErr(MeanErr<MeanErr(1)) = MeanErr(1);% no sensitivity can be smaller as stat only
                    SingleBarY(i,2:end) = sqrt(MeanErr(2:end).^2-MeanErr(1).^2);
                end
            end
            % SingleBarX = obj.GetRangeSingleSys('PlotRanges',StackBarX'); % Single Bar range values
            
            % MultiBar Plot (Stacked)
            if strcmp(HoldOn,'OFF')
                f55 = figure('Name','MultiBarPlot','Renderer','painters');
                set(f55, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
            end
            
            leg_str = {'Statistics', obj.SysEffectLeg{:}};
            
            [SingleBarY,ia] = sort(SingleBarY,'descend');
            leg_str = leg_str(ia);
            
            bsingle = cell(obj.nSys+1,1);
            t =cell(obj.nSys+1,1);
            SingleBarX = numel(SingleBarY):-1:1;      
            for i=1:numel(SingleBarY)
                bsingle{i}  = barh(SingleBarX(i),SingleBarY(i));
                hold on;
                bsingle{i}.FaceColor = obj.PlotColor{ia(i)}; bsingle{i}.LineStyle ='none';
                bsingle{i}.FaceAlpha=1;  
                bsingle{i}.BarWidth = 0.9;
                if i==1
                    bsingle{i}.FaceColor = rgb('LightGray');
                end
                
                if SingleBarY(i)<0.001
                    tstr = sprintf('<10^{-3}') ;
%                 elseif SingleBarY(i)<0.01 && SingleBarY(i)>0.001
%                     tstr = sprintf('%.0f\\cdot10^{-3}',SingleBarY(i)*1e3) ;
%                 elseif SingleBarY(i)>0.1
%                     tstr = sprintf('%.2f',SingleBarY(i)) ;
                else
                      tstr = sprintf('%.3f',SingleBarY(i)) ;
                     %tstr = sprintf('%.0f\\cdot10^{-3}',SingleBarY(i)*1e3) ;
                    % tstr = sprintf('%.1f\\cdot10^{-2}',SingleBarY(i)*1e2) ;
                end
                t{i}= text(1.4,SingleBarX(i),tstr,...%max(SingleBarY)*1.4
                    'HorizontalAlignment','right','FontSize',LocalFontSize-4,...
                    'Interpreter','tex');
            end
            
            % axis options
            if obj.ConfLevel==0
                cl = 0.683;
            else
                cl = obj.ConfLevel;
            end
            
            if Parameter==1
                if strcmp(obj.RunAnaObj.DataType,'Real')
                    xlabel(sprintf('Uncertainty {\\itm}^2_\\nu at %.0f%% C.L. (eV^{ 2})',cl*100));
                else
                    xlabel(sprintf('Sensitivity {\\itm}^2_\\nu at %.0f%% C.L. (eV^{ 2})',cl*100));
                end
            elseif Parameter==2
                xlabel(sprintf('Sensitivity E_{0_{fit}} at %.0f%% C.L. (eV)',cl*100));
            end
            
            ylim([min(SingleBarX)-1,max(SingleBarX)+1]);
            if ~isempty(xLimits)
                xmin = xLimits(1);
                xmax = xLimits(2);
            else
                xmin = 1e-03;
                xmax = 1.5;% max(SingleBarY)*2;
            end
            xlim([xmin, xmax]);
          
            PrettyFigureFormat('FontSize',LocalFontSize);
            % no y-ticks
            set(gca,'YMinorTick','off');
            h = gca;
            h.YRuler.TickLength = [0,0];
            set(gca,'XScale','log');
            xticks([0.0,0.0001,0.001,0.01,0.1,1])
            yticks(flip(SingleBarX))
            yticklabels(flip(leg_str));
            
            %save plot
            if strcmp(SavePlot,'ON') || strcmp(SavePlot,'pdf')
                savedir = [getenv('SamakPath'),'tritium-data/sensitivity/',obj.RunAnaObj.DataSet,'/plots/'];
                obj.CreateDir(savedir);
                ranges_str = strrep(num2str(Ranges),' ','');
                savefile = [savedir,sprintf('SensitivitySys%.0f_%s_%.0feV_%s_Par%.0f_%s%s_ranges%s_SysOnly.png',obj.RunAnaObj.SysBudget,obj.RunAnaObj.RunData.RunName,...
                    obj.GetRange,[obj.SysEffectsAll{:}],Parameter,obj.chi2sys,obj.GetFitLabel,ranges_str)];
                if strcmp(obj.RunAnaObj.AnaFlag,'Ring')
                    savefile = strrep(savefile,'.png',sprintf('_Ring%s.png',obj.RunAnaObj.RingMerge));
                end
                if strcmp(SavePlot,'ON')
                    obj.PlotWhiteSpace;
                    print(f55,savefile,'-dpng','-r100');
                else
                    export_fig(f55,strrep(savefile,'.png','.pdf'));
                end
                fprintf('Save plot to %s \n',savefile);
            end
        end
        function PlotScan(obj,varargin)
            p=inputParser;
            p.addParameter('SavePlot','OFF',@(x) ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            SavePlot = p.Results.SavePlot;
            
            if isempty(obj.ScanResults)%all(obj.ScanResults.chi2min(2:end))
                fprintf(2,'No Fit Results - Run NuMassScan first! \n');
                return
            end
            
            %plot chi2 curve
            if obj.ScanResults.mnuSq_i_Fit(end)>0 %scan of positive side
                mNuPlot = 0:0.01:obj.ScanResults.mnuSq_i_Fit(end);
            elseif obj.ScanResults.mnuSq_i_Fit(end)<0 %scan of negative side
                mNuPlot = 0:-0.01:obj.ScanResults.mnuSq_i_Fit(end);
            end
            %chi2func = interp1(obj.ScanResults.mnuSq_i_Fit(~isnan(obj.ScanResults.chi2min)),obj.ScanResults.chi2min(~isnan(obj.ScanResults.chi2min)),mNuPlot,'spline');
            p = polyfit(obj.ScanResults.mnuSq_i_Fit(~isnan(obj.ScanResults.chi2min)),obj.ScanResults.chi2min(~isnan(obj.ScanResults.chi2min)),2);
            chi2f = @(x) p(1).*x.^2+p(2).*x+p(3);
            chi2func = chi2f(mNuPlot);
            
            % Plot chi2 curve
            f11 = figure('Name','NuMassScanChi2Curve','Renderer','painters');
            set(f11, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
            x = [0,obj.Lpar(1)];
            y = [0,obj.ScanResults.DeltaChi2];
            
            pInterp = plot(mNuPlot,chi2func,'-','LineWidth',3,'Color',rgb('CadetBlue'));
            hold on;
            pLine1 = plot(x,obj.ScanResults.DeltaChi2.*ones(numel(x),1),'--','LineWidth',3,'Color',rgb('IndianRed'));
            pLine2 = plot(obj.Lpar(1).*ones(numel(y),1),y,'--','LineWidth',3,'Color',rgb('IndianRed'));
            pScan = plot(obj.ScanResults.mnuSq_i_Fit,obj.ScanResults.chi2min,'o','MarkerSize',10,'MarkerFaceColor',rgb('Navy'),'MarkerEdgeColor',rgb('Navy'));
            
            xlabel('m_{\beta}^2 (eV^2)')
            ylabel(['\chi^2 (',num2str(obj.ScanResults.dof(1)),' dof)']);
            
            if obj.ScanResults.ConfLevel==0
                cl = 0.683;
            else
                cl = obj.ScanResults.ConfLevel;
            end
            
            if obj.Lpar(1)>0
                sign = 1;
                title(sprintf('upper boundary (scan precision \\Delta\\chi2=%.2g)',obj.ScanResults.ScanPrcsn));
            elseif obj.Lpar(1)<0
                sign = -1;
                title(sprintf('lower boundary (scan precision \\Delta\\chi2=%.2g)',obj.ScanResults.ScanPrcsn'));
            end
            
            if obj.ScanResults.mNuSqmin>=0
                nuSign =1;
            elseif obj.ScanResults.mNuSqmin<0
                nuSign = -1;
            end
            
            scan_leg = ['m_{\nu}^2  = ',sprintf('%.3g \\pm %.3f eV^2 (%.0f%% C.L.)\n',obj.ScanResults.mNuSqmin,obj.Lpar(1),cl*100),...
                'm_{\nu}  = ',sprintf('%.3g \\pm %.3f eV   (%.0f%% C.L.)',nuSign*sqrt(obj.ScanResults.mNuSqmin),sign*sqrt(abs(obj.Lpar(1))),cl*100)];
            DeltaChi2 = obj.ConvertCl2Std('cl',obj.ScanResults.ConfLevel)^2;     % defines confidence level
            line_leg =sprintf('\\chi^2_{min}+%.2f',DeltaChi2);
            leg = legend([pInterp, pLine1],scan_leg,line_leg,'Location','northwest'); legend boxoff;
            
            if ~strcmp(obj.RunAnaObj.chi2,'chi2Stat')
                leg.Title.String = ['stat + ',obj.GetSysEffect];
            else
                leg.Title.String = 'stat';
            end
            PrettyFigureFormat;
            
            if obj.Lpar(1)>0      %positive side
                xlim([0 obj.Lpar(1)*1.5]);
            elseif obj.Lpar(1)<0  %negative side
                xlim([obj.Lpar(1)*1.5 0]);
                leg.Location = 'northeast';
            end
            
            ylim([min(chi2func),interp1(obj.ScanResults.mnuSq_i_Fit(~isnan(obj.ScanResults.chi2min)),obj.ScanResults.chi2min(~isnan(obj.ScanResults.chi2min)),obj.Lpar(1)*1.5,'spline')]);
            
            if strcmp(SavePlot,'ON')
                savedir = './plots/';
                obj.CreateDir(savedir);
                
                savename = sprintf('mNuScan_%s_%s_%s_%.0feV%s_%s.png',obj.RunAnaObj.RunData.RunName,...
                    obj.RunAnaObj.chi2,obj.GetSysEffect,obj.GetRange,obj.GetFitLabel,obj.ScanResults.ScanSide);
                print(f11,[savedir,savename],'-dpng','-r450');
                fprintf('save plot to %s \n',[savedir,savename]);
            end
        end
        function PlotCompareSensitivityTwinData(obj,varargin)
            p=inputParser;
            p.addParameter('SavePlot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Parameter',1,@(x)isfloat(x) & x>0);
            p.parse(varargin{:});
            SavePlot = p.Results.SavePlot;
            Parameter = p.Results.Parameter;
            %             if strcmp(obj.RunAnaObj.DataSet,'Knm1')
            %                 fprintf(2,'You are not allowed to see this! Wait for the unblinding ;) \n')
            %                 return
            %             end
            
            % Get breakdown for twins
            obj.RunAnaObj.DataType = 'Twin';
            obj.GetData;
            obj.PlotSysBreakdownBars('HoldOn','OFF','Ranges',obj.RunAnaObj.exclDataStart,'SavePlot','OFF','Parameter',Parameter,...
                'DispTitle','OFF','SysInfoBox','OFF');
            
            % Get breakdown for real data
            obj.RunAnaObj.DataType = 'Real';
            obj.GetData;
            obj.PlotSysBreakdownBars('HoldOn','ON','Ranges',obj.RunAnaObj.exclDataStart,'SavePlot','OFF','Parameter',Parameter,...
                'DispTitle','OFF','SysInfoBox','OFF');
            yticks([20,22])
            yticklabels({'Twin','Data'});
            ytickangle(0);
            % ylabel(sprintf('%.0f eV range',obj.GetRange));
            ylim([18.2 24.3]);
            switch Parameter
                case 1
                    xlim([0.76,1]);
                case 2
                    xlim([0.12,0.26]);
                    ylim([18.6 24.3]);
            end
            
            set(gca,'FontSize',22);
            %save plot
            if strcmp(SavePlot,'ON')
                savedir = [getenv('SamakPath'),'tritium-data/sensitivity/',obj.RunAnaObj.DataSet,'/plots/'];
                obj.CreateDir(savedir);
                savefile = [savedir,sprintf('SensitivitySysComparison_%s_%.0feV_%s_Par%.0f_%s%s.png',obj.RunAnaObj.RunData.RunName,...
                    obj.GetRange,[obj.SysEffectsAll{:}],1,obj.chi2sys,obj.GetFitLabel)];
                print(gcf,savefile,'-dpng','-r450');
            end
            
        end
        function PlotsUpperLimit(obj)
            if isempty(obj.NFitResults)
                fprintf('Neyman Confidence Interval: Fit Results empty \n')
                return
            end
            
            % 1. plots: plot all neutrino mass fit results
            mNuSq_t = obj.NFitResults.mNuSq_t; % "true" neutrino masses in simulation
            
            Colors = jet(3*numel(mNuSq_t));
            h = cell(numel(mNuSq_t),1);
            f1 = figure('Name','MultiBarPlot','Renderer','painters');
            set(f1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
            
            for i=1:numel(mNuSq_t)
                h{i} = histogram(obj.NFitResults.par(i,1,:));
                hold on;
                h{i}.FaceAlpha = 1-(i-1)/numel(mNuSq_t);
                h{i}.FaceColor = Colors(3*i,:);
                h{i}.EdgeColor = Colors(3*i,:);
                h{i}.LineWidth = 1.5;
            end
            PrettyFigureFormat;
            hold off
        end
        function PlotFCBelt(obj,varargin)
            p = inputParser;
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','png'}));
            p.addParameter('HoldOn','OFF',@(x)ismember(x,{'ON','OFF'})); % plot 90 and 95% at the same time
            p.addParameter('Lokov','OFF',@(x)ismember(x,{'ON','OFF'})); % plot Lokov belt instead
            p.addParameter('Sensitivity','OFF',@(x)ismember(x,{'ON','OFF'})); % show sensitivity instead of upper limit
            p.addParameter('Style','PRL',@(x)ismember(x,{'PRL','Pretty'}));
            p.addParameter('mNuSq_bf','',@(x)isfloat(x));
            p.addParameter('XLim',[-3,3],@(x)isfloat(x));
            p.parse(varargin{:});
            SavePlot = p.Results.SavePlot;
            HoldOn   = p.Results.HoldOn;
            Lokov    = p.Results.Lokov;
            Sensitivity = p.Results.Sensitivity;
            Style       = p.Results.Style;
            XLim = p.Results.XLim;
            mNuSq_bf = p.Results.mNuSq_bf;
            
            if strcmp(Style,'PRL')
                LocalFontSize = 30;
            else
                LocalFontSize = 28;
            end
            
            if isempty(obj.FC_x1)
                fprintf('No FC belt computed \n');
                return;
            end
            
            if strcmp(HoldOn,'OFF')
                f111 = figure('Renderer','painters');
                set(f111, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]); 
                LineArg = {'Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('SteelBlue'),'LineWidth',4,'MarkerSize',0.01};
                AreaArg = {'FaceColor',rgb('SteelBlue'),'FaceAlpha',0.1,'LineStyle','none'};
            else
                LineArg = {'Color',rgb('ForestGreen'),'MarkerFaceColor',rgb('PaleGreen'),'LineWidth',4,'MarkerSize',0.01};
                AreaArg = {'FaceColor',rgb('PaleGreen'),'FaceAlpha',0.1,'LineStyle','none'};
            end
            
            x1Index = find(~isnan(obj.FC_x1),1);
            plotX1 = linspace(min(obj.FC_x1(x1Index:end)),max(obj.FC_x1),1000)';
            plotY1 = interp1(obj.FC_x1(x1Index:end),obj.FC_mNuSqTrue(x1Index:end),plotX1,'spline');
            plotX2 = linspace(min(obj.FC_x2),max(obj.FC_x2),1000)';
            plotY2 = interp1(obj.FC_x2,obj.FC_mNuSqTrue,plotX2,'spline');
            
              switch Lokov
                case 'ON'
                 SensitivityLimit = interp1(obj.FC_x1(x1Index:end),obj.FC_mNuSqTrue(x1Index:end),0,'spline');  
                 xIndex = (plotX1<=0);
                 plotY1(xIndex) = SensitivityLimit;
                 mNuSqTrue  = obj.FC_mNuSqTrue;
                 mNuSqTrue(obj.FC_mNuSqTrue<SensitivityLimit) = SensitivityLimit;      
                 
                 % draw sharp line
                 yIndex1 = find(plotY2>=SensitivityLimit,1);
                 yIndex2 = numel(plotY2)-find(flipud(plotY2)<=SensitivityLimit,1);
                 plotY2(yIndex1:yIndex2) = SensitivityLimit;
                 
                 %add -3 sigma point
                 if min(plotX1)>-3
                    plotX1 =[-3;plotX1];
                    plotY1 = [SensitivityLimit;plotY1];
                 end
                  case 'OFF'
                      mNuSqTrue  = obj.FC_mNuSqTrue;
                      
              end
            
              a1 = area([plotX1;linspace(plotX1(end),3,1)],[plotY1;2],AreaArg{:});
              hold on
              a2 =area(plotX2,plotY2,0,'FaceColor',rgb('White'),'FaceAlpha',1,'LineStyle','none');
              
              if strcmp(Lokov,'ON')
                 
                  % plot continuous edge x1
                  IndexEdge = find(~xIndex,1);
                  % plot first part x2
                  EdgeIndex2 = max(find(obj.FC_x1<0));
                  EdgeX2 = obj.FC_x2(EdgeIndex2);
                  IndexEdgeX21= find(plotX2>=EdgeX2,1);
                  
                  % plot second part x2
                  EdgeX2 = obj.FC_x2(EdgeIndex2+1);
                  IndexEdgeX22 = find(plotX2>=EdgeX2,1);
                  % connection: lower part - horizontal line
                  x21_tmp = linspace(plotX2(IndexEdgeX21-10),plotX2(IndexEdgeX21+50),1000);
                  ploty2_edge = interp1(plotX2(IndexEdgeX21-10:IndexEdgeX21-2),plotY2(IndexEdgeX21-10:IndexEdgeX21-2),x21_tmp,'lin','extrap');   
                     
                  % connection to upper part
                   if strcmp(obj.RunAnaObj.DataSet,'Knm1')
                        x22_tmp = linspace(plotX2(IndexEdgeX22-50),plotX2(IndexEdgeX22+5),1000);
                 
                      ploty22_edge = interp1(plotX2(IndexEdgeX22:IndexEdgeX22+5),plotY2(IndexEdgeX22:IndexEdgeX22+5),x22_tmp,'lin','extrap');
                      xtmp = linspace(max(x21_tmp(ploty2_edge<=SensitivityLimit-0.01)),min(x22_tmp(ploty22_edge>=SensitivityLimit+0.008)),100);
                   else
                     x22_tmp = linspace(plotX2(IndexEdgeX22-52),plotX2(IndexEdgeX22+5),1000);
                     ploty22_edge = interp1(plotX2(IndexEdgeX22:IndexEdgeX22+4),plotY2(IndexEdgeX22:IndexEdgeX22+4),x22_tmp,'lin','extrap');
                     xtmp = linspace(max(x21_tmp(ploty2_edge<=SensitivityLimit)),min(x22_tmp(ploty22_edge>=SensitivityLimit)),100);
                   end
                  
                  % correct area
                  plotX2_new = [plotX2(1:IndexEdgeX21-5)',xtmp,plotX2(IndexEdgeX22:end)'];
                  ploty2_new = [plotY2(1:IndexEdgeX21-5)',SensitivityLimit.*ones(100,1)',plotY2(IndexEdgeX22:end)'];
                  a2.delete;
                  a2 =area(plotX2_new,ploty2_new,0,'FaceColor',rgb('White'),'FaceAlpha',1,'LineStyle','none');

                  %%% plot correct lines
                % plot first part x1
                  p1 =plot(plotX1(xIndex),plotY1(xIndex),'-',LineArg{:});%smooth(plotY1(xIndex),100)
                  % plot second part x1
                  plot(plotX1(~xIndex),plotY1(~xIndex),'-',LineArg{:});      
                  % plot continuous edge x1 
                  plot(plotX1(IndexEdge-1:IndexEdge+1),smooth(plotY1(IndexEdge-1:IndexEdge+1),100),'-',LineArg{:});
                  % plot first part x2
                  p2 =plot(plotX2(1:IndexEdgeX21-5),smooth(plotY2(1:IndexEdgeX21-5),100),'-',LineArg{:});
                  % plot second part x2
                  p2 =plot(plotX2(IndexEdgeX22:end),smooth(plotY2(IndexEdgeX22:end),100),'-',LineArg{:});
                 
                  % connection: lower part - horizontal line  
                 plot(x21_tmp(ploty2_edge<=SensitivityLimit+1e-04),ploty2_edge(ploty2_edge<=SensitivityLimit+1e-04),'-',LineArg{:});
                  % connection to upper part  
                plot(x22_tmp(ploty22_edge>=SensitivityLimit),ploty22_edge(ploty22_edge>=SensitivityLimit),'-',LineArg{:});
                 plot(xtmp,SensitivityLimit.*ones(100,1),'-',LineArg{:}) %horizontal line       
                % p4= plot(plotX2(IndexEdgeX21-1:IndexEdgeX21+1),smooth(plotY2(IndexEdgeX21-1:IndexEdgeX21+1),100),'-',LineArg{:})
                 
%                   % smooth edges
%                   xedge1_tmp = x21_tmp(ploty2_edge<=SensitivityLimit);
%                   yedge1_tmp = ploty2_edge(ploty2_edge<=SensitivityLimit);
%
%                   xSmoothEdge1 = [xedge1_tmp(end-1:end),xtmp(1:3)];
%                   YSmoothEdge1 = [yedge1_tmp(end-1:end),SensitivityLimit.*ones(1,3)];
%                   plot(xSmoothEdge1,smooth(YSmoothEdge1,100));
                  legStr = 'LT';
              else
                  legStr = 'FC';
                  p1 =plot(plotX1,plotY1,'-',LineArg{:});%smooth(plotY1,100)
                  p2 =plot(plotX2,plotY2,'-',LineArg{:});%smooth(plotY2,100)
              end
              
            if strcmp(SavePlot,'ON')
                xlabel(sprintf('Measured {\\itm}^2_\\nu (eV^{ 2})'));
                ylabel(sprintf('True {\\itm}^2_\\nu (eV^{ 2})'));
            else
                xlabel(sprintf('Measured {\\itm}^2_\\nu (eV^2)'));
                ylabel(sprintf('True {\\itm}^2_\\nu (eV^2)'));
            end
           
            if isempty(XLim)
                xmin = min(obj.FC_x1(obj.FC_x1>-10));
                xlim([xmin, max(obj.FC_x2)]);
            else
                xlim([min(XLim) max(XLim)]);
            end
            
            ylim([min(obj.FC_mNuSqTrue),max(obj.FC_mNuSqTrue)]);
            if strcmp(Sensitivity,'ON')
                mNuMeasured = 0;
                savestr = 'sensitivity_';
            else
                if isempty(mNuSq_bf)
                    if strcmp(obj.RunAnaObj.DataSet,'Knm1')
                        mNuMeasured = -0.98;
                    elseif strcmp(obj.RunAnaObj.DataSet,'Knm2')
                        mNuMeasured = 0.19;
                    else
                        mNuMeasured = 0;
                    end
                     savestr = '';
                else
                    mNuMeasured = mNuSq_bf;
                    savestr = sprintf('_bf%.2feV2',mNuSq_bf);
                end
               
            end
            x1 = obj.FC_x1(~isnan(obj.FC_x1));
            yticks(0:0.2:max(obj.FC_mNuSqTrue))
            mNuLimit = interp1(x1,obj.FC_mNuSqTrue(~isnan(obj.FC_x1)),mNuMeasured,'spline');
            if strcmp(Lokov,'ON') && mNuMeasured<0
                    mNuLimit = SensitivityLimit;
            end
                
            plimit = plot(mNuMeasured*[1,1],[0,mNuLimit],':','LineWidth',p1.LineWidth,'color',rgb('Orange'));
            plot([min(plotX1),mNuMeasured],[mNuLimit,mNuLimit],':','LineWidth',p1.LineWidth,'color',rgb('Orange'));
            if strcmp(obj.RunAnaObj.chi2,'chi2Stat')
                chi2str = '(stat. only)';
            else
                chi2str = '(stat. and sys.)';
            end
            
            if strcmp(Sensitivity,'ON')
                leg = legend([p2,plimit],[sprintf('%s confidence belt at %.4g%% C.L. ',legStr,obj.ConfLevel*100),chi2str],...
                    sprintf(' Sensitivity: {\\itm}^2_\\nu \\leq %.2f eV^2 \n              \\rightarrow {\\itm}_\\nu  \\leq %.2f eV',mNuLimit,sqrt(mNuLimit)),'Location','northwest');
            else
                leg = legend([p2,plimit],[sprintf('%s confidence belt at %.4g%% C.L. ',legStr,obj.ConfLevel*100),chi2str],...
                    sprintf('     {\\itm}^2_\\nu \\leq %.2f eV^{ 2} \n\\rightarrow {\\itm}_\\nu  \\leq %.2f eV',mNuLimit,sqrt(mNuLimit)),'Location','northwest');
            end
           %legend boxoff;
           PrettyLegendFormat(leg);
            %leg.EdgeColor = rgb('Silver');
            % axis style etc.
            % PrettyFigureFormat('FontSize',LocalFontSize);
            if strcmp(Style,'PRL')
                PRLFormat;
                set(gca,'FontSize',LocalFontSize)
                set(gca,'TickDir','out');
                % remove top and right ticks
                a = gca;
                set(a,'box','off','color','none')% set box property to off and remove background color
                b = axes('Position',a.Position,...
                    'box','on','xtick',[],'ytick',[],'LineWidth',1.5);% create new, empty axes with box but without ticks
                axes(a)% set original axes as active
                linkaxes([a b]) % link axes in case of zooming
            else
                PrettyFigureFormat('FontSize',LocalFontSize);
            end
            
            %%
            if ~strcmp(SavePlot,'OFF')
                savedir = [getenv('SamakPath'),'/tritium-data/FC/plots/'];
                obj.CreateDir(savedir);
                savefile = [savedir,sprintf('%s_FCbelt_%s%s_%.0fCL.pdf',obj.RunAnaObj.DataSet,savestr,obj.RunAnaObj.chi2,obj.ConfLevel*100)];
                if strcmp(Lokov,'ON')
                    leg.FontSize = get(gca,'FontSize')-6;
                    savefile = [savedir,sprintf('%s_Lokovbelt_%s%s_%.0fCL.pdf',obj.RunAnaObj.DataSet,savestr,obj.RunAnaObj.chi2,obj.ConfLevel*100)];
                end
                
                if strcmp(SavePlot,'png')
                    %obj.PlotWhiteSpace;
                    savefile = strrep(savefile,'.pdf','.png');
                    print(f111,savefile,'-dpng','-r100');
                else
                    export_fig(f111,savefile);
                end
                fprintf('Save plot to %s \n',savefile);
            end
            
        end
        function PlotFC_DeltaChi2(obj,varargin)
            p=inputParser;
            p.addParameter('mNuSq_t',0.6,@(x)isfloat(x));
            p.addParameter('mNuSq_bf',-0.98,@(x)isfloat(x)); 
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','pdf','png'}));
            p.addParameter('Mode','Chi2',@(x)ismember(x,{'Chi2','PDF','CompareBf'})); % plot PDF instead
            p.parse(varargin{:});
            mNuSq_t  = p.Results.mNuSq_t;
            SavePlot = p.Results.SavePlot;
            Mode      = p.Results.Mode;
            mNuSq_bf = p.Results.mNuSq_bf;
            
            LabelFontSize = 24; %
            
            f111 = figure('Renderer','painters');
            set(f111, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
            PlotArg = {'Color',rgb('DodgerBlue'),'LineWidth',4};
            
            Index = find(mNuSq_t==obj.FC_mNuSqTrue);
            if isempty(Index)
                fprintf('Delta Chi2 not available for this neutrino mass squared \n');
                return
            end
            
            x = obj.FC_mNuSqFit(Index,:);
           
            if strcmp(Mode,'Chi2')
                % plot chi2-profile
                y = obj.FC_DeltaChi2(Index,:);
                ystr = sprintf('\\Delta\\chi^2');
                xstr = sprintf('Measured {\\itm}^2_\\nu (eV^{2})');
                AuxLineMax1 = obj.FC_DeltaChi2C(Index);
                AuxLineMax2 = obj.FC_DeltaChi2C(Index);
            else
                % convert(original chi2, not delta chi2) into likelihood profile
                y = exp(-obj.FC_Chi2True./2);%exp(-obj.FC_DeltaChi2(Index,:)./2);%exp(-obj.FC_Chi2True./2);
%                 % tmp start
%                 Prob_tmp = y;
%                 mNuSq = obj.FC_mNuSqFit(Index,:);
%                 NormProb1 = 0.5.*simpsons(mNuSq(mNuSq<=mNuSq_t),Prob_tmp(mNuSq<=mNuSq_t));
%                 NormProb2 = 0.5.*simpsons(mNuSq(mNuSq>mNuSq_t),Prob_tmp(mNuSq>mNuSq_t));
%                 Prob_tmp(mNuSq<=mNuSq_t) = Prob_tmp(mNuSq<=mNuSq_t)./NormProb1;
%                 Prob_tmp(mNuSq>mNuSq_t) = Prob_tmp(mNuSq>mNuSq_t)./NormProb2;
%                 y = Prob_tmp;
%                 %end
              %  y = y./simpsons(x,y); 
                
                ystr = sprintf('Probability density (eV^{-2})');
                xstr = sprintf('Measured {\\itm}^2_\\nu (eV^{2})');  
                AuxLineMax1 =  interp1(x,y,obj.FC_x1(Index),'spline');
                AuxLineMax2 =  interp1(x,y,obj.FC_x2(Index),'spline');
            end
            
            % fill confidence area
            if ~strcmp(Mode,'CompareBf')
                if isnan(obj.FC_x1(Index))
                    FillLog = x<=obj.FC_x2(Index);
                    xmin = -(obj.FC_x2(Index)+(obj.FC_x2(Index)-mNuSq_t));
                else
                    FillLog =x >=obj.FC_x1(Index) & x<=obj.FC_x2(Index);
                    xmin = obj.FC_x1(Index)-(mNuSq_t-obj.FC_x1(Index));
                end
                p2 = area(x(FillLog),y(FillLog),0,...
                    'LineStyle','none','ShowBaseLine','OFF','FaceColor',rgb('Silver'));
                xlim([xmin,obj.FC_x2(Index)+(obj.FC_x2(Index)-mNuSq_t)]);
                hold on;
                if strcmp(Mode,'Chi2')
                    p2.BaseValue = obj.FC_DeltaChi2C(Index);
                    ytmp = ylim;
                    if mNuSq_t==0
                    ylim([-0.05,ytmp(2)]);
                    end
                end
            end
            
            % draw chi2 or likelihood profile
            p1 = plot(x,y,'-',PlotArg{:});
            hold on;
            
            if strcmp(Mode,'Chi2') || strcmp(Mode,'PDF')
                p3 =plot(obj.FC_x1(Index)*ones(10,1),linspace(0,AuxLineMax1,10),...
                    ':','LineWidth',3,'Color',rgb('Black'));
                plot(obj.FC_x2(Index)*ones(10,1),linspace(0,AuxLineMax2,10),...
                    ':','LineWidth',3,'Color',rgb('Black'));
            elseif strcmp(Mode,'CompareBf')
                p1.delete;
                a1 =area(x,y,'FaceColor',rgb('DodgerBlue'),'FaceAlpha',1,'EdgeColor',rgb('DodgerBlue'),...
                    'LineWidth',3);
                if mNuSq_bf<0
                    a2 =area(x(x<=mNuSq_bf),y(x<=mNuSq_bf),'FaceColor',rgb('DarkBlue'),'FaceAlpha',1,...
                        'EdgeColor',rgb('DarkBlue'),'LineWidth',3);
                else
                    a2 =area(x(x>=mNuSq_bf),y(x>=mNuSq_bf),'FaceColor',rgb('DarkBlue'),'FaceAlpha',1,...
                        'EdgeColor',rgb('DarkBlue'),'LineWidth',3);
                end
            end
            
            if ~strcmp(Mode,'Chi2')
                if strcmp(obj.RunAnaObj.DataSet,'Knm1')
                    xlim([mNuSq_t-5,mNuSq_t+5]);
                    xticks((-4:4))
                else
                    xlim([mNuSq_t-1.5,mNuSq_t+1.5]);
                    xticks((-1:0.5:1))
                end
            end
            
            xlabel(xstr);
            ylabel(ystr);
            
            PrettyFigureFormat;
            set(gca,'FontSize',LabelFontSize-2);
            set(get(gca,'XLabel'),'FontSize',LabelFontSize);
            set(get(gca,'YLabel'),'FontSize',LabelFontSize);
            set(gca,'FontSize',LabelFontSize-2);

            if mNuSq_t~=0
                plotx1 = sprintf('x_1 = %.2f eV^2  ,  ',obj.FC_x1(Index));
            else
                plotx1 = '';
            end
     
            if strcmp(Mode,'Chi2')
                leg = legend([p1,p2,p3],sprintf('\\chi^2 profile for true {\\itm}^2_\\nu = %.3g eV^2',mNuSq_t),...
                    sprintf('\\Delta\\chi^2_{crit.} = %.2f (%.0f%% C.L.)',obj.FC_DeltaChi2C(Index),obj.ConfLevel*100),...
                    sprintf('%sx_2 = %.2f eV^2',plotx1,obj.FC_x2(Index)));
                leg.Location = 'northwest';
            PrettyLegendFormat(leg,'alpha',0.7);
                savestr = '';
                savedir = [getenv('SamakPath'),'tritium-data/FC/plots/Chi2Profile/'];
            elseif strcmp(Mode,'PDF')
                leg = legend([p1,p2,p3],sprintf('Profile Likelihood for true {\\itm}^2_\\nu = %.3g eV^2',mNuSq_t),...
                    sprintf('Confidence region (%.0f%% C.L.)',obj.ConfLevel*100),...
                    sprintf('%sx_2 = %.2f eV^2',plotx1,obj.FC_x2(Index)));
                savedir = [getenv('SamakPath'),'tritium-data/FC/plots/Likelihood/'];
                savestr = 'PDF_1sigma';
                leg.Location = 'northwest';
                PrettyLegendFormat(leg,'alpha',0.7);
            elseif strcmp(Mode,'CompareBf')
                CumSum = GetCDF(x,y);
                ylim([0,max(y)*1.05])
                ax = gca;
                mNuSqFrac =100.*interp1(x,CumSum,mNuSq_bf,'spline');
                if mNuSq_bf>0
                    mNuSqFrac =100-mNuSqFrac;
                     t = text(mNuSq_bf+0.05,0.93,sprintf('{\\itP}({\\itm}_{measured}^2 \\geq {\\itm}_{bf}^2) = %.1f %s',mNuSqFrac,'%'),...
                        'FontSize',ax.FontSize,'FontWeight',ax.FontWeight,'Color',rgb('DarkBlue'));
                    
                    % make nice arrow
                    xStart = interp1(ax.XAxis.Limits,[0,1],x(find(x>=mNuSq_bf+0.2,1)));
                    yStart = interp1(ax.YAxis.Limits,[0,1],y(find(x>=mNuSq_bf+0.2,1)));
                    xEnd = interp1(ax.XAxis.Limits,[0,1],x(find(x>=mNuSq_bf+0.4,1)));  
                    arrow = annotation('textarrow',[xEnd,xStart],[0.8,yStart],'Color',rgb('DarkBlue'),'LineWidth',2);
                else
                    xStart = interp1(ax.XAxis.Limits,[0,1],x(find(x>=mNuSq_bf-0.1,1)));
                    yStart = interp1(ax.YAxis.Limits,[0,1],y(find(x>=mNuSq_bf-0.1,1)));
                    xEnd = interp1(ax.XAxis.Limits,[0,1],x(find(x>=mNuSq_bf,1)));
                    
                    t = text(min(xlim)+0.3,xEnd+0.02,sprintf('{\\itP}({\\itm}_{measured}^2 \\leq {\\itm}_{bf}^2) = %.1f %s',mNuSqFrac,'%'),...
                        'FontSize',ax.FontSize,'FontWeight',ax.FontWeight,'Color',rgb('DarkBlue'));
                    
                    % make nice arrow
                    annotation('textarrow',[xStart,xEnd],[0.8,yStart],'Color',rgb('DarkBlue'),'LineWidth',2);
                end
                
                savedir = [getenv('SamakPath'),'tritium-data/FC/plots/Likelihood/'];
                savestr = sprintf('PDF_CompareBf%.3geV2',mNuSq_bf);
          
            end
            
            
            %mylim = ylim;
            %  text(-3.85,max(mylim)*0.94,'d)','FontSize',get(gca,'FontSize')+6,'FontName',get(gca,'FontName'));
           
            obj.CreateDir(savedir);
            savefile = [savedir,sprintf('%s_DeltaChi2_%.3geV2_%s%s.pdf',obj.RunAnaObj.DataSet,mNuSq_t,obj.RunAnaObj.chi2,savestr)];
            
            if ismember(SavePlot,{'ON','pdf'})
                fprintf('Save plot to %s \n',savefile);
                export_fig(f111,savefile);
            elseif strcmp(SavePlot,'png')
             %   obj.PlotWhiteSpace;
                savefile = strrep(savefile,'.pdf','.png');
                print(f111,savefile,'-dpng','-r300');
                fprintf('Save plot to %s \n',savefile);
            end
        end
        function PlotFC_PDF(obj,varargin)
            p=inputParser;
            p.addParameter('mNuSq_t',0.6,@(x)isfloat(x));
            p.addParameter('nSamples',50000,@(x)isfloat(x));
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('KNM1Central','OFF',@(x)ismember(x,{'ON','OFF'})); % show central value of KNM1 instead of [x1,x2]
            p.parse(varargin{:});
            mNuSq_t = p.Results.mNuSq_t;
            nSamples = p.Results.nSamples;
            SavePlot = p.Results.SavePlot;
            KNM1Central = p.Results.KNM1Central;
            
            if strcmp(KNM1Central,'ON')
                mNuSq_t = 0;
                fprintf('Switching to mNuSqTrue = 0 eV^2 \n');
            end
            
            Index = find(mNuSq_t==obj.FC_mNuSqTrue);
            if isempty(Index)
                fprintf('PDF not available for this neutrino mass squared \n');
                return
            end
            
            f111 = figure('Renderer','painters');
            set(f111, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.7]);
            
            % Draw nSamples from this probability distribution
%               Prob_tmp = exp(-0.5.*obj.FC_Chi2True);
%               Prob = Prob_tmp./simpsons(obj.FC_mNuSqFit(1,:),Prob_tmp);
%                CumProb_tmp = cumsum(Prob);
%                CumProb = max(CumProb_tmp);
%                
            x_Samples = linspace(0,1,nSamples);
            [CumProb,ia,~] = unique(obj.FC_CumProb(Index,:));
            Prob_Samples = interp1(CumProb,obj.FC_mNuSqFit(Index,ia),x_Samples,'spline');
            Prob_Samples(isnan(Prob_Samples)) = [];
            
            % select confidence interval
            if strcmp(KNM1Central,'OFF')
                x2 = obj.FC_x2(Index);
                x2LineStyle = '--';
                x2color = rgb('Silver');
                plotx2 =  sprintf('x_2 = %.2f eV^2',obj.FC_x2(Index));
                if mNuSq_t~=0
                    ProbConfInt = Prob_Samples(Prob_Samples>=obj.FC_x1(Index) & Prob_Samples<=obj.FC_x2(Index));
                    plotx1 = sprintf('x_1 = %.2f eV^2  ,  ',obj.FC_x1(Index));
                elseif mNuSq_t==0
                    ProbConfInt = Prob_Samples(Prob_Samples<=obj.FC_x2(Index));
                    plotx1 = '';
                end
               
            elseif strcmp(KNM1Central,'ON')
                x2 = -0.98;
                x2color = rgb('Orange');
                x2LineStyle = '-';
                ProbConfInt = Prob_Samples(Prob_Samples<=x2);
                
                plotx2 = sprintf(' KATRIN best fit: {\\itm}^2_\\nu     =  - %.2g _{-  %.2f}^{+ %.2f} eV^2',...
                    0.98,1.06,0.89);
                %  plotx2 = sprintf('KATRIN best fit 
            end
            
            h1 =histogram(Prob_Samples,'FaceColor',rgb('DodgerBlue'),'FaceAlpha',0.9);
            hold on;
            h2 = histogram(ProbConfInt,'FaceColor',rgb('DarkBlue'),'FaceAlpha',0.9);
            h2.BinWidth = h1.BinWidth;
            
            if strcmp(KNM1Central,'OFF')
                histline =max(h1.Values)*0.8;
            elseif strcmp(KNM1Central,'ON')
                histline =max(h1.Values)*1.5;
            end
            
            plot(obj.FC_x1(Index)*ones(10,1),linspace(0,histline,10),...
                '--','LineWidth',4,'Color',rgb('Silver'));
            p3 =plot(x2.*ones(10,1),linspace(0,histline,10),...
                'LineStyle',x2LineStyle,'LineWidth',5,'Color',x2color);
            PrettyFigureFormat;
            ylim([0,max(h1.Values)*(1+1e-02)]);
            xlim([mNuSq_t-4.9,mNuSq_t+4.9])
            xlabel(sprintf('measured {\\itm}^2_\\nu (eV^2)'));
            ax = gca;
            ax.XLabel.FontSize = 23;
            if strcmp(KNM1Central,'ON')
                leg = legend([p3,h1,h2],plotx2,' Pseudo experiments',...
                    sprintf(' P( {\\itm}^2_\\nu \\leq -0.98 eV^2 | {\\itm}^2_\\nu = 0 eV^2 ) = %.1f %%',100*numel(ProbConfInt)/numel(Prob_Samples)));
                leg.Location = 'northeast';
            else
                leg = legend([h1,h2,p3],sprintf('pseudo experiments: true {\\itm}^2_\\nu = %.3g eV^2',mNuSq_t),...
                    sprintf('P( x_1 \\leq  m_{measured} \\leq x_2 | m_{true}) = %.1f %%',100*numel(ProbConfInt)/numel(Prob_Samples)),...
                    sprintf('%s%s',plotx1,plotx2));
                leg.Location = 'northwest';
            end
            
            legend boxoff;
           % leg.Location = 'northwest';
            if strcmp(SavePlot,'ON')
                obj.PlotWhiteSpace;
                savedir = [getenv('SamakPath'),'/tritium-data/FC/plots/'];
                obj.CreateDir(savedir);
                savefile = [savedir,sprintf('%s_PDFcovarage_%.3geV2_%s.pdf',obj.RunAnaObj.DataSet,mNuSq_t,obj.RunAnaObj.chi2)];
                if strcmp(KNM1Central,'ON')
                    savefile = strrep(savefile,'PDFcovarage','PDFprcntile');
                end
                %  print(f111,savefile,'-dpng','-r100');
                publish_figurePDF(f111,savefile);
                fprintf('Save plot to %s \n',savefile);
            end
           
        end
    end
    
    methods % small auxillary methods
        function out = LoadSaveSensitivity(obj,varargin)
            p = inputParser;
            p.addParameter('SysEffect','Stat',@(x)ischar(x));
            p.addParameter('Mode','load',@(x)ismember(x,{'load','save'}));
            p.parse(varargin{:});
            SysEffect_save = p.Results.SysEffect;
            SaveMode       = p.Results.Mode;
            
            % label
            savedir = [getenv('SamakPath'),'tritium-data/sensitivity/',obj.RunAnaObj.DataSet,'/'];
            MakeDir(savedir);
            switch obj.Mode
                case 'Scan'
                    clLabel = sprintf('_%.0fcl_%sSide',obj.ConfLevel*100,obj.ScanSide);
                case 'Asimov'
                    clLabel = '';
            end
            TDlabel = extractBefore(strrep(obj.RunAnaObj.RunData.RunName,'_',''),'E0');
            if isempty(TDlabel)
                TDlabel = obj.RunAnaObj.RunData.RunName;
            end
            StrCommon = [savedir,sprintf('RunSens%s_%s_%.0feV_%s_budget%.0f%s%s',obj.Mode,TDlabel,...
                obj.RunAnaObj.GetRange,obj.chi2sys,obj.RunAnaObj.SysBudget,clLabel,obj.GetFitLabel)];
            if strcmp(SysEffect_save,'Stat')
                savefile = [strrep(StrCommon,sprintf('%s_budget%.0f%',obj.chi2sys,obj.RunAnaObj.SysBudget),'chi2Stat'),'.mat'];
            else
                savefile = sprintf('%s_%s.mat',StrCommon,SysEffect_save);
            end
            
            if strcmp(obj.RunAnaObj.AnaFlag,'Ring')
                if ~contains(obj.RunAnaObj.fixPar,sprintf('fix %.0f',2*obj.RunAnaObj.nRings+10))
                    % free qU-Offsets
                     savefile = strrep(savefile,'.mat',sprintf('_Ring%s_qUOffsets.mat',obj.RunAnaObj.RingMerge));
                else
                    savefile = strrep(savefile,'.mat',sprintf('_Ring%s.mat',obj.RunAnaObj.RingMerge));
                end
            end
            % load or save file
            out = 1;
            switch SaveMode
                case 'load'
                    if ~exist(savefile,'file') || strcmp(obj.RecomputeFlag,'ON')
                        out = 0; % loading failed
                        %fprintf('NOT loading file %s \n',savefile);
                        return
                    end
                    d = importdata(savefile);
                    
%                     if d.nTrials<obj.GetnTrials(SysEffect_save)
%                           % if trials are not enough
%                         out = 0;
%                        return
%                     end
                    obj.MultiFitResults.(SysEffect_save) = [d.FitResult.par(1),d.FitResult.par(2)+obj.RunAnaObj.ModelObj.Q_i];            
                    obj.Lpar                         = d.Lpar;
                    obj.MultiLpar.(SysEffect_save)   = d.Lpar;
                    obj.CovMatFrac.(SysEffect_save)  = d.FitCMFrac;
                    obj.CovMatShape.(SysEffect_save) = d.FitCMShape;
                    if strcmp(obj.RunAnaObj.fitter,'minuit') && contains(obj.RunAnaObj.minuitOpt,'minos')
                        obj.LparNeg(1) = d.LparNeg(1);
                        obj.LparPos(1) = d.LparPos(1);
                        obj.MultiLparNeg.(SysEffect_save)   = d.LparNeg;
                        obj.MultiLparPos.(SysEffect_save)   = d.LparPos;
                        if isfield(d,'LparMinos') && strcmp(obj.AsymErr,'ON') % average asymmetric uncertainties
                            obj.Lpar                       = d.LparMinos;
                            obj.MultiLpar.(SysEffect_save) = d.LparMinos;                
                        end
                    end
                    fprintf('Load file %s \n',savefile);
                case 'save'
                    Lpar       = obj.Lpar;
                    LparNeg       = obj.LparNeg;
                    LparPos       = obj.LparPos;
                    LparMinos  = obj.LparMinos;
                    FitCMFrac  = obj.RunAnaObj.FitCMFrac;
                    FitCMShape = obj.RunAnaObj.FitCMShape;
                    nTrials    = obj.GetnTrials(SysEffect_save);
                    FitResult  = obj.RunAnaObj.FitResult;
                    save(savefile,'Lpar','FitCMFrac','FitCMShape','nTrials','LparMinos','LparNeg','LparPos','FitResult');
                    
                    obj.MultiLpar.(SysEffect_save)   = Lpar;
                    obj.CovMatFrac.(SysEffect_save)  = FitCMFrac;
                    obj.CovMatShape.(SysEffect_save) = FitCMShape;
                    
                    if strcmp(obj.RunAnaObj.fitter,'minuit') && contains(obj.RunAnaObj.minuitOpt,'minos')
                        if sum(obj.LparMinos)~=0 && strcmp(obj.AsymErr,'ON')% average asymmetric uncertainties
                            obj.Lpar                       = LparMinos;
                            obj.MultiLpar.(SysEffect_save) = LparMinos;
                            obj.MultiLparNeg.(SysEffect_save)   = LparNeg;
                            obj.MultiLparPos.(SysEffect_save)   = LparPos;
                        end
                    end
                    fprintf('Save file %s \n',savefile)
            end
        end
        function y = Convert2Stack(x,startValue)
            %x is array to be converted
            y = zeros(1,numel(x));
            y(1)= startValue;
            for i=2:numel(y)
                y(i)= x(i)-sum(y(1:i-1));
            end
        end
        function [SingleBarX] = GetRangeSingleSys(obj,varargin)
            % Get Range Vector for Single Contributions (MultiBarPlot)
            p = inputParser;
            p.addParameter('PlotRanges',[30, 45, 60],@(x)all(isfloat(x)));
            p.parse(varargin{:});
            PlotRanges = p.Results.PlotRanges;
            % Init Var
            
            Nranges = numel(PlotRanges);
            SingleBarX = zeros(Nranges,obj.nSys);
            % Set ranges (manually), found to be good
            for i=1:Nranges
                if Nranges==1
                    SingleBarX(i,:) = PlotRanges(i)-linspace(0.43,0.855,obj.nSys);
                elseif Nranges==2
                    SingleBarX(i,:) = PlotRanges(i)-linspace(5.39,10.85,obj.nSys);
                elseif Nranges==3
                    SingleBarX(i,:) = PlotRanges(i)-linspace(5.393,10.85,obj.nSys);
                else
                    SingleBarX(i,:) = PlotRanges(i)-linspace(5.393,10.85,obj.nSys);
                end
            end
        end
        function SysUncertainties = GetSysInfo(obj)
            % Systematic Budget Annotation
            [SysErr, ~] = GetSysErr(obj.RunAnaObj.SysBudget);
            
            if strcmp(SysErr.DataDriven,'ON')
                TASR_RelErr = mean(obj.RunAnaObj.Get_DataDriven_RelErr_TASR);
                [qUErr,~]=obj.RunAnaObj.Get_DataDriven_RelErr_qU;
                qU_AbsErr = mean(qUErr(obj.RunAnaObj.RunData.qU<18575));
            else
                TASR_RelErr = SysErr.WGTS_TASR_RelErr;
                qU_AbsErr = 0;
            end
            SysUncertainties = [sprintf('Response Function: '),...
                sprintf('\\Delta \\rhod\\sigma = %.2f%%, ',100*sqrt(SysErr.WGTS_CD_MolPerCm2_RelErr^2+SysErr.ISXsection_RelErr^2)),...
                sprintf('\\DeltaB_a = %.1f%%, \\DeltaB_{s} = %.1f%%, \\DeltaB_{max} = %.1f%%,  ',...
                100*SysErr.MACE_Ba_T_RelErr,...
                100*SysErr.WGTS_B_T_RelErr, 100*SysErr.MACE_Bmax_T_RelErr),...
                sprintf('\nTritium Activity Fluctuations %.1f%%,  Energy Loss: %s , Mean \\DeltaStacking = %.0fmV',100*TASR_RelErr,...
                obj.RunAnaObj.ELossFlag,qU_AbsErr*1e3),...
                sprintf('\nFinal State Distribution:'),...
                sprintf(' Normalization %.0f %%, ',SysErr.FSDNorm_RelErr*100),...
                sprintf('Bin-to-Bin uncorrelated %.0f %% (GS), %.0f %% (ES)',100*SysErr.FSDShapeGS_RelErr,...
                100*SysErr.FSDShapeES_RelErr)];
        end
        function SysEffect_str = GetSysEffect(obj)
            if strcmp(obj.RunAnaObj.chi2,'chi2Stat')
                SysEffect_str = '';
            else
                effects_logic = structfun(@(x)strcmp(x,'ON'),obj.RunAnaObj.FitCM_Obj.SysEffect);
                fields = fieldnames(obj.RunAnaObj.FitCM_Obj.SysEffect);
                SysEffect_str = '';
                for i=1:numel(fields)
                    %more readable labels:
                    if contains(fields{i},'_')
                        fields{i}=strrep(fields{i},'_','-');
                    end
                    if effects_logic(i)==1 && sum(effects_logic)>=2
                        SysEffect_str = strcat(SysEffect_str,' + ',fields(i));
                    elseif effects_logic(i)==1 && sum(effects_logic)==1
                        SysEffect_str = fields(i);
                    end
                end
                if contains(SysEffect_str, ' +RF-EL +RF-BF +RF-RX')
                    SysEffect_str{:} = strrep(SysEffect_str{:},' +RF-EL +RF-BF +RF-RX','Response Function');
                end
                if contains(SysEffect_str, '+TCoff-RAD +TCoff-OTHER')
                    SysEffect_str{:} = strrep(SysEffect_str{:},'+TCoff-RAD +TCoff-OTHER','+ Theoretical Correction ON/OFF');
                end
            end
        end
        function factor = ConvertCl2Std(obj,varargin)
            % factor to convert confidence level into sigma
            % when cl==0, factor =1
            p=inputParser;
            p.addParameter('cl',obj.ConfLevel,@(x)isfloat(x));
            p.parse(varargin{:});
            cl = p.Results.cl;
            
            if cl==0
                factor = 1;
            else
                if strcmp(obj.LimitFlag,'Central')
                    factor = erfinv(cl)*sqrt(2);
                elseif strcmp(obj.LimitFlag,'Up')
                    factor = norminv(cl);
                end
            end
            
        end
        function CreateDir(obj,dir)
            % check if directory exists, if not: create it
            if ~exist(dir,'dir')
                system(['mkdir -p ',dir]);
            end
        end
        function range = GetRange(obj)
            % data range below endpoint (eV)
            range = abs(round(obj.RunAnaObj.RunData.qU(obj.RunAnaObj.exclDataStart),0)-18575);
        end
        function FitLabel = GetFitLabel(obj)
            % label which fitter was used
            % for saving plots, fit results,...
            %special labels for non-default fitting
            if contains(obj.RunAnaObj.minuitOpt,'minos') && strcmp(obj.RunAnaObj.fitter,'minuit')
                FitLabel = '_MinuitMinosFit';
            elseif contains(obj.RunAnaObj.minuitOpt,'migrad') && strcmp(obj.RunAnaObj.fitter,'minuit')
                FitLabel =  '_MinuitMigradFit';
            elseif strcmp(obj.RunAnaObj.fitter,'matlab')
                FitLabel = '_MatlabFit';
            end
            
            FitLabel = [FitLabel,'_',obj.RunAnaObj.DataType];
        end
        function nTrials = GetnTrials(obj,SysEffect)
            if strcmp(obj.RunAnaObj.DataSet,'FirstTritium.katrin')
                nTrials = 1000;
            else
                if ischar(SysEffect)
                    nTrials = 5000;
                    if contains(SysEffect,'RF') %|| contains(SysEffect,'LongPlasma')
                        nTrials = 1000;
                    elseif contains(SysEffect,'LongPlasma') %&& strcmp(obj.RunAnaObj.AnaFlag,'Ring')
                        nTrials = 1000;
                    end
                else
                    if structfun(@(x) contains(x,'RF'),SysEffect)
                        nTrials = 1000;
                    else
                        nTrials = 5000;
                    end
                end
            end
        end
        function PlotWhiteSpace(obj)
            % remove white space around figure
            ax = gca;
            outerpos = ax.OuterPosition;
            ti = ax.TightInset;
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
            ax_height = outerpos(4)- ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];
        end
        function GetDefSysEffect(obj)
            
            switch obj.RunAnaObj.DataSet
                case 'Knm2'
                    %Default
                    obj.SysEffectsAll   = {'TCoff_OTHER','FSD','TASR','RF_EL','RF_BF','RF_RX','Stack','LongPlasma','FPDeff','NP','BkgPT','Bkg'}; %Bkg has to be last
                    obj.SysEffectLeg    = {'Theoretical corrections';'Final-state distribution';...
                        'Tritium activity fluctuations';'Energy-loss function';...
                        'Magnetic fields';'Source scattering';'HV fluctuations';...
                        'Long. source potential';...
                        'Detector efficiency';'Non-Poisson background';'Background PT';'Background slope'};
                case 'Knm1'
                    obj.SysEffectsAll      = {'TCoff_OTHER','FSD','TASR','RF_EL','RF_BF','RF_RX','Stack','FPDeff','NP','Bkg'}; %Bkg has to be last
                    obj.SysEffectLeg      = {'Theoretical corrections';'Final-state distribution';...
                        'Tritium activity fluctuations';'Energy-loss function';...
                        'Magnetic fields';'Source scattering';'HV fluctuations';...
                        'Detector efficiency';'Non-Poisson background';'Background slope'};
            end
        end
    end
end

