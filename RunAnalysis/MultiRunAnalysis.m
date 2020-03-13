% Analysis class for the KATRIN Experiment - Combining list of Runs
%-------------------------------------------------------------------------%
%  This class contains:
%  - Function to stack runs from KATRIN Tritum Data
%  - Computation of the model of the stacked runs
% 
% RunData = Structure with all information for run summaries
%-------------------------------------------------------------------------%
%  L.Schlueter        P. I. Morales Guzman       T. Lasserre
%     MPP/TUM             TUM / MPP                  CEA
%-------------------------------------------------------------------------%
%  Last Update:     (10/12/2018)

classdef MultiRunAnalysis < RunAnalysis & handle
    properties (Access=public)
        Debug;            % Debug Display Flag
        RunList;          % Vector of Runs
        nRuns;            % Number of Runs
        
        % Stack Runs
        StackTolerance;   % Tolerance in the qU to allow for stacking runs
        StackedRuns;      % Stores the numbers of the stacked runs
        NotStackedRuns;   % Stores the numbers of the not stacked runs
        StackFileName;    % Name of file created by stacking, includes all stacked runs
        StackqUCorrFlag;
        
        % Runwise analysis
        SingleRunData;    % Structure with all information for runsummaries of individual runs
        SingleRun_FitResults;            % structure with fit results from runwise fits
        
        % Taylor Expansion to improve stacking TBDIS
        StackTaylorFlag;              % Flag to Attempt to correct the StackSpectrum from being computed with the wrong qU's
        StackTaylorExpCorrection = 0; % Value for the corrections, per sub-runs
        TestRateCorr;                 % Test Correction Flag
        StackingCorr;                 % Stacking Correction (Vector):   TBDIS(corr) = TBDIS(uncorr).*StackingCorr
        % Ring
        CurrentRing;
       
        % Twin option for MultiRunAnalysis
        Twin_SameqUFlag;    %if flag is ON: compute twin runs with same qU as average qU in MultiRun, if OFF: take qU values from data, if a numer (eV) takes same plus gaussian distribution with this sigma
        Twin_SameqUfracFlag;
        Twin_SameCDFlag;
        Twin_SameIsotopFlag; % MolFrac_TT, HT,DT
        
        % Radiative Corrections Flag
        RadiativeFlag;
        
    end
    methods % Constructor
        function obj = MultiRunAnalysis(varargin)
            fprintf('-------------------Start MultiRunAnalysis Constructor----------------- \n')
            %--------------------------------- Parser Start----------------------------------------%
            p = inputParser;
            p.addParameter('DataType','Real',@(x)ismember(x,{'Real','Fake','Twin','FitriumTwin','KafitTwin'}));
            p.addParameter('StackTolerance',0.2,@(x)isfloat(x) && all(x)>0);%0.2
            p.addParameter('StackqUCorrFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('RunList',[],@(x)(isfloat(x) && all(x)>0 || ischar(x)));
            p.addParameter('StackTaylorFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Debug','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('TestRateCorr','OFF',@(x)ismember(x,{'OFF','ActivityScaling','EROS'}));
            p.addParameter('Twin_SameqUFlag','OFF',@(x)ismember(x,{'ON','OFF'}) || @(x)isfloat(x));
            p.addParameter('Twin_SameqUfracFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Twin_SameCDFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Twin_SameIsotopFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('RadiativeFlag','ON',@(x)ismember(x,{'ON','OFF'}));
           
            % Parse unmatched parameters to RunAnalysis.m
            p.KeepUnmatched=1;
            p.parse(varargin{:});
            if( isempty(fieldnames(p.Unmatched)) ); Unmatched={}; else
                Unmatched = reshape(...
                    [fieldnames(p.Unmatched),struct2cell(p.Unmatched)]',...
                    [1,length(fieldnames(p.Unmatched))*2]);
            end
            obj=obj@RunAnalysis(Unmatched{:},'DoRunAnalysisConstructor','OFF'); %Parse to superclass RunAnalysis.m
            
            obj.DataType            = p.Results.DataType;
            obj.RunList             = p.Results.RunList;
            obj.StackTolerance      = p.Results.StackTolerance;
            obj.StackTaylorFlag     = p.Results.StackTaylorFlag;   
            obj.StackqUCorrFlag     = p.Results.StackqUCorrFlag;
            obj.Debug               = p.Results.Debug;            
            obj.TestRateCorr        = p.Results.TestRateCorr;            
            obj.Twin_SameqUFlag     = p.Results.Twin_SameqUFlag;
            obj.Twin_SameqUfracFlag = p.Results.Twin_SameqUfracFlag;
            obj.Twin_SameCDFlag     = p.Results.Twin_SameCDFlag;
            obj.Twin_SameIsotopFlag = p.Results.Twin_SameIsotopFlag;
            obj.RadiativeFlag       = p.Results.RadiativeFlag;

            obj.SingleRun_FitResults = struct('chi2Stat','','chi2CMall','','chi2CMcorr','');
            %------------------------------- Parser End----------------------------------------%
            
            GetSamakPath; %sets current samak path as enviromental variable
            
            if isempty(obj.RunList)
                fprintf(2,'Provide a RUN LIST for MultiRunAnalysis! Exiting Constructor ... \n')
                return;
            end
            
            % StackName
            obj.StackFileName = obj.GetStackFileName; % For Labeling Files 
            
            % Run Lists 
            if ischar(obj.RunList)
                switch obj.DataType
                    case {'Real','Twin','FitriumTwin','KafitTwin'}
                        obj.RunList = obj.GetRunList;
                    case 'Fake'
                        obj.RunList = obj.GetFakeRunList;
                end
            end    

            % Data Set (Label)
            obj.DataSet       = GetDataSet(obj.RunList);
            obj.nRuns         = length(obj.RunList);
            
            obj.SetDefaultInit;
            
            if isempty(obj.PixList)
                obj.PixList = GetPixList(obj.DataSet);
            end
            switch obj.RingMerge
                case 'Default'
                    [obj.PixList,obj.RingPixList] = Ring2PixelDefCombi(obj.RingList,obj.PixList);
                    obj.nRings = 10;
                    obj.RingList = 1:10;
                case 'None'
                    [obj.PixList,obj.RingPixList] = Ring2Pixel(obj.RingList,obj.PixList);
                case 'Full'
                    [obj.PixList,obj.RingPixList] = Ring2PixelCombi(obj.RingList,obj.PixList);
                    obj.nRings = 4;
                    obj.RingList = 1:4;
                case 'Half'
                    [obj.PixList,obj.RingPixList] = Ring2PixelHalfCombi(obj.RingList,obj.PixList);
                    obj.nRings = 2;
                    obj.RingList = 1:2;
                case 'Azi'
                    [obj.PixList,obj.RingPixList] = AziPatch2PixelCombi(obj.RingList,obj.PixList);
                    obj.nRings = 5;
                    obj.RingList = 1:5;
            end

            % Init Analysis: Stack Data, Create Stack Model, Covariance Matrix
            switch obj.Debug
                case 'ON'
                    cprintf('blue','MultiRunAnalysis: StackRuns: CutOnSC=OFF CutOnFitSingleRuns=OFF\n');
            end
            obj.StackRuns('CutOnSC','OFF','SCsigma',3,'CutOnFitSingleRuns','OFF');
            if strcmp(obj.StackqUCorrFlag,'ON')
                obj.ComputeStackqU;
            end
            
            switch obj.Debug
                case 'ON'
                    cprintf('blue','MultiRunAnalysis: SimulateStackRuns\n');
            end
            obj.SetROI;
            obj.SimulateStackRuns();
            obj.InitFitPar;
            obj.SetNPfactor
            
            if ~strcmp(obj.chi2,'chi2Stat') && ~strcmp(obj.chi2,'chi2P')
                switch obj.Debug
                    case 'ON'
                        cprintf('blue','MultiRunAnalysis: ComputeCM\n');
                end
                % read/computes normal covariance matrix        
                % and stacking covariance matrix & sums both
                switch obj.DataType
                    case {'Real','Twin','FitriumTwin','KafitTwin'}
                obj.ComputeCM('RecomputeFlag','OFF');  
                    case 'Fake'
                obj.ComputeCM('RecomputeFlag','OFF'); 
                end
            end
        end
        
    end % Constructor
    
    methods % Data Import Methods begin
        
        function StackRuns(obj,varargin)
            % Build Stacked Data Object
            p = inputParser;
            p.addParameter('saveTD_DataBank','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('saveStackRun','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('CutOnSC','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SCsigma',3,@(x)isfloat(x) && x>0);
            p.addParameter('CutOnFitSingleRuns','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Fitsigma',3,@(x)isfloat(x) && x>0);
            p.parse(varargin{:});
            saveTD_DataBank    = p.Results.saveTD_DataBank;
            CutOnSC            = p.Results.CutOnSC; 
            SCsigma            = p.Results.SCsigma;
            CutOnFitSingleRuns = p.Results.CutOnFitSingleRuns;
            Fitsigma           = p.Results.Fitsigma;
            saveStackRun       = p.Results.saveStackRun;
            
            % Set Path
            switch obj.DataType
                case 'Real'
                    obj.RunData.matFilePath = [getenv('SamakPath'),'tritium-data/mat/',obj.DataSet,'/'];
                case 'Twin'
                    obj.RunData.matFilePath = [getenv('SamakPath'),'tritium-data/mat/','Twin',obj.DataSet,'/'];
                case 'Fake'
                    obj.RunData.matFilePath = [getenv('SamakPath'),'tritium-data/mat/',obj.FakeStudyName,'/'];   
                case 'FitriumTwin'
                    obj.RunData.matFilePath = [getenv('SamakPath'),'tritium-data/mat/Twin_Fitrium_',obj.DataSet,'/'];
                case 'KafitTwin'
                   obj.RunData.matFilePath = [getenv('SamakPath'),'tritium-data/mat/Twin_Kafit_',obj.DataSet,'/'];
                    
            end
            savename = [obj.RunData.matFilePath,obj.StackFileName,'.mat'];
            
            % Load data from single runs into struct
            obj.ReadSingleRunData(); 
            
            % Stacking Cuts: Choose which runs are stacked
           obj.ApplyStackingCuts('CutOnSC',CutOnSC,'sigmaSC',SCsigma,...
               'CutOnFitSingleRuns',CutOnFitSingleRuns,'sigmaFitSingleRuns',Fitsigma); % Apply Cuts -> obj.SingleRunData.Select
           obj.StackedRuns    = obj.SingleRunData.Runs(obj.SingleRunData.Select_all); 
           obj.NotStackedRuns = obj.SingleRunData.Runs(~obj.SingleRunData.Select_all);

            %stack time
            TimeSec                  = sum(obj.SingleRunData.TimeSec(:,obj.SingleRunData.Select_all),2);     
            TimeperSubRunperPixel    = squeeze(sum(obj.SingleRunData.TimeperSubRunperPixel(:,obj.SingleRunData.Select_all,:),2));
            qUfrac                   = TimeperSubRunperPixel./sum(TimeperSubRunperPixel,1);
            
            % not pixel and not subrun dependent variables
            WGTS_CD_MolPerCm2        = obj.StackWmean(obj.SingleRunData.WGTS_CD_MolPerCm2,obj.SingleRunData.TimeSec);
            WGTS_MolFrac_TT          = obj.StackWmean(obj.SingleRunData.WGTS_MolFrac_TT,obj.SingleRunData.TimeSec);
            WGTS_MolFrac_DT          = obj.StackWmean(obj.SingleRunData.WGTS_MolFrac_DT,obj.SingleRunData.TimeSec);
            WGTS_MolFrac_HT          = obj.StackWmean(obj.SingleRunData.WGTS_MolFrac_HT,obj.SingleRunData.TimeSec);
            ISXsection               = obj.StackWmean(obj.SingleRunData.ISXsection,obj.SingleRunData.TimeSec);
            if isfield(obj.SingleRunData,'RW_BiasVoltage')
                RW_BiasVoltage           = obj.StackWmean(obj.SingleRunData.RW_BiasVoltage,obj.SingleRunData.TimeSec);  
            end
            
            % pixel dependent, not subrun dependent
            SingleRunTimeperPixel = squeeze(mean(obj.SingleRunData.TimeperSubRunperPixel,1))'; % time per run per pixel
            MACE_Ba_T                = obj.StackWmean(obj.SingleRunData.MACE_Ba_T,SingleRunTimeperPixel);
            
            if isfield(obj.SingleRunData,'qU_RM') % rate monitor point
                qU_RM      = obj.StackWmean(obj.SingleRunData.qU_RM,SingleRunTimeperPixel);
                qUfrac_RM  = obj.StackWmean(obj.SingleRunData.qUfrac_RM,SingleRunTimeperPixel);
                TBDIS_RM   = squeeze(sum(obj.SingleRunData.TBDIS_RM(:,obj.SingleRunData.Select_all),2));
            end
            
            if isfield(obj.SingleRunData,'MACE_Bmax_T')
                MACE_Bmax_T = obj.StackWmean(obj.SingleRunData.MACE_Bmax_T,SingleRunTimeperPixel);
            end
            
            % not pixel, but subrun dependent variables
            SingleRunTimeperSubRun   = mean(obj.SingleRunData.TimeperSubRunperPixel(:,:,obj.PixList),3); % time per run per subrun
            WGTS_CD_MolPerCm2_SubRun = obj.StackWmean(obj.SingleRunData.WGTS_CD_MolPerCm2_SubRun,SingleRunTimeperSubRun)';
            WGTS_MolFrac_DT_SubRun   = obj.StackWmean(obj.SingleRunData.WGTS_MolFrac_DT_SubRun,SingleRunTimeperSubRun)';
            WGTS_MolFrac_HT_SubRun   = obj.StackWmean(obj.SingleRunData.WGTS_MolFrac_HT_SubRun,SingleRunTimeperSubRun)';
            WGTS_MolFrac_TT_SubRun   = obj.StackWmean(obj.SingleRunData.WGTS_MolFrac_TT_SubRun,SingleRunTimeperSubRun)';
            
            % pixel and subrun dependent
            EffCorr                  = obj.StackWmean(obj.SingleRunData.EffCorr,obj.SingleRunData.TimeperSubRunperPixel);
            TBDIS                    = squeeze(sum(obj.SingleRunData.TBDIS(:,obj.SingleRunData.Select_all,:),2));
            TBDISE                   = sqrt(TBDIS./EffCorr);            
            qU                       = obj.StackWmean(obj.SingleRunData.qU,obj.SingleRunData.TimeperSubRunperPixel);
            
            matFilePath       = obj.RunData.matFilePath;
            
            if isfield(obj.SingleRunData,'TBDIS14keV')
                TBDIS14keV      = squeeze(sum(obj.SingleRunData.TBDIS14keV(:,obj.SingleRunData.Select_all,:),2));
                TBDIS14keV_RM   = squeeze(sum(obj.SingleRunData.TBDIS14keV_RM(:,obj.SingleRunData.Select_all),2));
            end
            
            obj.RunData = struct(...
                'TBDIS',TBDIS,'TBDIS_Default',TBDIS,'TBDISE',TBDISE,'EffCorr',EffCorr,...
                'qU',qU,'TimeSec',TimeSec,'qUfrac',qUfrac,...
                'MACE_Ba_T',MACE_Ba_T,...
                'TimeperSubRunperPixel',TimeperSubRunperPixel,...
                'ISXsection',ISXsection,...
                'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2,'WGTS_CD_MolPerCm2_SubRun',WGTS_CD_MolPerCm2_SubRun,...
                'WGTS_MolFrac_TT',WGTS_MolFrac_TT,'WGTS_MolFrac_TT_SubRun',WGTS_MolFrac_TT_SubRun,...
                'WGTS_MolFrac_DT',WGTS_MolFrac_DT,'WGTS_MolFrac_DT_SubRun',WGTS_MolFrac_DT_SubRun,...
                'WGTS_MolFrac_HT',WGTS_MolFrac_HT,'WGTS_MolFrac_HT_SubRun',WGTS_MolFrac_HT_SubRun,...
                'matFilePath',matFilePath,'RunName',strrep(obj.StackFileName,'Twin',''));
            if isfield(obj.SingleRunData,'MACE_Bmax_T')
                obj.RunData.MACE_Bmax_T    = MACE_Bmax_T;
                obj.RunData.RW_BiasVoltage = RW_BiasVoltage;
            end
            if isfield(obj.SingleRunData,'qU_RM') %rate monitor
                obj.RunData.qU_RM            = qU_RM;
                obj.RunData.qUfrac_RM        = qUfrac_RM;
                obj.RunData.TBDIS_RM         = TBDIS_RM;
                obj.RunData.TBDIS_RM_Default = TBDIS_RM;
            end
            
            if isfield(obj.SingleRunData,'TBDIS14keV')
                obj.RunData.TBDIS14keV    = TBDIS14keV;
                obj.RunData.TBDIS14keV_RM = TBDIS14keV_RM; 
            end
            
            if strcmp(saveStackRun,'ON')
                save(savename,...
                    'TBDIS','TBDISE','EffCorr','qU','TimeSec','qUfrac',...
                    'MACE_Ba_T',...
                    'TimeperSubRunperPixel',...
                    'ISXsection',...
                    'WGTS_CD_MolPerCm2','WGTS_CD_MolPerCm2_SubRun',...
                    'WGTS_MolFrac_TT','WGTS_MolFrac_DT','WGTS_MolFrac_HT',...
                    'WGTS_MolFrac_DT_SubRun','WGTS_MolFrac_HT_SubRun','WGTS_MolFrac_TT_SubRun','matFilePath',...
                    '-v7.3','-nocompression');
                if isfield(obj.SingleRunData,'qU_RM') %rate monitor
                    save(savename,'qU_RM','-append')
                    save(savename,'qUfrac_RM','-append')
                    save(savename,'TBDIS_RM','-append')
                end
                if isfield(obj.SingleRunData,'MACE_Bmax_T') %KNM2 new variables
                    save(savename,'RW_BiasVoltage','-append')
                    save(savename,'MACE_Bmax_T','-append')
                end
                if isfield(obj.SingleRunData,'TBDIS14keV')
                     save(savename,'TBDIS14keV','-append');
                     save(savename,'TBDIS14keV_RM','-append'); 
                end
            end
            
            switch obj.Debug
                case 'ON'
                    fprintf(2,'Showing Stacked Data (data for Ring selected later) \n')
                    disp(obj.RunData);
            end
            
            obj.StackPixel;
            obj.StackPixelSingleRun;
            
            % Stacking Correction:
%             obj.ComputeStackingCorr('Mean_WGTS_MolFrac_TT_SubRun',WGTS_MolFrac_TT_SubRun',...
%                 'Mean_WGTS_MolFrac_DT_SubRun',WGTS_MolFrac_DT_SubRun',...
%                 'Mean_WGTS_MolFrac_HT_SubRun',WGTS_MolFrac_HT_SubRun',...
%                 'Mean_WGTS_CD_MolPerCm2_SubRun',WGTS_CD_MolPerCm2_SubRun',...
%                 'Mean_qUfrac',obj.RunData.qUfrac,...
%                 'qUfracFlag','OFF','CDFlag','OFF','IsoFlag','OFF');
            
            if any(obj.StackingCorr~=1)
                % Apply Stacking Correction: only apply to signal, not background
                % therefore: Bkg subtraction before correction
                qi=18573.7;
                BkgRate         = mean(obj.SingleRunData.TBDIS(obj.SingleRunData.qU(:,1)>qi,:)./...
                    (obj.SingleRunData.qUfrac(obj.SingleRunData.qU(:,1)>qi,:)...
                    .*obj.SingleRunData.TimeSec));
                
                Bkg             = BkgRate.*obj.SingleRunData.qUfrac.*obj.SingleRunData.TimeSec;
                TBDIScorr       = obj.SingleRunData.TBDIS-Bkg;
                TBDIScorr(obj.SingleRunData.qU(:,1)<qi,:) = TBDIScorr(obj.SingleRunData.qU(:,1)<qi,:).*obj.StackingCorr(obj.SingleRunData.qU(:,1)<qi);
                TBDIScorr       = TBDIScorr+Bkg;
                TBDIS = squeeze(sum(TBDIScorr(:,obj.SingleRunData.Select_all),2));
                fprintf(2,'Stacking Correction (Activity) Applied! \n');
            else
                TBDIS = squeeze(sum(obj.SingleRunData.TBDIS(:,obj.SingleRunData.Select_all,:),2));
            end
            
            obj.RunData.TBDIS           = TBDIS;
            
            switch obj.StackTaylorFlag
                case 'ON' %with Taylor Expansion Correction
                    % Taylor expansion around RunData.qU
                    if ismember(obj.AnaFlag,{'StackPixel'})
                        qUmean            = repmat(obj.RunData.qU,[1,numel(obj.StackedRuns)]);
                        qUsingleRun       = obj.SingleRunData.qU(:,obj.SingleRunData.Select_all,:);
                        Grad1Corr         = -diffxy(qUsingleRun,obj.SingleRunData.TBDIS(:,obj.SingleRunData.Select_all,:)).*(qUmean-qUsingleRun);
                        Grad2Corr         = -0.5*diffxy(qUsingleRun,Grad1Corr).*(qUmean-qUsingleRun).^2;
                        Grad3Corr         = -1/6*diffxy(qUsingleRun,Grad2Corr).*(qUmean-qUsingleRun).^3;
                        GradCorr          = Grad1Corr+Grad2Corr+Grad3Corr;
                        obj.RunData.TBDIS = obj.RunData.TBDIS + mean(GradCorr,2);
                        obj.StackTaylorExpCorrection   = obj.StackTaylorExpCorrection + GradCorr;
                    elseif ismember(obj.AnaFlag,{'MultiPixel','Ring'})
                        % doesnt work yet
                    end
                case 'OFF' %without Taylor Expansion Correction
            end
            
            switch saveTD_DataBank
                case 'ON'
                    switch obj.DataType
                        case 'Real'
                            TD = ['Run',obj.StackFileName]; RunTime = TimeSec; %ok
                            save([getenv('SamakPath'),'/simulation/katrinsetup/TD_DataBank/Run',...
                                obj.StackFileName,'.mat'],...
                                'qU','qUfrac','RunTime','TD','-v7.3','-nocompression')
                        case 'Fake'
                            TD = [obj.StackFileName]; RunTime = TimeSec; %ok
                            save([getenv('SamakPath'),'/simulation/katrinsetup/TD_DataBank/',...
                                obj.StackFileName,'.mat'],... %WARNING Thierry Feb 13 2019
                                'qU','qUfrac','RunTime','TD','-v7.3','-nocompression')
                    end
                case 'OFF'
                    % Don't save
            end
             
        end % StackRuns
        function ComputeStackingCorr(obj,varargin)
            % calculate correction to stacked data based on SC parameters
            % size: nqU
            
            % Build Stacked Data Object
            p = inputParser;
            p.addParameter('IsoFlag','ON',@(x)ismember(x,{'ON','OFF'}));                  
            p.addParameter('CDFlag','ON',@(x)ismember(x,{'ON','OFF'})); 
            p.addParameter('qUfracFlag','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Mean_WGTS_MolFrac_TT_SubRun',1,@(x)isfloat(x) && all(x>0));  
            p.addParameter('Mean_WGTS_MolFrac_DT_SubRun',1,@(x)isfloat(x) && all(x>0));
            p.addParameter('Mean_WGTS_MolFrac_HT_SubRun',1,@(x)isfloat(x) && all(x>0));
            p.addParameter('Mean_WGTS_CD_MolPerCm2_SubRun',1,@(x)isfloat(x) && all(x>0));
            p.addParameter('Mean_qUfrac',1,@(x)isfloat(x) && all(x>0));
            
            p.parse(varargin{:});
            
            IsoFlag                       = p.Results.IsoFlag;
            CDFlag                        = p.Results.CDFlag;
            qUfracFlag                    = p.Results.qUfracFlag;
            Mean_WGTS_MolFrac_TT_SubRun   = p.Results.Mean_WGTS_MolFrac_TT_SubRun;
            Mean_WGTS_MolFrac_DT_SubRun   = p.Results.Mean_WGTS_MolFrac_DT_SubRun;
            Mean_WGTS_MolFrac_HT_SubRun   = p.Results.Mean_WGTS_MolFrac_HT_SubRun;
            Mean_WGTS_CD_MolPerCm2_SubRun = p.Results.Mean_WGTS_CD_MolPerCm2_SubRun;
            Mean_qUfrac                   = p.Results.Mean_qUfrac;
            
            if isempty(obj.SingleRunData)
                obj.ReadSingleRunData;
            end
            
            if strcmp(obj.TestRateCorr,'OFF')
                obj.StackingCorr = ones(size(obj.SingleRunData.qU));
                return
            end
            
            if strcmp(IsoFlag,'ON')
                TTcorr = (Mean_WGTS_MolFrac_TT_SubRun./obj.SingleRunData.WGTS_MolFrac_TT_SubRun);
                DTcorr = (Mean_WGTS_MolFrac_DT_SubRun./obj.SingleRunData.WGTS_MolFrac_DT_SubRun);
                HTcorr = (Mean_WGTS_MolFrac_HT_SubRun./obj.SingleRunData.WGTS_MolFrac_HT_SubRun);
                if contains(obj.DataSet,'FirstTritium')
                    StackingCorrIso = Mean_WGTS_MolFrac_DT_SubRun.*DTcorr;
                else
                    StackingCorrIso = (Mean_WGTS_MolFrac_TT_SubRun.*TTcorr + Mean_WGTS_MolFrac_DT_SubRun.*DTcorr + Mean_WGTS_MolFrac_HT_SubRun.*HTcorr);
                end
            else
                StackingCorrIso = ones(size(obj.SingleRunData.qU));
            end
            
            
            if strcmp(CDFlag,'ON')
                StackingCorrCD  =  Mean_WGTS_CD_MolPerCm2_SubRun./obj.SingleRunData.WGTS_CD_MolPerCm2_SubRun;
            else
                StackingCorrCD = ones(size(obj.SingleRunData.qU));
            end
            
            if strcmp(qUfracFlag,'ON')
                StackingCorrqUfrac  =  Mean_qUfrac./obj.SingleRunData.qUfrac;
            else
                StackingCorrqUfrac = ones(size(obj.SingleRunData.qU));
            end
            
            obj.StackingCorr = StackingCorrIso.*StackingCorrCD.*StackingCorrqUfrac;
            
            %  rhodActivity = (obj.SingleRunData.WGTS_MolFrac_DT/2+obj.SingleRunData.WGTS_MolFrac_HT/2+obj.SingleRunData.WGTS_MolFrac_TT).*obj.SingleRunData.WGTS_CD_MolPerCm2;
            %            % obj.StackingCorr = (rhodActivity./mean(rhodActivity)).^-1;
            %                 case 'EROS'
%             obj.StackingCorr    = obj.R200RateErosCorrection();
%             end
            
        end
        function ApplyStackingCuts(obj,varargin)
            % Apply Cuts to define Golden Run Stacking List
            % 1) qU tolerance
            % 2) DT/TT/HT Nx sigmas
            % 3) Colmuns density Nx sigmas
            % 4) Fit Parameters Nx sigmas --> based on statistical fit only
            
            p = inputParser;
            p.addParameter('CutOnSC','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SigmaSC',3,@(x)isfloat(x) && x>0);
            p.addParameter('CutOnFitSingleRuns','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SigmaFitSingleRuns',3,@(x)isfloat(x) && x>0);
            p.parse(varargin{:});
            CutOnSC            = p.Results.CutOnSC;
            SigmaSC            = p.Results.SigmaSC;
            CutOnFitSingleRuns = p.Results.CutOnFitSingleRuns;
            SigmaFitSingleRuns = p.Results.SigmaFitSingleRuns;
            
            % Load Fit Results - Stat
            switch CutOnFitSingleRuns
                case 'ON'
                    save_name = sprintf('./results/Results_%s_%s_%s_%.0feVbelowE0_fixPar%s.mat',...
                obj.StackFileName,'chi2Stat','all',obj.ModelObj.Q_i-obj.ModelObj.qU(obj.exclDataStart),strrep(obj.fixPar,' ',''));
            if exist(save_name,'file')
                sprintf('FitResults Stat Only already computed - Loading from File \n')
                d = importdata(save_name);
                FitResults = d.FitResults;
                obj.SingleRun_FitResults.chi2Stat = FitResults;
            else
                fprintf('FitResults do not exist. Do you want to recompute? \n');
                return
            end
            end
                        
            % ------------------------Quality Cuts Start-------------------------:
            nNotStackedRuns = zeros(10,1); % Number of Not Stacked Runs

            % Prepare qU-Cut
            qUmeanPixel = obj.StackWmean(obj.SingleRunData.qU,obj.SingleRunData.TimeperSubRunperPixel);
            qUmean     = mean(qUmeanPixel(:,obj.PixList),2);
               
            % Prepare SC-Cut
            switch CutOnSC
                case 'ON'
            % Prepare DT-Cut
            DTmean     = obj.StackWmean(obj.SingleRunData.WGTS_MolFrac_DT,obj.SingleRunData.TimeSec);
            DTsigma    = std(obj.SingleRunData.WGTS_MolFrac_DT);
            % Prepare TT-Cut
            TTmean     = obj.StackWmean(obj.SingleRunData.WGTS_MolFrac_TT,obj.SingleRunData.TimeSec);
            TTsigma    = std(obj.SingleRunData.WGTS_MolFrac_TT);
            % Prepare HT-Cut
            HTmean     = obj.StackWmean(obj.SingleRunData.WGTS_MolFrac_HT,obj.SingleRunData.TimeSec);
            HTsigma    = std(obj.SingleRunData.WGTS_MolFrac_HT);
             % Prepare rhoD-Cit
            rhoDmean   = obj.StackWmean(obj.SingleRunData.WGTS_CD_MolPerCm2,obj.SingleRunData.TimeSec);
            rhoDsigma  = std(obj.SingleRunData.WGTS_CD_MolPerCm2);
            end
            
            % Prepare Fit-Cut
            switch CutOnFitSingleRuns
                case 'ON'
                    if ~isempty(obj.SingleRun_FitResults.chi2Stat)
                        E0mean  = wmean(obj.SingleRun_FitResults.chi2Stat.E0,1./obj.SingleRun_FitResults.chi2Stat.E0Err.^2);
                        E0sigma = std(obj.SingleRun_FitResults.chi2Stat.E0);
                        E0Errmean  = mean(obj.SingleRun_FitResults.chi2Stat.E0Err);
                        E0Errsigma = std(obj.SingleRun_FitResults.chi2Stat.E0Err);
                        
                        Bmean  = wmean(obj.SingleRun_FitResults.chi2Stat.B,1./obj.SingleRun_FitResults.chi2Stat.BErr.^2);
                        Bsigma = std(obj.SingleRun_FitResults.chi2Stat.B);
                        BErrmean  = mean(obj.SingleRun_FitResults.chi2Stat.BErr);
                        BErrsigma = std(obj.SingleRun_FitResults.chi2Stat.BErr);
                        
                        Nmean  = wmean(obj.SingleRun_FitResults.chi2Stat.N,1./obj.SingleRun_FitResults.chi2Stat.NErr.^2);
                        Nsigma = std(obj.SingleRun_FitResults.chi2Stat.N);
                        NErrmean  = mean(obj.SingleRun_FitResults.chi2Stat.NErr);
                        NErrsigma = std(obj.SingleRun_FitResults.chi2Stat.NErr);
                    end
            end
            
            % Apply Cuts
            for i=1:10 % iterations
                switch obj.AnaFlag
                    case {'StackPixel','Ring'}
                        % 1. qU-cut
                        obj.SingleRunData.Select_qU = squeeze(all( abs(...
                            mean(obj.SingleRunData.qU(obj.exclDataStart:end,:,obj.PixList),3) ...
                            - qUmean(obj.exclDataStart:end,:)) < obj.StackTolerance));
                        
                        switch CutOnSC
                            case 'ON'
                                % 2 DT-cut
                                obj.SingleRunData.Select_WGTS_MolFrac_DT = ...
                                    (abs(obj.SingleRunData.WGTS_MolFrac_DT - DTmean) < SigmaSC*DTsigma);
                                if ~strcmp(obj.DataSet,'FirstTritium.katrin')
                                    % 3 TT-cut
                                    obj.SingleRunData.Select_WGTS_MolFrac_TT = ...
                                        (abs(obj.SingleRunData.WGTS_MolFrac_TT - TTmean) < SigmaSC*TTsigma);
                                    % 4 HT-cut
                                    obj.SingleRunData.Select_WGTS_MolFrac_HT = ...
                                        (abs(obj.SingleRunData.WGTS_MolFrac_HT - HTmean) < SigmaSC*HTsigma);
                                end
                                % 5 rhod Cut
                                obj.SingleRunData.Select_WGTS_CD_MolPerCm2 = ...
                                    (abs(obj.SingleRunData.WGTS_CD_MolPerCm2 - rhoDmean) < SigmaSC*rhoDsigma);
                        end
                        
                        switch CutOnFitSingleRuns
                            case 'ON'
                                if ~isempty(obj.SingleRun_FitResults.chi2Stat)
                                % 6 E0 Fit Single Run (stat)
                                obj.SingleRunData.Select_E0 = ...
                                    (abs(obj.SingleRun_FitResults.chi2Stat.E0 - E0mean) < SigmaFitSingleRuns*E0sigma);
                                % 7 B Fit Single Run (stat)
                                obj.SingleRunData.Select_B = ...
                                    (abs(obj.SingleRun_FitResults.chi2Stat.B - Bmean) < SigmaFitSingleRuns*Bsigma);
                                % 8 N Fit Single Run (stat)
                                obj.SingleRunData.Select_N = ...
                                    (abs(obj.SingleRun_FitResults.chi2Stat.N - Nmean) < SigmaFitSingleRuns*Nsigma);
                                % 9 E0 Fit Single Run (stat)
                                obj.SingleRunData.Select_E0Err = ...
                                    (abs(obj.SingleRun_FitResults.chi2Stat.E0Err - E0Errmean) < SigmaFitSingleRuns*E0Errsigma);
                                % 10 B Fit Single Run (stat)
                                obj.SingleRunData.Select_BErr = ...
                                    (abs(obj.SingleRun_FitResults.chi2Stat.BErr - BErrmean) < SigmaFitSingleRuns*BErrsigma);
                                % 11 N Fit Single Run (stat)
                                obj.SingleRunData.Select_NErr = ...
                                    (abs(obj.SingleRun_FitResults.chi2Stat.NErr - NErrmean) < SigmaFitSingleRuns*NErrsigma);
                                end
                        end
                    case {'MultiPixel'}
                        % 1. qU-cut
                        obj.SingleRunData.Select_qU = ...
                            squeeze(all(all(abs(permute(obj.SingleRunData.qU - qUmean,[2,3,1])) < obj.StackTolerance)));
                end
                
                % Collect all cuts to "Selection_all" (so far only qU-cut)
                obj.SingleRunData.Select_all(~obj.SingleRunData.Select_qU)                       = 0;
                
                switch CutOnSC
                    case 'ON'
                        obj.SingleRunData.Select_all(~obj.SingleRunData.Select_WGTS_MolFrac_DT)  = 0;
                        obj.SingleRunData.Select_all(~obj.SingleRunData.Select_WGTS_MolFrac_TT)  = 0;
                        obj.SingleRunData.Select_all(~obj.SingleRunData.Select_WGTS_MolFrac_HT)  = 0;
                        obj.SingleRunData.Select_all(~obj.SingleRunData.Select_WGTS_CD_MolPerCm2)= 0;
                end
                switch CutOnFitSingleRuns
                    case 'ON'
                        if ~isempty(obj.SingleRun_FitResults.chi2Stat)
                            obj.SingleRunData.Select_all(~obj.SingleRunData.Select_E0)   = 0;
                            obj.SingleRunData.Select_all(~obj.SingleRunData.Select_B)    = 0;
                            obj.SingleRunData.Select_all(~obj.SingleRunData.Select_N)    = 0;
                            obj.SingleRunData.Select_all(~obj.SingleRunData.Select_E0Err)= 0;
                            obj.SingleRunData.Select_all(~obj.SingleRunData.Select_BErr) = 0;
                            obj.SingleRunData.Select_all(~obj.SingleRunData.Select_NErr) = 0;
                        end
                end
                
                % Compute new mean for next iteration
                qUmeanPixel = obj.StackWmean(obj.SingleRunData.qU,obj.SingleRunData.TimeperSubRunperPixel);
                qUmean     = mean(qUmeanPixel(:,obj.PixList),2);
                 
                switch CutOnSC
                    case 'ON'
                        DTmean = obj.StackWmean(obj.SingleRunData.WGTS_MolFrac_DT,obj.SingleRunData.TimeSec);
                        DTsigma = std(obj.SingleRunData.WGTS_MolFrac_DT(obj.SingleRunData.Select_all));
                        TTmean = obj.StackWmean(obj.SingleRunData.WGTS_MolFrac_TT,obj.SingleRunData.TimeSec);
                        TTsigma = std(obj.SingleRunData.WGTS_MolFrac_TT(obj.SingleRunData.Select_all));
                        HTmean = obj.StackWmean(obj.SingleRunData.WGTS_MolFrac_HT,obj.SingleRunData.TimeSec);
                        HTsigma = std(obj.SingleRunData.WGTS_MolFrac_HT(obj.SingleRunData.Select_all));
                        rhoDmean = obj.StackWmean(obj.SingleRunData.WGTS_CD_MolPerCm2,obj.SingleRunData.TimeSec);
                        rhoDsigma = std(obj.SingleRunData.WGTS_CD_MolPerCm2(obj.SingleRunData.Select_all));
                end
                
                switch CutOnFitSingleRuns
                    case 'ON'
                        if ~isempty(obj.SingleRun_FitResults.chi2Stat)
                            E0mean  = wmean(obj.SingleRun_FitResults.chi2Stat.E0(obj.SingleRunData.Select_all),1./obj.SingleRun_FitResults.chi2Stat.E0Err(obj.SingleRunData.Select_all).^2);
                            E0sigma = std(obj.SingleRun_FitResults.chi2Stat.E0(obj.SingleRunData.Select_all));
                            E0Errmean  = mean(obj.SingleRun_FitResults.chi2Stat.E0Err(obj.SingleRunData.Select_all));
                            E0Errsigma = std(obj.SingleRun_FitResults.chi2Stat.E0Err(obj.SingleRunData.Select_all));
                            
                            Bmean  = wmean(obj.SingleRun_FitResults.chi2Stat.B(obj.SingleRunData.Select_all),1./obj.SingleRun_FitResults.chi2Stat.BErr(obj.SingleRunData.Select_all).^2);
                            Bsigma = std(obj.SingleRun_FitResults.chi2Stat.B(obj.SingleRunData.Select_all));
                            BErrmean  = mean(obj.SingleRun_FitResults.chi2Stat.BErr(obj.SingleRunData.Select_all));
                            BErrsigma = std(obj.SingleRun_FitResults.chi2Stat.BErr(obj.SingleRunData.Select_all));
                            
                            Nmean  = wmean(obj.SingleRun_FitResults.chi2Stat.N(obj.SingleRunData.Select_all),1./obj.SingleRun_FitResults.chi2Stat.NErr(obj.SingleRunData.Select_all).^2);
                            Nsigma = std(obj.SingleRun_FitResults.chi2Stat.N(obj.SingleRunData.Select_all));
                            NErrmean  = mean(obj.SingleRun_FitResults.chi2Stat.NErr(obj.SingleRunData.Select_all));
                            NErrsigma = std(obj.SingleRun_FitResults.chi2Stat.NErr(obj.SingleRunData.Select_all));

                        end
                end
                
                % exit iteration when nothing changes
                nNotStackedRuns(i+1) = sum(~obj.SingleRunData.Select_all);
                if (nNotStackedRuns(i+1)-nNotStackedRuns(i))==0
                    fprintf(2,'ApplyStackingCuts: Exit Cut Routine after %.0f Iterations \n',i)
                    break
                end
            end
            
            
            % ----------------------Quality Cuts End------------------------
        end    
        function out = StackWmean(obj,par,time)
            %             nqU = size(obj.SingleRunData.qU,1);
            %   auxillary function: takes weighted mean over time
            %             DimPar = sum(size(par)~=1);
            %             if DimPar==1                                   % NOT pixel and NOT qU dependent, only run dependent
            %                 weight_m = obj.SingleRunData.TimeSec./time;
            %             elseif DimPar==2 && all(size(par)~=148)        % NOT pixel dependent, but qU and run dependent
            %                 TimeperSubRun = mean(obj.SingleRunData.TimeperSubRunperPixel(:,:,obj.PixList),3);
            %                 weight_m = squeeze(TimeperSubRun./time);
            %             elseif DimPar==2                               % NOT qU dependent, but pixel and run dependent
            %                  weight_m = squeeze(mean(obj.SingleRunData.TimeperSubRunperPixel./time))'; %
            %             elseif  DimPar==3                              % run, pixel and qU dependent
            %                  weight_m = permute(permute(obj.SingleRunData.TimeperSubRunperPixel,[1,3,2])./time,[1,3,2]);
            %             end
            out = squeeze(wmean(...
                par(:,obj.SingleRunData.Select_all,:),...
                time(:,obj.SingleRunData.Select_all,:),2));%weight_m(:,obj.SingleRunData.Select_all,:)
        end
        function [SingleRunRate,SingleRunqU] = ComputeStackqU(obj,varargin) 
            %% ----------------------------------------------------------------------------------------------------
            % calculate qU-vector for stacked spectrum according to shape of spectrum
            % method: 
            % 1. calculate nRun spectra, identical slow control parameter but different retarding potentials
            % 2. stack counts of nruns
            % 3. determine function qU(countrate) with a linear fit or interpolation
            % 4. determine qU_Stack= qU(countrate_stack)
            % 5. save results for this runlist
            %% ----------------------------------------------------------------------------------------------------
            
            % load if already computed
            savedir = [getenv('SamakPath'),'tritium-data/StackingCorrectionqU/'];
             qUStack_i = obj.RunData.qU;
            if ~exist(savedir,'dir')
              system(['mkdir ',savedir]);  
            end
            
            savename = [savedir,sprintf('StackqUCorr_%s_%.0fruns.mat',strrep(obj.StackFileName,'Twin',''),numel(obj.RunList))];
            
            if exist(savename,'file')
                load(savename,'SingleRunqU','SingleRunTBDIS','SingleRunTime','SingleRunRate','StackRate');
            else
                %% compute single run spectra
                %% add normalization and background from twins
                FitResultsPath = [getenv('SamakPath'),'tritium-data/fit/',obj.DataSet];
                fitfile = sprintf('%s/%s%s_%s_%.0fbE0_fixPar%s.mat',...
                    FitResultsPath,'Twin',...
                    strrep(obj.StackFileName,'Twin',''),...
                    'chi2Stat',...
                    39,...
                    strrep('1 2 5 6 7 8 9 10 11',' ',''));
                d = importdata(fitfile);
                
                % init
                SingleRunTBDIS = zeros(size(obj.SingleRunData.TBDIS));
                SingleRunqU = [repmat(qUStack_i(1),obj.nRuns,1)';...
                    obj.SingleRunData.qU(2:end-1,:);...
                    repmat(qUStack_i(end),obj.nRuns,1)'];
                
                % calculate spectra in vicinity of qU weighted avererage
                progressbar('calculte stack qU');
                for i=1:numel(obj.SingleRunData.Runs)
                    progressbar(i/numel(obj.SingleRunData.Runs))
                    obj.SimulateStackRuns('qU',SingleRunqU(:,i),...
                        'qUfrac',obj.SingleRunData.qUfrac(:,i),...
                        'TimeSec',obj.SingleRunData.TimeSec(i),...
                        'BKG_RateAllFPDSec',0);
                    
                    obj.ModelObj.ComputeTBDDS('N_bias',d.N,'B_bias',d.B); obj.ModelObj.ComputeTBDIS;
                    SingleRunTBDIS(:,i) = obj.ModelObj.TBDIS;
                end
                
                % save
                SingleRunqU = obj.SingleRunData.qU;
                SingleRunTime = obj.SingleRunData.qUfrac.*obj.SingleRunData.TimeSec;
                SingleRunRate = SingleRunTBDIS./SingleRunTime;
                StackRate = sum(SingleRunTBDIS,2)./sum(SingleRunTime,2);
                save(savename,'SingleRunqU','SingleRunTBDIS','SingleRunTime','SingleRunRate','StackRate');
            end
           %% get qUstack by either a linear fit or spline interpolation
            Mode = 'interpol'; % 'fit'
            FitDim = 1; % 1 = linear fit
            nqU = size(obj.RunData.qU,1);
            SubRunFitCoeff = zeros(nqU,FitDim+1);
            StackqU = qUStack_i;

            for i=2:nqU
                if obj.RunData.qU(i)>18575
                    continue % do nothing, just take normal average
                else
                    switch Mode
                        case 'fit'
                            SubRunFitCoeff(i,:) = polyfit(SingleRunRate(i,:),SingleRunqU(i,:),FitDim);
                            StackqU(i) = StackRate(i)*SubRunFitCoeff(i,1)+SubRunFitCoeff(i,2);
                        case 'interpol'
                            %delete redundant points
                            [qUUnique,iNsame,~] = unique(SingleRunqU(i,:),'stable');
                            Nsame = zeros(size(SingleRunRate(i,:)));
                            Nsame(iNsame) = 1;
                            RateUnique = SingleRunRate(i,logical(Nsame));
                            
                            % choose only points close to stackrate for interpolation
                            SelectPoints = abs(RateUnique-StackRate(i))<=0.1*std(RateUnique);
                            StackqU(i) = interp1(RateUnique(SelectPoints),qUUnique(SelectPoints),StackRate(i),'spline');
                           
                           % StackqU(i) = interp1(SingleRunRate(i,logical(Nsame)),qUUnique,StackRate(i),'spline');
                    end
                end
            end

            switch obj.Debug
                case 'ON'
                    f22 = figure('Renderer','opengl');
                    set(f22, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.5]);
                    
                    qu = 17;
                    Q_i = 18575;
                    pSingle = plot(SingleRunRate(qu,:),SingleRunqU(qu,:)-18575,'x','MarkerSize',9,'Color',rgb('DarkSlateGray'));
                    hold on;
                    pStackwMean = plot(StackRate(qu),wmean(SingleRunqU(qu,:),SingleRunTime(qu,:))-18575,...
                        'o','MarkerFaceColor',rgb('SlateGray'),'MarkerSize',9,'Color',rgb('SlateGray'));
                    pStack = plot(StackRate(qu),StackqU(qu)-18575,...
                        'o','MarkerFaceColor',rgb('FireBrick'),'MarkerSize',9,'Color',rgb('FireBrick'));
                    PrettyFigureFormat;
                    ylabel('qU below E0 (eV)');
                    xlabel('count rate (cps)');
                    legend('single runs', 'stack wmean',sprintf('stack %s',Mode)); legend boxoff 
                    grid on; 
                    set(gca,'FontSize',18); set(gca,'YScale','log');
                    hold off;
                    
                    savedir = [getenv('SamakPath'),'knm1ana/knm1_Stacking/plots/'];
                    MakeDir(savedir);
                    savename = [savedir,sprintf('StackqUCorr_%s_%.0f.png',obj.RunData.RunName,qu)];
                    print(f22,savename,'-dpng','-r450');
                    
                    xlim(StackRate(qu)+[-2, 2]*1e-03);
                    ylim(StackqU(qu)-18575+[-2,2]*1e-03);
                    savename = [savedir,sprintf('StackqUCorrZoom_%s_%.0f.png',obj.RunData.RunName,qu)];
                    print(f22,savename,'-dpng','-r450');        
            end
            
            obj.RunData.qU = StackqU;
            
            %qU = StackqU;
            %savename = [obj.RunData.matFilePath,obj.StackFileName,'.mat'];
            %save(savename,'qU','-append');
            
            obj.SimulateStackRuns;
        end
        function RunWiseErosFactor =  R200RateErosCorrection(obj)
            % Compute a Run Wise Correction Factor
            % R200 is the reference rate at -200 eV below E_0
            % Cancel Correlations between R200 and
            %  - Column Density
            %  - [TT]
            %  - [HT]
            %  - [DT]
            %  - qU200
            %
            % The corrected Rates is a linear combination written as:
            % Rcorrected = R -
            %              alpha . (rhoD - MeanrhoD)
            %              beta  . (TT - MeanTT)
            %              gamma . (HT - MeanHT)
            %              eta   . (DT - MeanDT)
            %              theta . (qU - MeanqU)
            %
            % All subruns from a single run undergo the same correction
            % Return a Run-wise Mutliplicative Correction factor
            %
            % Method Adapted from EROS2 experiment, T. Lasserre PhD, 2000
            % T. Lasserre, Last Modified 02/05/2019
            %
            
            % Reference Rate KNM1
%             try
%                 R200  = (sum(obj.SingleRunData.TBDIS(1,:,obj.PixList),3)./obj.SingleRunData.TimeperSubRun(1,:))';
%             catch
%                 R200  = (obj.SingleRunData.TBDIS(1,:)./obj.SingleRunData.TimeperSubRun(1,:))';
%             end
%             qU200 = obj.SingleRunData.qU(1,:,1);
            
            % Reference Rate KNM2
            R200  = sum(obj.SingleRunData.TBDIS_RM,1)'./(obj.SingleRunData.TimeSec.*mean(obj.SingleRunData.qUfrac_RM,1))';
            qU200 = mean(obj.SingleRunData.qU_RM,1);
            
            %% Definition of Parameters in use
            p1=obj.SingleRunData.WGTS_CD_MolPerCm2';
            p2=obj.SingleRunData.WGTS_MolFrac_TT';
            p3=obj.SingleRunData.WGTS_MolFrac_HT';
            p4=obj.SingleRunData.WGTS_MolFrac_DT';
            p5=qU200';
            
            TimeSec   = (obj.SingleRunData.TimeSec);
            m1        = obj.StackWmean(p1',TimeSec)';
            m2        = obj.StackWmean(p2',TimeSec)';
            m3        = obj.StackWmean(p3',TimeSec)';
            m4        = obj.StackWmean(p4',TimeSec)';
            m5        = obj.StackWmean(p5',TimeSec)';
            
            std1=std(p1);
            std2=std(p2);
            std3=std(p3);
            std4=std(p4);
            std5=std(p5);
            
            rhoR200p1=corr(R200,p1);
            rhoR200p2=corr(R200,p2);
            rhoR200p3=corr(R200,p3);
            rhoR200p4=corr(R200,p4);
            rhoR200p5=corr(R200,p5);
            
            rho12=corr(p1,p2);
            rho13=corr(p1,p3);
            rho23=corr(p2,p3);
            rho14=corr(p1,p4);
            rho24=corr(p2,p4);
            rho34=corr(p3,p4);
            rho15=corr(p1,p5);
            rho25=corr(p2,p5);
            rho35=corr(p3,p5);
            rho45=corr(p4,p5);
            
            %% Solving Equation for Cancelling Correlations Simultaneously
            M = [...
                std1         rho12*std2     rho13*std3    rho14*std4  rho15*std5      ; ...
                rho12*std1   std2           rho23*std3    rho24*std4  rho25*std5      ; ...
                rho13*std1   rho23*std2     std3          rho34*std4  rho35*std5      ; ...
                rho14*std1   rho24*std2     rho34*std3    std4        rho45*std5      ; ...
                rho15*std1   rho25*std2     rho35*std3    rho45*std4  std5            ; ...
                ];
            X = std(R200) * [rhoR200p1 rhoR200p2 rhoR200p3 rhoR200p4 rhoR200p5]';
            S = M\X;
            
            RunWiseErosFactor = ((R200  - S(1) * (p1-m1) - S(2) * (p2-m2) - S(3) * (p3-m3) - S(4) * (p4-m4) - S(5) * (p5-m5))./R200)';
            
        end
        
        function RunWiseErosFactor =  R200RateErosCorrectionNoqU(obj)
            % Compute a Run Wise Correction Factor
            % R200 is the reference rate at -200 eV below E_0
            % Cancel Correlations between R200 and
            %  - Column Density
            %  - [TT]
            %  - [HT]
            %  - [DT]
            %
            % The corrected Rates is a linear combination written as:
            % Rcorrected = R -
            %              alpha . (rhoD - MeanrhoD)
            %              beta  . (TT - MeanTT)
            %              gamma . (HT - MeanHT)
            %              eta   . (DT - MeanDT)
            %
            % All subruns from a single run undergo the same correction
            % Return a Run-wise Mutliplicative Correction factor
            %
            % Method Adapted from EROS2 experiment, T. Lasserre PhD, 2000
            % T. Lasserre, Last Modified 17/05/2019
            %
            
            % Reference Rate
            try
                R200  = (sum(obj.SingleRunData.TBDIS(1,:,obj.PixList),3)./obj.SingleRunData.TimeperSubRun(1,:))';
            catch
                R200  = (obj.SingleRunData.TBDIS(1,:)./obj.SingleRunData.TimeperSubRun(1,:))';
            end
            qU200 = obj.SingleRunData.qU(1,:,1);
            
            
            %% Definition of Parameters in use
            p1=obj.SingleRunData.WGTS_CD_MolPerCm2';
            p2=obj.SingleRunData.WGTS_MolFrac_TT';
            p3=obj.SingleRunData.WGTS_MolFrac_HT';
            p4=obj.SingleRunData.WGTS_MolFrac_DT';
            
            TimeSec   = sum(obj.SingleRunData.TimeSec(obj.SingleRunData.Select_all));
            m1        = obj.StackWmean(p1',TimeSec)';
            m2        = obj.StackWmean(p2',TimeSec)';
            m3        = obj.StackWmean(p3',TimeSec)';
            m4        = obj.StackWmean(p4',TimeSec)';
            
            std1=std(p1);
            std2=std(p2);
            std3=std(p3);
            std4=std(p4);
            
            rhoR200p1=corr(R200,p1);
            rhoR200p2=corr(R200,p2);
            rhoR200p3=corr(R200,p3);
            rhoR200p4=corr(R200,p4);
            
            rho12=corr(p1,p2);
            rho13=corr(p1,p3);
            rho23=corr(p2,p3);
            rho14=corr(p1,p4);
            rho24=corr(p2,p4);
            rho34=corr(p3,p4);
            
            %% Solving Equation for Cancelling Correlations Simultaneously
            M = [...
                std1         rho12*std2     rho13*std3    rho14*std4       ; ...
                rho12*std1   std2           rho23*std3    rho24*std4       ; ...
                rho13*std1   rho23*std2     std3          rho34*std4       ; ...
                rho14*std1   rho24*std2     rho34*std3    std4             ; ...
                ];
            X = std(R200) * [rhoR200p1 rhoR200p2 rhoR200p3 rhoR200p4]';
            S = M\X;
            
            RunWiseErosFactor = ((R200  - S(1) * (p1-m1) - S(2) * (p2-m2) - S(3) * (p3-m3) - S(4) * (p4-m4))./R200)';
            
        end
        
        function RunWiseErosFactor =  R200RateErosCorrectionLight(obj)
            % Compute a Run Wise Correction Factor
            % R200 is the reference rate at -200 eV below E_0
            % Cancel Correlations between R200 and
            %  - Column Density
            %  - qU200
            %
            % The corrected Rates is a linear combination written as:
            % Rcorrected = R -
            %              alpha . (rhoD - MeanrhoD)
            %              theta . (qU - MeanqU)
            %
            % All subruns from a single run undergo the same correction
            % Return a Run-wise Mutliplicative Correction factor
            %
            % Method Adapted from EROS2 experiment, T. Lasserre PhD, 2000
            % T. Lasserre, Last Modified 02/05/2019
            %
            
            % Reference Rate
            R200=(sum(obj.SingleRunData.TBDIS(1,:,obj.PixList),3)./obj.SingleRunData.TimeperSubRun(1,:))';
            % R200  = (obj.SingleRunData.TBDIS(1,:)./obj.SingleRunData.TimeperSubRun(1,:))';
            qU200 = obj.SingleRunData.qU(1,:,1);
            
            %% Definition of Parameters in use
            p1=obj.SingleRunData.WGTS_CD_MolPerCm2';
            p2=qU200';
            
            TimeSec   = sum(obj.SingleRunData.TimeSec(obj.SingleRunData.Select_all));
            m1        = obj.StackWmean(p1',TimeSec)';
            m2        = obj.StackWmean(p2',TimeSec)';
            
            std1=std(p1);
            std2=std(p2);
            
            rhoR200p1=corr(R200,p1);
            rhoR200p2=corr(R200,p2);
            
            rho12=corr(p1,p2);
            
            %% Solving Equation for Cancelling Correlations Simultaneously
            M = [...
                std1         rho12*std2    ; ...
                rho12*std1   std2          ; ...
                ];
            X = std(R200) * [rhoR200p1 rhoR200p2]';
            S = M\X;
            
            RunWiseErosFactor = ((R200  - S(1) * (p1-m1) - S(2) * (p2-m2))./R200)';
            
        end
        
        function RunWiseErosFactor =  R200RateErosCorrectionNoRhoD(obj)
            % Compute a Run Wise Correction Factor
            % R200 is the reference rate at -200 eV below E_0
            % Cancel Correlations between R200 and
            %  - qU
            %  - [TT]
            %  - [HT]
            %  - [DT]
            %
            % The corrected Rates is a linear combination written as:
            % Rcorrected = R -
            %              alpha . (qU - MeanqU)
            %              beta  . (TT - MeanTT)
            %              gamma . (HT - MeanHT)
            %              eta   . (DT - MeanDT)
            %
            % All subruns from a single run undergo the same correction
            % Return a Run-wise Mutliplicative Correction factor
            %
            % Method Adapted from EROS2 experiment, T. Lasserre PhD, 2000
            % T. Lasserre, Last Modified 4/7/2019
            %
            
            % Reference Rate KNM1
%             try
%                 R200  = (sum(obj.SingleRunData.TBDIS(1,:,obj.PixList),3)./obj.SingleRunData.TimeperSubRunperPixel(1,:,1))';
%             catch
%                 R200  = (obj.SingleRunData.TBDIS(1,:)./obj.SingleRunData.TimeperSubRunperPixel(1,:,1))';
%             end
%             qU200 = obj.SingleRunData.qU(1,:,1);
            % Reference Rate KNM2
            R200  = sum(obj.SingleRunData.TBDIS_RM,1)'./(obj.SingleRunData.TimeSec.*mean(obj.SingleRunData.qUfrac_RM,1))';
            qU200 = mean(obj.SingleRunData.qU_RM,1);
            
            %% Definition of Parameters in use
            p1=qU200';
            p2=obj.SingleRunData.WGTS_MolFrac_TT';
            p3=obj.SingleRunData.WGTS_MolFrac_HT';
            p4=obj.SingleRunData.WGTS_MolFrac_DT';
            
            TimeSec   = (obj.SingleRunData.TimeSec(obj.SingleRunData.Select_all));
            m1        = obj.StackWmean(p1',TimeSec)';
            m2        = obj.StackWmean(p2',TimeSec)';
            m3        = obj.StackWmean(p3',TimeSec)';
            m4        = obj.StackWmean(p4',TimeSec)';
            
            std1=std(p1);
            std2=std(p2);
            std3=std(p3);
            std4=std(p4);
            
            rhoR200p1=corr(R200,p1);
            rhoR200p2=corr(R200,p2);
            rhoR200p3=corr(R200,p3);
            rhoR200p4=corr(R200,p4);
            
            rho12=corr(p1,p2);
            rho13=corr(p1,p3);
            rho23=corr(p2,p3);
            rho14=corr(p1,p4);
            rho24=corr(p2,p4);
            rho34=corr(p3,p4);
            
            %% Solving Equation for Cancelling Correlations Simultaneously
            M = [...
                std1         rho12*std2     rho13*std3    rho14*std4       ; ...
                rho12*std1   std2           rho23*std3    rho24*std4       ; ...
                rho13*std1   rho23*std2     std3          rho34*std4       ; ...
                rho14*std1   rho24*std2     rho34*std3    std4             ; ...
                ];
            X = std(R200) * [rhoR200p1 rhoR200p2 rhoR200p3 rhoR200p4]';
            S = M\X;
            
            RunWiseErosFactor = ((R200  - S(1) * (p1-m1) - S(2) * (p2-m2) - S(3) * (p3-m3) - S(4) * (p4-m4))./R200)';
            
        end

        function RunWiseErosFactor =  R200RateErosCorrectionqUTpurity(obj)
            % Compute a Run Wise Correction Factor
            % R200 is the reference rate at -200 eV below E_0
            % Cancel Correlations between R200 and
            %  - Column Density
            %  - qU200
            %
            % The corrected Rates is a linear combination written as:
            % Rcorrected = R -
            %              alpha . (rhoD - MeanrhoD)
            %              theta . (qU - MeanqU)
            %
            % All subruns from a single run undergo the same correction
            % Return a Run-wise Mutliplicative Correction factor
            %
            % Method Adapted from EROS2 experiment, T. Lasserre PhD, 2000
            % T. Lasserre, Last Modified 02/05/2019
            %
            
            % % Reference Rate
            try
                R200  = (sum(obj.SingleRunData.TBDIS(1,:,obj.PixList),3)./obj.SingleRunData.TimeperSubRunperPixel(1,:,1))';
            catch
                R200  = (obj.SingleRunData.TBDIS(1,:)./obj.SingleRunData.TimeperSubRunperPixel(1,:,1))';
            end
            qU200 = obj.SingleRunData.qU(1,:,1);
            
            
            %% Definition of Parameters in use
            p1=obj.SingleRunData.WGTS_MolFrac_TT'+0.5*obj.SingleRunData.WGTS_MolFrac_HT'+0.5*obj.SingleRunData.WGTS_MolFrac_DT';
            p2=qU200';
            
            TimeSec   = (obj.SingleRunData.TimeSec(obj.SingleRunData.Select_all));
            m1        = obj.StackWmean(p1',TimeSec)';
            m2        = obj.StackWmean(p2',TimeSec)';
            
            std1=std(p1);
            std2=std(p2);
            
            rhoR200p1=corr(R200,p1);
            rhoR200p2=corr(R200,p2);
            
            rho12=corr(p1,p2);
            
            %% Solving Equation for Cancelling Correlations Simultaneously
            M = [...
                std1         rho12*std2    ; ...
                rho12*std1   std2          ; ...
                ];
            X = std(R200) * [rhoR200p1 rhoR200p2]';
            S = M\X;
            
            RunWiseErosFactor = ((R200  - S(1) * (p1-m1) - S(2) * (p2-m2))./R200)';
            
        end

        
        function RunWiseErosFactor =  RMRateErosCorrectionqUActivity(obj)
            % Compute a Run Wise Correction Factor
            % R200 is the reference rate at -200 eV below E_0
            % Cancel Correlations between R200 and
            %  - Purity x Column Density
            %  - qU200
            %
            % The corrected Rates is a linear combination written as:
            % Rcorrected = R -
            %              alpha . (Activity - MeanActivity)
            %              theta . (qU - MeanqU)
            %
            % All subruns from a single run undergo the same correction
            % Return a Run-wise Mutliplicative Correction factor
            %
            % Compatible with KNM2
            %
            %
            % Method Adapted from EROS2 experiment, T. Lasserre PhD, 2000
            % T. Lasserre, Last Modified 17/10/2019
            %
            
            % Reference Rate KNM2
            R200  = sum(obj.SingleRunData.TBDIS_RM,1)'./(obj.SingleRunData.TimeSec.*mean(obj.SingleRunData.qUfrac_RM,1))';
            qU200 = mean(obj.SingleRunData.qU_RM,1);
            
            
            %% Definition of Parameters in use
            p1=(obj.SingleRunData.WGTS_MolFrac_TT'+0.5*obj.SingleRunData.WGTS_MolFrac_HT'+0.5*obj.SingleRunData.WGTS_MolFrac_DT')./mean((obj.SingleRunData.WGTS_MolFrac_TT'+0.5*obj.SingleRunData.WGTS_MolFrac_HT'+0.5*obj.SingleRunData.WGTS_MolFrac_DT')).*obj.SingleRunData.WGTS_CD_MolPerCm2'./mean(obj.SingleRunData.WGTS_CD_MolPerCm2');
            p2=qU200'./mean(qU200);
            
            TimeSec   = (obj.SingleRunData.TimeSec);
            m1        = obj.StackWmean(p1',TimeSec)';
            m2        = obj.StackWmean(p2',TimeSec)';
            
            std1=std(p1);
            std2=std(p2);
            
            rhoR200p1=corr(R200,p1);
            rhoR200p2=corr(R200,p2);
            
            rho12=corr(p1,p2);
            
            %% Solving Equation for Cancelling Correlations Simultaneously
            M = [...
                std1         rho12*std2    ; ...
                rho12*std1   std2          ; ...
                ];
            X = std(R200) * [rhoR200p1 rhoR200p2]';
            S = M\X;
            
            RunWiseErosFactor = ((R200  - S(1) * (p1-m1) - S(2) * (p2-m2))./R200)';
            
                end
        
                
        function ROIPileUpDataCorrectionStackRuns(obj,run)
            %
            % Correction for
            % - ROI qU dependent Efficiency
            % - Pile-Up Pixel-rate dependent Efficiency
            %
            % +++++++++++++++++++++++++++++++++++++++++++++++
            % +             Analysis Warning                +
            % +++++++++++++++++++++++++++++++++++++++++++++++
            % DATA are corrected to match uncorrected Model
            % +++++++++++++++++++++++++++++++++++++++++++++++
            %
            % Thierry Lasserre
            % Last Modified: 2/8/2018
            %
            %
            switch obj.DataEffCorr
                case 'ROI'
                    RoiCorr = 1./obj.ModelObj.qUEfficiencyCorrectionFactor(obj.SingleRunData{run}.qU);
                    obj.SingleRunData{run}.TBDIS = obj.SingleRunData{run}.TBDIS .* RoiCorr;
                    %Plot(obj.SingleRunData.qU,RoiCorr);
                    fprintf(2,'Data: ROI Efficiency Correction Applied\n');
                case 'PileUp'
                    switch obj.ModelObj.FPD_Segmentation
                        case 'OFF'
                            nPix=124;
                        case {'SINGLEPIXEL','MULTIPIXEL'}
                            nPix=1;
                        case 'RING'
                            if obj.ModelObj.FPD_Ring == 1
                                nPix=4;
                            else
                                nPix=12;
                            end
                    end
                    PileUpCorr = 1./obj.ModelObj.PileUpEfficiencyCorrectionFactor(obj.SingleRunData{run}.TBDIS./obj.SingleRunData{run}.qUfrac./obj.SingleRunData{run}.TimeSec/nPix);
                    obj.SingleRunData{run}.TBDIS = obj.SingleRunData{run}.TBDIS .* PileUpCorr;
                    obj.SingleRunData{run}.TBDISE = sqrt(obj.SingleRunData{run}.TBDIS .* PileUpCorr);
                    %Plot(obj.SingleRunData.qU,PileUpCorr);
                    fprintf(2,'Data: Pile-Up Efficiency Correction Applied\n');
                case 'ROI+PileUp'
                    RoiCorr = 1./obj.ModelObj.qUEfficiencyCorrectionFactor(obj.SingleRunData{run}.qU);
                    switch obj.ModelObj.FPD_Segmentation
                        case 'OFF'
                            nPix=124;
                        case {'SINGLEPIXEL','MULTIPIXEL'}
                            nPix=1;
                        case 'RING'
                            if obj.ModelObj.FPD_Ring == 1
                                nPix=4;
                            else
                                nPix=12;
                            end
                    end
                    PileUpCorr = 1./obj.ModelObj.PileUpEfficiencyCorrectionFactor(obj.SingleRunData{run}.TBDIS./obj.SingleRunData{run}.qUfrac./obj.SingleRunData{run}.TimeSec/nPix);
                    obj.SingleRunData{run}.TBDIS  = obj.SingleRunData{run}.TBDIS .* RoiCorr .* PileUpCorr;
                    obj.SingleRunData{run}.TBDISE = sqrt(obj.SingleRunData{run}.TBDIS .* RoiCorr .* PileUpCorr);
                    %Plot(obj.SingleRunData.qU,RoiCorr.*PileUpCorr);
                    fprintf(2,'Data: ROI+Pile-Up Efficiency Corrections Applied\n');
            end
        end
        function rhoDLoadDataModel(obj)
            rhoDFileNames = {'668','763_764_765_766','926_927_928_929_930_931_932_933_934_935',...
                '794_795_796_797_798_799_800_801_802_803_804'};
            for mm = 1:length(obj.ModelObj.TBDCell)
                RunDataTemp = load(['Stack_',rhoDFileNames{mm},'ex2.mat']);
                obj.RunData.TBDIS(:,mm)  = RunDataTemp.TBDIS;
                obj.RunData.TBDISE(:,mm) = RunDataTemp.TBDISE;
                obj.RunData.qU(:,mm)     = RunDataTemp.qU;
                
            end
        end
        function RMCorrection(obj,varargin)
            %corrects the RM rate for activity and retarding potential
            p = inputParser;
            p.addParameter('ActivityCorrection','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('qUCorrection','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('pixlist','uniform',@(x)ismember(x,{'uniform','ring1','ring2','ring3','ring4'}));
            p.parse(varargin{:});
            ActivityCorrection = p.Results.ActivityCorrection;
            qUCorrection       = p.Results.ActivityCorrection;
            saveplot           = p.Results.saveplot;
            pixlist            = p.Results.pixlist;
            
            rhoD = obj.SingleRunData.WGTS_CD_MolPerCm2;
            T2 = obj.SingleRunData.WGTS_MolFrac_TT;
            DT = obj.SingleRunData.WGTS_MolFrac_DT;
            HT = obj.SingleRunData.WGTS_MolFrac_HT;
            Activity = rhoD.*(T2+0.5*DT+0.5*HT);
            meanActivity = mean(Activity)
            qU = obj.SingleRunData.qU_RM-mean(obj.SingleRunData.qU_RM);
            
            rate = obj.SingleRunData.TBDIS_RM./(obj.SingleRunData.qUfrac_RM.*obj.SingleRunData.TimeSec);
            correlation = corr(Activity(:),rate(:));
            fig88 = figure('Renderer','painters');
            set(fig88,'units','normalized','pos',[0.1, 0.1,1,0.6]);
            PlotStyle = { 'o','MarkerSize',8,'MarkerFaceColor',rgb('SkyBlue'),...%rgb('IndianRed'),...%
                'LineWidth',2,'Color',obj.PlotColor};
            e1 = errorbar(Activity,rate,repmat(0,[1,numel(Activity)]),PlotStyle{:});
            e1.CapSize = 0;
            e2 = errorbar(Activity,rate,repmat(std(rate),[1,numel(rate)]),PlotStyle{:});
            e2.CapSize = 0;
            xlabel(sprintf('Activity'));
            ylabel(sprintf('Rate (raw)'));
            leg = legend([e2],sprintf('Activity vs. rate'));
            title(leg,sprintf('Correlation factor = %.2f',correlation));
            leg.EdgeColor = rgb('Silver');
            leg.Location = 'best';
            PrettyFigureFormat;
            if strcmp(saveplot,'ON')
                SaveDir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/corrections/',obj.DataSet)];
                MakeDir(SaveDir);
                SaveName = sprintf('RateCorrection_raw_%s_%s.pdf',obj.RunData.RunName,pixlist);
                export_fig(fig88,[SaveDir,SaveName]);
            end
            
            if strcmp(ActivityCorrection,'ON')
                obj.SingleRunData.TBDIS_RM = obj.SingleRunData.TBDIS_RM.*(meanActivity./Activity);
                rate0 = rate;
                rate = obj.SingleRunData.TBDIS_RM./(obj.SingleRunData.qUfrac_RM.*obj.SingleRunData.TimeSec);
                DiffAct = abs(rate - rate0)/mean(rate0);
                fig88 = figure('Renderer','painters');
                set(fig88,'units','normalized','pos',[0.1, 0.1,1,0.6]);
                histogram(DiffAct)
                leg = legend(sprintf('\\sigma = %.2f x10^{-4}',std(DiffAct)*10.^4));
                title(leg,sprintf('Activity Correction percentage'));
                PrettyFigureFormat;
                if strcmp(saveplot,'ON')
                    SaveDir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/corrections/',obj.DataSet)];
                    MakeDir(SaveDir);
                    SaveName = sprintf('CorrectionPercentage_Activity_%s_%s.pdf',obj.RunData.RunName,pixlist);
                    export_fig(fig88,[SaveDir,SaveName]);
                end
                correlation = corr(Activity(:),rate(:));
                fig88 = figure('Renderer','painters');
                set(fig88,'units','normalized','pos',[0.1, 0.1,1,0.6]);
                PlotStyle = { 'o','MarkerSize',8,'MarkerFaceColor',rgb('SkyBlue'),...%rgb('IndianRed'),...%
                     'LineWidth',2,'Color',obj.PlotColor};
                e1 = errorbar(Activity,rate,repmat(std(Activity),[1,numel(Activity)]),PlotStyle{:});
                e1.CapSize = 0;
                e2 = errorbar(Activity,rate,repmat(std(rate),[1,numel(rate)]),PlotStyle{:});
                e2.CapSize = 0;
                xlabel(sprintf('Activity'));
                ylabel(sprintf('Rate (Activity corrected)'));
                leg = legend([e2],sprintf('Activity vs. rate'));
                title(leg,sprintf('Correlation factor = %.2f',correlation));
                leg.EdgeColor = rgb('Silver');
                leg.Location = 'best';
                PrettyFigureFormat;
                if strcmp(saveplot,'ON')
                    SaveDir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/corrections/',obj.DataSet)];
                    MakeDir(SaveDir);
                    SaveName = sprintf('RateCorrection_Activity_%s_%s.pdf',obj.RunData.RunName,pixlist);
                    export_fig(fig88,[SaveDir,SaveName]);
                end
                if strcmp(qUCorrection,'ON')
                    fitobject = fit(permute(qU,[2 1]),permute(rate,[2 1]),'poly1')
                    coeffs = coeffvalues(fitobject);
                    correlation = corr(qU(:),rate(:));
                    fig88 = figure('Renderer','painters');
                    set(fig88,'units','normalized','pos',[0.1, 0.1,1,0.6]);
                    plot(fitobject,qU,rate);
                    hold on;
                    PlotStyle = { 'o','MarkerSize',8,'MarkerFaceColor',rgb('SkyBlue'),...%rgb('IndianRed'),...%
                         'LineWidth',2,'Color',obj.PlotColor};
                    e1 = errorbar(qU,rate,repmat(std(qU),[1,numel(qU)]),PlotStyle{:});
                    e1.CapSize = 0;
                    e2 = errorbar(qU,rate,repmat(std(rate),[1,numel(rate)]),PlotStyle{:});
                    e2.CapSize = 0;
                    xlabel(sprintf('qU - \\langle qU \\rangle'));
                    ylabel(sprintf('Rate (Activity corrected)'));
                    leg = legend([e2],sprintf('qU vs. rate'));
                    title(leg,sprintf('Correlation factor = %.2f',correlation));
                    leg.EdgeColor = rgb('Silver');
                    leg.Location = 'best';
                    PrettyFigureFormat;
                    hold off;
                    if strcmp(saveplot,'ON')
                        SaveDir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/corrections/',obj.DataSet)];
                        MakeDir(SaveDir);
                        SaveName = sprintf('RateCorrection_qU_raw_%s_%s.pdf',obj.RunData.RunName,pixlist);
                        export_fig(fig88,[SaveDir,SaveName]);
                    end
                    slope = repmat(coeffs(1,1),[1,numel(obj.SingleRunData.TBDIS_RM)]).*qU;
                    rate1 = rate;
                    rate = rate-slope;
                    DiffqU = abs(rate - rate1)/mean(rate1);
                    fig88 = figure('Renderer','painters');
                    set(fig88,'units','normalized','pos',[0.1, 0.1,1,0.6]);
                    histogram(DiffqU)
                    leg = legend(sprintf('\\sigma = %.2f x10^{-4}',std(DiffqU)*10.^4));
                    title(leg,sprintf('qU Correction percentage'));
                    PrettyFigureFormat;
                    if strcmp(saveplot,'ON')
                        SaveDir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/corrections/',obj.DataSet)];
                        MakeDir(SaveDir);
                        SaveName = sprintf('CorrectionPercentage_qU_%s_%s.pdf',obj.RunData.RunName,pixlist);
                        export_fig(fig88,[SaveDir,SaveName]);
                    end
                    obj.SingleRunData.TBDIS_RM = rate.*(obj.SingleRunData.qUfrac_RM.*obj.SingleRunData.TimeSec);
                    correlation = corr(qU(:),rate(:));
                    fig88 = figure('Renderer','painters');
                    set(fig88,'units','normalized','pos',[0.1, 0.1,1,0.6]);
                    PlotStyle = { 'o','MarkerSize',8,'MarkerFaceColor',rgb('SkyBlue'),...%rgb('IndianRed'),...%
                         'LineWidth',2,'Color',obj.PlotColor};
                    e1 = errorbar(qU,rate,repmat(std(qU),[1,numel(qU)]),PlotStyle{:});
                    e1.CapSize = 0;
                    e2 = errorbar(qU,rate,repmat(std(rate),[1,numel(rate)]),PlotStyle{:});
                    e2.CapSize = 0;
                    xlabel(sprintf('qU - \\langle qU \\rangle'));
                    ylabel(sprintf('Rate (Activity & qU corrected)'));
                    leg = legend([e2],sprintf('qU vs. rate'));
                    title(leg,sprintf('Correlation factor = %.2f',correlation));
                    leg.EdgeColor = rgb('Silver');
                    leg.Location = 'best';
                    PrettyFigureFormat;
                    if strcmp(saveplot,'ON')
                        SaveDir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/corrections/',obj.DataSet)];
                        MakeDir(SaveDir);
                        SaveName = sprintf('RateCorrection_qU_%s_%s.pdf',obj.RunData.RunName,pixlist);
                        export_fig(fig88,[SaveDir,SaveName]);
                    end
                end
            end
        end
            
    end    % Data Import Methods end
    
    methods % Preparation for Fit Methods  
        function ReadSingleRunData(obj)
            %------------------------------------------------------------
            % Read the mat files from all runs in "RunList"
            % fill struct "SingleRunData"
            %------------------------------------------------------------
            % Init
            obj.SingleRunData = struct('Runs',obj.RunList);
            DataFiles = cell(obj.nRuns,1);
            
            % Load Files into temporary cell
            for i = 1:obj.nRuns
                if strcmp(obj.DataType,'Twin')
                    filename = [obj.RunData.matFilePath,'Twin',num2str(obj.RunList(i)),obj.SetTwinOrFakeFileName,'.mat'];
                elseif strcmp(obj.DataType,'Fake')
                    filename = [obj.RunData.matFilePath,obj.FakeStudyName,'_',num2str(obj.RunList(i)),'.mat'];
                elseif  strcmp(obj.DataType,'Real')
                    filename = [obj.RunData.matFilePath,num2str(obj.RunList(i)),'.mat'];
                elseif strcmp(obj.DataType,'FitriumTwin')
                    filename = [obj.RunData.matFilePath,'FitriumTwin',num2str(obj.RunList(i)),obj.SetTwinOrFakeFileName,'.mat'];   
                elseif strcmp(obj.DataType,'KafitTwin')
                    filename = [obj.RunData.matFilePath,'KafitTwin',num2str(obj.RunList(i)),obj.SetTwinOrFakeFileName,'.mat'];   
                end
                %disp(filename);
            
                if ~exist(filename,'file') && strcmp(obj.DataType,'Twin')                   
                    obj.ComputeTwinRun;
                elseif ~exist(filename,'file') && ismember(obj.DataType,{'FitriumTwin','KafitTwin'})
                   fprintf('%s Runs not avaiable / not found! \n',obj.DataType)
                elseif ~exist(filename,'file') && strcmp(obj.DataType,'Fake')  
                    fprintf('Fake Runs not avaiable / not found! \n')
                end
                
                DataFiles{i}   = load(filename);
               
            end
            RSvariables    = sort(fieldnames(DataFiles{1})); % variable names in the runsummary
             % put TBDIS_V to end of array
             if isfield(DataFiles{1},'TBDIS_V')
                 RSvariables{contains(RSvariables,'TBDIS_V')} = RSvariables{end-1};
                 tmp = RSvariables{end}; RSvariables{end-1} = tmp;
                 RSvariables{end} = 'TBDIS_V';
             end
               
             RSvariableSize =  cell2mat(cellfun(@(x) size(DataFiles{1}.(x)),RSvariables(~strcmp(RSvariables,'TBDIS_V')),'UniformOutput',0)); % size of every variable
             DimPar =  cell2mat(cellfun(@(x) sum(size(DataFiles{1}.(x))~=1),RSvariables(~strcmp(RSvariables,'TBDIS_V')),'UniformOutput',0)); % dimension of every variable
            % Fill structure with variables from mat files
            for i=1:numel(RSvariables) % Loop over variables
                if strcmp(RSvariables{i},'TBDIS_V')
                     tmpcell = cellfun(@(x) x.('TBDIS_V'),DataFiles,'UniFormOutput',0);
                     tmpSize = size(tmpcell{1});
                     partmp = zeros([tmpSize(1:2),obj.nRuns,tmpSize(3)]);
                     for r=1:obj.nRuns
                         partmp(:,:,r,:) = tmpcell{r};
                     end   
                     obj.SingleRunData.TBDIS_V = partmp;
                elseif any(RSvariableSize(i,:)==148) && DimPar(i)==2 %if multipixel variable and qU dependent
                    tmpcell = cellfun(@(x) x.(RSvariables{i}),DataFiles,'UniFormOutput',0); % cell with multipixel variable
                    partmp = zeros([RSvariableSize(i,1),obj.nRuns,RSvariableSize(i,2)]);
                    for r=1:obj.nRuns
                        partmp(:,r,:) = tmpcell{r};
                    end
                    obj.SingleRunData.(RSvariables{i}) = partmp;
                else % if variable is not pixel dependent
                    if ismember(RSvariables{i},{'matFileName','TBDarg'})
                         obj.SingleRunData.(RSvariables{i}) = ...
                            cellfun(@(x) x.(RSvariables{i}),DataFiles,'UniFormOutput',0);   
                    elseif ~contains(RSvariables{i},'RunTimeStart')
                        obj.SingleRunData.(RSvariables{i}) = ...
                            cell2mat(cellfun(@(x) x.(RSvariables{i})',DataFiles,'UniFormOutput',0))';
                    end
                end
                % add extra variable, which indicates if run passes criteria
                obj.SingleRunData.(['Select_',RSvariables{i}]) = true(1,obj.nRuns);
 
            end
            if numel(obj.PixList)>1
                obj.SingleRunData.TimeSec   = mean(obj.SingleRunData.TimeSec(obj.PixList,:));
            else
                obj.SingleRunData.TimeSec   = obj.SingleRunData.TimeSec(obj.PixList,:);
            end
            obj.SingleRunData.('Select_all') = true(1,obj.nRuns); % indicates if run is selected after passing all criterias
            if isfield(obj.SingleRunData,'StartTimeStamp')
                obj.SingleRunData.StartTimeStamp = datetime(obj.SingleRunData.StartTimeStamp, 'ConvertFrom', 'posixtime' );
            end
            obj.SingleRunData.TBDIS_Default = obj.SingleRunData.TBDIS;
            if isfield(obj.SingleRunData,'TBDIS_RM')
                obj.SingleRunData.TBDIS_RM_Default = obj.SingleRunData.TBDIS_RM;
            end
            
            obj.SetMosCorr;
        end
        function StackPixelSingleRun(obj)
            switch obj.AnaFlag
                case 'StackPixel'
                    obj.SingleRunData.qU          = mean(obj.SingleRunData.qU(:,:,obj.PixList),3);
                    obj.SingleRunData.qUfrac      = mean(obj.SingleRunData.qUfrac(:,:,obj.PixList),3);
                    obj.SingleRunData.EffCorr     = mean(obj.SingleRunData.EffCorr(:,:,obj.PixList),3);
                    obj.SingleRunData.TBDIS       = sum(obj.SingleRunData.TBDIS(:,:,obj.PixList),3);
                    obj.SingleRunData.TBDISE      = sqrt(obj.SingleRunData.TBDIS./obj.SingleRunData.EffCorr); % statstical uncertainty
                    obj.SingleRunData.MACE_Ba_T   = mean(obj.SingleRunData.MACE_Ba_T(obj.PixList,:));
                    obj.SingleRunData.MACE_Bmax_T = mean(obj.SingleRunData.MACE_Bmax_T(obj.PixList,:));
                    obj.SingleRunData.TBDIS_Default  = sum(obj.SingleRunData.TBDIS_Default(:,:,obj.PixList),3); % in case ROI is changed
                    
                    if isfield(obj.SingleRunData,'qU_RM')
                        obj.SingleRunData.qU_RM     = mean(obj.SingleRunData.qU_RM(obj.PixList,:));
                        obj.SingleRunData.qUfrac_RM = mean(obj.SingleRunData.qUfrac_RM(obj.PixList,:));
                        obj.SingleRunData.TBDIS_RM  = sum(obj.SingleRunData.TBDIS_RM(obj.PixList,:));
                    end
                    
                    if isfield(obj.SingleRunData,'qU200')
                        obj.SingleRunData.qU200     = mean(obj.SingleRunData.qU200(obj.PixList,:));
                        obj.SingleRunData.TBDIS200  = sum(obj.SingleRunData.TBDIS200(obj.PixList,:));
                    end
                    
                    if isfield(obj.SingleRunData,'TBDIS14keV')
                        obj.SingleRunData.TBDIS14keV    = sum(obj.SingleRunData.TBDIS14keV(:,:,obj.PixList),3);
                        obj.SingleRunData.TBDIS14keV_RM = sum(obj.SingleRunData.TBDIS14keV_RM(obj.PixList,:));
                        obj.SingleRunData.TBDIS_RM_Default  = sum(obj.SingleRunData.TBDIS_RM_Default(obj.PixList,:));
                    end
                    
                    if isfield(obj.RunData,'TBDIS_V')
                        obj.SinRegleRunData.TBDIS_V = squeeze(sum(obj.SingleRunData.TBDIS_V(:,:,:,obj.PixList),4));
                    end
                    
                case 'Ring'
                  dim = {size(obj.RunData.qU,1),numel(obj.SingleRunData.Runs),obj.nRings};
                   obj.SingleRunData.qU      = reshape(cell2mat(cellfun(@(x) mean(obj.SingleRunData.qU(:,:,x),3),obj.RingPixList,'UniformOutput',false)'),...
                       dim{:});
                   obj.SingleRunData.EffCorr = reshape(cell2mat(cellfun(@(x) mean(obj.SingleRunData.EffCorr(:,:,x),3),obj.RingPixList,'UniformOutput',false)'),...
                        dim{:});
                   obj.SingleRunData.TBDIS   = reshape(cell2mat(cellfun(@(x) sum(obj.SingleRunData.TBDIS(:,:,x),3),obj.RingPixList,'UniformOutput',false)'),...
                       dim{:});
                    obj.SingleRunData.TBDIS_Default   = reshape(cell2mat(cellfun(@(x) sum(obj.SingleRunData.TBDIS_Default(:,:,x),3),obj.RingPixList,'UniformOutput',false)'),...
                       dim{:});
                   obj.SingleRunData.TBDISE  = sqrt(obj.SingleRunData.TBDIS./obj.SingleRunData.EffCorr);
                   obj.SingleRunData.MACE_Ba_T = cell2mat(cellfun(@(x) mean(obj.SingleRunData.MACE_Ba_T(x,:)),...
                       obj.RingPixList,'UniformOutput',false))';
                   obj.SingleRunData.MACE_Bmax_T = cell2mat(cellfun(@(x) mean(obj.SingleRunData.MACE_Bmax_T(x,:)),...
                       obj.RingPixList,'UniformOutput',false))';
                   
                   if isfield(obj.SingleRunData,'TBDIS14keV')
                       obj.SingleRunData.TBDIS14keV  = reshape(cell2mat(cellfun(@(x) sum(obj.SingleRunData.TBDIS14keV(:,:,x),3),obj.RingPixList,'UniformOutput',false)'),...
                       dim{:});
                   end
                   
                   if strcmp(obj.RingMerge,'None')
                       % delete not used rings (otherwise problems with NaN)
                       obj.SingleRunData.qU(:,:,~ismember(1:13,obj.RingList)) = [];
                       obj.SingleRunData.EffCorr(:,:,~ismember(1:13,obj.RingList)) = [];
                       obj.SingleRunData.TBDIS(:,:,~ismember(1:13,obj.RingList)) = [];
                       obj.SingleRunData.TBDIS_Default(:,:,~ismember(1:13,obj.RingList)) = [];
                       obj.SingleRunData.TBDISE(:,:,~ismember(1:13,obj.RingList)) = [];
                       obj.SingleRunData.MACE_Ba_T(:,~ismember(1:13,obj.RingList)) = [];
                       obj.SingleRunData.MACE_Bmax_T(:,~ismember(1:13,obj.RingList)) = [];
                   end
            end
        end
        function LoadSingleRunObj(obj)
            % calculate 1 model for each run in RunList
            % fill models into SingleRunObj
            % replace StackedRuns by RunList - Thierry 13/2/19
            obj.SingleRunObj = cell(length(obj.RunList),1);

            tmpRunNr = obj.RunNr;       % temporary store
            tmpRunData = obj.RunData;

            [TTFSD,DTFSD,HTFSD] = obj.SetDefaultFSD;
            for r=1:length(obj.RunList)
                obj.RunNr = obj.RunList(r);
                obj.ReadData;
                obj.SingleRunObj{r} = ref_RunAnalysis(obj.RunData,...
                  'PixList',obj.PixList,'RingList',obj.RingList,...
                   'ISCS','Theory','recomputeRF','OFF','ELossFlag',obj.ELossFlag,'FPD_Segmentation','OFF',...
                   'DTFSD',DTFSD,'HTFSD',HTFSD,'TTFSD',TTFSD,'DopplerEffectFlag',obj.DopplerEffectFlag,'RadiativeFlag',obj.RadiativeFlag,'RingMerge',obj.RingMerge);

                obj.SingleRunObj{r}.ComputeTBDDS; obj.SingleRunObj{r}.ComputeTBDIS;    
            end
            
            obj.RunNr   = tmpRunNr;    % set back to initial value
            obj.RunData = tmpRunData;
        end      
        function SimulateStackRuns(obj,varargin)
            % calculate model for stacked data
            p=inputParser;
            p.addParameter('mNuSq_i',0,@(x)isfloat(x));
            p.addParameter('qU','',@(x)isfloat(x));
            p.addParameter('qUfrac','',@(x)isfloat(x));
            p.addParameter('TimeSec','',@(x)isfloat(x));
            p.addParameter('BKG_RateAllFPDSec','',@(x)isfloat(x));
            p.addParameter('WGTS_CD_MolPerCm2','',@(x)isfloat(x));
            p.addParameter('FSDflag',0,@(x)isfloat(x)); % 0 = nominal % 1-2 Thierry Test Systematics FSD
            p.addParameter('WGTS_B_T','',@(x)isfloat(x));
            p.addParameter('MACE_Bmax_T','',@(x)isfloat(x));
            p.addParameter('MACE_Ba_T','',@(x)isfloat(x));
            
            p.parse(varargin{:});
            mNuSq_i             = p.Results.mNuSq_i;
            qU                = p.Results.qU;
            qUfrac            = p.Results.qUfrac;
            TimeSec           = p.Results.TimeSec;
            BKG_RateAllFPDSec = p.Results.BKG_RateAllFPDSec;
            WGTS_CD_MolPerCm2 = p.Results.WGTS_CD_MolPerCm2;
            FSDflag           = p.Results.FSDflag;
            MACE_Bmax_T       = p.Results.MACE_Bmax_T;
            MACE_Ba_T         = p.Results.MACE_Ba_T;
            WGTS_B_T          = p.Results.WGTS_B_T;
            
            
            [TTFSD,DTFSD,HTFSD] = obj.SetDefaultFSD;
            
            % Test Thierry - begin
            if FSDflag==1
                DTFSD = 'HTFSD';
                HTFSD = 'SAENZ';
                TTFSD = 'SAENZ';
            end
            if FSDflag==2
                  DTFSD = 'HTFSD';
                  HTFSD = 'SAENZ';
                  TTFSD = 'DOSS';
              end
              % Test Thierry - end
             
            TBDarg = {obj.RunData,...
                'ISCS','Theory',...
                'recomputeRF','OFF',...
                'ELossFlag',obj.ELossFlag,...
                'PixList',obj.PixList,...
                'RingList',obj.RingList,...
                'DTFSD',DTFSD,'HTFSD',HTFSD,'TTFSD',TTFSD,...
                'DopplerEffectFlag',obj.DopplerEffectFlag,...
                'RadiativeFlag',obj.RadiativeFlag,...
                'RingMerge',obj.RingMerge...
                'mNuSq_i',mNuSq_i,...
                };
 
            if ~isempty(qU)
                TBDarg = {TBDarg{:},'qU',qU};
            end
            if ~isempty(qUfrac)
                TBDarg = {TBDarg{:},'qUfrac',qUfrac};
            end
            if ~isempty(TimeSec)
                TBDarg = {TBDarg{:},'TimeSec',TimeSec};
            end
            
            if ~isempty(BKG_RateAllFPDSec)
                TBDarg = {TBDarg{:},'BKG_RateAllFPDSec',BKG_RateAllFPDSec};
            end
            
            if ~isempty(WGTS_CD_MolPerCm2)
                TBDarg = {TBDarg{:},'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2};
            end
            if ~isempty(MACE_Bmax_T)
                TBDarg = {TBDarg{:},'MACE_Bmax_T',MACE_Bmax_T};
            end
            if ~isempty(MACE_Ba_T)
                TBDarg = {TBDarg{:},'MACE_Ba_T',MACE_Ba_T};
            end
            if ~isempty(WGTS_B_T)
                TBDarg = {TBDarg{:},'WGTS_B_T',WGTS_B_T};
            end
            
%             if strcmp(obj.MosCorrFlag,'ON') 
%                 % response function broadening
%                TBDarg = {TBDarg{:},'MACE_Sigma',std(obj.SingleRunData.qU,0,2)};
%             end
            
            switch obj.AnaFlag
                case 'StackPixel'
                    disp([obj.RunData.matFilePath,obj.StackFileName])
                    obj.ModelObj = ref_RunAnalysis(TBDarg{:},'FPD_Segmentation','OFF');
                case 'SinglePixel'
                    if length(obj.PixList) > 1; fprintf('Taking the first pixel in the pixel list \n'); end
                    obj.ModelObj = ref_RunAnalysis(TBDarg{:},'FPD_Segmentation','SINGLEPIXEL','nTeBinningFactor',20);
                case 'MultiPixel'
                    obj.ModelObj = ref_RunAnalysis(TBDarg{:},'FPD_Segmentation','MULTIPIXEL','nTeBinningFactor',20);
                case 'Ring'
                     obj.ModelObj = ref_RunAnalysis(TBDarg{:},'FPD_Segmentation','RING');
%                 case 'Ring'
%                     % Thierry - 5/4/19
%                     if isempty(obj.CurrentRing)
%                         obj.CurrentRing = 1;
%                     end
%                     obj.ModelObj = ref_RunAnalysis(obj.DataType,obj.RunData,[obj.RunData.matFilePath obj.StackFileName],obj.PixList,...
%                         'FPD_Segmentation','RING',...
%                         'FPD_Ring',obj.CurrentRing,...
%                         'nTeBinningFactor',30,...
%                         'ELossFlag',obj.ELossFlag);
% 
%                     if size(obj.RunData.TBDIS,2) > 1
%                         EffCorr = mean((obj.RunData.TBDISE(:,obj.ModelObj.ring{obj.CurrentRing})./ obj.RunData.TBDIS(:,obj.ModelObj.ring{obj.CurrentRing})).^2,2);
%                         obj.RunData.TBDIS = sum(obj.RunData.TBDIS(:,obj.ModelObj.ring{obj.CurrentRing}),2);
%                         obj.RunData.TBDISE = sqrt(obj.RunData.TBDIS.*EffCorr);
%                         obj.RunData.qU = mean(obj.RunData.qU(:,obj.ModelObj.ring{obj.CurrentRing}),2);
%                     elseif size(obj.RunData.qU,2) > 1
%                         obj.RunData.qU = mean(obj.RunData.qU(:,obj.ModelObj.ring{obj.CurrentRing}),2);
%                     end
            end
            obj.ModelObj.ComputeTBDDS;
            obj.ModelObj.ComputeTBDIS;
            
        end
    end % Preparation for Fit Methods end
    
    methods % Fit methods
        function FitAllRings(obj,varargin)
            % Ring Wise Fit of Data
            
            % Initilization and memory allocation
            if isempty (varargin)
                catsflag = 'OFF';
            else
                catsflag = varargin{1};
            end
            
            % Fit Results Storage
            obj.AllFitResults = cell(obj.nRings,1);
            par               = zeros(obj.nRings,12);
            err               = zeros(obj.nRings,12);
            chi2min           = zeros(obj.nRings,1);
            dof               = zeros(obj.nRings,1);
            errmat            = cell(obj.nRings,1);
            badrings          = false;

            % Loop on Rings
            for ri = 1:obj.nRings
                obj.RunData = load([obj.RunData.matFilePath,obj.StackFileName,'.mat']);
                obj.CurrentRing = obj.RingList(ri);
                obj.SimulateStackRuns();
                if ~strcmp(obj.chi2,'chi2Stat'); obj.ComputeCM(); end
                try
                    fprintf(2,'-------- MutliRunAnalysis: Fit Ring %0d ------------ \n',obj.CurrentRing)
                    obj.Fit('CATS',catsflag);
                    par(ri,:)     = obj.FitResult.par;
                    err(ri,:)     = obj.FitResult.err;
                    chi2min(ri,:) = obj.FitResult.chi2min;
                    errmat{ri}    = obj.FitResult.errmat;
                    dof(ri,:)     = obj.FitResult.dof;
                catch
                    par(ri,:)     = nan;
                    err(ri,:)     = nan;
                    chi2min(ri,:) = nan;
                    errmat{ri}    = nan;
                    dof           = nan;
                    badrings(obj.CurrentRing) = true;
                end
            end
            
            obj.RunData = load([obj.RunData.matFilePath,obj.StackFileName,'.mat']);
            obj.AllFitResults = {par,err,chi2min,errmat,dof};
            obj.RingList(badrings) = [];
        end % FitAllRings        
        function FitAllSinglePixels(obj)
            % Save PixList given by user to a variable, to be able to
            % modify "obj.PixList" for the fit loop
            PixListAllPixels = obj.PixList;
            nPixels = length(PixListAllPixels);
            
            % Initilization and memory allocation
            obj.AllFitResults = cell(nPixels,1);
            badpixels = false;
            par = zeros(nPixels,6);
            err = zeros(nPixels,6);t
            chi2min = zeros(nPixels,1);
            dof = zeros(nPixels,1);
            errmat = cell(nPixels,1);
            
            for p = 1:length(PixListAllPixels)
                obj.PixList = PixListAllPixels(p);
                fprintf('Fitting pixel %d \n',obj.PixList);
                obj.RunData = load([obj.RunData.matFilePath,obj.StackFileName,'mpix','.mat']);
                obj.RunData.TBDIS  = obj.RunData.TBDIS(:,obj.PixList);
                obj.RunData.TBDISE = obj.RunData.TBDISE(:,obj.PixList);
                obj.RunData.qU = obj.RunData.qU(:,obj.PixList);
                obj.SimulateStackRuns();
                if ~(strcmp(obj.chi2,'chi2Stat') || strcmp(obj.chi2,'chi2P')); obj.ComputeCM(); end
                
                %                 try
                obj.Fit();
                %                    obj.PlotFit();
                par(p,:) = obj.FitResult.par;
                err(p,:) = obj.FitResult.err;
                chi2min(p,:) = obj.FitResult.chi2min;
                errmat{p} = obj.FitResult.errmat;
                dof(p,:) = obj.FitResult.dof;
                %                 catch
                %                     par(p,:) = nan;
                %                     err(p,:) = nan;
                %                     chi2min(p,:) = nan;
                %                     errmat{p} = nan;
                %                     badpixels(obj.PixList) = true;
                %                 end
                
            end
            
            obj.AllFitResults = {par,err,chi2min,errmat,dof};
            PixListAllPixels(badpixels) = [];
            obj.PixList = PixListAllPixels;
            
        end    
    end
    
    
    methods % Run Wise Fit     
        function PlotSCdistributions(obj, varargin)
            % Plot Distribution of Slow control Parameters
            % Split : All/Stacked Runs
            p=inputParser;
            p.addParameter('SaveHisto','ON',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            SaveHisto   = p.Results.SaveHisto;
            myMainTitle=[sprintf('KATRIN Slow Control Parameters (%s)',obj.StackFileName)];
                    
            maintitle=myMainTitle;
            savefile=sprintf('plots/MRA_SC_%0.0fStackRuns.png',numel(obj.StackedRuns));
            fig1 = figure('Name','MRA Slow Control Parameter Distributions','NumberTitle','off','rend','painters','pos',[10 10 1400 2000]);
            a=annotation('textbox', [0 0.91 1 0.1], ...
                'String', maintitle, ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center');
            a.FontSize=24;a.FontWeight='bold';
            
            subplot(2,2,1)
            h1=histogram(obj.SingleRunData.WGTS_CD_MolPerCm2,10,'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.99,'LineWidth',2);
%             hold on
%             h2=histogram(obj.SingleRunData.WGTS_CD_MolPerCm2(obj.SingleRunData.Select_all>0),'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.5);
%             h2.BinWidth = h1.BinWidth;
%             hold off
            ylabel('runs');
            xlabel('\rho.d mol/cm^2');
%            leg = legend([h1,h2],sprintf('All: %0.0f',numel(obj.RunList)),sprintf('Stacked: %0.0f',numel(obj.StackedRuns)),'Location','northwest');
%             leg.Color = 'none'; legend boxoff;
            sp1 = sprintf('mean = %0.3g mol/cm^2 - \\sigma = %0.3g mol/cm^2',...
                mean(obj.SingleRunData.WGTS_CD_MolPerCm2),std(obj.SingleRunData.WGTS_CD_MolPerCm2));title(sp1)
            PrettyFigureFormat;set(gca,'Fontsize',18);
            subplot(2,2,2)
            h1=histogram(obj.SingleRunData.WGTS_MolFrac_TT,10,'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.99,'LineWidth',2);
%             hold on
%             h2=histogram(obj.SingleRunData.WGTS_MolFrac_TT(obj.SingleRunData.Select_all>0),'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.5);
%             h2.BinWidth = h1.BinWidth;
%             hold off
            ylabel('runs');
            xlabel('[TT]');
%             leg = legend([h1,h2],'All','Stacked','Location','northwest');
%             leg.Color = 'none'; legend boxoff;
            sp1 = sprintf('mean = %0.3g - \\sigma = %0.3g ',...
                mean(obj.SingleRunData.WGTS_MolFrac_TT),std(obj.SingleRunData.WGTS_MolFrac_TT));title(sp1)
            PrettyFigureFormat;set(gca,'Fontsize',18);
            subplot(2,2,3)
            h1=histogram(obj.SingleRunData.WGTS_MolFrac_HT,10,'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.99,'LineWidth',2);
%             hold on
%             h2=histogram(obj.SingleRunData.WGTS_MolFrac_HT(obj.SingleRunData.Select_all>0),'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.5);
%             h2.BinWidth = h1.BinWidth;
%             hold off
            ylabel('runs');
            xlabel('[HT]');
%                leg = legend([h1,h2],'All','Stacked','Location','northwest');
%             leg.Color = 'none'; legend boxoff;
            sp1 = sprintf('mean = %0.3g - \\sigma = %0.3g ',...
                mean(obj.SingleRunData.WGTS_MolFrac_HT),std(obj.SingleRunData.WGTS_MolFrac_HT));title(sp1)
            PrettyFigureFormat;set(gca,'Fontsize',18);
            subplot(2,2,4)
            h1=histogram(obj.SingleRunData.WGTS_MolFrac_DT,10,'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.99,'LineWidth',2);
%             hold on
%             h2=histogram(obj.SingleRunData.WGTS_MolFrac_DT(obj.SingleRunData.Select_all>0),'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.5);
%             h2.BinWidth = h1.BinWidth;
%             hold off
            ylabel('runs');
            xlabel('[DT]');
%                leg = legend([h1,h2],'All','Stacked','Location','northwest');
%             leg.Color = 'none'; legend boxoff;
            sp1 = sprintf('mean = %0.3g - \\sigma = %0.3g ',...
                mean(obj.SingleRunData.WGTS_MolFrac_DT),std(obj.SingleRunData.WGTS_MolFrac_DT));title(sp1)
            PrettyFigureFormat; set(gca,'Fontsize',18);
            switch SaveHisto
                case 'ON'
            %publish_figurePDF(gcf,savefile);           
            export_fig(gcf,savefile,'-q101','-m3');
            end

        end
        function FitResults = PlotFitRunList(obj,varargin)
            % work in progress, will be simplified later
            % Overview: Fit of all Runs
            % Overview of Fit Parameter + some slow control
            p=inputParser;
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'OFF','ON','pdf','png'}));
            p.addParameter('Parameter','E0',@(x)ismember(x,{'mNuSq','E0','N','B','pVal','RhoD','T2'}));
            p.addParameter('DisplayStyle','Abs',@(x)ismember(x,{'Rel','Abs'}));
            p.addParameter('YLim','',@(x)isfloat(x) || isempty(x));
            p.addParameter('DispHist','ON',@(x)ismember(x,{'OFF','ON','Separate'}));
            p.parse(varargin{:});
            saveplot      = p.Results.saveplot;
            Parameter     = p.Results.Parameter;
            DisplayStyle  = p.Results.DisplayStyle;
            YLim          = p.Results.YLim; % force ylimits
            DispHist      = p.Results.DispHist;
            
            LocalFontSize = 28;
            
            % Get Fit Results
            switch obj.chi2
                case 'chi2Stat'
                    if isempty(obj.SingleRun_FitResults.chi2Stat)
                        obj.FitRunList('Recompute','ON')
                    end
                    FitResults = obj.SingleRun_FitResults.chi2Stat;
                    chi2leg = 'statistics only';
                case {'chi2CM','chi2CMShape'}
                    if isempty(obj.SingleRun_FitResults.chi2CMall) || isempty(obj.SingleRun_fitResults.chi2CMcorr)
                        obj.FitRunList('Recompute','ON')
                    end
                    FitResults = obj.SingleRun_FitResults.chi2CMall;
                    chi2leg = 'stat. + syst.';
            end
            
            % prepare x-axis: time
            LiveTime = hours(obj.SingleRunData.StartTimeStamp-obj.SingleRunData.StartTimeStamp(1));
             % get rid of large jumps
             TimeBreakIndex = find(diff(LiveTime)>5*mean(diff(LiveTime)));
             x = zeros(numel(LiveTime),1);
             if isempty(TimeBreakIndex)
                 NxParts=1; TimeBreakIndex = numel(obj.RunList);
             else
                 NxParts = numel(TimeBreakIndex)+1;
             end
             for i=1:NxParts
                if i==1 %start
                    startI=1;
                    stopI = TimeBreakIndex(i);
                    t = 0;
                elseif i==NxParts %stop
                    startI= TimeBreakIndex(i-1)+1;
                    stopI = numel(LiveTime);
                    t = t+ abs(LiveTime(TimeBreakIndex(i-1))-LiveTime(TimeBreakIndex(i-1)+1));
                else
                    startI= TimeBreakIndex(i-1)+1;
                    stopI = TimeBreakIndex(i);
                  t = t + abs(LiveTime(TimeBreakIndex(i-1))-LiveTime(TimeBreakIndex(i-1)+1));
                end       
                x(startI:stopI) = LiveTime(startI:stopI)-0.75*t;
            end
            
            %% prepare y axis: depending on defined Parameter  
            switch Parameter
                case 'E0'
                    y = FitResults.E0;
                    yErr = FitResults.E0Err;
                    if strcmp(DisplayStyle,'Abs')
                        ystr = sprintf('{\\itE}_0^{fit} (eV)');
                    elseif strcmp(DisplayStyle,'Rel')
                        y = y-wmean(y,1./yErr.^2);
                        ystr = sprintf('{\\itE}_0^{fit} - \\langle{\\itE}_0^{fit}\\rangle (eV)');
                    end   
                case 'mNuSq'
                    y = FitResults.mnuSq;
                    yErr = FitResults.mnuSqErr;
                    if strcmp(DisplayStyle,'Abs')
                        ystr = sprintf('{\\itm}_{\\beta}^2 (eV^2)');
                    elseif strcmp(DisplayStyle,'Rel')
                        y = y-wmean(y,1./yErr.^2);
                        ystr = sprintf('{\\itm}_{\\beta}^2 - < m_{\\beta}^2 >  (eV^2)');
                    end
                case 'B'
                    y = FitResults.B.*1e3;
                    yErr = FitResults.BErr.*1e3;
                    if strcmp(DisplayStyle,'Abs')
                        ystr = sprintf('Background (mcps)');
                    elseif strcmp(DisplayStyle,'Rel')
                        y = y-wmean(y,1./yErr.^2);
                        ystr = sprintf('{\\itB} - \\langle{\\itB}\\rangle  (mcps)');
                    end
                case 'N'
                    y = FitResults.N+1;
                    yErr = FitResults.NErr;
                    if strcmp(DisplayStyle,'Abs')
                        ystr = sprintf('{\\itN}');
                    elseif strcmp(DisplayStyle,'Rel')
                        y = y-wmean(y,1./yErr.^2);
                        ystr = sprintf('{\\itN} - \\langle{\\itN}\\rangle');
                    end
                case 'pVal'
                    y = FitResults.pValue;
                    yErr = zeros(numel(y),1);
                    ystr = sprintf('p-value (%.0f dof)',FitResults.dof(1));
                case 'RhoD'
                    y = obj.SingleRunData.WGTS_CD_MolPerCm2;
                    yErr = std(obj.SingleRunData.WGTS_CD_MolPerCm2_SubRun);
                    if strcmp(DisplayStyle,'Abs')
                        ystr = sprintf('Column density (mol \\cdot cm^{-2})');
                    elseif strcmp(DisplayStyle,'Rel')
                        y = y-wmean(y,obj.SingleRunData.TimeSec);
                        ystr = sprintf('\\rhod - \\langle\\rhod\\rangle (mol \\cdot cm^{-2})');
                    end
                case 'T2'
                    y    = obj.SingleRunData.WGTS_MolFrac_TT;
                    yErr = zeros(numel(y),1);%std(obj.SingleRunData.WGTS_MolFrac_TT_SubRun);
                    ystr = sprintf('T_2 fraction');
            end
            
            %% plot
           fig88 = figure('Renderer','painters');
            if strcmp(DispHist,'ON')
                set(fig88,'units','normalized','pos',[0.1, 0.1,1.0,0.5]); 
                subplot(1,6,1:5);
            else
                 set(fig88,'units','normalized','pos',[0.1, 0.1,1.0,0.5]);
            end
            
            PlotStyle = { 'o','MarkerSize',8,'MarkerFaceColor',rgb('SkyBlue'),...%rgb('IndianRed'),...%
                'LineWidth',2,'Color',obj.PlotColor};
            % plot data
            xmin=min(x)-max(x)*0.02;
            xmax = max(x)*1.02;
            if strcmp(Parameter,'pVal')
                %l1 = plot(linspace(xmin,xmax,10),0.05*ones(10,1),'--','Color',rgb('Orange'),'LineWidth',4);
            elseif ismember(Parameter,{'T2'})
                l1 = plot(linspace(xmin,xmax,10),wmean(y(y>0),obj.SingleRunData.TimeSec(y>0))*ones(10,1),'-','Color',rgb('DimGray'),'LineWidth',2);
            elseif ismember(Parameter,{'RhoD'})
                if strcmp(DisplayStyle,'Abs')
                    l1 = plot(linspace(xmin,xmax,10),wmean(y,obj.SingleRunData.TimeSec)*ones(10,1),'-','Color',rgb('DimGray'),'LineWidth',2);
                else
                    l1 = plot(linspace(xmin,xmax,10),zeros(10,1),'-','Color',rgb('DimGray'),'LineWidth',2);
                end
            else
                l1 = plot(linspace(xmin,xmax,10),wmean(y,1./yErr.^2)*ones(10,1),'-','Color',rgb('DimGray'),'LineWidth',2);
            end
            hold on;
            e1 = errorbar(x,y,yErr,PlotStyle{:});
            e1.CapSize = 0;
            hold off;
            xlabel('Time (hours)');
            ylabel(ystr);
            PrettyFigureFormat('FontSize',LocalFontSize);
            
            % y-axis:
             ax = gca;
            if isempty(YLim)
                myylim=get(gca,'YLim'); % get current ylimits
                yrange = abs(ax.YTick(1)-ax.YTick(end))*0.04;
                ylim([myylim(1)-yrange,myylim(2)+yrange]);
                if strcmp(Parameter,'E0') && strcmp(DisplayStyle,'Abs')
                    yticklabels(string(round(get(gca,'YTick'),5)));
                elseif strcmp(Parameter,'pVal')
                    ylim([0,1.25])
                end
            else
                ylim([min(YLim),max(YLim)])
            end
            
            %x-axis: change xtick labels to live time, keep old ticks values 
             xlim([xmin,xmax])
            myxticks = xticks;
            myxtickslabels = cell(numel(myxticks),1);
            for i=1:numel(myxticks)
                Index = find(LiveTime>=myxticks(i),1);
                TimeDiff = LiveTime(Index)-x(Index);
               myxtickslabels{i} = round(myxticks(i)+TimeDiff,0);
            end
            xticklabels(myxtickslabels);
            % plot time break axis sybol
            myylim=get(gca,'YLim'); % get current ylimits
            if strcmp(Parameter,'B') && strcmp(DisplayStyle,'Abs')
                factor = 0.03;
            else
                factor = 0.08;
            end
            yOffset = abs(ax.YTick(1)-ax.YTick(2))*factor;
            for i=1:numel(TimeBreakIndex)
                if TimeBreakIndex(i)==numel(x) % do not plot for last run
                else
                    ht = text(x(TimeBreakIndex(i))+2,myylim(1)+yOffset,sprintf('\\mid\\mid'));
                    set(ht,'Rotation',-25)
                    set(ht,'FontSize',22)
                    set(ht,'FontWeight','bold');
                end
            end
            
            % title
            if ismember(Parameter,{'mNuSq','E0','N','B'})
                % calculate p-value with respect to straight line
                dof = numel(y)-1;
                chi2 = sum((y-wmean(y,1./yErr.^2)).^2./yErr.^2);
                pval = 1-chi2cdf(chi2,dof);
                %title(sprintf('\\chi2 = %.1f (%.0f dof) , p-value = %.2f',chi2,dof,pval))
            end
            
            % legend
            if strcmp(Parameter,'pVal')
                leg = legend([e1],sprintf('Scanwise fits (%s)',chi2leg));%,'p-value = 0.05'); %l1
            elseif ismember(Parameter,{'RhoD'})
                ytmp = obj.SingleRunData.WGTS_CD_MolPerCm2;
                leg = legend([e1,l1],'Scanwise',sprintf('Time weighted mean = %.2g mol\\cdotcm^{-2}',wmean(ytmp,obj.SingleRunData.TimeSec)));
            elseif ismember(Parameter,{'T2'})
                leg = legend([e1,l1],'Scanwise',sprintf('Time weighted mean = %.3f',wmean(y(y>0),obj.SingleRunData.TimeSec(y>0))));
            else
                leg = legend([e1,l1],sprintf('Scanwise fits (%s)',chi2leg),sprintf('Weighted mean, p-value = %.2f',pval));
            end
            legend boxoff;
            if ismember(Parameter,{'RhoD','T2'})
                leg.Location = 'northwest';
            else
                leg.Location = 'north';
                legpos = leg.Position;
                leg.Position = [legpos(1)-0.03,legpos(2)-0.02,legpos(3:4)];
                
            end
            leg.FontSize = get(gca,'FontSize');
            
            % remove white space around figure
            ax = gca;
            outerpos = ax.OuterPosition;
            ti = ax.TightInset;
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1);% - ti(3);
            ax_height = outerpos(4);% - ti(2) - ti(4);
            %ax.Position = [left bottom ax_width ax_height];
            %ax.Position = [left-0.02 bottom+0.01 ax_width+0.02 ax_height-0.2];
            ax.Position = [left-0.02 bottom+0.05 ax_width+0.02 ax_height-0.3];
            if strcmp(DispHist,'ON')
                myYlim = ylim;
                
               subplot(1,6,6);
                
                PlotStyle = { 'FaceColor',obj.PlotColor,...
                    'FaceAlpha',1};
               
                h1 = histogram(y,PlotStyle{:});
                h1.Orientation='horizontal';
                ylim(myYlim);
                if ismember(Parameter,{'T2'})
                    h1.BinWidth = 0.001;
                elseif ismember(Parameter,{'RhoD'})
                    if strcmp(DisplayStyle,'Rel')
                        h1.BinWidth = 0.3*1e15;
                    end
                elseif ~ismember(Parameter,{'T2','RhoD'})
                    hold on;
                    if strcmp(Parameter,'pVal')
                        h1.BinWidth = 0.04;
                        %                        x=linspace(min(y)*0.75,max(y)*1.2,100);
                        %                        yhist = chi2pdf(x,FitResults.dof(1));
                    else
                        if strcmp(DisplayStyle,'Rel')
                            x=linspace(min(y)*1.2,max(y)*1.2,100);
                        else
                            x=linspace(min(y)*0.8,max(y)*1.2,100);
                        end
                        yhist = gaussian(x,wmean(y,1./yErr.^2),std(y))./simpsons(x,gaussian(x,wmean(y,1./yErr.^2),std(y)));
                        p1 = plot(yhist.*h1.BinWidth.*numel(obj.RunList),x,...
                            '-','LineWidth',4,'Color',rgb('GoldenRod'));
                    end
                end
                PrettyFigureFormat;
                
                % get rid of box and x axis
                box off
                set(get(gca,'XAxis'),'Visible','off')
                set(get(gca,'YAxis'),'Visible','off')
                % adjust position
                ax2 = gca;
                %ax.Position = [left-0.02 bottom+0.01 ax_width+0.128 ax_height-0.2];
                ax.Position = [left-0.02 bottom+0.05 ax_width+0.1343 ax_height-0.3];
                ax2.Position(1) = 0.80+0.09;
                if strcmp(Parameter,'pVal')
                  ax2.Position(2)=0.158;
                 ax2.Position(4)=0.75;
                elseif strcmp(Parameter,'RhoD')
                     ax2.Position(2)=0.17;
                     ax2.Position(4) = 0.74;
                end
            end
           
            % save plot
            if ~strcmp(saveplot,'OFF')
  
                savedir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/',obj.DataSet)];
                if ~exist(savedir,'dir')
                    system(['mkdir -p ',savedir]);
                end
                savename = [savedir,sprintf('FitRunList_%s_%.0feVbE0_%s_%s%s_Hist%s',...
                    obj.RunData.RunName,obj.GetRange,obj.chi2,Parameter,DisplayStyle,DispHist)];
                switch saveplot
                    case {'ON','pdf'}
                      %  publish_figurePDF(fig88,[savename,'.pdf']);
                        export_fig(fig88,[savename,'.pdf']);
                    case 'png'
                        print(fig88,[savename,'.png'],'-dpng','-r100');
                end
                fprintf('save plot as  %s \n',savename)
            end
            
            %% histogram
            if ~strcmp(DispHist,'Separate')
                return;
            end
            
            if strcmp(Parameter,'pVal')
                y = FitResults.chi2min;
                ystr = sprintf('\\chi2');
            end
            
            fig89 = figure('Renderer','opengl');
            set(fig89,'units','normalized','pos',[0.1, 0.1,0.5,0.6]);
            PlotStyle = { 'FaceColor', rgb('DodgerBlue'),...%obj.PlotColor,...
                'FaceAlpha',0.9};
            h1 = histogram(y,PlotStyle{:});
            if ~strcmp(Parameter,{'T2','RhoD'})
                hold on;
                 if strcmp(DisplayStyle,'Rel') && strcmp(Parameter,'E0')        
                     h1.BinWidth = 0.05;
                 end
                if strcmp(Parameter,'pVal')
                    x=linspace(min(y)*0.75,max(y)*1.2,100);
                    yhist = chi2pdf(x,FitResults.dof(1));
                elseif strcmp(DisplayStyle,'Rel')
                    x=linspace(min(y)*1.2,max(y)*1.2,100);
                    dist = fitdist(y,'Normal');
                    yhist = gaussian(x,dist.mu,dist.std)./simpsons(x,gaussian(x,wmean(y,1./yErr.^2),std(y)));
                else
                    dist = fitdist(y,'Normal');
                    x=linspace(min(y),max(y),1000);
                    yhist = gaussian(x,dist.mu,dist.std)./simpsons(x,gaussian(x,wmean(y,1./yErr.^2),std(y)));
                end
                
                p1 = plot(x,yhist.*h1.BinWidth.*numel(obj.RunList),...
                    '-','LineWidth',4,'Color',rgb('GoldenRod'));
                xlim([min(x),max(x)]);
                
 
                if strcmp(Parameter,'pVal')
                    leg = legend([h1,p1],sprintf('Runwise fits (%s)',chi2leg),sprintf('\\chi2 (%.0f dof)',FitResults.dof(1)));
                elseif  strcmp(Parameter,'E0')
                     leg = legend([h1,p1],sprintf('Runwise fits (%s)',chi2leg),sprintf('Gaussian \\sigma = %.3g eV',std(y)));
                else
                    leg = legend([h1,p1],sprintf('Runwise fits (%s)',chi2leg),'Gaussian');
                end
                legend boxoff;
                leg.Location = 'northwest';
                
            end
            PrettyFigureFormat('FontSize',24);
            xlabel(ystr);
            ylabel('Occurrence');
            if ~strcmp(saveplot,'OFF')
                % remove white space around figure
                ax = gca;
                outerpos = ax.OuterPosition;
                ti = ax.TightInset;
                left = outerpos(1) + ti(1);
                bottom = outerpos(2) + ti(2);
                ax_width = outerpos(3) - ti(1) - ti(3);
                ax_height = outerpos(4) - ti(2) - ti(4);
                ax.Position = [left bottom ax_width ax_height];
                
                savename = [savedir,sprintf('FitRunListHist_%s_%.0feVbE0_%s_%s%s',...
                    obj.RunData.RunName,obj.GetRange,obj.chi2,Parameter,DisplayStyle)];
                switch saveplot
                    case {'ON','pdf'}
                        publish_figurePDF(fig89,[savename,'.pdf']);
                    case 'png'
                        print(fig89,[savename,'.png'],'-dpng','-r450');
                end
                fprintf('save plot as  %s \n',savename)
            end
        end
        function FitResults = PlotFitRunListCorr(obj, varargin)
            %correlations between Fit parameters and Slow Control
            p=inputParser;
            p.addParameter('Parametery','E0',@(x)ismember(x,{'mNuSq','E0','N','B','pVal','rate300','qU_RM'}));
            p.addParameter('Parameterx','RhoD',@(x)ismember(x,{'RhoD','Act','mNuSq','E0','N','B','pVal','qU_RM','time'}));
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'OFF','ON','pdf','png'}));
            p.addParameter('DisplayStyle','Abs',@(x)ismember(x,{'Rel','Abs'}));
            p.addParameter('YLim','',@(x)isfloat(x) || isempty(x));
            p.addParameter('Fit','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Detrend','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('pixlist','uniform',@(x)ismember(x,{'uniform','ring1','ring2','ring3','ring4'}));
            p.parse(varargin{:});
            Parametery    = p.Results.Parametery;
            Parameterx    = p.Results.Parameterx;            
            saveplot      = p.Results.saveplot;
            DisplayStyle  = p.Results.DisplayStyle;
            Fit           = p.Results.Fit;
            Detrend       = p.Results.Detrend;
            pixlist       = p.Results.pixlist;
            
            LocalFontSize = 24;
            
            % Get Fit Results if fit results are plotted
            if ismember(Parametery,{'mNuSq','E0','N','B','pVal'}) || ismember(Parameterx,{'mNuSq','E0','N','B','pVal'})
                switch obj.chi2
                    case 'chi2Stat'
                        if isempty(obj.SingleRun_FitResults.chi2Stat)
                            obj.FitRunList('Recompute','ON')
                        end
                        FitResults = obj.SingleRun_FitResults.chi2Stat;
                        chi2leg = 'fits (statistics only)';
                    case {'chi2CM','chi2CMShape'}
                        if isempty(obj.SingleRun_FitResults.chi2CMall) || isempty(obj.SingleRun_fitResults.chi2CMcorr)
                            obj.FitRunList('Recompute','ON')
                        end
                        FitResults = obj.SingleRun_FitResults.chi2CMall;
                        chi2leg = 'fits (stat. + syst.)';
                end
            else
                chi2leg = '';
            end
            %% prepare y axis: depending on defined Parametery  
            switch Parametery
                case 'E0'
                    UnitStr = 'eV';
                    y = FitResults.E0';
                    yErr = FitResults.E0Err';
                    if strcmp(DisplayStyle,'Abs')
                        ystr = sprintf('{\\itE}_0^{fit} (eV)');
                    elseif strcmp(DisplayStyle,'Rel')
                        y = y-wmean(y,1./yErr.^2);
                        ystr = sprintf('{\\itE}_0^{fit} - \\langle{\\itE}_0^{fit}\\rangle (eV)');
                    end   
                case 'mNuSq'
                    UnitStr = sprintf('eV^2');
                    y = FitResults.mnuSq;
                    yErr = FitResults.mnuSqErr;
                    if strcmp(DisplayStyle,'Abs')
                        ystr = sprintf('{\\itm}_{\\beta}^2 (eV^2)');
                    elseif strcmp(DisplayStyle,'Rel')
                        y = y-wmean(y,1./yErr.^2);
                        ystr = sprintf('{\\itm}_{\\beta}^2 - < m_{\\beta}^2 >  (eV^2)');
                    end
                case 'B'
                     UnitStr = 'cps';
                    y = FitResults.B.*1e3;
                    yErr = FitResults.BErr.*1e3;
                    if strcmp(DisplayStyle,'Abs')
                        ystr = sprintf('Background (mcps)');
                    elseif strcmp(DisplayStyle,'Rel')
                        y = y-wmean(y,1./yErr.^2);
                        ystr = sprintf('{\\itB} - \\langle{\\itB}\\rangle  (mcps)');
                    end
                case 'N'
                    UnitStr = '';
                    y = FitResults.N+1;
                    yErr = FitResults.NErr;
                    if strcmp(DisplayStyle,'Abs')
                        ystr = sprintf('{\\itN}');
                    elseif strcmp(DisplayStyle,'Rel')
                        y = y-wmean(y,1./yErr.^2);
                        ystr = sprintf('{\\itN} - \\langle{\\itN}\\rangle');
                    end
                case 'pVal'
                    UnitStr = '';
                    y = FitResults.pValue;
                    yErr = zeros(numel(y),1);
                    ystr = sprintf('p-value (%.0f dof)',FitResults.dof(1));
                case 'rate300'
                     UnitStr = 'cps';
                    y = obj.SingleRunData.TBDIS_RM./(obj.SingleRunData.qUfrac_RM.*obj.SingleRunData.TimeSec);
                    ystd = std(y);
                    yErr = repmat(ystd,[1,numel(y)]);
                    ystr = sprintf('Rate @ 300 eV below E0');
                case 'qU_RM'
                    UnitStr = 'eV';
                    y = obj.SingleRunData.qU_RM;
                    ystd = std(y);
                    yErr = repmat(ystd,[1,numel(y)]);
                    ystr = sprintf('qU @ 300 eV below E0');
            end
            
            %% prepare x axis: depending on defined Parameterx
            switch Parameterx
                case 'RhoD'
                    x = obj.SingleRunData.WGTS_CD_MolPerCm2;
                    %xErr = std(obj.SingleRunData.WGTS_CD_MolPerCm2_SubRun);
                    xErr = repmat(0,[1,numel(x)]);
                    if strcmp(DisplayStyle,'Abs')
                        xstr = sprintf('Column density (mol \\cdot cm^{-2})');
                    elseif strcmp(DisplayStyle,'Rel')
                        x = x-wmean(x,obj.SingleRunData.TimeSec);
                        xstr = sprintf('\\rhod - \\langle\\rhod\\rangle (mol \\cdot cm^{-2})');
                    end
                case 'Act'
                    x = obj.SingleRunData.WGTS_CD_MolPerCm2 .* (obj.SingleRunData.WGTS_MolFrac_TT + 0.5*obj.SingleRunData.WGTS_MolFrac_DT + 0.5*obj.SingleRunData.WGTS_MolFrac_HT);
                    xstd = std(x);
                    xErr = repmat(0,[1,numel(x)]);
                    xstr = sprintf('Activity');
                case 'E0'
                    x = FitResults.E0;
                    xErr = FitResults.E0Err;
                    if strcmp(DisplayStyle,'Abs')
                        xstr = sprintf('{\\itE}_0^{fit} (eV)');
                    elseif strcmp(DisplayStyle,'Rel')
                        x = x-wmean(y,1./xErr.^2);
                        xstr = sprintf('{\\itE}_0^{fit} - \\langle{\\itE}_0^{fit}\\rangle (eV)');
                    end   
                case 'mNuSq'
                    x = FitResults.mnuSq;
                    xErr = FitResults.mnuSqErr;
                    if strcmp(DisplayStyle,'Abs')
                        xstr = sprintf('{\\itm}_{\\beta}^2 (eV^2)');
                    elseif strcmp(DisplayStyle,'Rel')
                        x = x-wmean(x,1./xErr.^2);
                        xstr = sprintf('{\\itm}_{\\beta}^2 - < m_{\\beta}^2 >  (eV^2)');
                    end
                case 'B'
                    x = FitResults.B.*1e3;
                    xErr = FitResults.BErr.*1e3;
                    if strcmp(DisplayStyle,'Abs')
                        xstr = sprintf('Background (mcps)');
                    elseif strcmp(DisplayStyle,'Rel')
                        x = x-wmean(x,1./xErr.^2);
                        xstr = sprintf('{\\itB} - \\langle{\\itB}\\rangle  (mcps)');
                    end
                case 'N'
                    x = FitResults.N+1;
                    xErr = FitResults.NErr;
                    if strcmp(DisplayStyle,'Abs')
                        xstr = sprintf('{\\itN}');
                    elseif strcmp(DisplayStyle,'Rel')
                        x = x-wmean(x,1./xErr.^2);
                        xstr = sprintf('{\\itN} - \\langle{\\itN}\\rangle');
                    end
                case 'pVal'
                    x = FitResults.pValue;
                    xErr = zeros(numel(x),1);
                    xstr = sprintf('p-value (%.0f dof)',FitResults.dof(1));
                case 'qU_RM'
                    x = obj.SingleRunData.qU_RM;
                    xstd = std(x);
                    xErr = repmat(xstd,[1,numel(x)]);
                    xstr = sprintf('qU @ 300 eV below E0');
                case 'time'
                    LiveTime = hours(obj.SingleRunData.StartTimeStamp-obj.SingleRunData.StartTimeStamp(1));
                    % get rid of large jumps
                    TimeBreakIndex = find(diff(LiveTime)>5*mean(diff(LiveTime)));
                    x = zeros(numel(LiveTime),1);
                    if isempty(TimeBreakIndex)
                        NxParts=1; TimeBreakIndex = numel(obj.RunList);
                    else
                        NxParts = numel(TimeBreakIndex)+1;
                    end
                    for i=1:NxParts
                        if i==1 %start
                            startI=1;
                            stopI = TimeBreakIndex(i);
                            t = 0;
                        elseif i==NxParts %stop
                            startI= TimeBreakIndex(i-1)+1;
                            stopI = numel(LiveTime);
                            t = t+ abs(LiveTime(TimeBreakIndex(i-1))-LiveTime(TimeBreakIndex(i-1)+1));
                        else
                            startI= TimeBreakIndex(i-1)+1;
                            stopI = TimeBreakIndex(i);
                        t = t + abs(LiveTime(TimeBreakIndex(i-1))-LiveTime(TimeBreakIndex(i-1)+1));
                        end       
                        x(startI:stopI) = LiveTime(startI:stopI)-0.75*t;
                    end
                    xErr = zeros(numel(x),1);
                    xstr = sprintf('Time (hours)');
            end
            
            correlation = corr(x(:),y(:));
            meany = mean(y)
            
            %% Plot
            fig88 = figure('units','normalized','pos',[0.1, 0.1,1,0.6]);
            if strcmp(Fit,'ON')
                if strcmp(Parameterx,'time')
                    [fitobject,err,chi2min,dof] = linFit(x,permute(y,[2 1]),permute(yErr,[2 1]))
                elseif strcmp(Parameterx,'N')
                    [fitobject,err,chi2min,dof] = linFit(x,y,yErr)
                elseif strcmp(Parameterx,'Act') || strcmp(Parameterx,'RhoD')
                    [fitobject,err,chi2min,dof] = linFit(permute(x,[2 1]),permute(y,[2 1]),permute(yErr,[2 1]))
                else
                    [fitobject,err,chi2min,dof] = linFit(x,y,yErr)
                end
                    x1 = min(x):(max(x)-min(x))/100:max(x);
                    y1 = fitobject(1).*x1+fitobject(2);
                    pfit = plot(x1,y1,'LineWidth',2,'Color',rgb('GoldenRod'));
            end
            if(strcmp(Detrend,'ON'))
                y = y(:) - fitobject(1).*x;
                obj.SingleRunData.TBDIS_RM = permute(y,[2 1]).*(obj.SingleRunData.qUfrac_RM.*obj.SingleRunData.TimeSec);
            end
            hold on;
            % plot data
            PlotStyle = { 'o','MarkerSize',8,'MarkerFaceColor',rgb('SkyBlue'),...%rgb('IndianRed'),...%
                'LineWidth',2,'Color',obj.PlotColor};
            e1 = errorbar(x,y,xErr,'horizontal',PlotStyle{:});
            e1.CapSize = 0;
            e2 = errorbar(x,y,yErr,PlotStyle{:});
            e2.CapSize = 0;
            xlabel(xstr);
            ylabel(ystr);
            xmin = min(x);
            if strcmp(Parameterx,'time')
                xlim([min(x)-numel(x)*0.05 max(x)+numel(x)*0.05]);
            elseif xmin>0
                xlim([min(x)*0.99 max(x)*1.01]);
            elseif xmin<0
                xlim([min(x)*1.01 max(x)*1.01]);
            end
            if strcmp(Fit,'ON')
                if abs(fitobject(1))<0.01
                    fitStr = sprintf('%.2g\\pm%.2g m%s/h',fitobject(1)*1e3,err(1)*1e3,UnitStr);
                else
                    fitStr = sprintf('%.2g\\pm%.2g %s/h',fitobject(1),err(1),UnitStr);
                end
                leg = legend([e2,pfit],sprintf('Scanwise %s',chi2leg),sprintf('Fit slope: %s',fitStr));
            else
                leg = legend(e2,sprintf('Scanwise %s',chi2leg));
            end
            title(leg,sprintf('Correlation factor = %.2f',correlation));
            leg.EdgeColor = rgb('Silver');
            leg.Location = 'best';
            leg.Title.FontWeight = 'normal';
            PrettyFigureFormat('Fontsize',LocalFontSize);
            hold off;
            if strcmp(saveplot,'ON')
                if strcmp(obj.MosCorrFlag,'ON')
                    MosStr = '_MosCorr';
                else
                    MosStr = '';
                end
                if strcmp(obj.ROIFlag,'14keV')
                    RoiStr = '_14keV';
                else
                    RoiStr = '';
                end
                SaveDir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/correlations/',obj.DataSet)];
                MakeDir(SaveDir);
                SaveName = sprintf('Corr_%s_%s_%s_%s_%.0feV%s%s.pdf',...
                    Parameterx,Parametery,obj.RunData.RunName,pixlist,obj.GetRange,RoiStr,MosStr);
                export_fig(fig88,[SaveDir,SaveName]);
                fprintf('Save plot to %s \n',[SaveDir,SaveName]);
            end
        end
        function FitResults = FitRunList(obj,varargin)
            
            % Fit all Runs in (stacked or not) RunList independently
            % Save Results in struct: obj.SingleRun_FitResults
            p = inputParser;
            p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF','File'}));
            p.parse(varargin{:});
            RecomputeFlag  = p.Results.RecomputeFlag;

            %----------------- Load from memory or file, if RecomputeFlag OFF---------------------------------------
            % Fit all Runs from RunList and fill structure
            RunObj = cell(numel(obj.RunList),1);
            FitList = obj.RunList;
            FitListArg = {'exclDataStart',obj.exclDataStart,...
                'fixPar',obj.fixPar,...
                'AnaFlag',obj.AnaFlag,...
                'RingMerge',obj.RingMerge,...
                'DataType',obj.DataType,...
                'FitNBFlag',obj.FitNBFlag,...
                'FSDFlag',obj.FSDFlag,...
                'ELossFlag',obj.ELossFlag,...
                'NonPoissonScaleFactor',obj.NonPoissonScaleFactor,...
                'chi2',obj.chi2,...
                'NonPoissonScaleFactor',obj.NonPoissonScaleFactor,...
                'MosCorrFlag',obj.MosCorrFlag,...
                'ROIFlag',obj.ROIFlag};
            
            obj.InitFitPar;
            parAll     = zeros(numel(obj.RunList),obj.nPar);
            errAll     = zeros(numel(obj.RunList),obj.nPar);
            chi2minAll = zeros(numel(obj.RunList),1);
            dofAll     = zeros(numel(obj.RunList),1);
            
            % labeling
            savedir = [getenv('SamakPath'),'tritium-data/fit/',obj.DataSet,'/Uniform/'];
            MakeDir(savedir);
            switch obj.DataType
                case 'Twin'
                    DataTypeLabel = 'Twin_';
                case 'Real'
                    DataTypeLabel = '';
                case 'FitriumTwin'
                    DataTypeLabel = 'FitriumTwin_';
                case 'KafitTwin'
                    DataTypeLabel = 'KafitTwin_';
            end
            fixParstr = ConvertFixPar('freePar',obj.fixPar,'Mode','Reverse');
            
            switch obj.ROIFlag
                case 'Default'
                    RoiStr = '';
                case '14keV'
                    RoiStr = '_14keVROI';
            end
            
            if strcmp(obj.MosCorrFlag,'ON')
                MosStr = '_MosCorr';
            else
                MosStr = '';
            end
            
            savefile = arrayfun(@(x) ...
                sprintf('%sFit%s%.0f_%s_%.0fbE0_freePar%s%s%s.mat',...
                savedir,DataTypeLabel,...
                x,...
                obj.chi2,obj.ModelObj.Q_i-obj.ModelObj.qU(obj.exclDataStart),fixParstr,RoiStr,MosStr), ...
                obj.RunList,'UniformOutput',0);
           
            % load fit results, which are already computed
            fprintf('Attempt to load fit results from file %s \n',savefile{1});
           LoadFilesIndex             = cellfun(@(x) exist(x,'file')==2,savefile);            %logicals, indicate which runs are already fitted
           fprintf('Retrieve %.0f Fit Results from file\n',sum( LoadFilesIndex));
           LoadFiles                  = cellfun(@(x) importdata(x),savefile(LoadFilesIndex)); % import those runs
           parAll(LoadFilesIndex,:)   = cell2mat(arrayfun(@(x) x.FitResult.par,LoadFiles,'UniformOutput',false)); %asign to variables
           errAll(LoadFilesIndex,:)   = cell2mat(arrayfun(@(x) x.FitResult.err,LoadFiles,'UniformOutput',false)); 
           chi2minAll(LoadFilesIndex) = cell2mat(arrayfun(@(x) x.FitResult.chi2min,LoadFiles,'UniformOutput',false));
           dofAll(LoadFilesIndex)     = cell2mat(arrayfun(@(x) x.FitResult.dof,LoadFiles,'UniformOutput',false));
                
           % calculate missing fit results
           if sum(LoadFilesIndex)==numel(obj.RunList) % all runs are already fitted
               fprintf('Fit results from all runs loaded from file \n');
           else
               parfor i=1:numel(obj.RunList)
                   if ismember(i,find(LoadFilesIndex==0)) % if not already calculated
                       RunObj{i} = RunAnalysis('RunNr',FitList(i),FitListArg{:});
                       
                       if ~strcmp(RunObj{i}.chi2,'chi2Stat')
                           RunObj{i}.ComputeCM;
                       end
                       
                       RunObj{i}.Fit;
                       Run = FitList(i);
                       FitResult = RunObj{i}.FitResult;
                       parsave(savefile{i},FitResult,Run,FitListArg);
                       parAll(i,:)   = RunObj{i}.FitResult.par;
                       errAll(i,:)   = RunObj{i}.FitResult.err;
                       chi2minAll(i) = RunObj{i}.FitResult.chi2min;
                       dofAll(i)     = RunObj{i}.FitResult.dof;
                   end
               end
           end
           FitResults.mnuSq    = parAll(:,1);
           FitResults.mnuSqErr = errAll(:,1);
           FitResults.E0 = parAll(:,2)+obj.ModelObj.Q_i;
           FitResults.E0Err =  errAll(:,2);
           FitResults.B = parAll(:,3)+obj.ModelObj.BKG_RateSec_i;
           FitResults.BErr = errAll(:,3);
           FitResults.N = parAll(:,4);
           FitResults.NErr = errAll(:,4);
           FitResults.chi2min = chi2minAll;
           FitResults.dof =  dofAll;
           FitResults.pValue =  1-chi2cdf(chi2minAll,dofAll);
           
           switch obj.chi2
               case 'chi2Stat'
                   obj.SingleRun_FitResults.chi2Stat = FitResults;
               otherwise
                   obj.SingleRun_FitResults.chi2CM   = FitResults;
           end
        end    
        function PlotFitResultsDistributions(obj, varargin)
            % Plot Distribution of Fit Parameters Parameters
            % Stacked Runs
            p=inputParser;
            p.addParameter('SaveHisto','ON',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            SaveHisto   = p.Results.SaveHisto;
            myMainTitle=[sprintf('KATRIN: Run-wise Fit Parameter Distributions - [%.1f - %.1f] eV - %.0f subruns',...
                obj.RunData.qU(obj.exclDataStart),obj.RunData.qU(end), numel(obj.RunData.qU(obj.exclDataStart:end)))];
            
            maintitle=myMainTitle;
            savefile=sprintf('plots/runwise/MRA_FitSingleRuns_%0.0fStackRuns.png',numel(obj.RunList));
            fig1 = figure('Name','MRA Single Runs: Fit Parameter Distributions','NumberTitle','off','rend','painters','pos',[10 10 1400 2000]);
            a=annotation('textbox', [0 0.91 1 0.1], ...
                'String', maintitle, ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center');
            a.FontSize=20;a.FontWeight='bold';
            
            subplot(3,2,1)
            h1=histogram(obj.SingleRun_FitResults.chi2Stat.E0,'FaceColor',rgb('SteelBlue'),'EdgeAlpha',0.5,'LineWidth',2);
            %             hold on
            %             h2=histogram(obj.SingleRun_FitResults.chi2Stat.E0(obj.SingleRunData.Select_all>0),'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.5);
            %             h2.BinWidth = h1.BinWidth;
            %             hold off
            ylabel('runs');
            xlabel('E_0 (eV)');
            %             leg = legend([h1,h2],sprintf('All: %0.0f',numel(obj.RunList)),sprintf('Stacked: %0.0f',numel(obj.StackedRuns)),'Location','northwest');
            %             leg.Color = 'none'; legend boxoff;
            sp1 = sprintf('mean = %.03f eV - \\sigma = %.03f eV',...
                mean(obj.SingleRun_FitResults.chi2Stat.E0),std(obj.SingleRun_FitResults.chi2Stat.E0));title(sp1)
            PrettyFigureFormat;
            subplot(3,2,2)
            h1=histogram(obj.SingleRun_FitResults.chi2Stat.E0Err,'FaceColor',rgb('SteelBlue'),'EdgeAlpha',0.5,'LineWidth',2);
            %             hold on
            %             h2=histogram(obj.SingleRun_FitResults.chi2Stat.E0Err(obj.SingleRunData.Select_all>0),'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.5);
            %             h2.BinWidth = h1.BinWidth;
            %             hold off
            ylabel('runs');
            xlabel('E_0 Error (eV)');
            % leg = legend([h1,h2],'All','Stacked','Location','northeast');
            %leg.Color = 'none'; legend boxoff;
            sp1 = sprintf('mean = %.3f eV - \\sigma = %.3f eV',...
                mean(obj.SingleRun_FitResults.chi2Stat.E0Err),std(obj.SingleRun_FitResults.chi2Stat.E0Err));title(sp1)
            PrettyFigureFormat
            subplot(3,2,3)
            h1=histogram(obj.SingleRun_FitResults.chi2Stat.B*1e3,'FaceColor',rgb('SteelBlue'),'EdgeAlpha',0.5,'LineWidth',2);
            %             hold on
            %             h2=histogram(obj.SingleRun_FitResults.chi2Stat.B(obj.SingleRunData.Select_all>0)*1e3,'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.5);
            %             h2.BinWidth = h1.BinWidth;
            %             hold off
            ylabel('runs');
            xlabel('B (mcps)');
%            leg = legend([h1,h2],'All','Stacked','Location','northeast');
            %leg.Color = 'none'; legend boxoff;
            sp1 = sprintf('mean = %.3f mcps - \\sigma = %.3f mcps',...
                mean(obj.SingleRun_FitResults.chi2Stat.B)*1e3,std(obj.SingleRun_FitResults.chi2Stat.B)*1e3);title(sp1)
            PrettyFigureFormat;
            subplot(3,2,4)
            h1=histogram(obj.SingleRun_FitResults.chi2Stat.BErr*1e3,'FaceColor',rgb('SteelBlue'),'EdgeAlpha',0.5,'LineWidth',2);
%             hold on
%             h2=histogram(obj.SingleRun_FitResults.chi2Stat.BErr(obj.SingleRunData.Select_all>0)*1e3,'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.5);
%             h2.BinWidth = h1.BinWidth;
%             hold off
            ylabel('runs');
            xlabel('B Error (mcps)');
            %leg = legend([h1,h2],'All','Stacked','Location','northeast');
            %leg.Color = 'none'; legend boxoff;
            sp1 = sprintf('mean = %.3f mcps - \\sigma = %.3f mcps',...
                mean(obj.SingleRun_FitResults.chi2Stat.BErr)*1e3,std(obj.SingleRun_FitResults.chi2Stat.BErr)*1e3);title(sp1)
            PrettyFigureFormat
            subplot(3,2,5)
            h1=histogram(obj.SingleRun_FitResults.chi2Stat.N,'FaceColor',rgb('SteelBlue'),'EdgeAlpha',0.5,'LineWidth',2);
%             hold on
%             h2=histogram(obj.SingleRun_FitResults.chi2Stat.N(obj.SingleRunData.Select_all>0),'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.5);
%             h2.BinWidth = h1.BinWidth;
%             hold off
            ylabel('runs');
            xlabel('N ');
%           leg = legend([h1,h2],'All','Stacked','Location','northeast');
            %leg.Color = 'none'; legend boxoff;
            sp1 = sprintf('mean = %.3f  - \\sigma = %.3f ',...
                mean(obj.SingleRun_FitResults.chi2Stat.N),std(obj.SingleRun_FitResults.chi2Stat.N));title(sp1)
            PrettyFigureFormat;
            subplot(3,2,6)
            h1=histogram(obj.SingleRun_FitResults.chi2Stat.NErr,'FaceColor',rgb('SteelBlue'),'EdgeAlpha',0.5,'LineWidth',2);
%             hold on
%             h2=histogram(obj.SingleRun_FitResults.chi2Stat.NErr(obj.SingleRunData.Select_all>0),'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.5);
%             h2.BinWidth = h1.BinWidth;
%             hold off
            ylabel('runs');
            xlabel('N Error');
%             leg = legend([h1,h2],'All','Stacked','Location','northeast');
%             leg.Color = 'none'; legend boxoff;
            sp1 = sprintf('mean = %.3f - \\sigma = %.3f',...
                mean(obj.SingleRun_FitResults.chi2Stat.NErr),std(obj.SingleRun_FitResults.chi2Stat.NErr));title(sp1)
            PrettyFigureFormat
            switch SaveHisto
                case 'ON'
            %publish_figurePDF(gcf,savefile);
            export_fig(gcf,savefile);
            end
            
            if obj.DataType=='Real'
            % Print Table
            [MinE0 , i]    = min(obj.SingleRun_FitResults.chi2Stat.E0(obj.SingleRunData.Select_all>0));
            MinE0_Run      = obj.RunList(i);
            StartMinE0     = obj.SingleRunData.StartTimeStamp(i);
            [MaxE0 , i]    = max(obj.SingleRun_FitResults.chi2Stat.E0(obj.SingleRunData.Select_all>0));
            MaxE0_Run      = obj.RunList(i);
            StartMaxE0     = obj.SingleRunData.StartTimeStamp(i);
            MeanE0         = mean(obj.SingleRun_FitResults.chi2Stat.E0(obj.SingleRunData.Select_all>0));
            StdE0          = std(obj.SingleRun_FitResults.chi2Stat.E0(obj.SingleRunData.Select_all>0));
            
            [MinN , i]     = min(obj.SingleRun_FitResults.chi2Stat.N(obj.SingleRunData.Select_all>0));
            MinN_Run       = obj.RunList(i);
            StartMinN     = obj.SingleRunData.StartTimeStamp(i);
            [MaxN , i]     = max(obj.SingleRun_FitResults.chi2Stat.N(obj.SingleRunData.Select_all>0));
            MaxN_Run       = obj.RunList(i);
            StartMaxN     = obj.SingleRunData.StartTimeStamp(i);
            MeanN          = mean(obj.SingleRun_FitResults.chi2Stat.N(obj.SingleRunData.Select_all>0));
            StdN          = std(obj.SingleRun_FitResults.chi2Stat.N(obj.SingleRunData.Select_all>0));
            
            [MinB , i]     = min(obj.SingleRun_FitResults.chi2Stat.B(obj.SingleRunData.Select_all>0)*1e3);
            MinB_Run       = obj.RunList(i);
            StartMinB      = obj.SingleRunData.StartTimeStamp(i);
            [MaxB , i]     = max(obj.SingleRun_FitResults.chi2Stat.B(obj.SingleRunData.Select_all>0)*1e3);
            MaxB_Run       = obj.RunList(i);
            StartMaxB      = obj.SingleRunData.StartTimeStamp(i);
            MeanB         = mean(obj.SingleRun_FitResults.chi2Stat.B(obj.SingleRunData.Select_all>0))*1e3;
            StdB          = std(obj.SingleRun_FitResults.chi2Stat.B(obj.SingleRunData.Select_all>0))*1e3;
            
            [MinPV , i]    = min(obj.SingleRun_FitResults.chi2Stat.chi2min(obj.SingleRunData.Select_all>0))
            MinPV_Run      = obj.RunList(i);
            StartMinPV     = obj.SingleRunData.StartTimeStamp(i);
            [MaxPV , i]    = max(obj.SingleRun_FitResults.chi2Stat.chi2min(obj.SingleRunData.Select_all>0))
            MaxPV_Run      = obj.RunList(i);
            StartMaxPV     = obj.SingleRunData.StartTimeStamp(i);
            MeanPV         = mean(obj.SingleRun_FitResults.chi2Stat.chi2min(obj.SingleRunData.Select_all>0));
            StdPV          = std(obj.SingleRun_FitResults.chi2Stat.chi2min(obj.SingleRunData.Select_all>0));

            t = PrintTable('Run-wise Fit Results');
            t.addRow('','min','max','average','std');
            t.addRow('Endpoint (eV)',sprintf('%.2f',MinE0),sprintf('%.2f',MaxE0),sprintf('%.2f',MeanE0),sprintf('%.2f',StdE0));
            t.addRow('Run',MinE0_Run,MaxE0_Run,'','');
            %t.addRow('Start Time',datetime(StartMinE0),datetime(StartMaxE0),'','');
            t.addRow('','','','','');
            t.addRow('1 - Normalization',sprintf('%.2f',MinN),sprintf('%.2f',MaxN),sprintf('%.2f',MeanN),sprintf('%.2f',StdN));
            t.addRow('Run',MinN_Run,MaxN_Run,'','');
            %t.addRow('Start Time',StartMinN,StartMaxN,'','');
            t.addRow('','','','','');
            t.addRow('Background (cps)',sprintf('%.2f',MinB),sprintf('%.2f',MaxB),sprintf('%.2f',MeanB),sprintf('%.2f',StdB));
            t.addRow('Run',MinB_Run,MaxB_Run,'','');
            %t.addRow('Start Time',StartMinB,StartMaxB,'','');
            t.addRow('','','','','');
            t.addRow('chisquare ',sprintf('%.2f',MinPV),sprintf('%.2f',MaxPV),sprintf('%.2f',MeanPV),sprintf('%.2f',StdPV));
            t.addRow('Run',MinPV_Run,MaxPV_Run,'','');
            %t.addRow('Start Time',StartMinPV,StartMaxPV,'','');
            t.addRow('','','','','');
            t.display;
            %t.HasHeader = true;
            t.Format = 'tex';
            t.Caption = sprintf('Run-wise tritium beta decay spectra fits - %.0f runs',numel(obj.RunList)');
            t.print;
            end
            
%             figure(12345)
%             h2=histfit(obj.SingleRun_FitResults.chi2Stat.E0(obj.SingleRunData.Select_all>0)-mean((obj.SingleRun_FitResults.chi2Stat.E0(obj.SingleRunData.Select_all>0))));
%             set(h2(1),'facecolor',rgb('IndianRed'),'facealpha',.8,'edgecolor','none'); set(h2(2),'color',rgb('IndianRed'));
%             pd = fitdist(((obj.SingleRun_FitResults.chi2Stat.E0(obj.SingleRunData.Select_all>0))-mean((obj.SingleRun_FitResults.chi2Stat.E0(obj.SingleRunData.Select_all>0))))','Normal')
%             ylabel('runs');
%             xlabel('E_0 (eV)');
%             %leg = legend([h1,h2],sprintf('All: %0.0f',numel(obj.RunList)),sprintf('Stacked: %0.0f',numel(obj.StackedRuns)),'Location','northwest');
%             %leg.Color = 'none'; legend boxoff;
%             sp1 = sprintf('mean = %.3f eV - \\sigma = %.3f eV',...
%                 mean(obj.SingleRun_FitResults.chi2Stat.E0),std(obj.SingleRun_FitResults.chi2Stat.E0));title(sp1)
%             PrettyFigureFormat;
        end
        function SysCM = ComputeSysCMRunList(obj,varargin)
            % Runwise Endpoint (or in general fit parameter) Fit: 1x with stat, 1x stat + sys
            % Compute systematic uncertainty on Parameter for each run
            % Compute covariance matrix for fit to mean Parameter
            p = inputParser;
            p.addParameter('plotDist','ON');
            p.addParameter('Parameter','E0',@(x)ismember(x,{'E0','N','B'}));
            p.addParameter('ReFit','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('CorrCoeff',1,@(x)isfloat(x)); % Correlation Coefficient
            
            p.parse(varargin{:});
            plotDist  = p.Results.plotDist;
            Parameter = p.Results.Parameter;
            ReFit     = p.Results.ReFit;
            CorrCoeff = p.Results.CorrCoeff;
            
            chi2_prev = obj.chi2;
            
            %             if isempty(obj.SingleRun_FitResults.chi2CMall) || strcmp(ReFit,'ON') %compute sigma all
            %                 obj.chi2='chi2CM';
            %                 SysFlag = 'all';
            %                 file_name = sprintf('./results/FitRunListResults_%s_%s_%s_%.0feVbelowE0_fixPar%s.mat',...
            %                     obj.StackFileName,obj.chi2,SysFlag,obj.ModelObj.Q_i-obj.ModelObj.qU(obj.exclDataStart),strrep(obj.fixPar,' ',''));
            %                 try
            %                     d = importdata(file_name);
            %                     obj.SingleRun_FitResults.chi2CMall = d.FitResults;
            %                 catch
            %                     obj.FitRunList('displayHist','OFF','SysFlag','all');
            %                 end
            %             end
            %             if isempty(obj.SingleRun_FitResults.chi2CMcorr) || strcmp(ReFit,'ON') %compute sigma correlated
            %                 obj.chi2='chi2CM';
            %                 SysFlag = 'corr';
            %                 file_name = sprintf('./results/FitRunListResults_%s_%s_%s_%.0feVbelowE0_fixPar%s.mat',...
            %                     obj.StackFileName,obj.chi2,SysFlag,obj.ModelObj.Q_i-obj.ModelObj.qU(obj.exclDataStart),strrep(obj.fixPar,' ',''));
            %                 try
            %                     d = importdata(file_name);
            %                     obj.SingleRun_FitResults.chi2CMcorr = d.FitResults;
            %                 catch
            %                     obj.FitRunList('displayHist','OFF','SysFlag','corr');
            %                 end
            %             end
            %             if isempty(obj.SingleRun_FitResults.chi2Stat) || strcmp(ReFit,'ON') %compute sigma stat
            %                 obj.chi2='chi2Stat';
            %                 SysFlag = 'all';
            %                 file_name = sprintf('./results/FitRunListResults_%s_%s_%s_%.0feVbelowE0_fixPar%s.mat',...
            %                     obj.StackFileName,obj.chi2,SysFlag,obj.ModelObj.Q_i-obj.ModelObj.qU(obj.exclDataStart),strrep(obj.fixPar,' ',''));
            %                 try
            %                     d = importdata(file_name);
            %                     obj.SingleRun_FitResults.chi2Stat = d.FitResults;
            %                 catch
            %                     obj.FitRunList('displayHist','OFF');
            %                 end
            %             end
            
            %compute sigma stat
            obj.chi2='chi2Stat';
            obj.FitRunList('RecomputeFlag',ReFit)
            %compute sigma correlated
            obj.chi2='chi2CMShape';
            obj.FitRunList('RecomputeFlag',ReFit,'SysFlag','corr')
            %compute sigma all
            obj.FitRunList('RecomputeFlag',ReFit,'SysFlag','all')
            
            % sigma sys (corr + uncorr) for each run
            SigmaSysallRunWise = sqrt(obj.SingleRun_FitResults.chi2CMall.([Parameter,'Err']).^2-...
                obj.SingleRun_FitResults.chi2Stat.([Parameter,'Err']).^2);
            % sigma sys (only corr) for each run
            SigmaSysCorrRunWise = sqrt(obj.SingleRun_FitResults.chi2CMcorr.([Parameter,'Err']).^2-...
                obj.SingleRun_FitResults.chi2Stat.([Parameter,'Err']).^2);
            % build covariance matrix
            nRuns = numel(obj.SingleRun_FitResults.chi2CMcorr.([Parameter,'Err']));
            SysCM = CorrCoeff.*SigmaSysCorrRunWise.*ones(nRuns).*SigmaSysCorrRunWise';   % corr sys on off-diagonal
            SysCM = SysCM - diag(diag(SysCM));                                             % remove diagonal part
            SysCM = SysCM + diag(obj.SingleRun_FitResults.chi2CMall.([Parameter,'Err']).^2); % add uncorr + corr + stat on diagonal
            
            obj.chi2 = chi2_prev;
            
            switch plotDist  % Display Results
                case 'ON'
                    fig1 = figure(1);
                    set(fig1, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.9, 0.7]);
                    s1 = subplot(2,1,1);
                    x   = linspace(1,numel(obj.StackedRuns),numel(obj.StackedRuns));
                    plot(x,SigmaSysallRunWise,'o','color',rgb('CadetBlue'),'MarkerFaceColor',rgb('CadetBlue'))
                    %hold on;
                    %plot(x,repmat(SigmaSysStacked,numel(obj.StackedRuns),1),'-','LineWidth',2.5,'color',rgb('IndianRed'))
                    xlabel('run number');
                    ylabel(['\sigma_{sys} on',sprintf(' %s (eV)',Parameter)]);
                    leg = legend(['Runwise Fit <\sigma_{sys}> = ',sprintf('%.2f eV, std = %.2f eV',mean(SigmaSysallRunWise),std(SigmaSysallRunWise))], 'Location','northwest');
                    % ['Stacked Runs Fit \sigma_{sys} = ',sprintf('%.2f eV',SigmaSysStacked)],...
                    
                    leg.Color = 'none'; legend boxoff;
                    r = string(obj.StackedRuns);
                    xticks(x);
                    xticklabels(r);
                    xtickangle(45);
                    PrettyFigureFormat;
                    set(gca,'FontSize',18);
                    set(get(gca,'XAxis'),'FontSize',11);
                    set(get(gca,'XLabel'),'FontSize',18);
                    leg.FontSize = 12;
                    xlim([min(x) max(x)]);
                    title(sprintf('Samak %s Systematic Uncertainty - %.0f eV below Endpoint',Parameter,obj.ModelObj.Q_i-obj.ModelObj.qU(obj.exclDataStart)));
                    
                    s2 = subplot(2,1,2);
                    x   = linspace(1,numel(obj.StackedRuns),numel(obj.StackedRuns));
                    yyaxis left
                    yAx = get(gca,'YAxis');
                    set(yAx(1),'Color',rgb('CadetBlue'));
                    p1 = plot(x,obj.SingleRun_FitResults.chi2CMall.([Parameter,'Err']).^2,'o','MarkerFaceColor',rgb('CadetBlue'),'MarkerEdgeColor',rgb('CadetBlue'));
                    hold on;
                    ylabel(['\sigma_{all}^2 on ',sprintf('%s (eV)',Parameter)]);
                    xlabel('run number');
                    yyaxis right
                    set(yAx(2),'Color',rgb('CornflowerBlue'));
                    p2 = plot(x,obj.SingleRun_FitResults.chi2Stat.([Parameter,'Err']).^2,'o','MarkerFaceColor',rgb('CornflowerBlue'),'MarkerEdgeColor',rgb('CornflowerBlue'));
                    ylabel(['\sigma_{stat}^2 on ',sprintf('%s (eV)',Parameter)]);
                    r = string(obj.StackedRuns);
                    xticks(x);
                    xticklabels(r);
                    xtickangle(45);
                    PrettyFigureFormat;
                    set(gca,'FontSize',18);
                    set(get(gca,'XAxis'),'FontSize',11);
                    set(get(gca,'XLabel'),'FontSize',18);
                    xlim([min(x) max(x)]);
            end
        end
        function [PlotMean, err, chi2min, dof, SysCM, ResultCM] = FitRunList_AveragePar(obj,varargin)
            % fit to mean value of fit parameter over runs
            % input: fit Parameter of interest
            % ouput: mean and standard error of mean
            p=inputParser;
            p.addParameter('Parameter','E0',@(x)ismember(x,{'E0','N','B'}));
            p.addParameter('CorrCoeff',1,@(x)isfloat(x)); % Correlation Coefficient of systematics
            p.addParameter('ReFit','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF','png','fig','pdf'}));
            
            p.parse(varargin{:});
            CorrCoeff = p.Results.CorrCoeff;
            Parameter = p.Results.Parameter;
            ReFit     = p.Results.ReFit;
            saveplot  = p.Results.saveplot;
            
            SysCM = obj.ComputeSysCMRunList('plotDist','OFF','Parameter',Parameter,...
                'ReFit',ReFit,'CorrCoeff',CorrCoeff);
            fprintf('----------BEGIN FIT MINUIT-------------- \n');
            % Init Fit Parameter
            ResultCM = obj.SingleRun_FitResults.chi2CMall;
            parMean = mean(ResultCM.(Parameter));
            % Minuit Arguments
            tmparg = sprintf(['fix %s;set pri -10;'...
                'migrad ; hesse'],'');
            % Minuit Input: Init Fit Parameter, Data, Covariance Matrix
            Args = {parMean, {ResultCM.(Parameter), SysCM}, '-c', tmparg};
            [par, err, chi2min, errmat] = fminuit('chi2meanErr',Args{:});
            dof = numel(ResultCM.(Parameter))-1;
            results = {par, err, chi2min, errmat, dof};
            
            % display
            if strcmp(Parameter,'E0')
                PlotPar = ResultCM.(Parameter);
                PlotUnit = 'eV';
                PlotMean =par;
            elseif strcmp(Parameter,'B')
                PlotPar = ResultCM.(Parameter);%+obj.ModelObj.BKG_RateSec_i;
                PlotUnit = 'cps';
                PlotMean =par;%+obj.ModelObj.BKG_RateSec_i;
            else
                PlotPar = ResultCM.(Parameter);
                PlotUnit = '';
                PlotMean =par;
            end
            
            fprintf(2,'------------------------------------------------------------\n')
            fprintf(2,'Fit Results: \n')
            fprintf(2,'Average %s = %.2f +/- %.2f %s \n',Parameter,PlotMean,err,PlotUnit)
            fprintf(2,'chi2/dof   = %.2f/%.0f \n',chi2min,dof)
            fprintf(2,'------------------------------------------------------------\n')
            
            fig6 = figure(6);
            set(fig6, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.9, 0.7]);
            x = linspace(1,numel(PlotPar),numel(PlotPar));
            e1 = errorbar(x,PlotPar,ResultCM.([Parameter,'Err']),...
                'o','MarkerFaceColor',rgb('CadetBlue'),'Color',rgb('CadetBlue'),'MarkerSize',8);
            e1.LineWidth = 2;
            hold on
            p1 = plot(x,repmat(PlotMean,1,numel(x)),'--','LineWidth',2.5,'Color',rgb('IndianRed'));
            fit_leg = sprintf('Fitted Average \n%s= %.2f \\pm %.2f %s',...
                Parameter,PlotMean,err,PlotUnit);
            leg = legend([e1,p1],'Runwise Fit (stat + sys)',fit_leg);
            set(leg,'color','none'); legend boxoff;
            ylabel(sprintf('%s (%s)',Parameter,PlotUnit));
            xlabel('run');
            title(sprintf('Samak Fit to %.0f Runs \n%s Fit %.0f eV below Endpoint - %.1f Correlation Coefficient',...
                numel(obj.StackedRuns),Parameter,obj.ModelObj.Q_i-obj.ModelObj.qU(obj.exclDataStart),CorrCoeff));
            PrettyFigureFormat;
            xlim([1 max(x)]);
            %ylim([min(PlotPar-ResultCM.([Parameter,'Err']))-0.3*max(ResultCM.([Parameter,'Err'])) max(PlotPar+ResultCM.([Parameter,'Err']))+0.3*max(ResultCM.([Parameter,'Err']))]);
            set(gca,'FontSize',20);
            xticks(x);
            xticklabels(string(obj.StackedRuns));xtickangle(45);
            hold off;
            savename = sprintf('FitCMRunwiseMean%s_Runs%s_%.0feV_%.1fCorrCoeff',...
                Parameter,obj.StackFileName,obj.ModelObj.Q_i-obj.ModelObj.qU(obj.exclDataStart),CorrCoeff);
            if strcmp(saveplot,'ON') % save all formats
                export_fig(fig6,['./plots/',savename,'.','png']);
                savefig(fig6,['./plots/fig/',savename,'.','fig'],'compact');
                publish_figurePDF(fig6,['./plots/pdf/',savename,'.pdf']);
            elseif strcmp(saveplot,'png')
                export_fig(['./plots/',savename,'.png']);
            elseif strcmp(saveplot,'pdf')
                publish_figurePDF(fig6,['./plots/pdf/',savename,'.pdf']);
            elseif strcmp(saveplot,'fig')
                savefig(fig6,['./plots/fig/',savename,'.','fig'],'compact');
            end
        end
        function FitResults = RM_Fluctuations(obj,varargin)
            p = inputParser;
            p.addParameter('Dist','Normal',@(x)ismember(x,{'Normal','Poisson'}));
            p.parse(varargin{:});
            Dist = p.Results.Dist;
            
            rate300 = obj.SingleRunData.TBDIS_RM./(obj.SingleRunData.qUfrac_RM.*obj.SingleRunData.TimeSec);
            histogram(rate300)
            PrettyFigureFormat;
            hold on;
            pd = fitdist(rate300(:),Dist)
            x = min(rate300):1:max(rate300);
            y = 30*numel(rate300)*pdf(pd,x);
            plot(x,y)
            hold off;
        end
    end
    
    
    methods % Methods for Labeling and RunList
        function StackFileName = GetStackFileName(obj) %Labeling Files with stacked runs
            switch obj.DataType
                case {'Real','Twin','FitriumTwin','KafitTwin'}
                    if  ischar(obj.RunList)
                        StackFileName = obj.RunList; 
                    elseif isfloat(obj.RunList)
                        RunString           = strrep(num2str(obj.RunList),'  ','_');
                        StackFileName   = [strtok(obj.DataSet,'.'),'Stack',RunString];
                    end
                case 'Fake'
                 StackFileName = ['Stack_',obj.FakeStudyName,'_',num2str(numel(obj.RunList)),'Runs'];
            end
            if strcmp(obj.DataType,'Twin')
                StackFileName = ['Twin',StackFileName,obj.SetTwinOrFakeFileName];
            elseif strcmp(obj.DataType,'FitriumTwin')
                StackFileName = ['FitriumTwin',StackFileName,obj.SetTwinOrFakeFileName];
            elseif  strcmp(obj.DataType,'KafitTwin')
                StackFileName = ['KafitTwin',StackFileName,obj.SetTwinOrFakeFileName];
            end
        end
        function RunList = GetRunList(obj)
            % ------------------------------------------------------- %
            ListName = strrep(strrep(obj.StackFileName,'Twin',''),obj.SetTwinOrFakeFileName,'');
            
            if contains(ListName,'StackCD') || contains(ListName,'FT')
                % First Trititum runs
                RunList = obj.GetFTRunList(ListName);
            elseif contains(ListName,'KNM1')
                % KNM1
                RunList = obj.GetKNM1RunList(ListName);
            elseif contains(ListName,'KNM2')
                % KNM2
                RunList = obj.GetKNM2RunList(ListName);
            else
                % miscellaneous
                switch ListName
                    case 'KITNUSCAN'
                        % KIT NU SCAN - 2hMTD- March  6 - 15% RhoD - stable
                        RunList = sort([48872:48873,48879:48880,48920:48925]);
                        RunList=RunList(1:end);
                    case 'March132019'
                        % Nu MASS SCAN - 1hMTD - March 13 - 50%rhoD - unstable
                        RunList = sort([50380:50389]);RunList=RunList(1:end);
                    case 'March242019'
                        % Nu MASS SCAN - 1hMTD - March 24 - unstable
                        RunList = sort([50531:50577]);RunList=RunList(1:end);
                    case 'March2019'
                        % Nu MASS SCAN - 1hMTD - March 24 - unstable
                        RunList = sort([48872:48873,48879:48880,48920:48925,50380:50389,50531:50577]);RunList=RunList(1:end);
                    otherwise
                        fprintf('RunList Name unknown! \n');
                end
            end
        end
        function RunList = GetFTRunList(obj,ListName)
            % First tritium runs lists
            % (taken from BREW, http://katana.npl.washington.edu/~sanshiro/brew/FirstTritiumTest/)
            RunList_100up = [...% 100% column density - up scans (== from -1600eV to +40eV) with respect to E0
                40541, 40543, 40604, 40611, 40613, 40667, 40669, 40671, 40673, ...
                40675, 40677, 40679, 40681, 40683, 40685, 40687, 40689, 40691, 40693, ...
                40977, 40980, 40983, 40986, 40989, 40992, 40995, 41002, 41005, 41008, ...
                41011, 41014, 41019, 41022,  41026,  41029, 41032]; %
            RunList_100down = [...% 100% column density - down scans (== from +40eV to -1600eV) with respect to E0
                40531, 40540, 40542, 40603, 40668, 40670, 40672, 40674, ...
                40676, 40678, 40680, 40682, 40684, 40686, 40688, 40690, 40692, 40976, ...
                40979, 40982, 40985, 40988, 40991, 40994, 40997, 41007, 41010, 41013, ...
                41016, 41017, 41020, 41023, 41025,  41028, 41031];
            RunList_100random = [...% 100% column density - random scans
                41003,41006,41012, 41015,41018, 41021, 41024,  41027,  41030, 41033];%, 40610, 40612, 40539];
            RunList_FTextra = [40538,40998,41001,41004];
            
            switch ListName
                case 'StackCD100all'%all runs -  with 100% column density
                    RunList = sort([RunList_100up,RunList_100down]);
                case 'StackCD100up'
                    RunList = RunList_100up;
                case 'FTpaper'
                    RunList = importdata([getenv('Samak2.0'),'First-Tritium-Paper/RunList/RunListFT.mat']);
                    %RunList = sort([RunList_100up,RunList_100down,RunList_100random,RunList_FTextra]);
                case 'StackCD100down'
                    RunList = RunList_100down;
                case 'StackCD100random'
                    RunList = RunList_100random;
                case 'StackCD100_3hours'
                    RunList_CD100all = sort([RunList_100up,RunList_100down]);
                    start3h = find(RunList_CD100all==40667);
                    stop3h  = find(RunList_CD100all==40693);
                    RunList = RunList_CD100all(start3h:stop3h);
                case 'StackAnnaUp'
                    RunList = [41002 41005 41008 41011 41014 41017 41020 41023 41026 41029];
                case 'StackAnnaDown'
                    RunList = [41001 41004 41007 41010 41013 41016 41019 41022 41025 41028 41031];
                case 'StackAnnaRandom'
                    RunList = [41003 41006 41012 41015 41018 41021 41024 41027 41030];
                case 'FT_RD24'
                    RunList = [40771:40772,40794:40804];
                case 'FT_RD48'
                    RunList = [40763:40766];
                case 'FT_RD72'
                    RunList = [40926:40935];
            end
            
            % ------------------------------------------------------- %
            % FT Remove Runs
            %             RunList = RunList(RunList~=40531); % Large DT
            %             RunList = RunList(RunList~=40976); % Large DT
            %             RunList = RunList(RunList~=40977); % Large DT
            %             RunList = RunList(RunList~=40979); % Large DT
            %             RunList = RunList(RunList~=40980); % Large DT
            %             RunList = RunList(RunList~=40982); % Large DT
            %             RunList = RunList(RunList~=40983); % Large DT
            %             RunList = RunList(RunList~=41032); % Low DT
            %             RunList = RunList(RunList~=40995);  % Weird Time
            %             RunList = RunList(RunList~=40670); % Large RhoD
            %             RunList = RunList(RunList~=40671); % Large RhoD
            %             RunList = RunList(RunList~=40692); % Large RhoD
            %             RunList = RunList(RunList~=40693); % Large RhoD
            %             RunList = RunList(RunList~=40540); % Low RhoD
            %             RunList = RunList(RunList~=40541); % Low RhoD
            %             RunList = RunList(RunList~=40604); % Low RhoD
            %             RunList = RunList(RunList~=41019); % Low RhoD
            %             RunList = RunList(RunList~=41026); % Low RhoD
            % ------------------------------------------------------- %
            
        end
        function RunList = GetKNM1RunList(obj,ListName)
            Version = 'RunSummary-Durable3a-fpd00';
            % Excluded Runs: BREW SC
            knm1exclude_brew_R  = sort([51691 51710 51861 51896 51918]);
            % Excluded Runs: rhoD top-up
            knm1exclude_rhoD_R  = sort([51861 51877 51878]);
            knm1exclude_rhoD_Y  = sort([51445 51471 51578 51702 51918]);
            RunExcList          = unique([knm1exclude_brew_R knm1exclude_rhoD_Y knm1exclude_rhoD_R]);
            
            if    contains(ListName,'KNM1_m183mvRW')
                FirstRun = 51410; LastRun  = 51442;
            elseif contains(ListName,'KNM1_m149mvRW')
                FirstRun = 51443; LastRun  = 51937;
            elseif contains(ListName,'KNM1_part1_m149mvRW')
                FirstRun = 51531;  LastRun = 51531+20;
            elseif contains(ListName,'KNM1_part1')
                FirstRun = 51442;  LastRun = 51826;
                RunExcList = [51691,51641];
            elseif contains(ListName,'KNM1_part2_1')
                FirstRun = 51894;  LastRun = 1e9;
                RunExcList = [51861,51896,51897,51899,51918];
            elseif contains(ListName,'KNM1_part2')
                FirstRun = 51826;  LastRun = 1e9;
                RunExcList = [51861,51896,51897,51899,51918];
            elseif contains(ListName,'KNM1-RhoD')
                % Selection of rhod close to the mean +-1 sigma
                FirstRun = 51410;  LastRun = 51937;
                RunExcList = unique([RunExcList 51499  51500  51501  51502  51503  51516  51517  51521  51522  51547  51548  51549  51550  51551  51552  51553  51554  51555  51556  51557  51558  51559  51560  51561  51562  51563  51564  51565  51566  51579  51640  51641  51642  51643  51644  51645  51646  51647  51651  51652  51653  51654  51669  51670  51671  51672  51673  51674  51675  51676  51687  51688  51689  51690  51692  51693  51694  51695  51696  51701]);
            elseif contains(ListName,'KNM1-AntiRhoD')
                FirstRun = 51410;  LastRun = 51937;
                RunExcList = unique([RunExcList 51410  51411  51412  51413  51414  51415  51416  51417  51418  51419  51420  51421  51422  51423  51424  51425  51426  51441  51442  51443  51444  51446  51447  51448  51449  51450  51451  51452  51453  51454  51455  51456  51457  51458  51459  51460  51461  51462  51463  51464  51465  51466  51467  51468  51469  51470  51472  51473  51474  51475  51476  51477  51478  51479  51480  51481  51486  51487  51488  51489  51490  51491  51492  51493  51494  51495  51496  51497  51498  51523  51524  51525  51526  51527  51528  51529  51530  51531  51532  51533  51534  51535  51536  51537  51538  51539  51540  51541  51542  51543  51544  51545  51546  51580  51581  51582  51583  51584  51585  51586  51639  51655  51656  51657  51658  51659  51660  51664  51665  51677  51678  51679  51680  51681  51682  51683  51684  51685  51686  51703  51704  51705  51706  51707  51708  51709  51822  51823  51824  51825  51826  51827  51828  51829  51830  51831  51832  51833  51834  51835  51836  51837  51838  51839  51840  51841  51842  51843  51844  51845  51846  51847  51848  51849  51850  51851  51852  51853  51854  51855  51856  51857  51858  51859  51860  51870  51871  51872  51873  51874  51875  51876  51879  51880  51881  51882  51883  51884  51885  51886  51887  51888  51889  51890  51891  51892  51893  51894  51895  51898  51908  51909  51910  51911  51912  51913  51919  51920  51921  51922  51923  51924  51925  51926  51927  51928  51929  51930  51931  51932  51933  51934  51935  51936]);
            elseif contains(ListName,'KNM1-TpurityLow')
                FirstRun = 51410;  LastRun = 51937;
                RunExcList = unique([RunExcList 51556  51561  51562  51563  51564  51565  51566  51579  51580  51581  51582  51583  51584  51585  51586  51639  51640  51641  51642  51643  51644  51645  51646  51655  51656  51657  51658  51659  51660  51673  51674  51677  51679  51680  51681  51684  51685  51686  51687  51689  51690  51692  51693  51694  51696  51706  51707  51708  51709  51822  51823  51824  51825  51826  51827  51828  51829  51830  51831  51832  51833  51834  51835  51836  51837  51838  51839  51840  51841  51842  51843  51844  51845  51846  51847  51848  51849  51850  51851  51852  51853  51854  51855  51856  51857  51858  51859  51860  51870  51871  51872  51873  51874  51875  51876  51879  51880  51881  51882  51883  51884  51885  51886  51887  51888  51889  51890  51891  51892  51893  51894  51895  51898  51908  51909  51910  51911  51912  51913  51919  51920  51921  51922  51923  51924  51925  51926  51927  51928  51929  51930  51931  51932  51933  51934  51935  51936]);
            elseif contains(ListName,'KNM1-TpurityHigh')
                FirstRun = 51410;  LastRun = 51937;
                RunExcList = unique([RunExcList 51410  51411  51412  51413  51414  51415  51416  51417  51418  51419  51420  51421  51422  51423  51424  51425  51426  51441  51442  51443  51444  51446  51447  51448  51449  51450  51451  51452  51453  51454  51455  51456  51457  51458  51459  51460  51461  51462  51463  51464  51465  51466  51467  51468  51469  51470  51472  51473  51474  51475  51476  51477  51478  51479  51480  51481  51486  51487  51488  51489  51490  51491  51492  51493  51494  51495  51496  51497  51498  51499  51500  51501  51502  51503  51516  51517  51521  51522  51523  51524  51525  51526  51527  51528  51529  51530  51531  51532  51533  51534  51535  51536  51537  51538  51539  51540  51541  51542  51543  51544  51545  51546  51547  51548  51549  51550  51551  51552  51553  51554  51555  51557  51558  51559  51560  51647  51651  51652  51653  51654  51664  51665  51669  51670  51671  51672  51675  51676  51678  51682  51683  51688  51695  51701  51703  51704  51705]);
            elseif contains(ListName,'KNM1-FirstHalfTime')
                % Selection of 1st Half of the Runs
                FirstRun = 51410;  LastRun = 51536;
            elseif contains(ListName,'KNM1-MiddleHalfTime')
                % Selection of Middle Half of the Runs
                FirstRun = 51537;  LastRun = 51706;
            elseif contains(ListName,'KNM1-LastHalfTime')
                % Selection of Middle Half of the Runs
                FirstRun = 51707;  LastRun = 51937;
            elseif contains(ListName,'KNM1upScan')
                % Selection of Middle Half of the Runs
                FirstRun = 51410;  LastRun = 51937;
                RunExcList = unique([RunExcList  51411  51413  51415  51417  51419  51421  51423  51425  51442  51444  51446  51448  51450  51452  51454  51456  51458  51460  51462  51464  51466  51468  51470  51472  51474  51476  51478  51480  51487  51489  51491  51493  51495  51497  51499  51501  51503  51522  51524  51526  51528  51530  51532  51534  51536  51538  51540  51542  51544  51546  51548  51550  51552  51554  51556  51558  51560  51562  51564  51566  51579  51581  51583  51585  51640  51642  51644  51646  51652  51654  51656  51658  51660  51665  51670  51672  51674  51676  51678  51680  51682  51684  51686  51688  51690  51692  51694  51696  51704  51706  51708  51823  51825  51827  51829  51831  51833  51835  51837  51839  51841  51843  51845  51847  51849  51851  51853  51855  51857  51859  51871  51873  51875  51879  51881  51883  51885  51887  51889  51891  51893  51895  51909  51911  51913  51919  51921  51923  51925  51927  51929  51931  51933  51935 ]);
            elseif contains(ListName,'KNM1downScan')
                % Selection of Middle Half of the Runs
                FirstRun = 51410;  LastRun = 51937;
                RunExcList = unique([RunExcList  51410  51412  51414  51416  51418  51420  51422  51424  51426  51441  51443  51447  51449  51451  51453  51455  51457  51459  51461  51463  51465  51467  51469  51473  51475  51477  51479  51481  51486  51488  51490  51492  51494  51496  51498  51500  51502  51516  51517  51521  51523  51525  51527  51529  51531  51533  51535  51537  51539  51541  51543  51545  51547  51549  51551  51553  51555  51557  51559  51561  51563  51565  51580  51582  51584  51586  51639  51641  51643  51645  51647  51651  51653  51655  51657  51659  51664  51669  51671  51673  51675  51677  51679  51681  51683  51685  51687  51689  51693  51695  51701  51703  51705  51707  51709  51822  51824  51826  51828  51830  51832  51834  51836  51838  51840  51842  51844  51846  51848  51850  51852  51854  51856  51858  51860  51870  51872  51874  51876  51880  51882  51884  51886  51888  51890  51892  51894  51898  51908  51910  51912  51920  51922  51924  51926  51928  51930  51932  51934  51936]);
            elseif contains(ListName,'KNM1_Random')
                FirstRun = 51410;  LastRun = 51937;
                % exclude runs later
            elseif contains(ListName,'KNM1') % Must be the last
                FirstRun = 51410;  LastRun = 51937;
            elseif contains(ListName,'KNM1_300mvRW')
                RunList = sort(51022:51048);
                return
            elseif contains(ListName,'KNM1_175mvRW')
                RunList = sort([51279:51287,51289:51298,51300:51326,51327:51329,51392:51401]);
                return
            else
                fprintf('Run List not known \n');
                return
            end
            
            % Read All KNM1 HD5 File
            tmp = dir([getenv('SamakPath'), '/tritium-data/hdf5/',GetDataSet(FirstRun), '/*.h5']);
            h5list = arrayfun(@(x) x.name,tmp,'UniformOutput',0);
            h5list =  extractBefore(h5list,'.h5');
            h5list = str2double(extractAfter(h5list,Version));
            h5list(isnan(h5list)) = []; % %delete everything that has a different version
            
            % Truncate to exclude Runs from RunExcList:
            h5list(h5list<FirstRun)=[];
            h5list(h5list>LastRun)=[];
            h5list(ismember(h5list,RunExcList))=[];
            HDF5readallruns('h5runlist',h5list,'reConvert','OFF','DataSet',GetDataSet(FirstRun)); %looks for unconverted runs and converts if needed
            
            RunList=sort(h5list);
            
            if contains(obj.StackFileName,'KNM1_Random')
                r = randn(numel(RunList),1);
                r(r<=median(r)) = 0;
                r(r>median(r)) = 1;
                RunList = RunList(logical(r));
            end
            
        end
        function RunList = GetKNM2RunList(obj,ListName)
            % Run Summary version
            Version = 'RunSummary-Prompt4b-fpd00';
            
            % Excluded Runs
            RunExcListIE       = unique([56175,56185,56283,56318,56408,56410,56631, 56687, 56705, 57021, 57037]);      % Inner electrode power supply breakdown
            RunExclListNaN     = [56331,56630];                                                                      % NaN value in K35
            RunExclListT2      = [56348,56632:56638,56665:56668,56675:56683];    % No T2 value or LARA problem
            RunExclListFPD     = [56280, 56304]; % detector slow control
            RunExclListqU      = [56160:56277,56332];        % voltage spike or other HV problems
            RunExclListE0      = 56343;        % endpoint very low
            RunExclListOther   = [56415, 56604];        % other BREW comments
            RunExcList = [RunExcListIE,RunExclListNaN,RunExclListT2,RunExclListFPD,...
                RunExclListqU,RunExclListE0,RunExclListOther];

            if ismember(ListName,{'KNM2_Prompt','KNM2_RandHalf'})  % all runs 
                FirstRun = 56160;  LastRun = 57137;
            elseif strcmp(ListName,'KNM2_RW1')   % rear wall setting 1 (different names for back compatibility)
                FirstRun = 56160;  LastRun = 56479;
            elseif strcmp(ListName,'KNM2_RW2')   % rear wall setting 2
                FirstRun = 56560;  LastRun = 56713;
            elseif strcmp(ListName,'KNM2_RW3')   % rear wall setting 3
                FirstRun = 57015;  LastRun = 57137;
            elseif strcmp(ListName,'KNM2_RW12') % rear wall setting 1 + 2
                 FirstRun = 56160;  LastRun = 56713;
            elseif strcmp(ListName,'KNM2_RW23') % rear wall setting 2 + 3
                 FirstRun = 56560;  LastRun = 57137; 
            elseif strcmp(ListName,'KNM2_RW13') % rear wall setting 1 + 3
                 FirstRun = 56160;  LastRun = 57137; 
                 RunExcList = [RunExcList,56560:56713];
            else
                fprintf('RunList Name unknown! \n'); 
            end
            
            % Read All KNM2 HD5 File
            tmp = dir([getenv('SamakPath'), '/tritium-data/hdf5/',GetDataSet(FirstRun), '/*.h5']);
            h5list = arrayfun(@(x) x.name,tmp,'UniformOutput',0);
            h5list =  extractBefore(h5list,'.h5');
            h5list = str2double(extractAfter(h5list,Version));
            h5list(isnan(h5list)) = []; % %delete everything that has a different version
            
            % Truncate to exclude Runs from RunExcList:
            h5list(h5list<FirstRun)=[];
            h5list(h5list>LastRun)=[];
            h5list(ismember(h5list,RunExcList))=[];
            HDF5readallruns('h5runlist',h5list,'reConvert','OFF','DataSet',GetDataSet(FirstRun)); %looks for unconverted runs and converts if needed
            
            RunList=sort(h5list);
            
            if strcmp(ListName,'KNM2_RandHalf')              % random half
               RandIndex = randperm(numel(RunList));         % randomly permute runlist indices
               nRunsHalf = ceil(numel(RunList)/2);           % take only half of all runs
               RunList   = RunList(RandIndex(1:nRunsHalf));
            end
        end
        function RunList = GetFakeRunList(obj)
                
                switch obj.StackFileName
                    case [obj.FakeRunType 'KNM1_1d']
                        RunList = [1:1:12];
                    case [obj.FakeRunType 'KNM1_10d']
                        RunList = [1:1:120];
                    case [obj.FakeRunType 'KNM1_30d']
                        RunList = [1:1:360];
                    case [obj.FakeRunType 'KNM5y']
                        RunList = [1:1:90];
                    case [obj.FakeRunType 'KITNUSCAN']
                        RunList = [1:1:360];
                    otherwise
                        fprintf('Fake RunList Name unknown! \n');
                end
                
        end       
        end
        methods %Sanity Checks and Sanity Plots
        function qUDistribution(obj,varargin)
            % Plot qU Distribution for Stacked Runs or all runs
            p=inputParser;
            p.addParameter('saveplot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('legendFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Mode','Stacked',@(x)ismember(x,{'Stacked','All'}));
            p.parse(varargin{:});
            saveplot     = p.Results.saveplot;
            legendFlag   = p.Results.legendFlag;
            Mode         = p.Results.Mode;
            
            % Init
            legend_str = cell(obj.nRuns,1);
            %counter = 1; %counter
            figScatter = figure(2);
            set(figScatter, 'Units', 'normalized', 'Position', [0.9, 0.9, 1, 0.8]);
            x = linspace(obj.RunData.qU(obj.exclDataStart)-18700,obj.RunData.qU(end)-18400,10);
            plot(x,zeros(numel(x),1),'-','Color',rgb('DarkSlateGray'),'LineWidth',2);
            hold on;
            pt = plot(x,obj.StackTolerance.*ones(numel(x),1),'--','Color',rgb('DarkSlateGray'),'LineWidth',2);
            plot(x,-obj.StackTolerance.*ones(numel(x),1),'--','Color',rgb('DarkSlateGray'),'LineWidth',2);
            
            % Compute Actual qU mean of stacked runs
            qUmean      = obj.RunData.qU;
            %qUmean      =  obj.StackWmean(obj.SingleRunData.qU,obj.SingleRunData.TimeperSubRunperPixel);%(:,obj.SingleRunData.Select_all,:),2);
            
            for rr = 1:obj.nRuns
                c = colormap(jet(obj.nRuns));
                % Plot
                if(obj.SingleRunData.Select_all(rr))
                    p = scatter(obj.RunData.qU(obj.exclDataStart:end)-obj.ModelObj.Q_i,...
                        obj.SingleRunData.qU(obj.exclDataStart:end,rr)-qUmean(obj.exclDataStart:end),...
                        80,c(rr,:),'o','filled','MarkerEdgeColor',rgb('DimGray'));
                else
                    p = scatter(obj.RunData.qU(obj.exclDataStart:end)-obj.ModelObj.Q_i,...
                        obj.SingleRunData.qU(obj.exclDataStart:end,rr)-qUmean(obj.exclDataStart:end),...
                        80,c(rr,:),'x','filled','MarkerEdgeColor',rgb('DimGray'));
                end
                
                PrettyFigureFormat;
                hold on;
                set(gca,'FontSize',22);
                xlabel(sprintf('mean retarding potential < qU >  -  %.1f  (eV)',obj.ModelObj.Q_i),'FontSize',20);
                ylabel('qU  -  < qU > (eV)','FontSize',20);
                legend_str{rr} = sprintf('%u',obj.RunList(rr));
                %counter = counter +1;
            end
            
            title(sprintf('%.0f Runs Analyzed - %0.f Runs Stacked - Starting qU-E_0 = %g eV',obj.nRuns,numel(obj.StackedRuns),obj.RunData.qU(obj.exclDataStart)-obj.ModelObj.Q_i));
            
            if strcmp(legendFlag,'ON')
                myleg= legend(legend_str{:},'Location','bestoutside');
                myleg.NumColumns = 3;
                myleg.Title.String = sprintf('%.0f runs',obj.nRuns);
                else
                    myleg = legend([pt],'qU tolerance','Location','northwest');
                
            end
            legend('boxoff');
            myleg.FontSize = 22;
            grid on;
            xlim(1.05.*[min(obj.RunData.qU(obj.exclDataStart:end)-obj.ModelObj.Q_i) max(obj.RunData.qU(obj.exclDataStart:end)-obj.ModelObj.Q_i)]);
            hold off;
            
            if abs(min(ylim))<=obj.StackTolerance
                ylim([-max(ylim) ,max(ylim)]);
            elseif abs(min(ylim))<=abs(max(ylim))
                ylim([-max(ylim) ,max(ylim)]);
            elseif abs(max(ylim))<=abs(min(ylim))
                ylim([min(ylim) ,-min(ylim)]);
            end
            if strcmp(saveplot,'ON')
                 if exist('./plots','dir')~=7
                    mkdir plots
                end
%                 if exist('./plots/png','dir')~=7
%                     mkdir plots/png
%                 end
                save_name = sprintf('./plots/ScatterPlot-StackedRuns_%u-%u_excl%u.png',min(obj.StackedRuns),max(obj.StackedRuns),obj.exclDataStart);
                if contains(save_name,' ')
                    save_name=strrep(save_name,' ','_');
                end
                %publish_figurePDF(figScatter,save_name);
                print(figScatter,save_name,'-dpng','-r400');
            end
        end       
        function qUfracDistribution(obj,varargin)
            % Plot qU Distribution for Stacked Runs
            p=inputParser;
            p.addParameter('saveplot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('legendFlag','ON',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            saveplot = p.Results.saveplot;
            legendFlag   = p.Results.legendFlag;
            % Init
            legend_str = cell(obj.nRuns,1);
            %counter = 1; %counter
            figScatter = figure(2);
            set(figScatter, 'Units', 'normalized', 'Position', [0.9, 0.9, 1, 0.8]);
            x = linspace(obj.RunData.qU(obj.exclDataStart)-18600,obj.RunData.qU(end)-18500,10);
            plot(x,zeros(numel(x),1),'-','Color',rgb('DarkSlateGray'),'LineWidth',2);
            hold on;
            %pt = plot(x,obj.StackTolerance.*ones(numel(x),1),'--','Color',rgb('DarkSlateGray'),'LineWidth',2);
            %plot(x,-obj.StackTolerance.*ones(numel(x),1),'--','Color',rgb('DarkSlateGray'),'LineWidth',2);
            
            for rr = 1:obj.nRuns
                % Plot
                c = colormap(jet(obj.nRuns));
                p = scatter(obj.RunData.qU(obj.exclDataStart:end)-obj.ModelObj.Q_i,...
                    obj.SingleRunData.qUfrac(obj.exclDataStart:end,rr)-obj.RunData.qUfrac(obj.exclDataStart:end),...
                    80,c(rr,:),'o','filled','MarkerEdgeColor',rgb('DimGray'));
                PrettyFigureFormat;
                hold on;
                set(gca,'FontSize',20);
                    xlabel(sprintf('mean retarding potential < qU >  -  %.1f  (eV)',obj.ModelObj.Q_i),'FontSize',20);
                    ylabel('qUfrac  -  < qUfrac > ','FontSize',20);
                    legend_str{rr} = sprintf('%u',obj.RunList(rr));
                    %counter = counter +1;
            end
%             if strcmp(legendFlag,'ON')
 %               myleg= legend(legend_str{:},'Location','bestoutside');
 %                myleg.NumColumns = 3;
  %               myleg.Title.String = sprintf('%.0f runs',obj.nRuns);
%             else
%                % myleg = legend([pt],'qU tolerance');
%                 
%             end
    %        legend('boxoff');
    %        myleg.FontSize = 20;
            grid on;
            xlim(1.05.*[min(obj.RunData.qU(obj.exclDataStart:end)-obj.ModelObj.Q_i) max(obj.RunData.qU(obj.exclDataStart:end)-obj.ModelObj.Q_i)]);
            hold off;
            
%             if abs(min(ylim))<=obj.StackTolerance
%                 ylim([-max(ylim) ,max(ylim)]);
%             elseif abs(min(ylim))<=abs(max(ylim))
%                 ylim([-max(ylim) ,max(ylim)]);
%             elseif abs(max(ylim))<=abs(min(ylim))
%                 ylim([min(ylim) ,-min(ylim)]);
%             end
            if strcmp(saveplot,'ON')
                 if exist('./plots','dir')~=7
                    mkdir plots
                end
%                 if exist('./plots/png','dir')~=7
%                     mkdir plots/png
%                 end
                save_name = sprintf('./plots/ScatterPlot-StackedRuns_%u-%u_excl%u_qUfrac.png',min(obj.StackedRuns),max(obj.StackedRuns),obj.exclDataStart);
                if contains(save_name,' ')
                    save_name=strrep(save_name,' ','_');
                end
                %publish_figurePDF(figScatter,save_name);
                 print(figScatter,save_name,'-dpng','-r400');
            end
        end
        function PlotStackingData(obj,varargin)
            % DATA Plots to show effect of stacking Data
            % Only approximation here in Stacked Data: averaged qU
            % Plots:
            % 1. Data Stacked Count Rate & Single Run Count Rates
            % 2. Residuals
            %------------------------------------------------------------------

           
            % Plot
            fig123 = figure(123);
            set(fig123, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.8, 1.2]);
            s1 = subplot(2,1,1);
            for rr = 1:obj.nRuns
                %SingleRunData{rr}   = load([obj.RunData.matFilePath,num2str(obj.RunList(rr)),obj.mpix,obj.ringCutFlag,'.mat']);
                if ismember(obj.RunList(rr),obj.StackedRuns)
                    pruns = errorbar(obj.SingleRunData.qU(:,rr,1)-obj.ModelObj.Q_i,obj.SingleRunData.TBDIS(:,rr,1)./(obj.SingleRunData.TimeSec(rr)*obj.SingleRunData.qUfrac(:,rr,1)),...
                        sqrt(obj.SingleRunData.TBDIS(:,rr,1))./(obj.SingleRunData.TimeSec(rr)*obj.SingleRunData.qUfrac(:,rr,1)),...
                        'ko','MarkerSize',3,'MarkerEdgeColor' , [0 0 0]);
                    hold on
                elseif ismember(obj.RunList(rr),obj.NotStackedRuns)
                    fprintf(2,'MultiRunAnalysis:PlotStackingData: %g not stacked \n',obj.RunList(rr));
                end
            end
            pstack = errorbar(obj.RunData.qU-obj.ModelObj.Q_i, obj.RunData.TBDIS./(obj.RunData.TimeSec*obj.RunData.qUfrac),obj.RunData.TBDISE./(obj.RunData.TimeSec*obj.RunData.qUfrac),...
                'o','MarkerSize',4,'MarkerEdgeColor', rgb('Red'),'MarkerFaceColor', rgb('Red'),'Color',rgb('Red'));
            hold off;
            legSingleRuns = sprintf('%u Runs: %u - %u',numel(obj.StackedRuns),min(obj.StackedRuns),max(obj.StackedRuns));
            legStackRuns = sprintf('Stacked Runs');
            legend([pruns,pstack],legSingleRuns,legStackRuns);
            legend('boxoff');
            xlabel(sprintf('retarding potential - %.1f (eV)',obj.ModelObj.Q_i));
            ylabel('Count Rate (cps)');
            PrettyFigureFormat;
            title('Data Comparison: Stack - No Stack');
            set(gca,'FontSize',20);
            xlim([obj.RunData.qU(obj.exclDataStart)-obj.ModelObj.Q_i, obj.RunData.qU(end)-obj.ModelObj.Q_i]);
            
            %Residuals
            s2 = subplot(2,1,2);
            for r = 1:obj.nRuns
                if ismember(obj.RunList(r),obj.StackedRuns)
                    edata = errorbar(obj.SingleRunData.qU(:,r,1)-obj.ModelObj.Q_i,...
                        (obj.SingleRunData.TBDIS(:,r,1)./(obj.SingleRunData.TimeSec(r)*obj.SingleRunData.qUfrac(:,r,1))-obj.RunData.TBDIS./(obj.RunData.TimeSec*obj.RunData.qUfrac)),...
                        0*sqrt(obj.SingleRunData.TBDIS(:,r,1)./(obj.SingleRunData.TimeSec(r)*obj.SingleRunData.qUfrac(:,r,1))),...
                        'o','MarkerSize',3,'MarkerEdgeColor', [0 0 0]);
                    edata = boundedline(obj.SingleRunData.qU(:,r,1)-obj.ModelObj.Q_i,...
                        0.*obj.RunData.TBDISE./obj.RunData.TBDISE,...
                        obj.RunData.TBDISE./(obj.RunData.TimeSec*obj.RunData.qUfrac),'alpha','cmap','transparency', 0.2,rgb('CadetBlue'));
                    hold on;
                elseif ismember(obj.RunList(r),obj.NotStackedRuns)
                end
            end
            plot(obj.ModelObj.qU-obj.ModelObj.Q_i,zeros(obj.ModelObj.nqU,1),'--k');
            hold off;
            xlabel(sprintf('retarding potential - %.3f (eV)',obj.ModelObj.Q_i));
            ylabel('Residuals (cps)');
            legend(sprintf('%u Runs',numel(obj.StackedRuns)));
            PrettyFigureFormat;
            set(gca,'FontSize',20);
            xlim([obj.RunData.qU(obj.exclDataStart)-obj.ModelObj.Q_i, obj.RunData.qU(end)-obj.ModelObj.Q_i]);
            linkaxes([s1,s2],'x');
            publish_figurePDF(fig123,sprintf('./plots/StackingSpectrum_%s.pdf',obj.StackFileName));
                    
        end       
        function GetStackingData(obj,varargin)
            % Get Stacking Data - For test Only
            %------------------------------------------------------------------
            
            SingleRunData = cell(obj.nRuns,1);

            myqu       = [];
            mytbdis    = [];
            mytbdiserr = [];
            
            for rr = 1:obj.nRuns
                SingleRunData{rr}   = load([obj.RunData.matFilePath,num2str(obj.RunList(rr)),'.mat']);
                if ismember(obj.RunList(rr),obj.StackedRuns)
                    
                    myqu        = [myqu    SingleRunData{rr}.qU'];
                    mytbdis     = [mytbdis (SingleRunData{rr}.TBDIS./(SingleRunData{rr}.TimeSec*SingleRunData{rr}.qUfrac))'];
                    mytbdiserr  = [mytbdiserr (SingleRunData{rr}.TBDISE./(SingleRunData{rr}.TimeSec*SingleRunData{rr}.qUfrac))'];
                    
                elseif ismember(obj.RunList(rr),obj.NotStackedRuns)
                    fprintf(2,'MultiRunAnalysis:PlotStackingData: %g not stacked \n',obj.RunList(rr));
                end
            end
            
            figure(987)
            errorbar(myqu,mytbdis,mytbdiserr,'.','Color','red');
            set(gca,'yscale','log');
            
            MySubRun=6;
            mymin=(MySubRun-1)*numel(myqu)/numel(obj.ModelObj.qU);
            mymax=(MySubRun+1)*numel(myqu)/numel(obj.ModelObj.qU);
            p = polyfit(myqu(mymin:mymax),mytbdis(mymin:mymax),3);
            %y = @(xq) interp1(myqu(mymin:mymax),mytbdis(mymin:mymax),xq,'linear')
            x = linspace(myqu(mymin),myqu(mymax),1000);
            y = polyval(p,x);
            hold on
            %plot(x,y(x),'LineWidth',2);
            plot(x,y,'LineWidth',2);
            hold off
        end        
        function PlotStackingModel(obj,varargin)
            p=inputParser;
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            saveplot = p.Results.saveplot;
            % Model Comparison:
            % Stacked Run Model (ModelObj) and Sum of Single Run Models
            % Approximation in ModelObject:
            % - averaged qU
            % - averaged column density, DT-fraction
            % Plots:
            % 1. Data Stacked Count Rate & Single Run Count Rates
            % 2. Residuals (stat + stacking systematic)
            %-----------------------------------------------------------------
            
            obj.SimulateStackRuns; %stack model back to init values
               
            % Load Single Run Models (Stacked Only)
            if isempty(obj.SingleRunObj)
                obj.LoadSingleRunObj;
            end
           
            % Apply same correction as on stacked data
            %StackingCorr = obj.ComputeStackingCorr_FT(obj.ModelObj.WGTS_MolFrac_DT);            
            StackingCorr  = 1;
            %StackingCorr    = [0.99308     0.99336      0.9933     0.99371     0.99356     0.99463     0.99435      0.9941     0.99423     0.99384      0.9938     0.99366     0.99503     0.99436     0.99477     0.99462     0.99526     0.99649     0.99525     0.99461     0.99535     0.99549     0.99513     0.99529     0.99548     0.99636     0.99653     0.99687     0.99684     0.99668     0.99719     0.99683     0.99655     0.99745     0.99601       0.997     0.99765     0.99749     0.99737     0.99766     0.99789     0.99829     0.99827     0.99861     0.99879     0.99965     0.99963      1.0003     0.99945     0.99964     0.99956     0.99961     0.99964     0.99953     0.99912           1      1.0004      1.0009      1.0045      1.0042      1.0048      1.0049      1.0053      1.0059      1.0061       1.007       1.008      1.0087      1.0083      1.0099      1.0095      1.0107      1.0116      1.0108      1.0119      1.0119       1.012      1.0122      0.9919     0.99115     0.99122     0.99239     0.99257     0.99305     0.99323     0.99367     0.99355     0.99493     0.99496     0.99574     0.99613     0.99668     0.99681     0.99753     0.99866      1.0011      1.0014      1.0013      1.0018      1.0031       1.004      1.0039      1.0043      1.0047      1.0057      1.0056      1.0068       1.007      1.0077      1.0089      1.0092      1.0103      1.0118       1.011      1.0118      1.0135      1.0143      1.0142      1.0167      1.0184      1.0197      1.0204      1.0216      1.0223     0.98823     0.99111     0.99292      0.9948     0.99523     0.99609     0.99756      0.9985      1.0041      1.0062      1.0081      1.0091      1.0103      1.0112      1.0122      1.0133     0.98633     0.98774       0.989     0.99059     0.99126     0.99256     0.99355     0.99451     0.99494     0.99795      1.0018     0.97723     0.97904      0.9803     0.98227     0.98434      0.9853     0.98748     0.98895     0.99016     0.99139     0.99255     0.99511     0.99705     0.99833     0.99945      1.0015      1.0035      1.0048      1.0055       1.008      1.0094      1.0113      1.0142      1.0158      1.0171      1.0189      1.0195     0.98731      0.9881     0.99046     0.99217      0.9942     0.99654     0.99739     0.99871      1.0001];

            BkgRate = mean(cell2mat(cellfun(@(x) x.TBDIS(x.qU>obj.ModelObj.Q_i)./(x.qUfrac(x.qU>obj.ModelObj.Q_i).*x.TimeSec),obj.SingleRunObj,'UniformOutput',0)'));
            Bkg = BkgRate.*cell2mat(cellfun(@(x) x.qUfrac.*x.TimeSec,obj.SingleRunObj,'UniformOutput',0)');
            TBDIS_tmp = cell2mat(cellfun(@(x) x.TBDIS,obj.SingleRunObj,'UniformOutput',0)')-Bkg;
            TBDIS_tmp(obj.ModelObj.qU<obj.ModelObj.Q_i,:) = TBDIS_tmp(obj.ModelObj.qU<obj.ModelObj.Q_i,:).*StackingCorr; %WARNING
            TBDIS_tmp = TBDIS_tmp + Bkg;

            qUfracCorr   = obj.ModelObj.qUfrac./cell2mat(cellfun(@(x) x.qUfrac,obj.SingleRunObj,'UniformOutput',0)');
            TBDIS_tmp = TBDIS_tmp.*qUfracCorr; 
            TBDIS_Stack = sum(TBDIS_tmp,2);

            %Plot
            fig72 = figure('Renderer','opengl');
            set(fig72, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.7, 1]);
            subplot(2,1,1); % Overlay of Spectra
            plot(obj.ModelObj.qU-obj.ModelObj.Q_i,TBDIS_Stack,'--','Color', rgb('Orange'),'LineWidth',4);
            hold on;
            plot(obj.ModelObj.qU-obj.ModelObj.Q_i, obj.ModelObj.TBDIS,'-','Color',rgb('SteelBlue'),'LineWidth',2);
            PrettyFigureFormat;
            leg = legend('\Sigma single run models','stack model','Location','north'); legend boxoff
                    xlabel(sprintf('mean retarding potential -  %.1f  (eV)',obj.ModelObj.Q_i),'FontSize',20);
            ylabel('counts');
            xlim([obj.RunData.qU(obj.exclDataStart)-obj.ModelObj.Q_i, obj.RunData.qU(end)-obj.ModelObj.Q_i]);
            %title(obj.GetRunTitle)
            set(gca,'FontSize',20);
            leg.FontSize = 20;
            
            subplot(2,1,2); % Residuals
            obj.FitCM_Obj.ComputeCM('SysEffect',struct('Stack','ON'))
            Sigma = sqrt(TBDIS_Stack + (diag(obj.FitCM_Obj.CovMat)));
            
            plot(obj.ModelObj.qU-obj.ModelObj.Q_i,(TBDIS_Stack-obj.ModelObj.TBDIS)./Sigma,...
                '-','Color',rgb('SteelBlue'),'LineWidth',4);
            hold on;
            plot(obj.ModelObj.qU-obj.ModelObj.Q_i, zeros(obj.ModelObj.nqU,1),'k--','LineWidth',2);
            PrettyFigureFormat;
            leg = legend('\Sigma single run models - stack model'); legend boxoff
            leg.Location = 'north';
                    xlabel(sprintf('mean retarding potential -  %.1f  (eV)',obj.ModelObj.Q_i),'FontSize',20);
            ylabel('norm. residuals');
            xlim([obj.RunData.qU(obj.exclDataStart)-obj.ModelObj.Q_i, obj.RunData.qU(end)-obj.ModelObj.Q_i]);
            set(gca,'FontSize',20);
            leg.FontSize = 20;
            if strcmp(saveplot,'ON')
                plotname = sprintf('./plots/StackingModel_%s.png',obj.GetRunTitle);
                print(fig72,plotname,'-dpng','-r400');
            end
        end
        function Runtitle = GetRunTitle(obj)
            %Nice title for plots for many runs
            if obj.nRuns == 1
                Runtitle = ['Run',num2str(obj.RunList)];
            elseif obj.nRuns >1
                Runtitle = [num2str(length(obj.StackedRuns)),'Runs-',num2str(obj.StackedRuns(1)),'-',num2str(obj.StackedRuns(end))];
            end
        end
        
        function PlotSCdata_RhoD(obj,varargin)
            % Plot Slow Control Data & Stack Value
            % 1) RhoD Versus Run
            % 2) RhoD Distribution - Stacked Value - Cut Values
            %------------------------------------------------------------------
            
            % Plot
            fig10100 = figure(10100);
            set(fig10100, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.95, 0.7]);
            
            % Recover Data
            RhoD_StackedRuns    = obj.SingleRunData.WGTS_CD_MolPerCm2(obj.SingleRunData.Select_all);
            RhoD_StackedRunsErr = std(obj.SingleRunData.WGTS_CD_MolPerCm2_SubRun(:,obj.SingleRunData.Select_all));
            RhoD_NotStackedRuns = obj.SingleRunData.WGTS_CD_MolPerCm2(:,obj.SingleRunData.Select_all==0);

            % Plot Data
            subplot(1,4,[1 3])
            hstack=plot(obj.SingleRunData.StartTimeStamp(obj.SingleRunData.Select_all)',...
                RhoD_StackedRuns,'ks-.','MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
            hold on                
            hold off

            xlabel('Time');
            ylabel('\rho d (mol/cm^2)');
            legend([hstack],sprintf('<\\rho d> = %.3g mol/cm^2',obj.RunData.WGTS_CD_MolPerCm2));
            legend boxoff;
            title(sprintf('KATRIN WGTS Column Density (%s)',obj.StackFileName));
            %datetick('x', 'dd-mm-yyyy');
            xtickangle(45);
            PrettyFigureFormat;set(gca,'FontSize',20);
            
            subplot(1,4,4)
%             fiths = histfit(RhoD_StackedRuns,10);
%             fiths(1).FaceColor = rgb('IndianRed');
%             fiths(1).LineStyle = '-';
            
            fiths = histogram(RhoD_StackedRuns,10,'FaceColor',rgb('IndianRed'),'EdgeColor',rgb('Black'),'LineWidth',2);
            
            hold on
            if ~isempty(RhoD_NotStackedRuns)
                hs = histogram(RhoD_NotStackedRuns);hns(1).FaceColor = rgb('Black') ;
            end
            hold off

            xlabel('\rho d');            
            xtickangle(45);
            ylabel('Runs');
%             if ~isempty(RhoD_NotStackedRuns)
%             legend([fiths hns],'Stacked','Not Stacked','Location','NorthWest'); legend boxoff;
%             else
%             legend([fiths],'Stacked','Location','NorthWest');legend boxoff;
%             end   
            legend([fiths],sprintf('%.0f runs',numel(obj.RunList)),'Location','NorthWest'); legend boxoff;
            PrettyFigureFormat;set(gca,'FontSize',20);
            
            % Save
            savefile=sprintf('plots/MRA_SC%0.fdata_RhoD.png',numel(obj.RunList));
            export_fig(gcf,savefile,'-q101','-m3');
            
            % Print Table
            Mymin  = min(RhoD_StackedRuns);
            Mymax  = max(RhoD_StackedRuns);
            Myav   = mean(RhoD_StackedRuns);
            Myst   = std(RhoD_StackedRuns);
            
            t = PrintTable('Column Density (mol/cm^2)');
            t.addRow('minimum',Mymin,'mol/cm^2');
            t.addRow('maximum',Mymax,'mol/cm^2');
            t.addRow('std',Myst,'mol/cm^2');
            t.addRow('average',Myav,'mol/cm^2');
            t.addRow('stacked',obj.RunData.WGTS_CD_MolPerCm2,'mol/cm^2');
            t.display;
            %t.HasHeader = true;
            t.Format = 'tex';
            t.Caption = sprintf('Column Density Statistics - %.0f runs',numel(obj.RunList)');
            t.print;
            
        end    
        function PlotSCdata_RhoDError(obj,varargin)
            % Plot Slow Control Data & Stack Value
            % 1) RhoD Standard Deviation over all subruns
            %------------------------------------------------------------------
            
            % Plot
            fig10100 = figure(10100);
            set(fig10100, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.95, 0.7]);
            
            % Recover Data
            RhoD_StackedRuns    = std(obj.SingleRunData.WGTS_CD_MolPerCm2_SubRun,1);
           
            % Plot Data
            subplot(1,4,[1 3])
            hstack=plot(obj.SingleRunData.StartTimeStamp(obj.SingleRunData.Select_all)',...
                RhoD_StackedRuns,'ks-.','MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
            hold on                
            hold off

            xlabel('Time');
            ylabel('\rho d subrun-wise std (mol/cm^2)');
            legend([hstack],sprintf('stacked \\rho d = %.3g mol/cm^2',obj.RunData.WGTS_CD_MolPerCm2));
            legend boxoff;
            title(sprintf('KATRIN WGTS Column Density subrun-wise standard deviation (%s)',obj.StackFileName));

            %datetick('x', 'dd-mm-yyyy');
            xtickangle(45);
            PrettyFigureFormat;set(gca,'FontSize',20);
            
            subplot(1,4,4)
%             fiths = histfit(RhoD_StackedRuns,10);
%             fiths(1).FaceColor = rgb('IndianRed');
%             fiths(1).LineStyle = '-';
             fiths = histogram(RhoD_StackedRuns,10,'FaceColor',rgb('IndianRed'),'EdgeColor',rgb('Black'),'LineWidth',2);

            
            xlabel('\rho d std');            
            xtickangle(45);
            ylabel('Runs');
            legend([fiths],sprintf('%.0f runs',numel(obj.RunList)),'Location','northeast'); legend boxoff;
            PrettyFigureFormat;set(gca,'FontSize',20);
            
            % Save
            savefile=sprintf('plots/MRA_SC%0.fdata_RhoDstd.png',numel(obj.RunList));
            export_fig(gcf,savefile,'-q101','-m3');
            
            % Print Table
            Mymin  = min(RhoD_StackedRuns);
            Mymax  = max(RhoD_StackedRuns);
            Myav   = mean(RhoD_StackedRuns);
            Myst   = std(RhoD_StackedRuns);
            
            t = PrintTable('Column Density (mol/cm^2)');
            t.addRow('minimum',Mymin,'mol/cm^2');
            t.addRow('maximum',Mymax,'mol/cm^2');
            t.addRow('std',Myst,'mol/cm^2');
            t.addRow('average',Myav,'mol/cm^2');
            t.addRow('stacked',obj.RunData.WGTS_CD_MolPerCm2,'mol/cm^2');
            t.display;
            %t.HasHeader = true;
            t.Format = 'tex';
            t.Caption = sprintf('Column Density Statistics - %.0f runs',numel(obj.RunList)');
            t.print;
        end    
        function PlotSCdata_Isotopologue(obj,varargin)
            % Plot Slow Control Data & Stack Value
            % 1) [TT] Versus Time
            % 2) [HT] Versus Time
            % 3) [DT] Versus Time
            %------------------------------------------------------------------
            
            p = inputParser;
            p.addParameter('Molecule','TT',@(x)ismember(x,{'TT','HT','DT'}));
            p.parse(varargin{:});
            Molecule        = p.Results.Molecule;

            % Plot
            fig10100 = figure(10100);
            set(fig10100, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.95, 0.7]);
            
            % Recover Data
            switch Molecule
                case 'TT'
                    XX_StackedRuns    = obj.SingleRunData.WGTS_MolFrac_TT(obj.SingleRunData.Select_all);
                    XX_NotStackedRuns = obj.SingleRunData.WGTS_MolFrac_TT(:,obj.SingleRunData.Select_all==0);
                case 'HT'
                    XX_StackedRuns    = obj.SingleRunData.WGTS_MolFrac_HT(obj.SingleRunData.Select_all);
                    XX_NotStackedRuns = obj.SingleRunData.WGTS_MolFrac_HT(:,obj.SingleRunData.Select_all==0);
                case 'DT'
                    XX_StackedRuns    = obj.SingleRunData.WGTS_MolFrac_DT(obj.SingleRunData.Select_all);
                    XX_NotStackedRuns = obj.SingleRunData.WGTS_MolFrac_DT(:,obj.SingleRunData.Select_all==0);
            end
            
          
            % Plot Data
            subplot(1,4,[1 3])
            hstack=plot(obj.SingleRunData.StartTimeStamp(obj.SingleRunData.Select_all)',...
                XX_StackedRuns,'ks-.','MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
            hold on                
            hold off

            xlabel('Time');
            ylabel(sprintf('%s fraction',Molecule));
            legend([hstack],sprintf('<%s> = %.3g ',Molecule,mean(XX_StackedRuns)),'Location','Best');
            legend boxoff;
            title(sprintf('KATRIN %s fraction in WGTS (%s)',Molecule,obj.StackFileName));

            xtickangle(45);
            PrettyFigureFormat;set(gca,'FontSize',20);
            
            subplot(1,4,4)
%             fiths = histfit(XX_StackedRuns,10);
%             fiths(1).FaceColor = rgb('IndianRed');
%             fiths(1).LineStyle = '-';
                         fiths = histogram(XX_StackedRuns,10,'FaceColor',rgb('IndianRed'),'EdgeColor',rgb('Black'),'LineWidth',2);

            hold on
            if ~isempty(XX_NotStackedRuns)
                hs = histogram(XX_NotStackedRuns);hns(1).FaceColor = rgb('Black') ;
            end
            hold off
            xlabel(sprintf('%s',Molecule));            
            xtickangle(45);
            ylabel('Runs');
            legend([fiths],sprintf('%.0f runs',numel(obj.RunList)),'Location','northoutside'); legend boxoff;
            PrettyFigureFormat;set(gca,'FontSize',20);
            
            % Print Table
            Mymin  = min(XX_StackedRuns);
            Mymax  = max(XX_StackedRuns);
            Myav   = mean(XX_StackedRuns);
            Myst   = std(XX_StackedRuns);
            
            t = PrintTable(sprintf('%s fraction',Molecule));
            t.addRow('minimum',Mymin,'(fraction)');
            t.addRow('maximum',Mymax,'(fraction)');
            t.addRow('std',Myst,'(fraction)');
            t.addRow('average',Myav,'(fraction)');
            switch Molecule
                case 'TT'
                    t.addRow('stacked',obj.RunData.WGTS_MolFrac_TT,'(fraction)');
                case 'HT'
                    t.addRow('stacked',obj.RunData.WGTS_MolFrac_HT,'(fraction)');
                case 'DT'
                    t.addRow('stacked',obj.RunData.WGTS_MolFrac_DT,'(fraction)');
            end
            t.display;
            %t.HasHeader = true;
            t.Format = 'tex';
            t.Caption = sprintf('%s fraction Statistics - %.0f runs',Molecule,numel(obj.RunList)');
            t.print;
            
            
            % Save
            savefile=sprintf('plots/MRA_SC%0.fdata_%s.png',numel(obj.RunList),Molecule);
            export_fig(gcf,savefile,'-q101','-m3');
        end     
        function PlotSCdata_TritiumPurity(obj,varargin)
            % Plot Slow Control Data & Stack Value
            %Tritium Purity
            %------------------------------------------------------------------
            
            p = inputParser;
            p.parse(varargin{:});

            % Plot
            fig10100 = figure(10100);
            set(fig10100, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.95, 0.7]);
            
            % Recover Data
            XX_StackedRuns    = obj.SingleRunData.WGTS_MolFrac_TT(obj.SingleRunData.Select_all) + ...
                0.5 * obj.SingleRunData.WGTS_MolFrac_HT(obj.SingleRunData.Select_all) + ...
                0.5 * obj.SingleRunData.WGTS_MolFrac_DT(obj.SingleRunData.Select_all);
            XX_NotStackedRuns = obj.SingleRunData.WGTS_MolFrac_TT(:,obj.SingleRunData.Select_all==0) + ...
                0.5 * obj.SingleRunData.WGTS_MolFrac_HT(:,obj.SingleRunData.Select_all==0) + ...
                0.5 * obj.SingleRunData.WGTS_MolFrac_DT(:,obj.SingleRunData.Select_all==0);
                      
            % Plot Data
            subplot(1,4,[1 3])
            hstack=plot(obj.SingleRunData.StartTimeStamp(obj.SingleRunData.Select_all)',...
                XX_StackedRuns,'ks-.','MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
            hold on                
            hold off

            xlabel('Time');
            ylabel(sprintf('Tritium Purity'));
            legend([hstack],sprintf('<%s> = %.3g ','\epsilon _T',mean(XX_StackedRuns)),'Location','Best');
            legend boxoff;
                                    title(sprintf('KATRIN Tritium Purity in WGTS (%s)',obj.StackFileName));

            xtickangle(45);
            PrettyFigureFormat;set(gca,'FontSize',20);
            
            subplot(1,4,4)
%             fiths = histfit(XX_StackedRuns,10);
%             fiths(1).FaceColor = rgb('IndianRed');
%             fiths(1).LineStyle = '-';
                                     fiths = histogram(XX_StackedRuns,10,'FaceColor',rgb('IndianRed'),'EdgeColor',rgb('Black'),'LineWidth',2);

            hold on
            if ~isempty(XX_NotStackedRuns)
                hs = histogram(XX_NotStackedRuns);hns(1).FaceColor = rgb('Black') ;
            end
            hold off
            xlabel(sprintf('Tritium Purity'));            
            xtickangle(45);
            ylabel('Runs');
            legend([fiths],sprintf('%.0f runs',numel(obj.RunList)),'Location','NorthWest'); legend boxoff;
            PrettyFigureFormat;set(gca,'FontSize',20);
            
            % Print Table
            Mymin  = min(XX_StackedRuns);
            Mymax  = max(XX_StackedRuns);
            Myav   = mean(XX_StackedRuns);
            Myst   = std(XX_StackedRuns);
            
            t = PrintTable(sprintf('Tritium Purity'));
            t.addRow('minimum',Mymin,'(fraction)');
            t.addRow('maximum',Mymax,'(fraction)');
            t.addRow('std',Myst,'(fraction)');
            t.addRow('average',Myav,'(fraction)');
            t.addRow('stacked',obj.RunData.WGTS_MolFrac_TT+0.5*obj.RunData.WGTS_MolFrac_DT+0.5*obj.RunData.WGTS_MolFrac_HT,'(fraction)');
               
            t.display;
            %t.HasHeader = true;
            t.Format = 'tex';
            t.Caption = sprintf('Tritium Purity Statistics - %.0f runs',numel(obj.RunList)');
            t.print;
            
            
            % Save
            savefile=sprintf('plots/MRA_SC%0.fdata_%s.png',numel(obj.RunList),'TritiumPurity');
            export_fig(gcf,savefile,'-q101','-m3');
        end
        function PlotSCdata_TritiumActivity(obj,varargin)
            % Plot Slow Control Data & Stack Value
            % Tritium Activity
            %------------------------------------------------------------------
            
            p = inputParser;
            p.parse(varargin{:});

            % Plot
            fig10100 = figure(10100);
            set(fig10100, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.95, 0.7]);
            
            % Recover Data
            XX_StackedRuns    = obj.SingleRunData.WGTS_MolFrac_TT(obj.SingleRunData.Select_all) + ...
                0.5 * obj.SingleRunData.WGTS_MolFrac_HT(obj.SingleRunData.Select_all) + ...
                0.5 * obj.SingleRunData.WGTS_MolFrac_DT(obj.SingleRunData.Select_all);
            XX_StackedRuns    = XX_StackedRuns .* obj.SingleRunData.WGTS_CD_MolPerCm2(obj.SingleRunData.Select_all);
            XX_StackedRuns    = XX_StackedRuns .* obj.ModelObj.ISXsection * 1e4;        
            XX_StackedRuns    = XX_StackedRuns/mean(XX_StackedRuns);
            
            XX_NotStackedRuns = obj.SingleRunData.WGTS_MolFrac_TT(:,obj.SingleRunData.Select_all==0) + ...
                0.5 * obj.SingleRunData.WGTS_MolFrac_HT(:,obj.SingleRunData.Select_all==0) + ...
                0.5 * obj.SingleRunData.WGTS_MolFrac_DT(:,obj.SingleRunData.Select_all==0);
            XX_NotStackedRuns    = XX_NotStackedRuns .* obj.SingleRunData.WGTS_CD_MolPerCm2(:,obj.SingleRunData.Select_all==0);
            XX_NotStackedRuns    = XX_NotStackedRuns .* obj.ModelObj.ISXsection * 1e4;
            XX_NotStackedRuns    = XX_NotStackedRuns/mean(XX_NotStackedRuns);
            
            % Plot Data
            subplot(1,4,[1 3])
            hstack=plot(obj.SingleRunData.StartTimeStamp(obj.SingleRunData.Select_all)',...
                XX_StackedRuns,'ks-.','MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
            hold on                
            hold off

            xlabel('Time');
            ylabel(sprintf('Tritium Activity'));
            legend([hstack],sprintf('<%s> = %.3g ','Activity Factor',mean(XX_StackedRuns)),'Location','Best');
            legend boxoff;
            title(sprintf('KATRIN Tritium Activity in WGTS (%s)',obj.StackFileName));

            xtickangle(45);
            PrettyFigureFormat;set(gca,'FontSize',20);
            
            subplot(1,4,4)
%             fiths = histfit(XX_StackedRuns,10);
%             fiths(1).FaceColor = rgb('IndianRed');
%             fiths(1).LineStyle = '-';
                                     fiths = histogram(XX_StackedRuns,10,'FaceColor',rgb('IndianRed'),'EdgeColor',rgb('Black'),'LineWidth',2);

            hold on
            if ~isempty(XX_NotStackedRuns)
                hs = histogram(XX_NotStackedRuns);hns(1).FaceColor = rgb('Black') ;
            end
            hold off
            xlabel(sprintf('Tritium Activity'));            
            xtickangle(45);
            ylabel('Runs');
            legend([fiths],sprintf('%.0f runs',numel(obj.RunList)),'Location','NorthWest'); legend boxoff;
            PrettyFigureFormat;set(gca,'FontSize',20);
            
            % Print Table
            Mymin  = min(XX_StackedRuns);
            Mymax  = max(XX_StackedRuns);
            Myav   = mean(XX_StackedRuns);
            Myst   = std(XX_StackedRuns);
            
            t = PrintTable(sprintf('Tritium Activity'));
            t.addRow('minimum',Mymin,'(fraction)');
            t.addRow('maximum',Mymax,'(fraction)');
            t.addRow('std',Myst,'(fraction)');
            t.addRow('average',Myav,'(fraction)');
               
            t.display;
            %t.HasHeader = true;
            t.Format = 'tex';
            t.Caption = sprintf('Tritium Activity Statistics - %.0f runs',numel(obj.RunList)');
            t.print;
            
            
            % Save
            savefile=sprintf('plots/MRA_SC%0.fdata_%s.png',numel(obj.RunList),'TritiumActivity');
            export_fig(gcf,savefile,'-q101','-m3');
        end
        function PlotSCdata_IsotopologueError(obj,varargin)
            % Plot Slow Control Data & Stack Value
            % 1) [TT] Error Versus Time
            % 2) [HT] Error Versus Time
            % 3) [DT] v Versus Time
            %------------------------------------------------------------------
            
            p = inputParser;
            p.addParameter('Molecule','TT',@(x)ismember(x,{'TT','HT','DT'}));
            p.parse(varargin{:});
            Molecule        = p.Results.Molecule;

            % Plot
            fig10100 = figure(10100);
            set(fig10100, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.95, 0.7]);
            
            % Recover Data
            switch Molecule
                case 'TT'
                    XX_StackedRuns    = obj.SingleRunData.WGTS_MolFrac_TT_SubRun_error_mean(obj.SingleRunData.Select_all)./...
                        obj.SingleRunData.WGTS_MolFrac_TT(obj.SingleRunData.Select_all)*100;
                    XX_NotStackedRuns = obj.SingleRunData.WGTS_MolFrac_TT_SubRun_error_mean(:,obj.SingleRunData.Select_all==0)./...
                        obj.SingleRunData.WGTS_MolFrac_TT(:,obj.SingleRunData.Select_all==0)*100;
                case 'HT'
                    XX_StackedRuns    = obj.SingleRunData.WGTS_MolFrac_HT_SubRun_error_mean(obj.SingleRunData.Select_all)./...
                        obj.SingleRunData.WGTS_MolFrac_HT(obj.SingleRunData.Select_all)*100;
                    XX_NotStackedRuns = obj.SingleRunData.WGTS_MolFrac_HT_SubRun_error_mean(:,obj.SingleRunData.Select_all==0)./...
                        obj.SingleRunData.WGTS_MolFrac_HT(:,obj.SingleRunData.Select_all==0)*100;
                case 'DT'
                    XX_StackedRuns    = obj.SingleRunData.WGTS_MolFrac_DT_SubRun_error_mean(obj.SingleRunData.Select_all)./...
                        obj.SingleRunData.WGTS_MolFrac_DT(obj.SingleRunData.Select_all)*100;
                    XX_NotStackedRuns = obj.SingleRunData.WGTS_MolFrac_DT_SubRun_error_mean(:,obj.SingleRunData.Select_all==0)./...
                        obj.SingleRunData.WGTS_MolFrac_DT(:,obj.SingleRunData.Select_all==0)*100;
            end
            
          
            % Plot Data
            subplot(1,4,[1 3])
            hstack=plot(obj.SingleRunData.StartTimeStamp(obj.SingleRunData.Select_all)',...
                XX_StackedRuns,'ks-.','MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
            hold on                
            hold off

            xlabel('Time');
            ylabel(sprintf('%s fraction uncertainty (%%)',Molecule));
            legend([hstack],sprintf('<%s> = %.3g percent',Molecule,nanmean(XX_StackedRuns)),'Location','Best');
            legend boxoff;
            title(sprintf('%s fraction uncertainty in relative (%%)',Molecule));
            title(sprintf('KATRIN WGTS %s fraction: relative uncertainty (%s)',Molecule, obj.StackFileName));
            xtickangle(45);
            PrettyFigureFormat;set(gca,'FontSize',20);
            
            subplot(1,4,4)
%             fiths = histfit(XX_StackedRuns,10);
%             fiths(1).FaceColor = rgb('IndianRed');
%             fiths(1).LineStyle = '-';
                                     fiths = histogram(XX_StackedRuns,10,'FaceColor',rgb('IndianRed'),'EdgeColor',rgb('Black'),'LineWidth',2);

            hold on
            if ~isempty(XX_NotStackedRuns)
                hs = histogram(XX_NotStackedRuns);hns(1).FaceColor = rgb('Black') ;
            end
            hold off
            xlabel(sprintf('\\sigma _{%s} (%%)',Molecule));            
            xtickangle(45);
            ylabel('Runs');
            legend([fiths],sprintf('%.0f runs',numel(obj.RunList)),'Location','NorthWest'); legend boxoff;
            PrettyFigureFormat;set(gca,'FontSize',20);
            
            % Print Table
            Mymin  = min(XX_StackedRuns);
            Mymax  = max(XX_StackedRuns);
            Myav   = nanmean(XX_StackedRuns);
            Myst   = std(XX_StackedRuns);
            MyNan  = obj.RunList(find(isnan(XX_StackedRuns)));
            
            t = PrintTable(sprintf('%s fraction uncertainty',Molecule));
            t.addRow('NaN',MyNan,'(Run Number)');
            t.addRow('minimum',Mymin,'percent');
            t.addRow('maximum',Mymax,'percent');
            t.addRow('std',Myst,'percent');
            t.addRow('average',Myav,'percent');
            t.display;
            %t.HasHeader = true;
            t.Format = 'tex';
            t.Caption = sprintf('%s fraction uncertainty statistics in percent - %.0f runs',Molecule,numel(obj.RunList)');
            t.print;
            
            
            % Save
            savefile=sprintf('plots/MRA_SC%0.fdata_%s_uncertainty.png',numel(obj.RunList),Molecule);
            export_fig(gcf,savefile,'-q101','-m3');
        end  
        function Display_AllInfos(obj)
            
            if exist('./info','dir')~=7
                mkdir ./info
            end
            obj.ModelObj.DisplayTDBInfo('Output','file','filename','TBDinfo.m');
            obj.ModelObj.DisplayWGTSMACEInfo('Output','file','filename','WGTSMACEinfo.m');
            obj.FitCM_Obj.DisplayCMInfo('Output','file','filename','COVMATinfo.m');
            options.format    = 'pdf';
            options.outputDir = './info';
            options.evalCode  = false;
            options.catchError = false;
            
            publish('./TBDinfo.m',options);
            publish('./WGTSMACEinfo.m',options);
            publish('./COVMATinfo.m',options);
            
            !rm ./TBDinfo.m
            !rm ./WGTSMACEinfo.m
            !rm ./COVMATinfo.m
            
        end    
        function ColumnDensityDrift(obj,varargin)
            % Compute the velocity drift of the column density
            %
            p=inputParser;
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('RunsSmoothed',1);
            p.addParameter('FirstRun',min(obj.RunList));
            p.addParameter('LastRun',max(obj.RunList));
            p.parse(varargin{:});
            saveplot     = p.Results.saveplot;
            RunsSmoothed = p.Results.RunsSmoothed;
            FirstRun     = p.Results.FirstRun;
            LastRun      = p.Results.LastRun;
            
            % Build Run List Considered
            FirstRunIndex = find(obj.SingleRunData.Runs==FirstRun);
            LastRunIndex  = find(obj.SingleRunData.Runs==LastRun);
            WGTS_CD_MolPerCm2 = obj.SingleRunData.WGTS_CD_MolPerCm2(FirstRunIndex:LastRunIndex);
            StartTimeStamp    = obj.SingleRunData.StartTimeStamp(FirstRunIndex:LastRunIndex);
            % Convert StartTimeStamp
            StartTimeStampF = datetime(StartTimeStamp,...
                'ConvertFrom','posixTime',...
                'TimeZone','Europe/Zurich',...
                'Format','dd-MMM-yyyy HH:mm:ss');
            
            % Compute Drift Velocity
            driftVelocityPerDay = zeros(1,numel(WGTS_CD_MolPerCm2));
            %SmoothedRhoD=smooth(WGTS_CD_MolPerCm2,1);
            SmoothedRhoD= WGTS_CD_MolPerCm2;
            for i=1:numel(WGTS_CD_MolPerCm2)
                if i>RunsSmoothed
                    % Compute Drift Velocity
                    deltaRhoD = (SmoothedRhoD(i)-SmoothedRhoD(i-RunsSmoothed));
                    deltaT    = double((obj.SingleRunData.StartTimeStamp(i)-obj.SingleRunData.StartTimeStamp(i-RunsSmoothed)));
                    driftVelocityPerDay(i) = deltaRhoD./deltaT*86400;
                    %driftVelocityPerDayPerCent(i) = driftVelocityPerDay(i)./ SmoothedRhoD(i) * 100;
                end
            end
            % Select Top-up periods
            driftVelocityPerDay(driftVelocityPerDay>1e14)=NaN;
            
            % Plot
            fig = figure('Renderer','opengl');
            set(fig,'units','normalized','pos',[0.1, 0.1,1,0.8]);
            s1=subplot(2,1,1)
            plot(StartTimeStampF,WGTS_CD_MolPerCm2,...
                'ks','MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
            xlabel('Start Run Date and time');
            ylabel('Column Density (mol/cm^2)');
            grid on
            title('Column Density Drift Velocity');
            PrettyFigureFormat;   set(gca,'FontSize',20);
            s2=subplot(2,1,2)
            plot(StartTimeStampF(1:end),smooth(driftVelocityPerDay,RunsSmoothed),...
                'ks-.','MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
            %
            %            idx = isnan(driftVelocityPerDay);
            %            P = polyfit(StartTimeStamp(~idx),driftVelocityPerDay(~idx),1);
            %            yfit = P(1)*StartTimeStamp+P(2);
            %            hold on;
            %            plot(StartTimeStampF,yfit,'r-.','LineWidth',2,'Color',rgb('SteelBlue'));
            %            hold off
            %
            ylim([min(driftVelocityPerDay) 2e14]); grid on;
            xlabel('Start Run Date and time');
            ylabel('Drift Velocity (mol/cm^2/day)');
            PrettyFigureFormat;  set(gca,'FontSize',20);
            linkaxes([s1,s2],'x');
            
            return;
            
            fig = figure('Renderer','opengl');
            set(fig,'units','normalized','pos',[0.2, 0.2,1,0.6]);
            ss1=subplot(2,1,1)
            y = detrend(smooth(driftVelocityPerDay));
            %y = sin(2*pi.*obj.SingleRunData.StartTimeStamp/86400/2);
            plot(StartTimeStampF,y,...
                'ks-.','MarkerSize',12,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
            xlabel('Start Run Date and time');
            ylabel('Detrend Signal');
            PrettyFigureFormat ; set(gca,'FontSize',20);
            ss2=subplot(2,1,2)
            [f,P,prob] = lomb(StartTimeStamp',y',4,1);
            plot(1./f./86400,P,...
                'ks-.','MarkerSize',12,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
            xlabel('Period (d^{-1})');
            ylabel('Strength');
            [Pmax,jmax] = max(P);
            [Psort Pindex] = sort(P);
            
            for i=1:numel(Pindex)
                if prob(Pindex(i))<0.1
                    disp(['Most significant period is ',num2str(1/f(Pindex(i))/86400),...
                        ' days with FAP of ',num2str(prob(Pindex(i)))]);
                end
            end
            PrettyFigureFormat ; set(gca,'FontSize',20);
            xlim([0 3]);
        end    
        function PlotPISdistributions(obj,varargin)
            % Plot Inelastic Scattering Probability Distributions
            % For All Runs from actual run list
            %------------------------------------------------------------------
            
            p=inputParser;
            p.addParameter('saveplot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('LatexOutput','ON',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            saveplot     = p.Results.saveplot;
            LatexOutput  = p.Results.LatexOutput;
            
            if isempty(obj.SingleRunObj)
                obj.LoadSingleRunObj;
            end
            
            %% Loop over all KNM1 runs
            % Compute IS probabilities
            % Store Results in is Arrays
            save_name = sprintf('%s/inputs/WGTSMACE/WGTS_ISProb/RunWise_ISprobabilities_%.d_%s.mat',...
                getenv('SamakPath'),obj.ModelObj.NIS,obj.StackFileName);
            if exist(save_name,'file')
                sprintf('IS probabilities already computed - Loading from File \n')
                load(save_name);
            else
                fprintf('Compute Run-Wise IS probabilities for % Run List (%d scatterings) \n',...
                    obj.StackFileName, obj.ModelObj.NIS);
                for i=1:1:numel(obj.RunList)
                    fprintf('MultiRunAnalysis: PlotPISdistributions: Compute IS probability Run %d\n',...
                        obj.SingleRunData.Runs(i));
                    runwise_pis(i,:) = obj.SingleRunObj{i}.ComputeISProb;
                end
                save(save_name,'runwise_pis');
            end
            
            %% Get Mean IS probabilities from Run List
            mean_is = mean(runwise_pis,1);
            
            %% Get Stacked IS probabilities
            stacked_is = obj.ModelObj.ComputeISProb;
            
            %% Plot Distributions
            titlename = sprintf('%s Run-wise Inelastic Scattering Probability Distributions',...
                obj.StackFileName);
            figure('Units', 'pixels','Position', [0 0 1200 1200]);
            a=annotation('textbox', [0 0.91 1 0.1], ...
                'String', titlename, ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center');
            a.FontSize=20;a.FontWeight='bold';
            
            for i=1:8
                s(i) = subplot(4,2,i);
                h(i)=histogram(runwise_pis(:,i),'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.5);
                title(sprintf('%d-fold scattering',i-1));
                xlabel('%');
                ylabel('runs');
                
                hold on
                m(i)=line([mean_is(i) mean_is(i)],[0 max(h(i).Values)],'Color',rgb('Red'),'LineWidth',2);
                l(i)=line([stacked_is(i) stacked_is(i)],[0 max(h(i).Values)],'Color',rgb('SteelBlue'),'LineWidth',2);
                hold off
                
                leg = legend([h(i) l(i) m(i)],sprintf('\\sigma / mean: %.2g %% \n min= %.3g %% \n max= %.3g %%',...
                    100*std(runwise_pis(:,i)/mean(runwise_pis(:,i))),min(runwise_pis(:,i)),max(runwise_pis(:,i))),'stacked','mean','Location','northwest');
                leg.Color = 'none'; legend boxoff; leg.FontSize = 10;
                
                PrettyFigureFormat
            end
            
            switch saveplot
                case 'ON'
                    if exist('./plots','dir')~=7
                        mkdir plots
                    end
                    export_fig(gcf,sprintf('plots/%s_RunWiseISprobabilities_Mean_Stacked',obj.StackFileName),'-m3');
            end
            
            %% Print Latex Table
            switch LatexOutput
                case 'ON'
                    t = PrintTable('Inelastic Scattering Probabilities');
                    t.addRow('Number of Scattering','Mean Value','Stacked Value','Relative Difference (%)')
                    for i=1:8
                        t.addRow(sprintf('%d',i),...
                            sprintf('%5.6f %',mean_is(i)),sprintf('%5.6f %',stacked_is(i)),...
                            (mean_is(i)-stacked_is(i))/((mean_is(i)+stacked_is(i))/2)*100);
                    end
                    t.display;
                    t.Format = 'tex';
                    t.Caption = sprintf('%s Inelastic Scattering Probabilities - %.0f runs',obj.StackFileName,numel(obj.RunList)');
                    t.print;
            end
            
        end
        function PrintDataStatistics(obj,varargin)
            % Print Main Statistics of the Run List
            % being considered for the analysis
            % T. Lasserre
            % Last Modified: 10/07/2019
            
            switch obj.StackFileName
                case 'KNM1'
            MyqU1 = 2;
            MyqU2 = 14;
            MyqU3 = 35;
                case {'KNM2','Others','KNM2_after2ndFix'}
            MyqU1 = 1;
            MyqU2 = 11;
            MyqU3 = 33;
            end
            
            %%
t = PrintTable(sprintf('%s Statistics'),obj.StackFileName);
t.addRow('Parameter',sprintf('Run List : %s',obj.StackFileName));
t.addRow('Column Density',sprintf('%.2g \\rm{ \\, mol/cm}^2 (\\sigma=%.2g \\%%)',obj.RunData.WGTS_CD_MolPerCm2,std(obj.SingleRunData.WGTS_CD_MolPerCm2)/obj.RunData.WGTS_CD_MolPerCm2*100 ));
Tpuritylocal   = (obj.RunData.WGTS_MolFrac_TT+1/2*obj.RunData.WGTS_MolFrac_HT+1/2*obj.RunData.WGTS_MolFrac_DT) * 100;
TpuritylocalRW = (obj.SingleRunData.WGTS_MolFrac_TT+1/2*obj.SingleRunData.WGTS_MolFrac_HT+1/2*obj.SingleRunData.WGTS_MolFrac_DT) * 100;
t.addRow('tritium Purity',sprintf('%.2g \\%% (\\sigma=%.2g \\%%)',Tpuritylocal,std(TpuritylocalRW)/Tpuritylocal*100));
t.addRow('Number of Runs',sprintf('%.0f',numel(obj.RunList)));
t.addRow('Run Duration',sprintf('%.1f hours',hours(seconds(obj.RunData.TimeSec/numel(obj.RunList)))));
t.addRow('First Run',sprintf('%0.f - %s',obj.RunList(1),obj.SingleRunData.StartTimeStamp(1)));
t.addRow('Last Run',sprintf('%0.f - %s',obj.RunList(end),obj.SingleRunData.StartTimeStamp(end)));
localTimetotal=hours(obj.SingleRunData.StartTimeStamp(end)-obj.SingleRunData.StartTimeStamp(1));
t.addRow('Total Time',sprintf('%.1f hours',localTimetotal));
t.addRow('$\beta$-scan Time',sprintf('%.1f \\rm{ \\, hours} - (%.1f \\%%)',hours(seconds(obj.RunData.TimeSec)),hours(seconds(obj.RunData.TimeSec))/localTimetotal*100));

t.addRow('Signal+Background','');
locaTime=hours(seconds(sum(obj.RunData.TimeperSubRunperPixel(MyqU1:end,1),1)));
t.addRow(sprintf('At qU$>$%.1f eV',obj.ModelObj.qU(MyqU1)),...
         sprintf('%.1f \\rm{ \\, hours} - (%.1f \\%%) - %.2g \\rm{ \\, electrons}',locaTime, locaTime/localTimetotal*100,sum(obj.RunData.TBDIS(MyqU1:end))));
locaTime=hours(seconds(sum(obj.RunData.TimeperSubRunperPixel(MyqU2:end,1),1)));
t.addRow(sprintf('At qU$>$%.1f eV',obj.ModelObj.qU(MyqU2)),...
         sprintf('%.1f \\rm{ \\, hours} - (%.1f \\%%) - %.2g \\rm{ \\, electrons}',locaTime, locaTime/localTimetotal*100,sum(obj.RunData.TBDIS(MyqU2:end))));

     t.addRow('Only Background','');
locaTime=hours(seconds(sum(obj.RunData.TimeperSubRunperPixel(MyqU3:end,1),1)));
t.addRow(sprintf('At qU$>$%.1f eV',obj.ModelObj.qU(MyqU3)),...
         sprintf('%.1f \\rm{ \\, hours} - (%.1f \\%%) - %.2g \\rm{ \\, electrons}',locaTime, locaTime/localTimetotal*100,sum(obj.RunData.TBDIS(MyqU3:end))));

t.display;
t.HasHeader = true;
t.Format = 'tex';
t.print;                        
          end
        
    end
    
end