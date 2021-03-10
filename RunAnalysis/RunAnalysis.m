% Analysis class for the KATRIN Experiment (Run Specific)
% Perform Fit and Various Analysis on a single Run
%-------------------------------------------------------------------------%
% This class contains:
% -Import and conversion of KATRIN Tritium Data
% -Fit (different options) and display of fit results
% -Computation of Covariance Matrices
% -More specific analysis methods
%
% Mandatory input:
% RunNr
%-------------------------------------------------------------------------%
%       L.Schlueter        P. I. Morales Guzman      T. Lasserre
%         MPP/TUM             TUM / MPP                  CEA
%-------------------------------------------------------------------------%

classdef RunAnalysis < handle & matlab.mixin.Copyable
     properties (Access=public,NonCopyable=true)
         ModelObj;
     end
    properties (Access=public)
        DataType;       % Switch for analysing Real Data of Fake Data
        FakeRunType;    % Prefix for Fake run to be considered
        RunNr;          % Run Number
        DataSet;        % Information about Measurement period, filled automatically
        RunData;        % Structure with all information for runsummaries
    %    ModelObj;       % TBD Object
        SingleRunObj;   % cell with TBD Objects of singl runs
        AnaFlag;        % SinglePixel, MultiPixel or StackPixel, Ring
        ELossFlag;
        KTFFlag;
        DopplerEffectFlag; 
        FSDFlag         % final state distributions: Sibille, Sibille0p5eV, BlindingKNM1, OFF ...
        FSD_Sigma;
        SynchrotronFlag;
        AngularTFFlag;
        RadiativeFlag;
        ROIFlag;        % region of interest
        MosCorrFlag;    % correct data qU by monitor spectrometer drift 
        ISCSFlag;       % inel. scattering cross section flag
        BKG_PtSlope;
        
        %Covariance Matrices
        FitCM_Obj;      % Fit Covariance Matrix Object
        FitCM;          % Current (Combined Systematics) Fit Covariance Matrices
        FitCMFrac;      % Current (Combined Systematics) Fit Fractional Covariance Matrix
        FitCMShape;     % Current (Combined Systematics) Fit Shape Only Covariance Matrix
        FitCMFracShape; % Current (Combined Systematics) Fit Fractional Shape Only Covariance Matrix
        
        FitCMNorm;      % Norm Only Covariance Matrix
        SysBudget;      % choose among predefined set of systematic uncertainties
        
        % Fit Option
        chi2;           % chi2Stat, chi2CM, chi2CMFrac,chi2CMShape, chi2P, chi2CMStack
        minuitOpt;      % Minuit fitter String (ex: min ; minors)
        fitter;         % minuit, matlab
        fixPar;         % string of fixed Parameters: e.g. '1 4'
        nPar;           % number of available parameter
        exclDataStart;  % Start Fit at bin i (i=1 -> qU(min)=-1.7keV, i=9 ->qU(min)=200eV)
        exclDataStop;   % Stop Fit at bin i e.g. exclDataStart=5, exclDataStop=15 %  -> used Data points [5:15] 
        FitResult;      % struct of Fit results
        PixList;        % list of pixels to be analyzed e.g. [3:7,9:136] 
        pulls;          % pulls
        pullFlag;       % 1 = nu mass pull, 2 = nu mass + FSD pull, 3 = FSD pull
        NonPoissonScaleFactor; % Non poissonian background fluctuation factor - Horizontal Array if ring wise
        BkgFluct;       % Type of background fluctuations: poissonian or not
        i_mnu;          % Fit Parameter Init: neutrino mass squared
        i_Q;            % Fit Parameter Init: endpoint
        i_B;            % Fit Parameter Init: background
        i_N;            % Fit Parameter Init: normalization
        i_DTGS          % Fit Parameter Init: FSD GS
        i_DTES          % Fit Parameter Init: FSD ES
        i_HTGS          % Fit Parameter Init: FSD GS
        i_HTES          % Fit Parameter Init: FSD ES
        i_TTGS          % Fit Parameter Init: FSD GS
        i_TTES          % Fit Parameter Init: FSD ES
        i_qUOffset      % Fit Parameter Init: qU Offset per ring
        i_mTSq          % Fit Parameter Init: tachyonic neutrino mass (ringwise)
        CatsResult;     % Results for CATS diagnostics
        RingMerge;      % How you choose to merge the rings

        % Plot Option
        PlotColor;      % Color of Plots, different for Twins & Data
        PlotColorLight; % Color of Plots, different for Twins & Data
        ErrorBarScaling;% Scaling of Error Bars for Plots only 
        
        % rhoD scan
        CDmin          % Column density at minimum chi2 after performing the rhoD scan
        rhoDlowerUnc;  % Lower bound of uncertainty on column density
        rhoDupperUnc;  % Upper bound of uncertainty on column density
        
        % Pixel/Ring-wise analysis
        RingList;      % list of rings to be analyzed iteratively
        AllFitResults; % Saves all fit results from rings or pixels (not an input)
        nRings;        % number of rings to analyze (not an input)
        RingPixList;   % cell with pixels per ring  (not an input)
        
        % Data Efficiency Correction - WARNING: Correction Applied to DATA
        DataEffCorr;   %  ROI+Pile-Up qU/Rate dependent efficiency
                       % 'OFF'                         = Take uncorrected Data
                       % 'ROI', 'PileUp','ROI+PileUp'  = Take uncorrected Data and correct it in Samak
                       % 'RunSummary'                  = Take corrected Data directly from RunSummary      

        % Class handling options
        DoRunAnalysisConstructor;
        
       % Monte Carlo Twin options
        TwinBias_WGTS_CD_MolPerCm2; % relative (%) with respet to real data
        TwinBias_WGTS_MolFrac_TT;   % relative (%)
        TwinBias_WGTS_MolFrac_HT,   % relative (%)
        TwinBias_WGTS_MolFrac_DT;   % relative (%)
        TwinBias_qU;                % absolute (eV) shift
        TwinBias_qUfrac;
        TwinBias_Time;              % absolute (s)
        TwinBias_Bkg;               % relative (%) can be ringwise or scalar
        TwinBias_mnuSq;             % absolute value for neutrino mass
        TwinBias_Q;                 % absolute value for endpoint or 'Fit' -> take fit value      
        FitNBFlag;                  % use normlization and background from fit
        TwinBias_FSDSigma;          % broadening of fsd in eV
        TwinBias_BKG_PtSlope;      % subrun-wise background slope from penning trap
        TwinFakeLabel;              % for labeling twin or fake runs with extra info: e.g. qU-bias,...
         
        %Fake MC option
        FakeInitFile %name of study -> Init file
        % 
    end
    methods % Constructor
        function obj = RunAnalysis(varargin)
            cprintf('blue','---------------------Start RunAnalysis Constructor------------------- \n')
            p = inputParser;
            p.addParameter('DataType','Real',@(x)ismember(x,{'Real','Twin','Fake','FitriumTwin','KafitTwin'}));
            p.addParameter('FakeRunType','Fake1');
            p.addParameter('RunNr',[],@(x)(isfloat(x) && x>0));
            p.addParameter('AnaFlag','StackPixel',@(x)ismember(x,{'StackPixel', 'SinglePixel', 'MultiPixel', 'Ring'}));
            p.addParameter('ELossFlag','',@(x)ismember(x,{'Aseev','Abdurashitov','CW_GLT','CW_G2LT','KatrinD2','KatrinT2','KatrinT2A20'}));%default given later
            p.addParameter('FSDFlag','Sibille0p5eV',@(x)ismember(x,{'SAENZ','BlindingKNM1','Sibille','Sibille0p5eV','OFF','SibilleFull','BlindingKNM2','KNM2','KNM2_0p5eV','KNM2_0p1eV'}));
            p.addParameter('FSD_Sigma',0,@(x)isfloat(x));
            p.addParameter('DopplerEffectFlag','OFF',@(x)ismember(x,{'OFF','FSD','FSD_Knm1'}));%default given later
            p.addParameter('ROIFlag','Default',@(x)ismember(x,{'Default','14keV'})); % default->default counts in RS, 14kev->[14,32]keV ROI
            p.addParameter('MosCorrFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('KTFFlag','WGTSMACE',@(x)ismember(x,{'WGTSMACE','MACE','WGTSMACE_NIS1'}));
            p.addParameter('SynchrotronFlag','ON',@(x)ismember(x,{'OFF','ON'}));
            p.addParameter('AngularTFFlag','ON',@(x)ismember(x,{'OFF','ON'}));
            p.addParameter('RadiativeFlag','ON',@(x)ismember(x,{'ON','OFF'}));
          
            p.addParameter('ISCSFlag','Edep',@(x)ismember(x,{'Aseev','Theory','Edep'}));
            % Fit Options
            p.addParameter('chi2','chi2Stat',@(x)ismember(x,{'chi2Stat', 'chi2CM', 'chi2CMFrac','chi2CMShape', 'chi2P','chi2Nfix'}));
            p.addParameter('fitter','minuit',@(x)ismember(x,{'minuit','matlab'}));
            p.addParameter('minuitOpt','min;minos',@(x)ischar(x));
            p.addParameter('exclDataStart',1,@(x)isfloat(x));
            p.addParameter('exclDataStop',9999,@(x)isfloat(x)); %9999== use all points from exclDataStart to above endpoint 
            p.addParameter('PixList',[],@(x)isfloat(x) && all(x)>0);
            p.addParameter('RingList',1:12,@(x)isfloat(x) && all(x>0));
            p.addParameter('fixPar','',@(x)ischar(x)); %default given in constructor: FSD and qUOffset fixed
            p.addParameter('pulls',[],@(x)isfloat(x) && all(x>0));
            p.addParameter('pullFlag',99);%@(x)ismember(x,{1,2,3}) 
            p.addParameter('RingMerge','None',@(x)ismember(x,{'Default','None','Full','Half','Azi','AziHalfNS','AziHalfEW'}));
            p.addParameter('NonPoissonScaleFactor',[],@(x)isfloat(x) && all(x>0)); 
            p.addParameter('BKG_PtSlope',0,@(x)isfloat(x));     
            p.addParameter('SysBudget','',@(x)isfloat(x)); %if none given-> defined according to data set

            % initialization of fitted parameters
            p.addParameter('i_mnu',0,@(x)isfloat(x));
            p.addParameter('i_Q',0,@(x)isfloat(x));
            p.addParameter('i_B',[],@(x)isfloat(x));
            p.addParameter('i_N',[],@(x)isfloat(x));
            p.addParameter('i_DTGS',0,@(x)isfloat(x));
            p.addParameter('i_DTES',0,@(x)isfloat(x));
            p.addParameter('i_HTGS',0,@(x)isfloat(x));
            p.addParameter('i_HTES',0,@(x)isfloat(x));
            p.addParameter('i_TTGS',0,@(x)isfloat(x));
            p.addParameter('i_TTES',0,@(x)isfloat(x));
            p.addParameter('i_qUOffset',[],@(x)isfloat(x));
            p.addParameter('i_mTSq',[],@(x)isfloat(x));
             
            % Monte Carlo Twin options
            p.addParameter('TwinBias_WGTS_CD_MolPerCm2',1,@(x)isfloat(x)); % relative (%)
            p.addParameter('TwinBias_WGTS_MolFrac_TT',1,@(x)isfloat(x));       % relative (%)
            p.addParameter('TwinBias_WGTS_MolFrac_HT',1,@(x)isfloat(x));       % relative (%)
            p.addParameter('TwinBias_WGTS_MolFrac_DT',1,@(x)isfloat(x));       % relative (%)
            p.addParameter('TwinBias_qU','',@(x)isfloat(x));                   % absolute shift (eV)
            p.addParameter('TwinBias_qUfrac','',@(x)isfloat(x));               % absolute (eV)
            p.addParameter('TwinBias_Time','',@(x)all(isfloat(x)));           % absolute 
            p.addParameter('TwinBias_Bkg',1,@(x)isfloat(x));                   % absolute (eV)
            p.addParameter('TwinBias_mnuSq',0,@(x)isfloat(x));                 % absolute (eV^2)
            p.addParameter('TwinBias_Q',18573.73,@(x)isfloat(x) || ismember(x,{'Fit'}));  % absolute (eV)
            p.addParameter('FitNBFlag','ON',@(x)ismember(x,{'ON','OFF','NormOnly'}));
            p.addParameter('TwinBias_FSDSigma',0,@(x)isfloat(x));                 % absolute (eV)
            p.addParameter('TwinBias_BKG_PtSlope',0,@(x)isfloat(x));
           
            %fake MC option
            p.addParameter('FakeInitFile','',@(x)isa(x,'function_handle') || isempty(x));
            % Efficiency Correction Applied to Data 
            p.addParameter('DataEffCorr','RunSummary',@(x)ismember(x,{'OFF', 'ROI', 'PileUp','ROI+PileUp','RunSummary'}));
 
            % Class handling options
            p.addParameter('DoRunAnalysisConstructor','ON',@(x)ismember(x,{'ON','OFF'}));
            
            p.parse(varargin{:});
     
            obj.DataType          = p.Results.DataType;
            obj.FakeRunType       = p.Results.FakeRunType;
            obj.RunNr             = p.Results.RunNr;
            obj.AnaFlag           = p.Results.AnaFlag;
            obj.chi2              = p.Results.chi2;
            obj.ELossFlag         = p.Results.ELossFlag;
            obj.FSDFlag           = p.Results.FSDFlag;
            obj.FSD_Sigma         = p.Results.FSD_Sigma;
            obj.DopplerEffectFlag = p.Results.DopplerEffectFlag;
            obj.ROIFlag           = p.Results.ROIFlag;
            obj.MosCorrFlag       = p.Results.MosCorrFlag;
            obj.KTFFlag           = p.Results.KTFFlag;
            obj.SynchrotronFlag   = p.Results.SynchrotronFlag;
            obj.AngularTFFlag     = p.Results.AngularTFFlag;
            obj.RadiativeFlag       = p.Results.RadiativeFlag;
            obj.ISCSFlag          = p.Results.ISCSFlag;
            obj.fitter            = p.Results.fitter;
            obj.minuitOpt         = p.Results.minuitOpt;
            obj.exclDataStart     = p.Results.exclDataStart;
            obj.exclDataStop      = p.Results.exclDataStop;
            obj.PixList           = p.Results.PixList;
            obj.RingList          = p.Results.RingList;
            obj.fixPar            = p.Results.fixPar;
            obj.pulls             = p.Results.pulls;
            obj.pullFlag          = p.Results.pullFlag;
            obj.NonPoissonScaleFactor = p.Results.NonPoissonScaleFactor;
            obj.i_mnu         = p.Results.i_mnu;
            obj.i_Q           = p.Results.i_Q;
            obj.i_B           = p.Results.i_B;
            obj.i_N           = p.Results.i_N;
            obj.i_DTGS        = p.Results.i_DTGS;
            obj.i_DTES        = p.Results.i_DTES;
            obj.i_HTGS        = p.Results.i_HTGS;
            obj.i_HTES        = p.Results.i_HTES;
            obj.i_TTGS        = p.Results.i_TTGS;
            obj.i_TTES        = p.Results.i_TTES;
            obj.i_qUOffset    = p.Results.i_qUOffset;
            obj.i_mTSq        = p.Results.i_mTSq;
            obj.DataEffCorr   = p.Results.DataEffCorr;
            obj.ELossFlag     = p.Results.ELossFlag;
            obj.RingMerge     = p.Results.RingMerge;
            obj.SysBudget     = p.Results.SysBudget;
            obj.BKG_PtSlope   = p.Results.BKG_PtSlope;
            
            % Monte Carlo Twin options
            obj.TwinBias_WGTS_CD_MolPerCm2  = p.Results.TwinBias_WGTS_CD_MolPerCm2; % relative (%) with respect to real data
            obj.TwinBias_WGTS_MolFrac_TT    = p.Results.TwinBias_WGTS_MolFrac_TT;   % relative (%)
            obj.TwinBias_WGTS_MolFrac_HT    = p.Results.TwinBias_WGTS_MolFrac_HT;   % relative (%)
            obj.TwinBias_WGTS_MolFrac_DT    = p.Results.TwinBias_WGTS_MolFrac_DT;   % relative (%)
            obj.TwinBias_qU                 = p.Results.TwinBias_qU;                % absolute (eV)
            obj.TwinBias_qUfrac             = p.Results.TwinBias_qUfrac;            % absolute (eV)
            obj.TwinBias_Time               = p.Results.TwinBias_Time;              % absolute (s)
            obj.TwinBias_Bkg                = p.Results.TwinBias_Bkg;               % relative (%)
            obj.FitNBFlag                   = p.Results.FitNBFlag;
            obj.TwinBias_mnuSq              = p.Results.TwinBias_mnuSq;
            obj.TwinBias_Q                  = p.Results.TwinBias_Q;
            obj.TwinBias_FSDSigma           = p.Results.TwinBias_FSDSigma;
            obj.TwinBias_BKG_PtSlope        = p.Results.TwinBias_BKG_PtSlope;
            
           if isempty(obj.TwinBias_qU) && strcmp(obj.AnaFlag,'Ring')
                obj.TwinBias_qU = zeros(1,numel(obj.RingList));
                obj.TwinBias_qUfrac = zeros(1,numel(obj.RingList));
            elseif isempty(obj.TwinBias_qU) && strcmp(obj.AnaFlag,'StackPixel')
                obj.TwinBias_qU     = 0; 
                obj.TwinBias_qUfrac = 0;
            end
                        
            obj.FakeInitFile             = p.Results.FakeInitFile;
            
            obj.DoRunAnalysisConstructor = p.Results.DoRunAnalysisConstructor;
            
            GetSamakPath; %sets current samak path as enviromental variable
            switch obj.DataType
                case {'Real','Twin','FitriumTwin','KafitTwin'}
                    obj.DataSet = GetDataSet(obj.RunNr);
                case 'Fake'
                    obj.DataSet = GetFakeDataSet(func2str(obj.FakeInitFile));
                %case {'FitriumTwin','KafitTwin'}
                %    obj.DataSet = 'Knm1';
            end
            
            obj.SetDefaultInit;
            
            % select pixels according to data set (if not specified otherwise)
            if isempty(obj.PixList) && ~isempty(obj.DataSet)
                obj.PixList = GetPixList(obj.DataSet);
            end
            
            % select pixels according to selected rings
            if strcmp(obj.AnaFlag,'Ring')
                switch obj.RingMerge
                    case 'Default'
                        [obj.PixList,obj.RingPixList] = Ring2PixelDefCombi(obj.RingList,obj.PixList);
                        obj.RingList = 1:10;
                    case 'None'
                        [obj.PixList,obj.RingPixList] = Ring2Pixel(obj.RingList,obj.PixList);
                        obj.RingList = 1:12;
                    case 'Full'
                        [obj.PixList,obj.RingPixList] = Ring2PixelCombi(obj.RingList,obj.PixList);
                        obj.RingList = 1:4;%1:numel();
                    case 'Half'
                        [obj.PixList,obj.RingPixList] = Ring2PixelHalfCombi(obj.RingList,obj.PixList);
                        obj.RingList = 1:2;
                    case 'Azi'
                        [obj.PixList,obj.RingPixList] = AziPatch2PixelCombi(obj.RingList,obj.PixList);
                        obj.RingList = 1:5;
                    case {'AziHalfNS','AziHalfEW'}
                        [obj.PixList,obj.RingPixList] = AziHalfPatch2PixelCombi(obj.RingList,obj.PixList,obj.RingMerge);
                        obj.RingList = 1:2;
                end
            end
            
            obj.nRings = numel(obj.RingPixList);
            
            %Init Analysis: Get Data, Model, Covariance Matrix
            switch obj.DoRunAnalysisConstructor
                case 'ON'
                    obj.ReadData;
                    obj.SetROI; 
                    obj.SimulateRun;
                    obj.InitFitPar;
                    obj.SetNPfactor;
                    
                    if all(~ismember(obj.DataEffCorr,{'OFF','RunSummary'}))
                        obj.ROIPileUpDataCorrection();
                    end
                    if all(~ismember(obj.chi2,{'chi2Stat','chi2P'}))
                        obj.ComputeCM('RecomputeFlag','OFF');
                    end
                    if strcmp(obj.chi2,'chi2Stat')
                        [StatCM, StatCMFrac] = obj.ComputeCM_StatPNP(varargin);
                        obj.FitCM = StatCM;            obj.FitCMShape = StatCM;
                        obj.FitCMFrac = StatCMFrac;    obj.FitCMFracShape = StatCMFrac;
                    end
                    obj.GetPlotColor;
                case 'OFF'
                    % Don't do it
            end
            
            cprintf('blue','---------------------End RunAnalysis Constructor------------------- \n');
        end
    end % Constructor
    
    methods % Data Import Methods
        function ReadData(obj) % Read Tritium Data From .mat File in Samak2.0/tritium-data/mat/
            %try
            switch obj.DataType
                case 'Real'
                    matFilePath = [getenv('SamakPath'),'tritium-data/mat/',obj.DataSet,'/'];
                    filename = [matFilePath,num2str(obj.RunNr),'.mat'];
                    if ~exist(filename,'file')
                        HDF5readallruns('h5runlist',obj.RunNr,'DataSet',GetDataSet(obj.RunNr));
                    end
                    obj.RunData = load(filename);
                    
                    if strcmp(obj.ROIFlag,'14keV')
                        obj.RunData.TBDIS ;
                    end
                case 'Twin'
                    matFilePath = [getenv('SamakPath'),'tritium-data/mat/','Twin',obj.DataSet,'/'];
                    matFileName = ['Twin',num2str(obj.RunNr),obj.SetTwinOrFakeFileName,'.mat'];
                    filename = [matFilePath,matFileName];
                    if ~exist(filename,'file')
                        obj.RunData. matFilePath = matFilePath;
                        obj.ComputeTwinRun;
                    end
                    obj.RunData = load(filename);
                case {'FitriumTwin','KafitTwin'}
                    MCdata = extractBefore(obj.DataType,'Twin');
                    if obj.TwinBias_mnuSq~=0
                        ExtraStrmNu = 'nonzero_';
                    else
                        ExtraStrmNu = '';
                    end
                    if obj.TwinBias_Time~=0
                        ExtraStrTime = sprintf('_%.0ftimes',obj.TwinBias_Time);
                    else
                        ExtraStrTime = '';
                    end
                    
                    matFilePath = [getenv('SamakPath'),sprintf('tritium-data/mat/%s_Twin%s/',obj.DataSet,MCdata)];
                    obj.RunData.matFileName = sprintf('Twin%s_%s%.0f%s.mat',MCdata,ExtraStrmNu,obj.RunNr,ExtraStrTime);
              
                    filename = [matFilePath,obj.RunData.matFileName];
                    obj.RunData = load(filename);
                case 'Fake'
                    if contains(func2str(obj.FakeInitFile),'KNM1')
                        matFilePath = [getenv('SamakPath'),'tritium-data/mat/FakeKnm1/'];
                    elseif contains(func2str(obj.FakeInitFile),'KNM2')
                        matFilePath = [getenv('SamakPath'),'tritium-data/mat/FakeKnm2/'];
                    else
                        matFilePath = [getenv('SamakPath'),'tritium-data/mat/Fake/'];
                    end
                    matFileName = [extractAfter(func2str(obj.FakeInitFile),'ref_'),obj.SetTwinOrFakeFileName,...
                          sprintf('_Run%.0f',obj.RunNr),'.mat'];
                    filename = [matFilePath,matFileName];
                    
                    if ~exist(filename,'file')
                        obj.RunData.matFilePath = matFilePath;
                        obj.ComputeTwinRun;
                    end
                    obj.RunData = load(filename);
            end
            
            obj.RunData.matFilePath   = matFilePath;
            obj.RunData.TimeSec       = mean(obj.RunData.TimeSec(obj.PixList));
            obj.RunData.TBDIS_Default = obj.RunData.TBDIS;
            obj.StackPixel;
            if isfield(obj.RunData,'TBDIS_RM')
                obj.RunData.TBDIS_RM_Default = obj.RunData.TBDIS_RM;
            end
            obj.RunData.RunName          = num2str(obj.RunNr);

            try 
                obj.RunData.StartTimeStamp = datetime(obj.RunData.StartTimeStamp, 'ConvertFrom', 'posixtime' );
            catch
            end

            obj.SetMosCorr;
            fprintf(2,'RunAnalysis:ReadData: Reading Tritium Data from File:\n')
            disp(obj.RunData);
            %             catch ME
            %                 if contains(ME.identifier,'couldNotReadFile')
            %                     fprintf('Error %s \n',ME.message);
            %                     fprintf(2,'RunAnalysis:ReadData: WARNING - mat file not existing! Reading HDF5 file instead \n')
            %                     try
            %                         if contains(obj.DataSet,'FirstTritium.katrin')
            %                             HDF5Reader('RunNr',obj.RunNr,'version','4a-fpd00');
            %                         elseif contains(obj.DataSet,'Knm1')
            %                             HDF5Reader('RunNr',obj.RunNr,'version','RunSummary-Durable2a-fpd00');
            %                         elseif contains(obj.DataSet,'Knm2')
            %                                 HDF5Reader('RunNr',obj.RunNr,'version','RunSummary-Prompt4b-fpd00');
            %                         end
            %                         obj.ReadData();
            %                     catch ME
            %                         error('ERROR:   HDF5 file not existing!')
            %                     end
            %                 else
            %                     fprintf('Error %s \n',ME.message);
            %                     return
            %                 end
            %end
        end
        function StackPixel(obj)
                % Stack Pixels
                switch obj.AnaFlag
                    case {'StackPixel'}
                        % Average Rest over Pixel
                        obj.RunData.qU            = mean(obj.RunData.qU(:,obj.PixList),2);
                        obj.RunData.qUfrac        = mean(obj.RunData.qUfrac(:,obj.PixList),2);
                        obj.RunData.EffCorr       = mean(obj.RunData.EffCorr(:,obj.PixList),2);
                        obj.RunData.TBDIS         = sum(obj.RunData.TBDIS(:,obj.PixList),2);
                        obj.RunData.TBDIS_Default = sum(obj.RunData.TBDIS_Default(:,obj.PixList),2);
                        obj.RunData.TBDISE        = sqrt(obj.RunData.TBDIS./obj.RunData.EffCorr); % statstical uncertainty
                        obj.RunData.MACE_Ba_T     = mean(obj.RunData.MACE_Ba_T(obj.PixList));
                        obj.RunData.MACE_Bmax_T   = mean(obj.RunData.MACE_Bmax_T(obj.PixList));
                        
                        if isfield(obj.RunData,'TBDIS_V')
                            obj.RunData.TBDIS_V = sum(obj.RunData.TBDIS_V(:,obj.PixList),2);
                        end
                        
                        if isfield(obj.RunData,'TBDIS_RM')
                            obj.RunData.TBDIS_RM   = mean(obj.RunData.TBDIS_RM(obj.PixList));
                            obj.RunData.TBDIS_RM_Default   = obj.RunData.TBDIS_RM;
                            obj.RunData.qU_RM      = mean(obj.RunData.qU_RM(obj.PixList));
                            obj.RunData.qUfrac_RM  = mean(obj.RunData.qUfrac_RM(obj.PixList));
                        end
                        
                        if isfield(obj.RunData,'TBDIS14keV')
                            obj.RunData.TBDIS14keV      = sum(obj.RunData.TBDIS14keV(:,obj.PixList),2);
                            obj.RunData.TBDIS14keV_RM   = mean(obj.RunData.TBDIS14keV_RM(obj.PixList));
                        end
                    case {'MultiPixel','SinglePixel'}
                        %Select pixellist
                        obj.RunData.TimeSec      = mean(obj.RunData.TimeSec(obj.PixList));
                        obj.RunData.qUfrac       = obj.RunData.qUfrac(:,obj.PixList);
                        obj.RunData.qU           = obj.RunData.qU(:,obj.PixList);
                        obj.RunData.EffCorr      = obj.RunData.EffCorr(:,obj.PixList);
                        obj.RunData.TBDIS        = obj.RunData.TBDIS(:,obj.PixList);
                        obj.RunData.TBDIS_Default= obj.RunData.TBDIS_Default(:,obj.PixList);
                        obj.RunData.TBDISE       = obj.RunData.TBDISE(:,obj.PixList);
                        obj.RunData.MACE_Ba_T    = obj.RunData.MACE_Ba_T(obj.PixList);
                        obj.RunData.MACE_Bmax_T  = obj.RunData.MACE_Bmax_T(obj.PixList);
                        
                        if isfield(obj.RunData,'TBDIS_V')
                            obj.RunData.TBDIS_V = obj.RunData.TBDIS_V(:,obj.PixList,:);
                        end
                        if isfield(obj.RunData,'TBDIS_RM')
                            obj.RunData.TBDIS_RM   = obj.RunData.TBDIS_RM(obj.PixList);
                            obj.RunData.qU_RM      = obj.RunData.qU_RM(obj.PixList);
                            obj.RunData.qUfrac_RM  = obj.RunData.qUfrac_RM(obj.PixList);
                        end
                        
                    case {'Ring'}
                        %Set time to average per ring
                        obj.RunData.TimeperSubRun = cell2mat(cellfun(@(x) mean(obj.RunData.TimeperSubRunperPixel(:,x),2),obj.RingPixList,'UniformOutput',false)');
                        if strcmp(obj.RingMerge,'None')
                            %obj.RunData.TimeperSubRun(:,~ismember(1:13,obj.RingList)) = []; %truncate
                        end
                        obj.RunData.TimeSec       = sum(obj.RunData.TimeperSubRun);
                        obj.RunData.qUfrac        = obj.RunData.TimeperSubRun./obj.RunData.TimeSec;% warning -> should not be 1, because of RM point
                        
                        % average or sum rest ringwise
                        obj.RunData.qU      = cell2mat(cellfun(@(x) mean(obj.RunData.qU(:,x),2),obj.RingPixList,'UniformOutput',false)');
                        obj.RunData.EffCorr = cell2mat(cellfun(@(x) mean(obj.RunData.EffCorr(:,x),2),obj.RingPixList,'UniformOutput',false)');
                        obj.RunData.TBDIS   = cell2mat(cellfun(@(x) sum(obj.RunData.TBDIS(:,x),2),obj.RingPixList,'UniformOutput',false)');
                        obj.RunData.TBDIS_Default   = cell2mat(cellfun(@(x) sum(obj.RunData.TBDIS_Default(:,x),2),obj.RingPixList,'UniformOutput',false)');
                        obj.RunData.TBDISE  = sqrt(obj.RunData.TBDIS./obj.RunData.EffCorr); % statstical uncertainty
                        obj.RunData.MACE_Ba_T = cell2mat(cellfun(@(x) mean(obj.RunData.MACE_Ba_T(x)),obj.RingPixList,'UniformOutput',false)');
                        obj.RunData.MACE_Bmax_T =  mean(obj.RunData.MACE_Bmax_T(obj.PixList));%cell2mat(cellfun(@(x) mean(obj.RunData.MACE_Bmax_T(x)),obj.RingPixList,'UniformOutput',false)');
                        
                        % Correction Thierry 1/4/2020
                        if isfield(obj.RunData,'TBDIS1ISCS4keV')
                            obj.RunData.TBDIS14keV      = cell2mat(cellfun(@(x) sum(obj.RunData.TBDIS14keV(:,x),2),obj.RingPixList,'UniformOutput',false)');
                        end
                        
                        % delete not used rings (otherwise problems with NaN)
                        if strcmp(obj.RingMerge,'None')
%                             obj.RunData.qU(:,~ismember(1:13,obj.RingList)) = [];
%                             obj.RunData.EffCorr(:,~ismember(1:13,obj.RingList)) = [];
%                             obj.RunData.TBDIS(:,~ismember(1:13,obj.RingList)) = [];
%                             obj.RunData.TBDISE(:,~ismember(1:13,obj.RingList)) = [];
%                             obj.RunData.MACE_Ba_T(~ismember(1:13,obj.RingList)) = [];
%                             obj.RunData.MACE_Bmax_T(~ismember(1:13,obj.RingList)) = [];
                        end
                        
                        if isfield(obj.RunData,'TBDIS_V')
                            obj.RunData.TBDIS_V = cell2mat(cellfun(@(x) sum(obj.RunData.TBDIS_V(:,:,x),3),obj.RingPixList,'UniformOutput',false)');
                        end
                        
                        if isfield(obj.RunData,'TBDIS_RM')
                            obj.RunData.qU_RM      = cell2mat(cellfun(@(x) mean(obj.RunData.qU_RM(x)),obj.RingPixList,'UniformOutput',false)');
                            obj.RunData.qUfrac_RM  = cell2mat(cellfun(@(x) mean(obj.RunData.qUfrac_RM(x)),obj.RingPixList,'UniformOutput',false)');
                            obj.RunData.TBDIS_RM   = cell2mat(cellfun(@(x) mean(obj.RunData.TBDIS_RM(x)),obj.RingPixList,'UniformOutput',false)');
                            % Correction Thierry 1/4/2020
                            obj.RunData.TBDIS14keV_RM   = cell2mat(cellfun(@(x) mean(obj.RunData.TBDIS14keV_RM(x)),obj.RingPixList,'UniformOutput',false)');
                        end
                end
        end
        function ReadFakeData(obj)
            % Read Tritium Data From .mat File
            % Data in Samak2.0/fake-data/mat/
            try
                [obj.RunData.matFilePath] = deal([getenv('SamakPath'),'/tritium-fakedata/mat/']);
                
                switch obj.AnaFlag
                    case 'StackPixel'
                        ringCutFlag = 'ex2'; %% TEMPORARY
                        file = [obj.RunData.matFilePath,obj.FakeRunType,'Run',num2str(obj.RunNr),ringCutFlag,'.mat'];
                        obj.RunData = load(file);
                    case 'SinglePixel'
                        return;
                    case {'MultiPixel','Ring'}
                        return;
                end
                
                fprintf(2,'RunAnalysis:ReadFakeData: Reading Fake Tritium Data from File:\n')
                disp(obj.RunData);
            catch
                error('ERROR: Fake MC data file not existing!')
            end
        end     
        function ROIPileUpDataCorrection(obj)
            % Correct Data for ROI/Pile-up ineficiencies
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
                    RoiCorr = 1./obj.ModelObj.qUEfficiencyCorrectionFactor(obj.RunData.qU);
                    obj.RunData.TBDIS = obj.RunData.TBDIS .* RoiCorr;
                    %Plot(obj.RunData.qU,RoiCorr);
                    fprintf(2,'RunAnalysis:ROIPileUpDataCorrection: Data - ROI Efficiency Correction Applied\n');
                case 'PileUp'
                    switch obj.ModelObj.FPD_Segmentation
                        case 'OFF'
                            nPix=numel(obj.RunData.PixelList);
                        case {'SINGLEPIXEL','MULTIPIXEL'}
                            nPix=1;
                        case 'RING'
                            if obj.ModelObj.FPD_Ring == 1
                                nPix=4;
                            else
                                nPix=12;
                            end
                    end
                    PileUpCorr = 1./obj.ModelObj.PileUpEfficiencyCorrectionFactor(obj.RunData.TBDIS./obj.RunData.qUfrac./obj.RunData.TimeSec/nPix);
                    obj.RunData.TBDIS = obj.RunData.TBDIS .* PileUpCorr;
                    %Plot(obj.RunData.qU,PileUpCorr);
                    fprintf(2,'RunAnalysis:ROIPileUpDataCorrection: Data - Pile-Up Efficiency Correction Applied\n');
                case 'ROI+PileUp'
                    RoiCorr = 1./obj.ModelObj.qUEfficiencyCorrectionFactor(obj.RunData.qU);
                    switch obj.ModelObj.FPD_Segmentation
                        case 'OFF'
                            nPix=numel(obj.RunData.PixelList);
                        case {'SINGLEPIXEL','MULTIPIXEL'}
                            nPix=1;
                        case 'RING'
                            if obj.ModelObj.FPD_Ring == 1
                                nPix=4;
                            else
                                nPix=12;
                            end
                    end
                    PileUpCorr = 1./obj.ModelObj.PileUpEfficiencyCorrectionFactor(obj.RunData.TBDIS./obj.RunData.qUfrac./obj.RunData.TimeSec/nPix);
                    obj.RunData.TBDIS = obj.RunData.TBDIS .* RoiCorr .* PileUpCorr;
                    %Plot(obj.RunData.qU,RoiCorr.*PileUpCorr);
                    fprintf(2,'RunAnalysis:ROIPileUpDataCorrection: Data - ROI+Pile-Up Efficiency Corrections Applied\n');
            end
        end
        
    end % Data Import Methods
    methods % Twin Monte Carlo Generator Methods
        function ComputeTwinRun(obj,varargin)
            p = inputParser;
            p.addParameter('StatFluct','ON',@(x)ismember(x,{'ON','OFF'}));   
            p.parse(varargin{:});
            
            StatFluct              = p.Results.StatFluct;

            MC = McRunGenerator('RunObj',obj);
            switch obj.DataType
                case 'Twin'
                    MC.ComputeTwinRun;
                case 'Fake'
                    MC.InitFile = obj.FakeInitFile;
                    MC.ComputeFakeRun;
            end
        end
        
        function [FitPar, FitErr, FitChi2min, dof,TBDIS]  = FitTwin(obj,varargin)
            p = inputParser;
            p.addParameter('nSamples',10,@(x)isfloat(x) & x>0);
            p.addParameter('CATS','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            nSamples = p.Results.nSamples;
            CATS     = p.Results.CATS;
        
           if strcmp(obj.DataType,'Real')
               fprintf('Data Type real - must be Twin \n')
               return
           end
           
           % randomize twins (add fluctiations)
           if strcmp(obj.chi2,'chi2Stat')
              [StatCM, StatCMFrac]         =  obj.ComputeCM_StatPNP;
               obj.FitCM     = StatCM;       obj.FitCMShape = StatCM;
               obj.FitCMFrac = StatCMFrac;   obj.FitCMFracShape = StatCMFrac;
               TBDIS = mvnrnd(obj.RunData.TBDIS',obj.FitCM,nSamples)';
           else
              if isempty(obj.FitCMShape)
                  obj.ComputeCM;
              end
               % Modif Thierry - 30/3/2020
               % TBDIS = mvnrnd(obj.RunData.TBDIS',obj.FitCMShape,nSamples)';
               TBDIS = mvnrnd(obj.RunData.TBDIS(:)',obj.FitCMShape,nSamples)';
               % End Modif Thierry - 30/3/2020
           end
           
            % init
            FitPar     = zeros(obj.nPar,nSamples);
            FitErr     = zeros(obj.nPar,nSamples);
            FitChi2min = zeros(nSamples,1);
            dof        = 0;
           
            progressbar('Fit twins with stat. fluctuations')
            
            for i=1:nSamples
                 progressbar(i/nSamples)
                 
                obj.RunData.TBDIS = TBDIS(:,i); 
                obj.RunData.TBDISE = sqrt(TBDIS(:,i));   
                
                if nSamples<2
                obj.Fit('CATS',CATS);
                obj.PlotFit();
                else
                obj.Fit;
                end                    
                FitPar(:,i)   = obj.FitResult.par;
                FitErr(:,i)   = obj.FitResult.err;
                FitChi2min(i) = obj.FitResult.chi2min;
                dof           = obj.FitResult.dof;
                
                obj.ModelObj.SetFitBias(0);
            end
            
            % reset to asimov data
            if isa(obj,'MultiRunAnalysis')
                obj.StackRuns;
            else
                obj.ReadData;
            end
        end
        
    end
    
    methods % Preparation for Fit Methods
        function SimulateRun(obj,varargin)
            % Compute Model Associated to a given Run
            
            p=inputParser;
            p.addParameter('qU','',@(x)isfloat(x));
            p.addParameter('qUfrac','',@(x)isfloat(x));
            p.addParameter('TimeSec','',@(x)isfloat(x));
            p.addParameter('BKG_RateAllFPDSec','',@(x)isfloat(x));
            p.addParameter('WGTS_CD_MolPerCm2','',@(x)isfloat(x));
            
            p.parse(varargin{:});
            
            qU                = p.Results.qU;
            qUfrac            = p.Results.qUfrac;
            TimeSec           = p.Results.TimeSec;
            BKG_RateAllFPDSec = p.Results.BKG_RateAllFPDSec;
            WGTS_CD_MolPerCm2 = p.Results.WGTS_CD_MolPerCm2;
            
            [TTFSD,DTFSD,HTFSD] = obj.SetDefaultFSD;
            
            if strcmp(obj.KTFFlag,'WGTSMACE_NIS1')
                obj.KTFFlag = 'WGTSMACE';
                NIS = 1;
            else
                NIS = 7;
            end
            
            TBDarg  = {obj.RunData,...
                'ISCS',obj.ISCSFlag,...
                'recomputeRF','OFF',...
                'ELossFlag',obj.ELossFlag,...
                'PixList',obj.PixList,...
                'RingList',obj.RingList,...
                'DTFSD',DTFSD,...
                'HTFSD',HTFSD,...
                'TTFSD',TTFSD,...
                'DopplerEffectFlag',obj.DopplerEffectFlag,...
                'RadiativeFlag',obj.RadiativeFlag,...
                'RingMerge',obj.RingMerge...
                'KTFFlag',obj.KTFFlag,...
                'NIS',NIS,...
                'SynchrotronFlag',obj.SynchrotronFlag,...
                'AngularTFFlag',obj.AngularTFFlag,...
                'FSD_Sigma',obj.FSD_Sigma,...
                'BKG_PtSlope',obj.BKG_PtSlope};
            
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

            ref_func = @ref_RunAnalysis;
            
            switch obj.AnaFlag
                case 'StackPixel'
                    obj.ModelObj =  ref_func(TBDarg{:},'FPD_Segmentation','OFF');
                case 'SinglePixel'
                    obj.ModelObj = ref_func(TBDarg{:},'FPD_Segmentation','SINGLEPIXEL');
                case 'MultiPixel'
                    obj.ModelObj = ref_func(TBDarg{:},'FPD_Segmentation','MULTIPIXEL','nTeBinningFactor',5);
                case 'Ring'
                    obj.ModelObj =  ref_func(TBDarg{:},'FPD_Segmentation','RING');
                    %                 case 'Ring'
                    %                     obj.ModelObj = ref_func(TBDarg{:},'FPD_Segmentation','RING','FPD_Ring',obj.RingList,...
                    %                         'nTeBinningFactor',20);
                    %                     if size(obj.RunData.TBDIS,2) > 1
                    %                         obj.RunData.TBDIS = sum(obj.RunData.TBDIS(:,obj.ModelObj.ring{obj.ring}),2);
                    %                         obj.RunData.qU = mean(obj.RunData.qU(:,obj.ModelObj.ring{obj.ring}),2);
                    %                     end
            end
            
            obj.ModelObj.ComputeTBDDS;
            obj.ModelObj.ComputeTBDIS;
        end
        function PlotDataModel(obj,varargin)
            % Plot overlaying Data and Model - Dedicated to FT
            % Includes Error Band depending on the Fit option
            % Model = TBDIS
            % (not a Fit)
            p=inputParser;
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF','png','pdf','eps'}));
           p.addParameter('Style','Reg',@(x)ismember(x,{'Reg','PRD'}));
            p.parse(varargin{:});
            saveplot = p.Results.saveplot;
            Style    = p.Results.Style;
            
            switch obj.chi2
                case {'chi2Stat','chi2P'}
                    PlotCM = diag(obj.RunData.TBDISE(:).^2);
                case {'chi2CM','chi2CMFrac'}
                    PlotCM = obj.FitCM;
                case 'chi2CMShape'
                    PlotCM = obj.FitCMShape;
            end
            
            % Model / Data overlay
            fig6 = figure(6);
            set(fig6, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
            
            h1=errorbar(obj.RunData.qU-18575,...
                obj.RunData.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                obj.RunData.TBDISE./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                's','MarkerSize',3,'MarkerEdgeColor' , rgb('DarkSlateGray'));
            h1.MarkerSize = 3; h1.CapSize = 0;h1.LineStyle= 'none';h1.LineWidth= 2;legend hide
            hold on
            hfit1 = boundedline(obj.ModelObj.qU-18575,...
                obj.ModelObj.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                diag(sqrt(PlotCM))./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                'alpha','cmap',rgb('CadetBlue'));
            [ll la]= boundedline(obj.ModelObj.qU-18575,...
                obj.ModelObj.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                sqrt(diag(PlotCM))./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                'alpha','cmap',rgb('CadetBlue'));
            ll.LineStyle='--'; legend hide
            hold off
            
            set(gca, 'YScale', 'lin');
            xlabel(sprintf('Retarding energy - %.0f (eV)',18574),'FontSize',16);
            ylabel('Rate (cps)','FontSize',16);
            set(gca,'FontSize',24);
            l1 = sprintf('Data (run %s, %.0f minutes)',num2str(obj.RunNr),obj.ModelObj.TimeSec/60);
            l2 = sprintf('model and uncorrelated error band (not a fit)');
            axis([min(obj.ModelObj.qU-obj.ModelObj.Q) max(obj.ModelObj.qU-obj.ModelObj.Q) -1000 max(obj.RunData.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec)*1.1 ]);
            title('');
            PrettyFigureFormat
            set(gca,'FontSize',24);
            
            % zoomPlot to highlight a portion of the major plot
            [~,~] = zoomPlotError(obj.RunData.qU-obj.ModelObj.Q,...
                obj.RunData.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                obj.RunData.TBDISE./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                [-220 31],[0.4 0.45 0.45 0.3],[3 4]);
            hold on
            set(gca, 'YScale', 'log');
            [ll , ~]=boundedline(obj.ModelObj.qU-obj.ModelObj.Q,...
                obj.ModelObj.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                sqrt(diag(PlotCM))./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                'alpha','cmap',rgb('CadetBlue'));
            ll.LineStyle='-';ll.LineWidth=1; legend hide
            hold off
            a = legend([h1 la],l1,l2,'Location','NorthEast'); a.String=a.String(1:2);
            legend('boxoff');
            set(gca,'yscale','log');
            %             if obj.RunNr > 1
            %                 title(sprintf('%s Data - Samak Fit (%s) \n(18575 eV - qUmin) = %.0f eV',obj.StackFileName,obj.chi2,obj.ModelObj.Q_i-obj.RunData.qU(obj.exclDataStart)),'Interpreter', 'none');
            %             else
            %                 title(sprintf('Run %u Data - Samak Fit (%s) \n(18575 eV - qU_{min}) = %.0f eV',obj.RunNr,obj.chi2,obj.ModelObj.Q_i-obj.RunData.qU(obj.exclDataStart)));
            %             end
            legend('boxoff')
            PrettyFigureFormat;
            set(gca,'YMinorTick','off');
            set(gca,'TickLength',[0.01 0.01]);
            set(gca,'FontSize',24);
            set(gca,'YTick',[1 100]);
            
            if strcmp(saveplot,'ON') || strcmp(saveplot,'pdf')
                save_name = sprintf('./plots/DataModel%s_%s-excl%u.pdf',num2str(obj.RunNr),obj.chi2,obj.exclDataStart);
                export_fig(fig6,save_name);
            elseif ~strcmp(saveplot,'OFF')
                save_name = sprintf('./plots/DataModel%s_%s-excl%u.%s',num2str(obj.RunNr),obj.chi2,obj.exclDataStart,saveplot);
                export_fig(fig6,save_name);
            end
        end
        function PlotDataModel_KNM1(obj,varargin)
            % Plot overlaying Data and Model - Dedicated to KNM1
            % Includes Error Band depending on the Fit option
            % Model = TBDIS
            % (not a Fit)
            p=inputParser;
            p.addParameter('saveplot','ON',@(x)ismember(x,{'ON','OFF','png','pdf','eps'}));
            p.addParameter('scale','log',@(x)ismember(x,{'lin','log'}));
            p.addParameter('RunLabel',numel(obj.RunList));
            p.addParameter('ErrorBarScaling',50,@(x)isfloat(x)); % limits for norm. residuals
            p.addParameter('Style','PRD',@(x)ismember(x,{'Reg','PRD'}));
            p.addParameter('TickDir','Out',@(x)ismember(x,{'In','Out'}));

            p.parse(varargin{:});
            saveplot            = p.Results.saveplot;
            scale               = p.Results.scale;
            RunLabel            = p.Results.RunLabel;
            Style               = p.Results.Style;
            TickDir             = p.Results.TickDir;
            obj.ErrorBarScaling = p.Results.ErrorBarScaling;

            switch obj.chi2
                case {'chi2Stat','chi2P'}
                    PlotCM = diag(obj.RunData.TBDISE(:).^2);
                case {'chi2CM','chi2CMFrac'}
                    PlotCM = obj.FitCM;
                case 'chi2CMShape'
                    PlotCM = obj.FitCMShape;
            end
            
            % Model / Data overlay
            fig6 = figure(6);
            set(fig6, 'Units', 'normalized', 'Position', [0., 0.1, 0.75 , 0.8]);
           
            AxisLabelFontSize = 35;
            LocalFontSize = 27;
            
           % best fit
            [lmodel, ~]= boundedline(obj.ModelObj.qU(obj.exclDataStart:end),...
                obj.ModelObj.TBDIS(obj.exclDataStart:end)./obj.ModelObj.qUfrac(obj.exclDataStart:end)./obj.ModelObj.TimeSec,...
                sqrt(diag(PlotCM((obj.exclDataStart:end),(obj.exclDataStart:end))))...
                ./obj.ModelObj.qUfrac(obj.exclDataStart:end)./obj.ModelObj.TimeSec,...
                'alpha','cmap',rgb('DodgerBlue'));
            lmodel.LineWidth= 3;
            lmodel.LineStyle='-'; %legend hide
            hold on;

             % Background Line
            [lbkg, ~]= boundedline(obj.ModelObj.qU(obj.exclDataStart:end),...
                obj.ModelObj.BKG_RateSec./obj.ModelObj.qUfrac(obj.exclDataStart:end).*obj.ModelObj.qUfrac(obj.exclDataStart:end),...
                obj.FitResult.err(3)./obj.ModelObj.qUfrac(obj.exclDataStart:end).*obj.ModelObj.qUfrac(obj.exclDataStart:end),...
                'alpha','cmap',rgb('IndianRed')); lbkg.LineWidth= 3;
            lbkg.LineStyle='-.'; legend hide
      
            % Tritium
            [lsignal, ~]= boundedline(obj.ModelObj.qU(obj.exclDataStart:end),...
                obj.ModelObj.TBDIS(obj.exclDataStart:end)./obj.ModelObj.qUfrac(obj.exclDataStart:end)./obj.ModelObj.TimeSec-...
                obj.ModelObj.BKG_RateSec./obj.ModelObj.qUfrac(obj.exclDataStart:end).*obj.ModelObj.qUfrac(obj.exclDataStart:end),...
                0./obj.ModelObj.qUfrac(obj.exclDataStart:end).*obj.ModelObj.qUfrac(obj.exclDataStart:end),...
                'alpha','cmap',rgb('Orange')); lsignal.LineWidth= 3;
            lsignal.LineStyle=':'; legend hide
            
            pdata =errorbar(obj.RunData.qU(obj.exclDataStart:end),...
                obj.RunData.TBDIS(obj.exclDataStart:end)./obj.ModelObj.qUfrac(obj.exclDataStart:end)./obj.ModelObj.TimeSec,...
                obj.RunData.TBDISE(obj.exclDataStart:end)./obj.ModelObj.qUfrac(obj.exclDataStart:end)./obj.ModelObj.TimeSec,...
                '.','MarkerSize',27,'MarkerEdgeColor' , rgb('Black'),'MarkerFaceColor' , rgb('Black'));
            pdata.Color =rgb('Black') ; pdata.CapSize = 0;pdata.LineStyle= 'none';pdata.LineWidth= 3;legend hide

            if ~strcmp(Style,'PRD')
                if obj.RunNr
                    title(sprintf('KATRIN First m_^2(\\nu_e) Mass Run (KNM1) - run %s, %.0f hours, %.0f e^-',num2str(obj.RunNr),obj.ModelObj.TimeSec/3600,sum(obj.ModelObj.TBDIS(obj.exclDataStart:end))));
                else
                    %title(sprintf('KATRIN - m_{\\beta} Test Scan - March 2019 - %0.f run stacked, %.0f hours, %.0f e^-',RunLabel,obj.ModelObj.TimeSec/3600,sum(obj.ModelObj.TBDIS(obj.exclDataStart))));
                    title(sprintf('KATRIN First m^2(\\nu_e) Mass Run (KNM1) - %.0f hours, %.0g e^- in [E_0-%.0f ; E_0+%.0f] eV',obj.ModelObj.TimeSec/3600,sum(obj.ModelObj.TBDIS(obj.exclDataStart:end)),...
                        obj.ModelObj.Q_i-obj.ModelObj.qU(obj.exclDataStart),obj.ModelObj.qU(end)-obj.ModelObj.Q_i));
                end
            end
            
            l1 = sprintf(' KATRIN data with 1\\sigma error bars');
            switch obj.chi2
                case 'chi2Stat'
                    l2 = sprintf(' Fit result with 1\\sigma uncertaintes (stat. only)');
                otherwise
                    l2 = sprintf(' Fit result with 1\\sigma uncertaintes (stat. and syst.)');
            end
             set(gca, 'YScale', scale);
            
            PRLFormat;
            xlabel(sprintf('Retarding energy (eV)'));
            ylabel('Rate (cps)');
            set(gca,'FontSize',LocalFontSize);
            set(get(gca,'XLabel'),'FontSize',AxisLabelFontSize);
            set(get(gca,'YLabel'),'FontSize',AxisLabelFontSize)
            %xlim([-100,+51]); ylim([0.1,1000]);
            legend([pdata lmodel lbkg lsignal],l1,l2,' Background',' Tritium signal','Location','northwest','Box','off'); %a.String=a.String(1:2);
            % set(a,'Color',rgb('White'),'Box', 'on');
          
            grid off;
            hold off;
            
            if strcmp(TickDir,'Out')
                set(gca,'TickDir','out');
%               remove top and right ticks
                a = gca;
                set(a,'box','off','color','none')% set box property to off and remove background color
                b = axes('Position',a.Position,...
                    'box','on','xtick',[],'ytick',[],'LineWidth',1.5);% create new, empty axes with box but without ticks
                axes(a)% set original axes as active
                linkaxes([a b]) % link axes in case of zooming
            end
             axis([0.9999*obj.ModelObj.qU(obj.exclDataStart) 1.0001*max(obj.ModelObj.qU)  obj.ModelObj.BKG_RateSec/2 max(obj.RunData.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec)*1.1 ]);
            ax = gca;
            ax.XAxis.Exponent = 0;
                 
            %xlim([-40 50]);
            ylim([0.2 100]) %ylim([0.2 30]);
            % zoomPlot to highlight a portion of the major plot
            [~,~] = zoomPlotError(obj.RunData.qU,...
                obj.RunData.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                obj.RunData.TBDISE./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
               [obj.ModelObj.Q_i-25 obj.ModelObj.Q_i+5],[0.4 0.33 0.45 0.45],[]);%[-25 5],[0.5 0.35 0.35 0.5],[]);
          
            hold on;
            % best fit
            [lzoomBestFit,~] = boundedline(obj.ModelObj.qU,...
                obj.ModelObj.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                diag(sqrt(PlotCM))./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec.*1,...
                'alpha','cmap',rgb('DodgerBlue'));
            lzoomBestFit.LineWidth = lmodel.LineWidth;
            hold on
            % data
            zoomData=errorbar(obj.RunData.qU,...
                obj.RunData.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                obj.RunData.TBDISE./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec*obj.ErrorBarScaling,...
                '.','MarkerSize',pdata.MarkerSize,'Color',pdata.Color);
            zoomData.CapSize = pdata.CapSize; zoomData.LineStyle= pdata.LineStyle; zoomData.LineWidth=pdata.LineWidth; legend hide
            legZoom = legend(zoomData,sprintf(' Zoom on {\\itm}_\\nu^2 ROI \n KATRIN data with 1 \\sigma error bars \\times %.0f',obj.ErrorBarScaling));
            hold off
            set(gca, 'YScale', 'lin');
           % a = legend([pdata la lbkg lsignal],l1,l2,'background    ','tritium    ','Location','SouthWest','Box','off'); a.String=a.String(1:4);
           % a.NumColumns=4;
          %  set(a,'Color',rgb('White'),'EdgeColor',rgb('White'),'Box', 'on');
            legend('boxoff');
            set(gca,'yscale','lin');
            PRLFormat
            
            set(gca,'TickLength',[0.01 0.01]);
            set(gca,'YTick',[0.1, 0.5, 1, 2, 3, 4, 5]);
            set(gca,'YAxisLocation','right');
            %aZ.Title.String = ' Zoom on m_{\beta} ROI';
            set(gca,'FontSize',LocalFontSize);
            if strcmp(TickDir,'Out')
                set(gca,'TickDir','out');
%                remove top and right ticks
                a = gca;
                set(a,'box','off','color','none')% set box property to off and remove background color
                b = axes('Position',a.Position,...
                    'box','on','xtick',[],'ytick',[],'LineWidth',1.5);% create new, empty axes with box but without ticks
                axes(a)% set original axes as active
               % linkaxes([a b]) % link axes in case of zooming
            end
             ax = gca;
            ax.XAxis.Exponent = 0;
            if strcmp(saveplot,'ON')
                save_name = sprintf('./plots/KNM1_DataModel%s_%s-excl%u.png',num2str(obj.RunNr),obj.chi2,obj.exclDataStart);
                %export_fig(fig6,save_name,'-q101','-m3');
                export_fig(fig6,strrep(save_name,'.png','.pdf'));
                fprintf('save to %s \n',save_name);
            end
          end
        function PlotDataModel_KSN1(obj,varargin)
            % Plot overlaying Data and Model - Dedicated to KSN1
            % Includes Error Band depending on the Fit option
            % Model = TBDIS
            % (not a Fit)
            p=inputParser;
            p.addParameter('saveplot','ON',@(x)ismember(x,{'ON','OFF','png','pdf','eps'}));
            p.addParameter('scale','log',@(x)ismember(x,{'lin','log'}));
            p.addParameter('RunLabel',numel(obj.RunList));
            p.addParameter('ErrorBarScaling',50,@(x)isfloat(x)); % limits for norm. residuals            
            
            p.parse(varargin{:});
            saveplot = p.Results.saveplot;
            scale = p.Results.scale;
            RunLabel = p.Results.RunLabel;
            obj.ErrorBarScaling= p.Results.ErrorBarScaling;

            switch obj.chi2
                case {'chi2Stat','chi2P'}
                    PlotCM = diag(obj.RunData.TBDISE(:).^2);
                case {'chi2CM','chi2CMFrac'}
                    PlotCM = obj.FitCM;
                case 'chi2CMShape'
                    PlotCM = obj.FitCMShape;
            end
            
            % Model / Data overlay
            fig6 = figure(6);
            set(fig6, 'Units', 'normalized', 'Position', [0., 0.1, 0.9 , 0.8]);
            
            h1=errorbar(obj.RunData.qU-obj.ModelObj.Q_i,...
                obj.RunData.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                obj.RunData.TBDISE./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec*0,...
                'o','MarkerSize',10,'MarkerEdgeColor' , rgb('DarkSlateGray'),'MarkerFaceColor' , rgb('DarkSlateGray'));
             h1.Color =rgb('DarkSlateGray') ; h1.MarkerSize = 6; h1.CapSize = 0;h1.LineStyle= 'none';h1.LineWidth= 3;legend hide
            hold on
            hfit1 = boundedline(obj.ModelObj.qU-obj.ModelObj.Q_i,...
                obj.ModelObj.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                diag(sqrt(PlotCM))./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                'alpha','cmap',rgb('IndianRed'));
            [ll la]= boundedline(obj.ModelObj.qU-obj.ModelObj.Q_i,...
                obj.ModelObj.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                sqrt(diag(PlotCM))./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                'alpha','cmap',rgb('IndianRed'));
            ll.LineWidth= 3
            ll.LineStyle='--'; legend hide
            
            % Background Line
            [lll lla]= boundedline(obj.ModelObj.qU-obj.ModelObj.Q_i,...
                obj.ModelObj.BKG_RateSec./obj.ModelObj.qUfrac.*obj.ModelObj.qUfrac,...
                obj.FitResult.err(3)./obj.ModelObj.qUfrac.*obj.ModelObj.qUfrac,...
                'alpha','cmap',rgb('SeaGreen')); lll.LineWidth= 3;
            lll.LineStyle='--'; legend hide
            % Tritium
            [llll llla]= boundedline(obj.ModelObj.qU-obj.ModelObj.Q_i,...
                obj.ModelObj.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec-obj.ModelObj.BKG_RateSec./obj.ModelObj.qUfrac.*obj.ModelObj.qUfrac,...
                0./obj.ModelObj.qUfrac.*obj.ModelObj.qUfrac,...
                'alpha','cmap',rgb('RoyalBlue')); llll.LineWidth= 3;
            llll.LineStyle='--'; legend hide
            
           % data again, for superimposition 
           h1=errorbar(obj.RunData.qU-obj.ModelObj.Q_i,...
                obj.RunData.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,...
                obj.RunData.TBDISE./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec*0,...
                'o','MarkerSize',10,'MarkerEdgeColor' , rgb('DarkSlateGray'),'MarkerFaceColor' , rgb('DarkSlateGray'));
            h1.Color =rgb('DarkSlateGray') ;  h1.MarkerSize = 6; h1.CapSize = 0;h1.LineStyle= 'none';h1.LineWidth= 3;legend hide
            
            hold off
            axis([1.1*min(obj.ModelObj.qU-obj.ModelObj.Q) 1.1*max(obj.ModelObj.qU-obj.ModelObj.Q)  obj.ModelObj.BKG_RateSec/2 max(obj.RunData.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec)*1.1 ]);
            set(gca, 'YScale', scale);
            xlabel(sprintf('retarding energy - %.1f (eV)',obj.ModelObj.Q_i),'FontSize',16);
            ylabel('rate (cps)','FontSize',16);
            set(gca,'FontSize',24);
            if obj.RunNr
            title(sprintf('KATRIN First eV-Sterile Neutrino Run Run (KSN1) - run %s, %.0f hours, %.0f e^-',num2str(obj.RunNr),obj.ModelObj.TimeSec/3600,sum(obj.ModelObj.TBDIS(obj.exclDataStart:end))));
            else
            %title(sprintf('KATRIN - m_{\\beta} Test Scan - March 2019 - %0.f run stacked, %.0f hours, %.0f e^-',RunLabel,obj.ModelObj.TimeSec/3600,sum(obj.ModelObj.TBDIS(obj.exclDataStart))));
            title(sprintf('KATRIN First eV-Sterile Neutrino Run (KSN1) - %.0f hours, %.0g e^- in [E_0-%.0f ; E_0+%.0f] eV',obj.ModelObj.TimeSec/3600,sum(obj.ModelObj.TBDIS(obj.exclDataStart:end)),...
                obj.ModelObj.Q_i-obj.ModelObj.qU(obj.exclDataStart),obj.ModelObj.qU(end)-obj.ModelObj.Q_i));
            end
            l1 = 'data    ';
            l2 = sprintf('model + uncertainty    ');
            PrettyFigureFormat
            set(gca,'FontSize',22);
            %xlim([-100,+51]); ylim([0.1,1000]);
            a = legend([h1 la lll llll],l1,l2,'background','tritium','Location','NorthEast','Box','off'); %a.String=a.String(1:2);
            set(a,'Color',rgb('White'),'Box', 'on');
            xlim([-90 50]);
            ylim([0.2 400]);
            grid on;
            
 
            if strcmp(saveplot,'ON') 
                save_name = sprintf('./plots/KSN1s_DataModel%s_%s-excl%u.png',num2str(obj.RunNr),obj.chi2,obj.exclDataStart);
                export_fig(fig6,save_name,'-q101','-m3');
            end
           
          end
          
          

    end
    methods % Covariance Matrix methods
        function [StatCM, StatCMFrac] = ComputeCM_StatPNP(obj,varargin)
            % Compute Statistical Covariance Matrix including P/NP fluctuations
            
            % Seperate Signal and Background
            Signal = obj.ModelObj.TBDIS - obj.ModelObj.TimeSec .* obj.ModelObj.qUfrac .* obj.ModelObj.BKG_RateSec;
            Signal(obj.RunData.qU>obj.ModelObj.Q) = 0;
            signalE         = sqrt(Signal);
%            backgroundE     = obj.NonPoissonScaleFactor * sqrt(obj.ModelObj.TimeSec .* obj.ModelObj.qUfrac .* obj.ModelObj.BKG_RateSec);%obj.ModelObj.BKG_RateSec);
             backgroundE     = sqrt(obj.ModelObj.TimeSec .* obj.ModelObj.qUfrac .* obj.ModelObj.BKG_RateSec .*obj.NonPoissonScaleFactor.^2);%obj.ModelObj.BKG_RateSec);
          
            % Reshape if Multiring
            if strcmp(obj.AnaFlag,'Ring')
                signalE = reshape(signalE,[obj.ModelObj.nqU*obj.ModelObj.nRings,1]);
                backgroundE = reshape(backgroundE,[obj.ModelObj.nqU*obj.ModelObj.nRings,1]);
            end
            
            % Build statsitical covariance matrix
            StatCM      = diag( (signalE.^2 + backgroundE.^2) );
            StatCMFrac  = diag(1./(signalE.^2 + backgroundE.^2));      
        end
        function ComputeCM(obj,varargin)
            % Compute Covariance Matrix
             defaultEffects =  obj.GetDefaultEffects;
            
            [SysErr,~] = GetSysErr(obj.SysBudget);
            
            p = inputParser;
            %RunAnalysis settings
            p.addParameter('InitNormFit','ON',@(x)ismember(x,{'ON','OFF'})); % Init Model Normalization + Background with Fit
            p.addParameter('BkgCM','ON',@(x)ismember(x,{'ON','OFF'}));       % Use Background CovMat ON/OFF
            p.addParameter('BkgPtCM','ON',@(x)ismember(x,{'ON','OFF'}));      % Use Background CovMat ON/OFF
            % CovMat settings
            p.addParameter('SysEffects',defaultEffects,@(x)isstruct(x));
            p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('nTrials',1000,@(x)isfloat(x));
            p.addParameter('DataDriven',SysErr.DataDriven,@(x)ismember(x,{'ON','OFF'}));
%           % default systematic uncertainties
            p.addParameter('WGTS_TASR_RelErr',SysErr.WGTS_TASR_RelErr,@(x)all(isfloat(x)));
            p.addParameter('Stack_qUErr',0,@(x)all(isfloat(x)));
            p.addParameter('Stack_qUfracRelErr',0,@(x)all(isfloat(x)));
            p.addParameter('nStack','',@(x)isfloat(x));
            p.addParameter('FSDNorm_RelErr',SysErr.FSDNorm_RelErr,@(x)isfloat(x) && x>=0);
            p.addParameter('FSDShapeGS_RelErr',SysErr.FSDShapeGS_RelErr,@(x)isfloat(x) && x>=0);
            p.addParameter('FSDShapeES_RelErr',SysErr.FSDShapeES_RelErr,@(x)isfloat(x) && x>=0);
            p.addParameter('MACE_Ba_T_RelErr',SysErr.MACE_Ba_T_RelErr,@(x)isfloat(x) && x>=0);
            p.addParameter('MACE_Bmax_T_RelErr',SysErr.MACE_Bmax_T_RelErr,@(x)isfloat(x) && x>=0);
            p.addParameter('WGTS_B_T_RelErr',SysErr.WGTS_B_T_RelErr,@(x)isfloat(x) && x>=0);
            p.addParameter('WGTS_CD_MolPerCm2_RelErr',SysErr.WGTS_CD_MolPerCm2_RelErr,@(x)isfloat(x) && x>=0);
            p.addParameter('ISXsection_RelErr',SysErr.ISXsection_RelErr,@(x)isfloat(x) && x>=0);
            p.addParameter('FPDeff_RelErr',SysErr.FPDeff_RelErr,@(x)isfloat(x) && x>=0); 
            p.addParameter('PlotSaveCM','OFF',@(x)ismember(x,{'ON','OFF'})); %plot individual covmats and saves plot
            p.addParameter('MaxSlopeCpsPereV',SysErr.MaxSlopeCpsPereV,@(x)isfloat(x));
            p.addParameter('BkgRange',-5,@(x)isfloat(x));
            p.addParameter('RW_SigmaErr',0.1,@(x)isfloat(x) );
            p.addParameter('RW_MultiPosErr',0,@(x)isfloat(x));
            p.addParameter('MACE_VarErr',SysErr.MACE_VarErr,@(x)isfloat(x));
            p.addParameter('is_EOffsetErr',SysErr.is_EOffsetErr,@(x)isfloat(x));
            p.addParameter('BkgRingCorrCoeff',0,@(x)isfloat(x));
            p.addParameter('BkgScalingOpt',2,@(x)isfloat(x));
            p.addParameter('BkgMode','Gauss',@(x)ismember(x,{'SlopeFit','Gauss'}));
            p.addParameter('BKG_PtSlopeErr',SysErr.BKG_PtSlopeErr,@(x)isfloat(x));
            p.parse(varargin{:});
            
            InitNormFit              = p.Results.InitNormFit;
            BkgCM                    = p.Results.BkgCM;
            BkgPtCM                  = p.Results.BkgPtCM;
            SysEffects               = p.Results.SysEffects;
            RecomputeFlag            = p.Results.RecomputeFlag;
            nTrials                  = p.Results.nTrials;
            DataDriven               = p.Results.DataDriven;
            MACE_Ba_T_RelErr         = p.Results.MACE_Ba_T_RelErr;
            MACE_Bmax_T_RelErr       = p.Results.MACE_Bmax_T_RelErr;
            WGTS_B_T_RelErr          = p.Results.WGTS_B_T_RelErr;
            WGTS_CD_MolPerCm2_RelErr = p.Results.WGTS_CD_MolPerCm2_RelErr;
            ISXsection_RelErr        = p.Results.ISXsection_RelErr;
            WGTS_TASR_RelErr         = p.Results.WGTS_TASR_RelErr;
            FSDNorm_RelErr           = p.Results.FSDNorm_RelErr;
            FSDShapeGS_RelErr        = p.Results.FSDShapeGS_RelErr;
            FSDShapeES_RelErr        = p.Results.FSDShapeES_RelErr;
            FPDeff_RelErr            = p.Results.FPDeff_RelErr;
            Stack_qUErr              = p.Results.Stack_qUErr;
            Stack_qUfracRelErr       = p.Results.Stack_qUfracRelErr;
            PlotSaveCM               = p.Results.PlotSaveCM;
            nStack                   = p.Results.nStack;
            MaxSlopeCpsPereV         = p.Results.MaxSlopeCpsPereV;
            BkgRange                 = p.Results.BkgRange;
            RW_SigmaErr              = p.Results.RW_SigmaErr;
            RW_MultiPosErr           = p.Results.RW_MultiPosErr;
            MACE_VarErr              = p.Results.MACE_VarErr; %longitudinal plasma
            is_EOffsetErr            = p.Results.is_EOffsetErr; %longitudinal plasma
            BkgRingCorrCoeff         = p.Results.BkgRingCorrCoeff; % ring to ring correlation coefficient
            BkgScalingOpt            = p.Results.BkgScalingOpt;
            BkgMode                  = p.Results.BkgMode;
            BKG_PtSlopeErr           = p.Results.BKG_PtSlopeErr;
            % --------------------- END PARSER ---------------------------------%
            
            % --------------------  Initialize Covariance Matrix----------------------------%
            if isempty(nStack) && isa(obj,'MultiRunAnalysis')
                nStack =numel(obj.StackedRuns);
            else
                nStack=1;
            end
            
            % init
            obj.FitCM          = zeros(obj.ModelObj.nqU*obj.ModelObj.nPixels,obj.ModelObj.nqU*obj.ModelObj.nPixels);
            obj.FitCMFrac      = zeros(obj.ModelObj.nqU*obj.ModelObj.nPixels,obj.ModelObj.nqU*obj.ModelObj.nPixels);
            obj.FitCMShape     = zeros(obj.ModelObj.nqU*obj.ModelObj.nPixels,obj.ModelObj.nqU*obj.ModelObj.nPixels);
            obj.FitCMNorm      = zeros(obj.ModelObj.nqU*obj.ModelObj.nPixels,obj.ModelObj.nqU*obj.ModelObj.nPixels);
            obj.FitCMFracShape = zeros(obj.ModelObj.nqU*obj.ModelObj.nPixels,obj.ModelObj.nqU*obj.ModelObj.nPixels);
            
            %Initialize Normalization and Background with a stat. Fit
            if strcmp(InitNormFit,'ON')
                % 40 eV range, stat only, free parameters: E0, Bkg, Norm
                obj.InitModelObj_Norm_BKG('RecomputeFlag','OFF');
            end
            
            obj.FitCM_Obj = CovarianceMatrix('StudyObject',obj.ModelObj, 'nTrials',nTrials,...
                'SysEffect',SysEffects,'RecomputeFlag',RecomputeFlag,'SanityPlots','OFF',...
                'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr,...
                'MACE_Bmax_T_RelErr',MACE_Bmax_T_RelErr,...
                'MACE_Ba_T_RelErr',MACE_Ba_T_RelErr,...
                'WGTS_B_T_RelErr',WGTS_B_T_RelErr,...
                'ISXsection_RelErr',ISXsection_RelErr,...
                'WGTS_TASR_RelErr',WGTS_TASR_RelErr,...
                'FPDeff_RelErr',FPDeff_RelErr,...
                'FSDNorm_RelErr',FSDNorm_RelErr,'FSDShapeGS_RelErr',FSDShapeGS_RelErr,'FSDShapeES_RelErr',FSDShapeES_RelErr,...
                'ShapeCMnormIndex',obj.exclDataStart,'nStack',nStack,...
                'Stack_qUfracRelErr',Stack_qUfracRelErr,'Stack_qUErr',Stack_qUErr,...
                'RW_SigmaErr',RW_SigmaErr,...
                'RW_MultiPosErr',RW_MultiPosErr,...
                'NonPoissonScaleFactor',obj.NonPoissonScaleFactor,...
                'MACE_VarErr',MACE_VarErr,...
                'is_EOffsetErr',is_EOffsetErr,...
                'BKG_PtSlopeErr',BKG_PtSlopeErr);
            
            if strcmp(DataDriven,'ON')
                obj.Get_DataDriven_RelErr_TASR; % replaces default values with Data Driven values
                obj.Get_DataDriven_RelErr_qU;
            end
            
            % -------------------- Computation of Covariance Matrix -------------------------------%
            % RecomputeFlag=='ON' : Computation of Covariance Matrix
            % RecomputeFlag=='OFF': Read adequard Covariance Matrix if available, otherwise: computation
            
            if any(structfun(@(x) contains(x,'ON'),SysEffects))
                obj.FitCM_Obj.ComputeCM('PlotSaveCM',PlotSaveCM);
                
                
                % Move to Fit CM and do renormalization with current statistics
                obj.FitCMFrac = obj.FitCM_Obj.CovMatFrac;
                TBDIS_NoBKG   = obj.ModelObj.TBDIS-(obj.ModelObj.BKG_RateSec.*obj.ModelObj.TimeSec.*obj.ModelObj.qUfrac);
                if strcmp(obj.AnaFlag,'Ring')
                    TBDIS_NoBKG = reshape(TBDIS_NoBKG,[obj.ModelObj.nqU*obj.ModelObj.nRings,1]);
                end
                obj.FitCM     = TBDIS_NoBKG.*obj.FitCMFrac.*TBDIS_NoBKG';
                
                % Compute Shape only covariance matrix for desired fit range
                
                try
                    if sum(sum(obj.FitCMFrac))~=0
                        [obj.FitCMShape,obj.FitCMFracShape] = obj.FitCM_Obj.DecomposeCM('CovMatFrac',obj.FitCMFrac,'exclDataStart',obj.exclDataStart);
                    else
                        obj.FitCMShape     = zeros(obj.ModelObj.nqU,obj.ModelObj.nqU);
                        obj.FitCMFracShape = zeros(obj.ModelObj.nqU,obj.ModelObj.nqU);
                    end
                catch
                    fprintf('Decomposition doesnt work - TASR covariance matrix probably too small (Knm2) \n')
                    fprintf('Temporary fix: take regular covariance matrix instead - no shape only \n')
                    obj.FitCMShape = obj.FitCM;
                    obj.FitCMFracShape = obj.FitCMFrac;
                end
                
            end
            % Background Covariance Matrix: Compute and Add to Signal Covariance Matrix
            if strcmp(BkgCM,'ON')
                obj.FitCM_Obj.ComputeCM_Background('Display',PlotSaveCM,...
                    'MaxSlopeCpsPereV',MaxSlopeCpsPereV,'BkgRange',BkgRange,...
                    'RingCorrCoeff',BkgRingCorrCoeff,'ScalingOpt',BkgScalingOpt,...
                    'Mode',BkgMode);
                
                obj.FitCM           = obj.FitCM          + obj.FitCM_Obj.CovMat;     % regular covmat:    add background covmat to signal covmat 
                obj.FitCMFrac       = obj.FitCMFrac      + obj.FitCM_Obj.CovMatFrac; % fractional covmat: add background covmat to signal covmat 
                
                % shape only:
                [BkgCMShape,BkgCMFracShape] = obj.FitCM_Obj.DecomposeCM('CovMatFrac',obj.FitCM_Obj.CovMatFrac,...
                    'exclDataStart',obj.exclDataStart,'BkgCM','ON'); 
                obj.FitCMShape      = obj.FitCMShape     + BkgCMShape;     % shape only covmat:            add background covmat to signal covmat 
                obj.FitCMFracShape  = obj.FitCMFracShape + BkgCMFracShape; % fractional shape only covmat: add background covmat to signal covmat 
            else
                cprintf('blue','RunAnalysis:ComputeCM: Background Covariance Matrix = OFF  \n')
            end
            
            % Penning trap covariance matrix
            if BKG_PtSlopeErr~=0 && strcmp(BkgPtCM,'ON')
                obj.ModelObj.nRuns = obj.nRuns;
                obj.FitCM_Obj.ComputeCM_BackgroundPT('Display',PlotSaveCM,'nTrials_loc',1e4);
                
                % add to previous cm
                obj.FitCM           = obj.FitCM          + obj.FitCM_Obj.CovMat;     % regular covmat:    add background covmat to signal covmat
                obj.FitCMFrac       = obj.FitCMFrac      + obj.FitCM_Obj.CovMatFrac; % fractional covmat: add background covmat to signal covmat
               
                % shape only:
                [BkgCMPtShape,BkgCMPtFracShape] = obj.FitCM_Obj.DecomposeCM('CovMatFrac',obj.FitCM_Obj.CovMatFrac,...
                    'exclDataStart',obj.exclDataStart,'BkgCM','PT');
                obj.FitCMShape      = obj.FitCMShape     + BkgCMPtShape;     % shape only covmat:            add background covmat to signal covmat
                obj.FitCMFracShape  = obj.FitCMFracShape + BkgCMPtFracShape; % fractional shape only covmat: add background covmat to signal covmat
            end
            
            % Compute Statistical Uncertainties,including P/NP fluctuations
            [StatCM, StatCMFrac] = obj.ComputeCM_StatPNP(varargin);
            
            %Add Statistical uncertainties, 
            switch obj.chi2
                case {'chi2Stat'}
                    obj.FitCM          =  StatCM;
                    obj.FitCMFrac      =  StatCMFrac;
                    obj.FitCMFracShape =  StatCMFrac;
                    obj.FitCMShape     =  StatCM;
                case {'chi2CMShape','chi2CM'}
                    obj.FitCM          = obj.FitCM      + StatCM;
                    obj.FitCMFrac      = obj.FitCMFrac  + StatCMFrac;
                    obj.FitCMFracShape = obj.FitCMFracShape + StatCMFrac;
                    obj.FitCMShape     = obj.FitCMShape + StatCM;
            end
            
            % Test for PositiveSemiDefinite
            obj.TestCM_PositiveSemiDefinite;
        end    
        
        function [WGTS_TASR_RelErr,SubRunActivity,TASR_CorrMat] = Get_DataDriven_RelErr_TASR(obj,varargin)
            p = inputParser;
            p.addParameter('SanityPlot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            SanityPlot = p.Results.SanityPlot;
            
            % Tritium Activity Relative Fluctuations
            % Read tritium activity from RunSummaries
            % Calculate relative error
            
            % Select qU values relevant for TASR calculation
            % Could Have A Cut OFF at E0-5 eV
             TritiumqUIndexes = (obj.ModelObj.qU(:,1)<=obj.ModelObj.Q_i);
             
            if isa(obj,'MultiRunAnalysis') %multiple runs: MultiRunAnalysis
                % take standard error of the mean
                switch obj.DataSet
                    case 'FirstTritium.katrin'
                        WGTS_TASR_RelErr = (std(obj.SingleRunData.WGTS_MolFrac_DT_SubRun(:,obj.SingleRunData.Select_all)')./...
                            (sqrt(numel(obj.StackedRuns)).*obj.ModelObj.WGTS_MolFrac_DT_SubRun))';
                        
                        % Set Fluctuations to ZERO above enpoint & Save
                        obj.FitCM_Obj.WGTS_TASR_RelErr = WGTS_TASR_RelErr .* TritiumqUIndexes;
                        
                        % only used in >= KNM1
                        SubRunActivity = NaN; 
                        TASR_CorrMat   = NaN;
                    case {'Knm1','Knm2'}
                        % compute tritium activity per run and per subrun
                        SubRunActivity = ...
                            (obj.SingleRunData.WGTS_MolFrac_TT_SubRun(:,obj.SingleRunData.Select_all)        + ...
                            0.5 * obj.SingleRunData.WGTS_MolFrac_HT_SubRun(:,obj.SingleRunData.Select_all)   + ...
                            0.5 * obj.SingleRunData.WGTS_MolFrac_DT_SubRun(:,obj.SingleRunData.Select_all)) .* ...
                            obj.SingleRunData.WGTS_CD_MolPerCm2_SubRun(:,obj.SingleRunData.Select_all)      .* ...
                            obj.ModelObj.ISXsection(obj.ModelObj.Q_i)*1e4;
                        
                        % renormalize to mean activity
                        MeanActivityStackedRun = mean(nanmean(SubRunActivity,2));
                        MeanActivitySingleRuns = nanmean(SubRunActivity,1);
                        ActiCorrFactor = repmat(MeanActivityStackedRun./MeanActivitySingleRuns,size(SubRunActivity,1),1);
                        
                        SubRunActivity(~isnan(SubRunActivity))=SubRunActivity(~isnan(SubRunActivity)).*ActiCorrFactor(~isnan(SubRunActivity));
                         
                        %take standard error of the mean
                        WGTS_TASR_AbsErr = std(SubRunActivity,0,2)./sqrt(numel(obj.StackedRuns));
                        
                        % relative uncertainty -> fractional covariance
                        WGTS_TASR_RelErr = WGTS_TASR_AbsErr./mean(SubRunActivity,2);
                        
                        % Set Fluctuations to ZERO above enpoint
                        WGTS_TASR_RelErr = WGTS_TASR_RelErr .* TritiumqUIndexes;
                      
                        if strcmp(obj.DataSet,'Knm1')
                            % pass to covariance matrix object
                            obj.FitCM_Obj.WGTS_TASR_RelErr =  diag(WGTS_TASR_RelErr.^2);
                        elseif strcmp(obj.DataSet,'Knm2')
                            % calculate correlation matrix
                            %SubRunActivity(isnan(SubRunActivity))=1;
                            TASR_CorrMat = corrcoef(SubRunActivity','Rows','pairwise');
                            %TASR_CorrMat(isnan(TASR_CorrMat)) = (nanmean(TASR_CorrMat(TASR_CorrMat>0)));
                            obj.FitCM_Obj.WGTS_TASR_RelErr = TASR_CorrMat.*WGTS_TASR_RelErr.*WGTS_TASR_RelErr';
                            WGTS_TASR_RelErr = obj.FitCM_Obj.WGTS_TASR_RelErr;
                            
                        end
                end
                
                if strcmp(SanityPlot,'ON') && strcmp(obj.DataSet,'Knm2')
                    x  = repmat(obj.RunData.qU-18574,obj.nRuns,1);
                    y = reshape(SubRunActivity,size(obj.RunData.qU,1).*obj.nRuns,1);
                    
                    f1 = figure('Units','normalized','Position',[0.1,0.1,0.9,0.5]);
                    s1 = subplot(1,2,1);
                    pref = plot(linspace(min(x)-5,max(x),10),MeanActivityStackedRun.*ones(10,1),'k-','LineWidth',2);
                    hold on;
                    p1 = scatter(x,y,'o','filled',...
                        'MarkerFaceColor',rgb('DodgerBlue'),'MarkerFaceAlpha',0.2);
                    xlabel('Retarding potential - 18574 (eV)');
                    ylabel(sprintf('Activity (molecules)'));
                    PrettyFigureFormat('FontSize',24)
                    xlim([min(x)-5,5]);
                    leg = legend([pref,p1],sprintf('Mean activity = %.2f molecules',MeanActivityStackedRun),...
                        sprintf('Subrun activity %.0f runs',obj.nRuns));
                    leg.EdgeColor = rgb('Silver');
                    leg.Location='northwest';
                    title('Single run activity normalized to mean activity','FontWeight','normal');
                     
                     %% subplot 2
                    s2 = subplot(1,2,2);
                    y2 = reshape(SubRunActivity./ActiCorrFactor,size(obj.RunData.qU,1).*obj.nRuns,1);
                    
                    pref = plot(linspace(min(x)-5,max(x),10),MeanActivityStackedRun.*ones(10,1),'k-','LineWidth',2);
                    hold on;
                    p1 = scatter(x,y2,'o','filled',...
                        'MarkerFaceColor',rgb('Orange'),'MarkerFaceAlpha',0.2);
                    xlabel('Retarding potential - 18574 (eV)');
                    ylabel(sprintf('Activity (molecules)'));
                    PrettyFigureFormat('FontSize',24)
                    xlim([min(x)-5,5]);
                    leg = legend([pref,p1],sprintf('Mean activity = %.2f molecules',MeanActivityStackedRun),...
                        sprintf('Subrun activity %.0f runs',obj.nRuns));
                    leg.EdgeColor = rgb('Silver');
                    leg.Location='northwest';
                    ylim([min(0.995*[y;y2]),1.005*max([y;y2])]);
                    linkaxes([s1,s2],'y')
                    title(sprintf('Single activity not normalized to mean activity'),'FontWeight','normal');
                end
                
            else % single run: RunAnalysis
                switch obj.DataSet
                    case 'FirstTritium.katrin'
                        WGTS_TASR_RelErr = std(obj.RunData.WGTS_MolFrac_DT_SubRun)./...
                            (sqrt(numel(obj.RunData.qU)).*obj.RunData.WGTS_MolFrac_DT_SubRun)';
                        obj.FitCM_Obj.WGTS_TASR_RelErr = WGTS_TASR_RelErr;
                    case 'Knm1'
                        % to do
                        WGTS_TASR_RelErr = 0.001; %temporary value
                end
            end
        end
        function [Stack_qUErr, Stack_qUfracRelErr] = Get_DataDriven_RelErr_qU(obj,varargin)
            % read qU from RunSummaries
            % calculate rel. error
            
            p=inputParser;
            p.addParameter('Debug','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            Debug = p.Results.Debug;

            if ~isa(obj,'MultiRunAnalysis') 
                %this makes only sense for stacked spectra
                fprintf('WARNING: single run --> set StackCM to OFF \n');
                return
            end
            
            % take standard deviation
            Stack_qUErr        = std(obj.SingleRunData.qU(:,obj.SingleRunData.Select_all)')';%./sqrt(numel(obj.StackedRuns));     % absolute error!
           % Stack_qUfracRelErr = std(obj.SingleRunData.qUfrac(:,obj.SingleRunData.Select_all)')'...
            %    ./mean(obj.SingleRunData.qUfrac(:,obj.SingleRunData.Select_all),2);%./sqrt(numel(obj.StackedRuns)); % relative error!
            
            % set to zero above endpoint
            TritiumqUIndexes   = (obj.ModelObj.qU<=obj.ModelObj.Q_i); 
            Stack_qUErr        = Stack_qUErr.* TritiumqUIndexes;
            Stack_qUfracRelErr =  0;%Stack_qUfracRelErr.*TritiumqUIndexes;
            
            % pass to covariance matrix
            obj.FitCM_Obj.Stack_qUErr        =  Stack_qUErr;
            obj.FitCM_Obj.Stack_qUfracRelErr =  Stack_qUfracRelErr;
            
             %sanity plot (distribution of qU-values)
            if strcmp(Debug,'ON')
                fig123 = figure('Renderer','painters');
                set(fig123, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6,1.3]);
                qu = 1:1:obj.ModelObj.nqU-2;
                for i=1:numel(qu)
                    subplot(ceil(numel(qu)/6),6,i);
                    h = histogram(obj.SingleRunData.qU(qu(i),obj.SingleRunData.Select_all)-18574);
                    h.FaceColor = rgb('LightSkyBlue');
                    h.FaceAlpha = 1;
                    PrettyFigureFormat; 
                    xticks(mean(obj.SingleRunData.qU(qu(i),obj.SingleRunData.Select_all)-18574))
                    xticklabels(sprintf('%.0f eV',mean(obj.SingleRunData.qU(qu(i),obj.SingleRunData.Select_all)-18574)))
                    ylabel(''); yticklabels('');
                    set(gca,'FontSize',18);
                end
                
                tstr= sprintf('HV fluctuations within subrun (%s stacked data) \n \\sigma_{min} = %.0f meV, \\sigma_{max} = %.0f meV, \\langle\\sigma\\rangle = %.0f meV',...
                    extractBefore(strrep(obj.RunData.RunName,'_',' '),' E0'), ...
                    min(std(obj.SingleRunData.qU,0,2))*1e3,max(std(obj.SingleRunData.qU,0,2))*1e3,1e3*mean(std(obj.SingleRunData.qU,0,2)));
                t =sgtitle(tstr);
                t.FontSize = get(gca,'FontSize');
                savedir = './plots/';
                MakeDir(savedir);
                savename = sprintf('%sDataDriven_%s_qUErr.pdf',savedir,obj.RunData.RunName);
                export_fig(fig123,savename);
                fprintf('Save plot to %s \n',savename)
            end
        end
        function  defaultEffects = GetDefaultEffects(obj)
            switch obj.DataSet
                case 'Knm1'
                    defaultEffects = struct(...
                        'RF_EL','ON',...   % Response Function(RF) EnergyLoss
                        'RF_BF','ON',...   % RF B-Fields
                        'RF_RX','ON',...   % Column Density, inel cross ection
                        'FSD','ON',...
                        'TASR','ON',...
                        'TCoff_RAD','OFF',...
                        'TCoff_OTHER','ON',...
                        'DOPoff','OFF',...
                        'Stack','ON',...
                        'FPDeff','ON');
                case 'Knm2'
                     defaultEffects = struct(...
                        'RF_EL','ON',...   % Response Function(RF) EnergyLoss
                        'RF_BF','ON',...   % RF B-Fields
                        'RF_RX','ON',...   % Column Density, inel cross ection
                        'FSD','ON',...
                        'TASR','ON',...
                        'TCoff_RAD','OFF',...
                        'TCoff_OTHER','ON',...
                        'DOPoff','OFF',...
                        'Stack','ON',...
                        'FPDeff','ON',...
                        'LongPlasma','ON');
                case 'FirstTritium.katrin'
                    defaultEffects = struct(...
                        'RF_EL','ON',...   % Response Function(RF) EnergyLoss
                        'RF_BF','ON',...   % RF B-Fields
                        'RF_RX','ON',...   % Column Density, inel cross ection
                        'FSD','ON',...
                        'TASR','ON',...
                        'TCoff_RAD','ON',...
                        'TCoff_OTHER','ON',...
                        'DOPoff','OFF',...
                        'Stack','ON',...
                        'FPDeff','ON');
            end
        end
        function InitModelObj_Norm_BKG(obj,varargin)
            % Init Model (Normalization and Background) with Data
            % Should be done before computing Stacking Covariance Matrix
            % -------------------------------------------------------------%
            p=inputParser;
            p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            RecomputeFlag = p.Results.RecomputeFlag;
            
            % label
            savedir  = [getenv('SamakPath'),sprintf('tritium-data/fit/%s/InitFit/',obj.DataSet)];
            MakeDir(savedir);
            fixPar_init = 'E0 Norm Bkg';
            
            if strcmp(obj.DataSet,'Knm2') && strcmp(obj.ROIFlag,'14keV')
                RoiStr = '_14keV';
            else
                RoiStr = '';
            end
            
            savename = [savedir,sprintf('FitResult_InitModelObj_%s_%s_%s_%.0fpixels%s.mat',...
                strrep(fixPar_init,' ',''),obj.RunData.RunName,obj.DataType,numel(obj.PixList),RoiStr)];
            if strcmp(obj.AnaFlag,'Ring')
                savename = strrep(savename,'.mat',sprintf('Ring%s.mat',obj.RingMerge));
            end
            
            if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
                load(savename,'FitResults','BKG_i');
                nPixels = numel(obj.RunData.MACE_Ba_T); %always gives the correct number of pseudo rings
                obj.ModelObj.BKG_RateSec_i = BKG_i;
                obj.ModelObj.ComputeTBDDS(...
                    'E0_bias',FitResults.par(2),...
                    'B_bias',FitResults.par(3:2+nPixels),...
                    'N_bias',FitResults.par(3+nPixels:3+2*nPixels-1));
                obj.ModelObj.ComputeTBDIS;
                  cprintf('blue','RunAnalysis: Loading (%s) for Initialization from file %s \n',fixPar_init,savename)
            else
                cprintf('blue','RunAnalysis: Initialize Model with fit %s \n',fixPar_init)
                chi2_prev     = obj.chi2;
                fixPar_prev   = obj.fixPar;
                exclDataStart_prev = obj.exclDataStart;
                exclDataStop_prev  = obj.exclDataStop;
                obj.chi2='chi2Stat';
                
                obj.fixPar = ConvertFixPar('freePar',fixPar_init,'nPar',obj.nPar,'nPixels',numel(obj.RunData.MACE_Ba_T));
                
                if strcmp(obj.DataSet,'FirstTritium.katrin')
                    obj.exclDataStart=12;
                    obj.exclDataStop = obj.ModelObj.nqU;
                else
                    % 40 eV range
                    obj.exclDataStart = find((obj.ModelObj.qU)>=18574-40,1);
                    obj.exclDataStop = obj.ModelObj.nqU;
%                    % 90 eV range
%                    obj.exclDataStart = find((obj.ModelObj.qU)>=18574-90,1);

                end
                
               % stat only: Model is not initialized to data statistics yet -> use stat. uncertainty from data 
                TBDIS_Data         = reshape(obj.RunData.TBDIS,[obj.ModelObj.nqU.*numel(obj.ModelObj.MACE_Ba_T),1]);
                obj.FitCM          = diag(TBDIS_Data);
                obj.FitCMFrac      = diag(1./TBDIS_Data);
                obj.FitCMShape     = diag(TBDIS_Data);
                obj.FitCMFracShape = diag(1./TBDIS_Data);
                obj.Fit('InitNB','ON');
                %obj.Fit('InitNB','ON');
               
                FitResults = obj.FitResult;
                BKG_i = obj.ModelObj.BKG_RateSec_i;
                obj.chi2=chi2_prev;
                obj.fixPar = fixPar_prev;
                obj.exclDataStart = exclDataStart_prev;
                obj.exclDataStop  = exclDataStop_prev;
                % save
                MakeDir(savedir);
                save(savename,'FitResults','BKG_i');
            end
        end
        function CM = TestCM_PositiveSemiDefinite(obj)
            % Test of Covariance Matrices are positive semi definite!
            % If not conversion with nearestSPD
            CM = {obj.FitCM; obj.FitCMFrac; obj.FitCMShape; obj.FitCMFracShape};
            if strcmp(obj.AnaFlag,'Ring')
                TBDIS = reshape(obj.ModelObj.TBDIS,[obj.ModelObj.nqU*obj.ModelObj.nRings,1]);
            else
                TBDIS = obj.ModelObj.TBDIS;
            end
            for i=1:numel(CM)
                try mvnrnd(TBDIS,CM{i},1);
                catch
                    fprintf(2,'RunAnalysis:TestCM_PositiveSemiDefinite: WARNING - FitCM not positive semi definite! Apply nearestSPD \n');
                    CM{i} = nearestSPD(CM{i});
                end
            end
        end
    end
    methods % Fit
        
        function Fit(obj,varargin)
            % Fit of the Data/Sim
            % -------------------------------------------------------------%

            p = inputParser;
            p.addParameter('CATS','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('InitNB','OFF', @(x) ismember(x,{'ON','OFF'}));
            p.addParameter('SaveFit','OFF', @(x) ismember(x,{'ON','OFF'}));
            
            p.parse(varargin{:});
            
            CATS       = p.Results.CATS;
            InitNB     = p.Results.InitNB;
            SaveFit    = p.Results.SaveFit;
            
            if isempty(obj.pulls)
                obj.pulls = [inf,inf(1,obj.ModelObj.nPixels*2+7)];
            end
            Data = [obj.RunData.qU, obj.RunData.TBDIS, obj.RunData.TBDISE];
            
            if isempty(obj.i_qUOffset)
                obj.i_qUOffset =zeros(1,obj.ModelObj.nPixels);
            end
            if isempty(obj.i_mTSq)
                obj.i_mTSq =zeros(1,obj.ModelObj.nPixels);
            end
            
            if strcmp(obj.chi2,'chi2Stat') && strcmp(InitNB,'OFF')
                %                 if strcmp(obj.AnaFlag,'StackPixel')
                %                     obj.InitModelObj_Norm_BKG('RecomputeFlag','ON');
                %                 else
                if ~contains(obj.fixPar,'fix 3 ;')
                    obj.InitModelObj_Norm_BKG('RecomputeFlag','OFF');
                end
                
               [StatCM, StatCMFrac] = obj.ComputeCM_StatPNP(varargin);
               obj.FitCM = StatCM;         
               obj.FitCMShape = StatCM;
               obj.FitCMFrac = StatCMFrac; 
               obj.FitCMFracShape = StatCMFrac;
            end
            
            F = FITC('SO',obj.ModelObj,'DATA',Data,'fitter',obj.fitter,...
                'chi2name',obj.chi2,'minuitOpt',obj.minuitOpt,...
                'COVMAT', real(obj.FitCM),'COVMATFrac', real(obj.FitCMFrac),...
                'COVMATShape', real(obj.FitCMShape),'COVMATNorm',obj.FitCMNorm,...
                'COVMATFracShape',real(obj.FitCMFracShape),...
                'pulls',obj.pulls,...
                'pullFlag',obj.pullFlag,...
                'fixPar',obj.fixPar,...
                'exclDataStart',obj.exclDataStart,...
                'exclDataStop',obj.exclDataStop,...
                'i_mnu',obj.i_mnu,...
                'i_Q',obj.i_Q,...
                'i_B',obj.i_B,...
                'i_N',obj.i_N,...
                'i_DTGS',obj.i_DTGS,...
                'i_DTES',obj.i_DTES,...
                'i_HTGS',obj.i_HTGS,...
                'i_HTES',obj.i_HTES,...
                'i_TTGS',obj.i_TTGS,...
                'i_TTES',obj.i_TTES,...
                'i_qUOffset',obj.i_qUOffset,...
                'i_mTSq',obj.i_mTSq,...
                'i_FracTm',0);
            
            obj.FitResult = struct(...
                'par',F.RESULTS{1},....
                'err',F.RESULTS{2},....
                'chi2min',F.RESULTS{3},...
                'errmat',F.RESULTS{4},...
                'dof',F.RESULTS{5});
            
            % Asymmetric uncertainty
            if numel(F.RESULTS)>5 && strcmp(obj.fitter,'minuit') && contains(obj.minuitOpt,'minos')
                  obj.FitResult.errNeg = F.RESULTS{6}';
                  obj.FitResult.errPos = F.RESULTS{7}';
            end
            
            % Calling CATS
            switch CATS
                case 'ON'
                    F.catss; 
                    F.Samakcats_StandResidual_Leverage;
                    F.Samakcats_StandResidualLeverage2D;
                    F.Samakcats_Dffits;
                    F.Samakcats_Dfbetas_Rel;
                    F.Samakcats_Dfbetas_Abs;
                    F.Samakcats_DMSE;
                    F.Samakcats_Delete1Variance;
                    %F.CATSplot; % Old Version
            end
            
            if strcmp(SaveFit,'ON')
                % label
                
                switch obj.AnaFlag
                    case 'StackPixel'
                        savedir = [getenv('SamakPath'),'tritium-data/fit/',obj.DataSet,'/Uniform/'];
                    case 'Ring'
                        savedir = [getenv('SamakPath'),'tritium-data/fit/',obj.DataSet,'/MultiRing/'];
                end
                MakeDir(savedir);
                switch obj.ROIFlag
                    case 'Default'
                        RoiStr = '';
                    case '14keV'
                        RoiStr = '_14keVROI';
                end
                freePar = ConvertFixPar('freePar',obj.fixPar,'nPar',obj.nPar,'nPixels',numel(obj.RunData.MACE_Ba_T),'Mode','Reverse');
                savename = [savedir,sprintf('Fit%s_%s_%s_%.0fpix%s.mat',...
                    obj.RunData.RunName,obj.DataType,freePar,numel(obj.PixList),RoiStr)];
                FitResult = obj.FitResult;
                save(savename,'FitResult');
            end
        end
        function FitBackground(obj,varargin)
            % Fit Background Only
            % Consider all sub-runs with qU >= E0
            % -------------------------------------------------------------%
            
            p = inputParser;
            p.addParameter('UncMPixFlag',true, @(x) islogical(x));
            p.parse(varargin{:});
            UncMPixFlag     = p.Results.UncMPixFlag;
            
            if isempty(obj.pulls)
                obj.pulls = [inf,inf(1,obj.ModelObj.nPixels*2+7)];
            end
            BKGindex = find((obj.ModelObj.qU>=obj.ModelObj.Q_i),1);BKGindex=BKGindex-1;
            Data = [obj.RunData.qU, obj.RunData.TBDIS, obj.RunData.TBDISE];
            
            % Thierry - WARNING - 22/2/2019
            %if strcmp(obj.chi2,'chi2Stat')
            %    obj.ComputeCM_StatPNP(varargin);
            %end
            % End
            fixPar_bkg = ConvertFixPar('freePar','Bkg','nPar',obj.nPar,'nPixels',numel(A.RunData.MACE_Ba_T)); % only fit background
            
            F=FITC('SO',obj.ModelObj,'DATA',Data,'fitter',obj.fitter,...
                'chi2name',obj.chi2,'minuitOpt',obj.minuitOpt,...
                'COVMAT', real(obj.FitCM),'COVMATFrac', real(obj.FitCMFrac),...
                'COVMATShape', real(obj.FitCMShape),'COVMATNorm',obj.FitCMNorm,...
                'pulls',obj.pulls,...
                'pullFlag',obj.pullFlag,...
                'fixPar',fixPar_bkg,...
                'UncMPixFlag',UncMPixFlag,...
                'i_mnu',obj.i_mnu,...
                'i_Q',obj.i_Q,...
                'i_B',obj.i_B,...
                'i_N',obj.i_N,...
                'i_DTGS',obj.i_DTGS,...
                'i_DTES',obj.i_DTES,...
                'i_HTGS',obj.i_HTGS,...
                'i_HTES',obj.i_HTES,...
                'i_TTGS',obj.i_TTGS,...
                'i_TTES',obj.i_TTES,...
                'i_qUOffset',obj.i_qUOffset,...
                'exclDataStart',BKGindex);
            
            
            cprintf('blue','RunAnalysis: FitBackground: qU>%g - index=%g - %g data points\n',...
                obj.RunData.qU(BKGindex),BKGindex,numel(obj.RunData.qU(BKGindex:end)));
            
            obj.FitResult = struct(...
                'par',F.RESULTS{1},....
                'err',F.RESULTS{2},....
                'chi2min',F.RESULTS{3},...
                'errmat',F.RESULTS{4},...
                'dof',F.RESULTS{5});
        end
        function FitMCTwin(obj)
            % work in progress, not finished yet
            if ~strcmp(obj.DataType,'Twin') | ~ismember(obj.DataType,{'FitriumTwin','KafitTwin'})
                fprintf(2,'Data Type has to be Twin');
                return
            end
            % temporary store asimov spectrum
            TBDIS_asimov = obj.RunData.TBDIS;
            nRuns = size(obj.RunData.TBDIS_V); nRuns = nRuns(2);
            nRuns = 10;
            FitResultsPar  = zeros(nRuns,10);
            FitResultsChi2 = zeros(nRuns,1);
            for i=1:nRuns
                obj.RunData.TBDIS = obj.RunData.TBDIS_V(:,i);
                obj.Fit;
                FitResultsPar(i,:) = obj.FitResult.par;
                FitResultsChi2(i)  = obj.FitResult.chi2min;
            end
            obj.RunData.TBDIS = TBDIS_asimov;
            %display
            nhist(FitResultsPar(:,1));
        end
        function ScanResults = GetAsymFitError(obj,varargin)
            % Get asymetric fit uncertainties for 1 fit parameter
            % Method: Scan through Chi^2 profile and calculate Chi^2+1
            % To be run after the regular fit
            p = inputParser;
            p.addParameter('Parameter','mNu',@(x)ismember(x,{'mNu','E0'}));
            p.addParameter('ScanPrcsn',0.01,@(x)isfloat(x)); % start value
            p.addParameter('nFitMax',20,@(x)isfloat(x));   % max number of fits
            p.addParameter('ParScanMax',1,@(x)isfloat(x)); % maximal deviation
            p.addParameter('Mode','Smart',@(x)ismember(x,{'Smart','Uniform'})); %scanning strategy
            p.addParameter('SanityPlot','OFF',@(x)ismember(x,{'ON','OFF'}));
            
            p.parse(varargin{:});
            Parameter  = p.Results.Parameter;
            nFitMax    = p.Results.nFitMax;
            ScanPrcsn  = p.Results.ScanPrcsn;
            ParScanMax = p.Results.ParScanMax;
            Mode       = p.Results.Mode;
            SanityPlot = p.Results.SanityPlot;
            fixPar_prev = obj.fixPar; % save previous setting
            FitResult_prev = obj.FitResult;
            
            % init
            ScanResults = struct('ParScan',zeros(nFitMax,2),...% scan vector
                'par',zeros(obj.nPar,nFitMax,2),... % Fit Parameter
                'err',zeros(obj.nPar,nFitMax,2),... % Error on Fit Parameter
                'chi2min',zeros(nFitMax,2),...
                'dof',zeros(nFitMax,2),...
                'AsymErr',zeros(1,2));
            
            switch Parameter
                case 'mNu'
                    ParScanIndex = 1;
                    ModelPar_i = 'mnuSq_i';
                    xstr = sprintf('{\\itm}_\\nu^2');
                case 'E0'
                    ParScanIndex = 2;
                    ModelPar_i = 'Q_i';
            end
            obj.fixPar = [sprintf('fix %.0f ;',ParScanIndex),fixPar_prev];
            
            if strcmp(Mode,'Smart')
                % start scan: calculate chi2 fit result from regular fit routine
                ScanResults.ParScan(1,1)   = obj.FitResult.par(ParScanIndex);
                ScanResults.ParScan(2,1)   = obj.FitResult.par(ParScanIndex)+1.2*obj.FitResult.err(ParScanIndex);
                ScanResults.ParScan(1,2)   = obj.FitResult.par(ParScanIndex);
                ScanResults.ParScan(2,2)   = obj.FitResult.par(ParScanIndex)-1.2*obj.FitResult.err(ParScanIndex);
            elseif strcmp(Mode,'Uniform')
                % scan in equidistant through chi2 profile
                ScanResults.ParScan(:,1) = linspace(obj.FitResult.par(ParScanIndex),obj.FitResult.par(ParScanIndex)+ParScanMax,nFitMax);
                ScanResults.ParScan(:,2) = sort(linspace(obj.FitResult.par(ParScanIndex)-ParScanMax,obj.FitResult.par(ParScanIndex),nFitMax),'descend');
            end
            
            for d=1:2 % direction: upper and lower uncertainty
                if d==1
                    progressbar('Calculate uppper fit error with scan');
                else
                    progressbar('Calculate lower fit error with scan');
                end
                
                for i=1:nFitMax
                    progressbar(i/nFitMax);
                    obj.ModelObj.(ModelPar_i) = ScanResults.ParScan(i,d);
                    obj.Fit;
                    ScanResults.par(:,i,d) = obj.FitResult.par;
                    ScanResults.err(:,i,d) = obj.FitResult.err;
                    ScanResults.chi2min(i,d) = obj.FitResult.chi2min;
                    ScanResults.dof(i,d)     = obj.FitResult.dof;
                    if strcmp(Mode,'Smart')
                        if i==1
                            chi2min = obj.FitResult.chi2min;
                            continue
                        elseif i==2 % test ParScanMax. if chi2 less than 1 larger than minimal chi2 -> test range not wide enough
                            if ScanResults.chi2min(i,d)-chi2min<1
                                fprintf('ParScanMax of %.2g is too close to best fit \n',ParScanMax);
                                break
                            end
                        end
                        DeltaChi2  = ScanResults.chi2min(i,d)-chi2min;
                        if  abs(DeltaChi2-1)<ScanPrcsn % exit condition
                            fprintf('1 sigma uncertainty found - exit scan \n');
                            ScanResults.ParScan(i+1:end,d)=NaN; % set all not fitted to NaN
                            ScanResults.chi2min(i+1:end,d) = NaN;
                            ScanResults.dof(i+1:end,d)     = NaN;
                            ScanResults.par(:,i+1:end,d)   = NaN;
                            ScanResults.err(:,i+1:end,d)   = NaN;
                            break
                        elseif (DeltaChi2<1 && d==1) || (DeltaChi2>1 && d==2)
                            ParScanLarger= min(ScanResults.ParScan(ScanResults.ParScan(:,d)>ScanResults.ParScan(i,d),d)); % find next higher ParScan
                            %ScanResults.ParScan(i+1,d) =  0.5*(ScanResults.ParScan(i,d)+ParScanLarger); % find middle between current and next higher one
                            ScanResults.ParScan(i+1,d) =  ScanResults.ParScan(i,d) + abs(diff([ScanResults.ParScan(i,d),ParScanLarger])/2);
                        elseif (DeltaChi2>1 && d==1) || (DeltaChi2<1 && d==2)
                            ParScanSmaller= max(ScanResults.ParScan(ScanResults.ParScan(:,d)<ScanResults.ParScan(i,d),d)); % find next smaller mNuSq
                            ScanResults.ParScan(i+1,d) =  ScanResults.ParScan(i,d) - abs(diff([ScanResults.ParScan(i,d),ParScanSmaller])/2);  % find middle between current and next smaller one
                        end
                    elseif strcmp(Mode,'Uniform') && i==1 && d==1
                        chi2min = obj.FitResult.chi2min;
                    end
                end
            end
            
            % find parameters (positive and negative), for which chi2 = chi2min + 1
            ParScanUp = interp1(ScanResults.chi2min(~isnan(ScanResults.chi2min(:,1)),1),...
                ScanResults.ParScan(~isnan(ScanResults.chi2min(:,1)),1),chi2min+1,'spline');
            ParScanDown = interp1(ScanResults.chi2min(~isnan(ScanResults.chi2min(:,2)),2),...
                ScanResults.ParScan(~isnan(ScanResults.chi2min(:,2)),2),chi2min+1,'spline');
            % get asymmetric uncertainties
            ScanResults.AsymErr(1) = diff([ScanResults.ParScan(1,1),ParScanUp]);
            ScanResults.AsymErr(2) = diff([ScanResults.ParScan(1,2),ParScanDown]);
            
            % set everything back to previous setting and fit results
            obj.fixPar = fixPar_prev;
            obj.FitResult = FitResult_prev;
            nPixels = numel(obj.RunData.MACE_Ba_T);
            obj.ModelObj.ComputeTBDDS(...
                'mSq_bias',obj.FitResult.par(1),...
                'E0_bias',obj.FitResult.par(2),...
                'B_bias',obj.FitResult.par(3:2+nPixels),...
                'N_bias',obj.FitResult.par(3+nPixels:3+2*nPixels-1),...
                'DTGS_bias',obj.FitResult.par(2*nPixels+3),...
                'DTES_bias',obj.FitResult.par(2*nPixels+4),...
                'HTGS_bias',obj.FitResult.par(2*nPixels+5),...
                'HTES_bias',obj.FitResult.par(2*nPixels+6),...
                'TTGS_bias',obj.FitResult.par(2*nPixels+7),...
                'TTES_bias',obj.FitResult.par(2*nPixels+8),...
                'qUOffset_bias',obj.FitResult.par((2*nPixels+9):(3*nPixels+8)),...
                'BSlope_bias',obj.FitResult.par(3*nPixels+9),...
                'mTSq_bias',obj.FitResult.par(3*nPixels+10:4*nPixels+9));
            obj.ModelObj.ComputeTBDIS;
            
            if strcmp(SanityPlot,'ON')
                obj.PlotChi2Curve;
            end 
        end
        function out = PlotChi2Curve(obj,varargin)
            p=inputParser;
            p.addParameter('Parameter','mNu',@(x)ismember(x,{'mNu','E0'}));
            p.addParameter('ScanResult','',@(x)isstruct(x));
            p.addParameter('FitResult','',@(x)isstruct(x));
            p.addParameter('HoldOn','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            Parameter   = p.Results.Parameter;
            ScanResult  = p.Results.ScanResult;
            FitResult   = p.Results.FitResult;
            HoldOn      = p.Results.HoldOn;
            
            if isempty(ScanResult)
                fprintf('ScanResults are missing \n');
                return
            end
            
            switch Parameter
                case 'mNu'
                    xstr = sprintf('{{\\it m}_\\nu}^2');
                    xUnit = sprintf('eV^{ 2}');
                    ParIndex = 1;
                case 'E0'
                    xstr = sprintf('{\\it E}_0');
                    xUnit = sprintf('eV');
                    ParIndex = 2;
            end
            nFitMax = size(ScanResult.chi2min,1);
            if strcmp(HoldOn,'OFF')
                f4 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
            end
            mypar  = reshape(ScanResult.ParScan,[nFitMax*2,1]);
            mypar = mypar(~isnan(mypar));
            [mypar,sortI] = sort(mypar);
            mychi2min = reshape(ScanResult.chi2min,[nFitMax*2,1]);
            mychi2min = mychi2min(~isnan(mychi2min));
            mychi2min = mychi2min(sortI);
            if ~strcmp(obj.chi2,'chi2Stat')
                PlotStyle = {'-','LineWidth',3,'Color',rgb('DodgerBlue')};
            else
                PlotStyle = {'-.','LineWidth',3,'Color',rgb('SkyBlue')};
            end
            
            if strcmp(obj.DataType,'Real')
                plot(FitResult.par(ParIndex).*ones(100,1),linspace(0,1e2,1e2),':','LineWidth',2.5,'Color',rgb('Gray'))
                hold on;
            end
            pchi2 = plot(mypar,mychi2min,PlotStyle{:});
            hold on;
            p1 =plot(FitResult.par(ParIndex)+ScanResult.AsymErr(1).*ones(100,1),linspace(0,100,100),...
                ':','LineWidth',2,'Color',rgb('Gray'));
            p2 =plot(FitResult.par(ParIndex)+ScanResult.AsymErr(2).*ones(100,1),linspace(0,100,100),...
                ':','LineWidth',2,'Color',rgb('Gray'));
         %   p3 = plot(linspace(ScanResult.AsymErr(2)+FitResult.par(ParIndex),...
         %       ScanResult.AsymErr(1)+FitResult.par(ParIndex),100),FitResult.chi2min+ones(1,100),...
         %       ':','LineWidth',2,'Color',rgb('Gray'));
           
            PrettyFigureFormat('FontSize',24);
            xlabel(sprintf(' %s (%s)',xstr,xUnit));
            ylabel(sprintf('\\chi^2 (%.0f dof)',ScanResult.dof(1,1) - 1 ));
            
           if abs(FitResult.par(ParIndex))>0.05
               parStr = sprintf('%.3f',FitResult.par(ParIndex));
           else
               parStr = sprintf('%.3f',FitResult.par(ParIndex));
           end
            leg = legend(sprintf('%s = %s (%.3f +%.3f) %s',...
                xstr,parStr,...
                ScanResult.AsymErr(2),ScanResult.AsymErr(1),xUnit));
            leg.EdgeColor = rgb('Silver');
            leg.Location = 'northwest';
             xlim([min(ScanResult.ParScan(:,2)),max(ScanResult.ParScan(:,1))]);
            if strcmp(obj.DataType,'Real')
                ylim([FitResult.chi2min-1 max(max(ScanResult.chi2min))])
            end
            out = {pchi2,p1,p2};
        end
        function GetPlotColor(obj)
            % Real / Twin Color Flag
            if strcmp(obj.DataSet,'FirstTritium.katrin')
                obj.PlotColor = rgb('CadetBlue');
                obj.PlotColorLight = rgb('PowderBlue');
            else
                switch obj.DataType
                    case 'Real'
                        obj.PlotColor = rgb('DodgerBlue');
                        obj.PlotColorLight = rgb('PowderBlue');
                    case {'Twin','FitriumTwin','KafitTwin'}
                        obj.PlotColor = [0.8,0,0];%rgb('Crimson');
                        obj.PlotColorLight =  rgb('LightCoral');
                    case 'Fake'
                        obj.PlotColor = rgb('ForestGreen');
                        obj.PlotColorLight =  rgb('LimeGreen');
                end
            end
        end
        function PlotFit(obj,varargin)
            % Plot Fit Results
            % Call After Fit, otherwise crash...
            % -------------------------------------------------------------%
            p=inputParser;
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF','pdf','png','eps'}));
            p.addParameter('ResidualsFlag','Norm',@(x)ismember(x,{'ON','Norm'}));
            p.addParameter('Mode','Rate',@(x)ismember(x,{'Count','Rate'}));
            p.addParameter('FitResultsFlag','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('LabelFlag','data',@(x)ismember(x,{'data','simulation','FinalKNM1'}));
            p.addParameter('DisplayStyle','Default',@(x)ismember(x,{'Default','PRL'}));
            p.addParameter('DisplayMTD','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('YLimRes','',@(x)isfloat(x)); % limits for norm. residuals
            p.addParameter('XLims','',@(x)isfloat(x));  % x-axis limits
            p.addParameter('ErrorBarScaling',50,@(x)isfloat(x)); % scale error bars in fit
            p.addParameter('Colors','RGB',@(x)ismember(x,{'RGB','BW'})); % color or back and white
            p.addParameter('qUDisp','Rel',@(x)ismember(x,{'Rel','Abs'}));
            p.addParameter('ring',1,@(x)isfloat(x));
            p.addParameter('MaxBkgRange',40,@(x)isfloat(x)); %(eV) maximum range of bkg points shown
            p.addParameter('TickDir','In',@(x)ismember(x,{'In','Out'}));

            p.parse(varargin{:});
            saveplot       = p.Results.saveplot;
            ResidualsFlag  = p.Results.ResidualsFlag;
            Mode           = p.Results.Mode;
            FitResultsFlag = p.Results.FitResultsFlag;
            LabelFlag      = p.Results.LabelFlag;
            YLimRes        = p.Results.YLimRes;
            XLims          = p.Results.XLims;
            Colors         = p.Results.Colors;
            DisplayStyle   = p.Results.DisplayStyle;
            qUDisp         = p.Results.qUDisp;
            ring           = p.Results.ring;
            DisplayMTD     = p.Results.DisplayMTD;
            MaxBkgRange    = p.Results.MaxBkgRange;
            TickDir        = p.Results.TickDir;
            obj.ErrorBarScaling= p.Results.ErrorBarScaling;
            
           BkgEnd = find(obj.RunData.qU(:,1)>=obj.ModelObj.Q+MaxBkgRange,1);
            
            % PlotStyle
            if strcmp(Colors,'RGB')
                obj.GetPlotColor;
            else
                obj.PlotColor = rgb('Silver');
            end
            
            if strcmp(obj.AnaFlag,'StackPixel')
                ring = 1;
            end
            %  else
            MarkerSize = 4;
            LocalFontSize = 17;
            LocalLineWidth = 3;
            %  end
            if contains(obj.DataSet,'FirstTritium')
                LocalFontSize = 9;
                LocalLineWidth = 1;
                MarkerSize = 1.5;
            end
            FitStyleArg = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',MarkerSize,'MarkerEdgeColor',rgb('Black')};
            ResStyleArg = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',MarkerSize,'Color',rgb('Black')};
            
            % Chi2 Flag
            [StatCM, ~] = obj.ComputeCM_StatPNP;
            CMdim = ((ring-1)*obj.ModelObj.nqU+1):(ring*obj.ModelObj.nqU); %ring dimension in covariance matrix
            StatErr = diag(StatCM(CMdim,CMdim));
            %StatErr = obj.RunData.TBDISE(:,ring).^2;
            switch obj.chi2
                case {'chi2Stat','chi2P'}
                    PlotErr = StatErr;
                case {'chi2CM','chi2CMFrac'}
                    PlotErr = diag(obj.FitCM(CMdim,CMdim));
                case 'chi2CMShape'
                    PlotErr = diag(obj.FitCMShape(CMdim,CMdim));
            end
            
            % Spectrum + Fit with Residuals
            fig5 = figure('Renderer','painters');
            if contains(obj.DataSet,'FirstTritium')
                  set(fig5, 'Units', 'centimeters', 'Position', [0.001, 0.001,8.4 ,8.4]);
            else
                 set(fig5, 'Units', 'normalized', 'Position', [0.001, 0.001,0.45, 0.7]);
            end
              
            if strcmp(DisplayStyle,'PRL') || strcmp(DisplayMTD,'ON')
                s1= subplot(4,1,[1 2]);
            else
                s1= subplot(3,1,[1 2]);
            end
            
            switch Mode
                case 'Rate'
                    yfit= obj.ModelObj.TBDIS((obj.exclDataStart:BkgEnd),:)./obj.ModelObj.qUfrac(obj.exclDataStart:BkgEnd,:)./obj.ModelObj.TimeSec;
                    ydata = obj.RunData.TBDIS(obj.exclDataStart:BkgEnd,:)./obj.ModelObj.qUfrac(obj.exclDataStart:BkgEnd,:)./obj.ModelObj.TimeSec;
                    yerrdata = obj.RunData.TBDISE(obj.exclDataStart:BkgEnd,:)./obj.ModelObj.qUfrac(obj.exclDataStart:BkgEnd,:)./obj.ModelObj.TimeSec*obj.ErrorBarScaling;
                    ystr = 'Count rate (cps)';
                case 'Count'
                    yfit = obj.ModelObj.TBDIS(obj.exclDataStart:BkgEnd);
                    ydata = obj.RunData.TBDIS(obj.exclDataStart:BkgEnd);
                    yerrdata = obj.RunData.TBDISE(obj.exclDataStart:BkgEnd)*obj.ErrorBarScaling;
                    ystr = 'Counts';
            end

            if strcmp(qUDisp,'Rel')
                qU = obj.ModelObj.qU(obj.exclDataStart:BkgEnd,:)-obj.ModelObj.Q_i;
                xstr = sprintf('Retarding energy - %.0f (eV)',obj.ModelObj.Q_i);
                textx = -38;
            elseif strcmp(qUDisp,'Abs')
                qU = obj.ModelObj.qU(obj.exclDataStart:BkgEnd,:);
                xstr = sprintf('Retarding energy (eV)');
                myxticks = (round(min(qU),0):20:max(qU));
                textx =min(qU)+0.5;
            end
            
             if strcmp(obj.AnaFlag,'Ring') % show only 1 ring
                yfit     = yfit(:,ring);
                ydata    = ydata(:,ring);
                yerrdata = yerrdata(:,ring);
                qU = qU(:,ring);
             end
                 
            lfit = plot(qU,yfit,'Color',obj.PlotColor,'LineWidth',LocalLineWidth);
            hold on;
            hebar = errorbar(qU,ydata,yerrdata,FitStyleArg{:},'CapSize',0);
            hdata = plot(qU,ydata,'o','MarkerSize',MarkerSize,'MarkerFaceColor',...
                rgb('Black'),'MarkerEdgeColor',rgb('Black'),'LineWidth',LocalLineWidth);
           
            if ismember(DisplayStyle,'PRL')
                PRLFormat;
            elseif contains(obj.DataSet,'FirstTritium')
                FTpaperFormat;
            else
                PrettyFigureFormat('FontSize',LocalFontSize);
            end
            %% subplot1: legend
            d=obj.GetPlotDescription(LabelFlag,ring);
            if strcmp(obj.DataType,'Real')
                if obj.ErrorBarScaling==1
                    datalabel = sprintf(' KATRIN data with %.0f\\sigma error bars',obj.ErrorBarScaling);
                else
                    datalabel = sprintf(' KATRIN data with 1 \\sigma error bars \\times %.0f',obj.ErrorBarScaling);
                end
            else
                if obj.ErrorBarScaling==1
                    datalabel = sprintf(' KATRIN MC data with %.0f\\sigma error bars',obj.ErrorBarScaling);
                else
                    datalabel = sprintf(' KATRIN MC data with 1 \\sigma error bars \\times %.0f',obj.ErrorBarScaling);
                end
            end
                     
             pnone = plot(qU,NaN.*ones(size(qU)),'Color',rgb('White'));
             hold off
             
             if strcmp(obj.AnaFlag,'Ring')
                 myleg = legend([hebar(1) lfit(1)],datalabel,sprintf(' Fit result (pseudo-ring %.0f out of %.0f)',ring,numel(obj.RingPixList)) ,'Location','Northeast');
             else
                 switch FitResultsFlag
                     case 'ON' % show Fit Results in plot
                         myleg = legend([hebar(1) lfit,pnone],datalabel,' Fit result',d.fitleg,'Location','Northeast') ;
                     case 'OFF' %hide Fit Results in plot
                         myleg = legend([hebar(1) lfit],datalabel,' Fit result' ,'Location','Northeast');
                 end
             end
             
             if strcmp(obj.DataSet,'FirstTritium.katrin') && obj.exclDataStart<=7
                 myleg.Location = 'Southwest';
             end
             legend boxoff
            
             %% subplot1: labels and style
             ylabel(ystr);
             set(gca,'yscale','log');%grid on;
             switch Mode
                 case 'Rate'
                     ylim([0.9*min(ydata-obj.ErrorBarScaling*yerrdata),...
                         1.1*max(ydata+obj.ErrorBarScaling*yerrdata)]);
             end
             if contains(obj.DataSet,'FirstTritium')
                 ylim([0.18 1.35*max(obj.ModelObj.TBDIS(obj.exclDataStart:BkgEnd)./obj.ModelObj.qUfrac(obj.exclDataStart:BkgEnd,:)./obj.ModelObj.TimeSec)]);
             end
             
             if strcmp(qUDisp,'Rel')
                 xlim([min(qU)*1.04 max(max(qU))*1.04]);
             elseif strcmp(qUDisp,'Abs')
                 xlim([min(min(qU))*0.9999 max(max(qU))*1.0001]);
             end
             
             ax = gca;
             mypos = ax.Position;
             ax.Position = [mypos(1)+0.05 mypos(2) mypos(3:4)];
             if strcmp(DisplayStyle,'PRL')
                 if strcmp(obj.AnaFlag,'Ring')
                     ylim([0.02 2*max(max(obj.ModelObj.TBDIS(obj.exclDataStart:BkgEnd,ring)./obj.ModelObj.qUfrac(obj.exclDataStart:BkgEnd,ring)./obj.ModelObj.TimeSec(ring)))]);
                 else
                     ylim([0.1 2*max(max(obj.ModelObj.TBDIS(obj.exclDataStart:BkgEnd,ring)./obj.ModelObj.qUfrac(obj.exclDataStart:BkgEnd,ring)./obj.ModelObj.TimeSec(ring)))]);
                     
                 end
                 PRLFormat;
                 myleg.FontSize = get(gca,'FontSize')+4;
                 mylim = ylim;
                 %    text(-57,log(mean(mylim)),'a)','FontSize',get(gca,'FontSize')+4,'FontName',get(gca,'FontName'));
                 text(textx(1),max(mylim)*0.7,'a)','FontSize',get(gca,'FontSize')+4,'FontName',get(gca,'FontName'));
             elseif contains(obj.DataSet,'FirstTritium')
                 FTpaperFormat;
                 set(gca,'FontSize',LocalFontSize);
                 set(get(gca,'YLabel'),'FontSize',LocalFontSize);
                 xlim([min(qU)-6 max(qU)+6]);
                 myleg.FontSize = LocalFontSize;
             else
                 PrettyFigureFormat;
                 set(gca,'FontSize',LocalFontSize);
                 set(gca,'YMinorTick','off');
                 set(gca,'TickLength',[0.01 0.01]);
                 set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);
                 myleg.FontSize = LocalFontSize;
             end
             
             if strcmp(qUDisp,'Abs')
                 xticks(myxticks);
                 ax = gca;
                 ax.XAxis.Exponent = 0;
             end
             
             ax1 = gca;
             
             if strcmp(TickDir,'Out')
                 set(gca,'TickDir','out');
                 ax1.Position = [ax1.Position(1) ax1.Position(2)+0.01, ax1.Position(3:4)];
                 % remove top and right ticks
                 a = gca;
                 set(a,'box','off','color','none')% set box property to off and remove background color
                 b = axes('Position',[ax1.Position(1) ax1.Position(2)+0.01, ax1.Position(3:4)],...
                     'box','on','xtick',[],'ytick',[],'LineWidth',1.5);% create new, empty axes with box but without ticks
                 axes(a)% set original axes as active
                 % linkaxes([a b]) % link axes in case of zooming
             end
             %% residuals
             switch ResidualsFlag
                 case 'Norm' %normalized Residuals
                     if strcmp(DisplayStyle,'PRL') || strcmp(DisplayMTD,'ON')
                         s2= subplot(4,1,3);
                     else
                         s2= subplot(3,1,3);
                     end
                     
                     DataStat = [qU,zeros(size(qU)),sqrt(StatErr(obj.exclDataStart:BkgEnd)./PlotErr(obj.exclDataStart:BkgEnd))];
                     
                     if strcmp(obj.chi2,'chi2Stat')
                         [lstat , pstat]  = boundedline(DataStat(:,1),DataStat(:,2),DataStat(:,3));
                         lstat.LineStyle= '--'; lstat.Color = rgb('DarkSlateGray'); lstat.LineWidth=LocalLineWidth;
                         pstat.FaceColor = obj.PlotColorLight;
                         
                     elseif ~strcmp(obj.chi2,'chi2Stat') && numel(hdata)==1
                         DataSys = [qU,zeros(numel(qU),1), ones(numel(qU),1)];
                         [l,p] = boundedline(DataSys(:,1),DataSys(:,2),DataSys(:,3),...
                             '-b*',DataStat(:,1),DataStat(:,2),DataStat(:,3),'--ro');
                         lsys = l(1);  lstat = l(2);
                         psys = p(1);  pstat = p(2);
                         if strcmp(Colors,'RGB')
                             psys.FaceColor  = obj.PlotColor; %psys.FaceAlpha=0.3;
                             pstat.FaceColor = obj.PlotColorLight;
                             lstat.Color     = rgb('Silver');
                         else
                             pstat.FaceColor = rgb('Black')';%obj.PlotColor;
                             psys.FaceColor = rgb('Black'); psys.FaceAlpha=0.4;
                             lstat.Color = rgb('Black');
                         end
                         lsys.LineStyle= 'none';
                         lstat.LineStyle= '--';  lstat.LineWidth=LocalLineWidth;
                         lstat.Marker = 'none'; lsys.Marker = 'none';
                     end
                     
                     yres = (obj.RunData.TBDIS(obj.exclDataStart:BkgEnd,ring)-obj.ModelObj.TBDIS(obj.exclDataStart:BkgEnd,ring))...
                         ./sqrt(PlotErr(obj.exclDataStart:BkgEnd));
                     
                     hold on;
                     pRes = errorbar(qU,yres,...
                         zeros(numel(qU),1),...
                         ResStyleArg{:});
                     pRes.CapSize = 0;
                    if ~strcmp(obj.chi2,'chi2Stat')
                        leg = legend([pstat psys],'Stat.','Stat. and syst.','Location','Northeast'); %hsyst
                    elseif strcmp(obj.chi2,'chi2Stat')
                        leg = legend(pstat,'stat','Location','Northeast'); %hsyst
                    end
                    legend('boxoff');
                    
                    leg.NumColumns = 2;
                    hold off;
                    % xlabel(sprintf('retarding energy - %.1f (eV)',obj.ModelObj.Q_i),'FontSize',LocalFontSize);
                    if strcmp(DisplayStyle,'PRL') || strcmp(DisplayMTD,'ON')
                        % xlabel(sprintf('retarding energy - %.0f (eV)',obj.ModelObj.Q_i));
                        if ismember(DisplayStyle,'PRL')
                            PRLFormat;
                        else
                            PrettyFigureFormat
                        end
                        leg.FontSize = get(gca,'FontSize')+4;
                        %pstat.delete; psys.delete;
                        lstat.Color = rgb('DarkGray');
                    elseif contains(obj.DataSet,'FirstTritium')
                        xlabel(xstr,'FontSize',LocalFontSize);
                        leg.FontSize = LocalFontSize;
                    else
                        xlabel(xstr,'FontSize',LocalFontSize);
                        leg.FontSize = LocalFontSize+4;
                    end
                    
                    ylabel(sprintf('Residuals (\\sigma)'));
                    if strcmp(qUDisp,'Rel')
                        xlim([floor(min(qU)) max(qU)*1.04]);
                    elseif strcmp(qUDisp,'Abs')
                        xlim([min(qU)*0.9999 max(qU)*1.0001]);
                    end
                    ymin = 1.1*min(yres);
                    ymax = 1.1*max(yres);
                   if ymin<=-1 || ymax<=1
                    ylim([ymin,ymax]);%ylim([-2 2])
                   else
                       ylim([ymin,ymax]);
                   end
                   
                   if ismember(DisplayStyle,'PRL')
                       PRLFormat;
                       xlim([-40 max(qU)+1]);
                   elseif contains(obj.DataSet,'FirstTritium')
                       FTpaperFormat;
                       set(gca,'FontSize',LocalFontSize);
                       set(get(gca,'YLabel'),'FontSize',LocalFontSize);
                       set(get(gca,'XLabel'),'FontSize',LocalFontSize);
                   else
                       PrettyFigureFormat; set(gca,'FontSize',LocalFontSize);
                       set(gca,'YMinorTick','off');
                       set(gca,'XMinorTick','off');
                       set(gca,'TickLength',[0.01 0.01]);
                       set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);
                       set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
                   end
                    if numel(pRes)>1
                        for i=1:numel(pRes)
                            pRes(i).Color = rgb('Black'); % Warning Cut/Past Bug before
                        end
                    end
                    if ~isempty(YLimRes)
                        ylim([min(YLimRes) max(YLimRes)]);
                    end
                case 'ON' % not normalized Residuals
                    s2=subplot(2,1,2);
                    Residuals = (obj.RunData.TBDIS-obj.ModelObj.TBDIS);
                    if ~strcmp(obj.chi2,'chi2Stat') && numel(hdata)==1
                        [lsys,psys] = boundedline(qU,Residuals,sqrt(PlotErr-obj.RunData.TBDISE.^2),...
                            'alpha','transparency',0.4,'cmap',obj.PlotColor);
                        lsys.LineStyle= '--';
                    elseif strcmp(obj.chi2,'chi2Stat')
                        plot(qU,zeros(obj.ModelObj.nqU,1),'--','Color',obj.PlotColor);
                    end
                    hold on;
                    pResi = errorbar(qU,Residuals,obj.RunData.TBDISE,...
                        'o','Color',rgb('DimGray'));
                    if ~strcmp(obj.chi2,'chi2Stat')
                        leg = legend([pResi psys],'stat','sys (CM)','Location','Southeast'); %hsyst
                        leg = legend(pResi,'stat','Location','Southeast'); %hsyst
                    end
                    PrettyFigureFormat;set(gca,'FontSize',LocalFontSize);
                    legend('boxoff');
                    leg.FontSize = LocalFontSize;
                    hold off;
                    xlabel(xstr,'FontSize',LocalFontSize);
                    ylabel('residuals');
                    if strcmp(qUDisp,'Rel')
                        xlim([min(qU)*1.04 max(qU)*1.04]);
                    elseif strcmp(qUDisp,'Abs')
                        xlim([min(qU)*0.9999 max(qU)*1.0001]);
                    end
                    
                   %ylim([1.4*min((obj.RunData.TBDIS-obj.ModelObj.TBDIS)./sqrt(diag(PlotCM)))  1.4*max((obj.RunData.TBDIS-obj.ModelObj.TBDIS)./sqrt(diag(PlotCM)))]);%2*max((Data(:,2)-obj.ModelObj.TBDIS)./sqrt(diag(FitCM)))])
                    set(gca,'YMinorTick','off');
                    set(gca,'XMinorTick','off');
                    set(gca,'TickLength',[0.01 0.01]);
                    set(gca,'FontSize',LocalFontSize);
             end
             
             if strcmp(qUDisp,'Abs')
                 xticks(myxticks);
                 ax = gca;
                 ax.XAxis.Exponent = 0;
             end
            
             
             if strcmp(qUDisp,'Rel')
                 xlim([min(qU)*1.04 max(qU)*1.04]);
             elseif strcmp(qUDisp,'Abs')
                 xticks(myxticks);
                 ax = gca;
                 ax.XAxis.Exponent = 0;
                 xlim([min(qU)*0.9999 max(qU)*1.0001]);
             end
             tmp = get(ax1,'YLabel');
             tmp2 = get(gca,'YLabel');
             set(get(gca,'YLabel'),'Position',[tmp.Position(1),tmp2.Position(2),tmp2.Position(3)]);
             ax = gca; ax2 = gca;
             mypos = ax.Position;
             ax.Position = [mypos(1)+0.05 mypos(2)+0.015 mypos(3:4)];
             
              if strcmp(TickDir,'Out')
                 set(gca,'TickDir','out');
                 % remove top and right ticks
                 a = gca;
                 set(a,'box','off','color','none')% set box property to off and remove background color
                 b = axes('Position',ax.Position,...
                     'box','on','xtick',[],'ytick',[],'LineWidth',1.5);% create new, empty axes with box but without ticks
                 axes(a)% set original axes as active
              end
             
             linkaxes([s1,s2],'x');
             if ~isempty(XLims)
                    xlim([min(XLims),max(XLims)])
                end
             if strcmp(DisplayStyle,'PRL')  || strcmp(DisplayMTD,'ON')
                 mylim = ylim;
                 %  text(-57,mean(mylim),'b)','FontSize',get(gca,'FontSize')+4,'FontName',get(gca,'FontName'));
                text(textx(1),max(mylim)*0.65,'b)','FontSize',get(gca,'FontSize')+4,'FontName',get(gca,'FontName'));
                   
                 s3= subplot(4,1,4);
                 b1 = bar(qU,obj.RunData.qUfrac(obj.exclDataStart:BkgEnd,ring).*obj.RunData.TimeSec(ring)./(60*60));
                 b1.FaceColor = obj.PlotColor;
                 b1.EdgeColor = obj.PlotColor;
                 xlabel(xstr);
                 ylabel('Time (h)')
                 ylim([0,max(obj.RunData.qUfrac(obj.exclDataStart:obj.exclDataStop,ring).*obj.RunData.TimeSec(ring))/(60*60)*1.05]);
                 if strcmp(DisplayStyle,'PRL')
                     PRLFormat;
                 else
                     PrettyFigureFormat;
                 end
                 tmp = get(ax1,'YLabel');
                 tmp2 = get(gca,'YLabel');
                 set(get(gca,'YLabel'),'Position',[tmp.Position(1),tmp2.Position(2),tmp2.Position(3)]);
                 b1.BarWidth=0.6;
                
                if strcmp(qUDisp,'Abs')
                    xticks(myxticks);
                    ax = gca;
                    ax.XAxis.Exponent = 0;
                end
                linkaxes([s1,s2,s3],'x');
              
                ax = gca;
                mypos = ax.Position;
                ax.Position = [mypos(1)+0.05 mypos(2)+0.01 mypos(3:4)];

                if ~isempty(XLims)
                    xlim([min(XLims),max(XLims)])
                end
                % linkaxes([s1,s2,s3],'x');
                mylim = ylim;
                %  text(-57,mean(mylim),'c)','FontSize',get(gca,'FontSize')+4,'FontName',get(gca,'FontName'));
                text( textx(1),max(mylim)*0.8,'c)','FontSize',get(gca,'FontSize')+4,...
                    'FontName',get(gca,'FontName'),'FontWeight',get(gca,'FontWeight'));
                if strcmp(TickDir,'Out')
                    set(gca,'TickDir','out');
                    
                    ax1.Position = [ax1.Position(1) ax1.Position(2)+0.01, ax1.Position(3:4)];
                    % remove top and right ticks
                    a = gca;
                    set(a,'box','off','color','none')% set box property to off and remove background color
                    b = axes('Position',a.Position,...
                        'box','on','xtick',[],'ytick',[],'LineWidth',1.5);% create new, empty axes with box but without ticks
                    axes(a)% set original axes as active
                end
             else
                 xlabel(xstr);
             end
             
             if contains(obj.DataSet,'FirstTritium')
                mypos = leg.Position;
                leg.Position = [mypos(1)+0.017,mypos(2)+0.016,mypos(3:4)];
             end
             
            % if obj.exclDataStart==1
            %     set(get(gca,'YLabel'),'Position',[tmp.Position(1),tmp2.Position(2),tmp2.Position(3)]);
            %     set(get(ax2,'YLabel'),'Position',[tmp.Position(1),tmp2.Position(2)-YLimRes(2)/2,tmp2.Position(3)]);
            % end
             
             if ~strcmp(saveplot,'OFF')
                 if strcmp(Colors,'BW')
                     savename = [d.savename,'_BW'];
                 else
                     savename = d.savename;
                 end
                 if strcmp(saveplot,'ON') || strcmp(saveplot,'png')
                     % save_name = sprintf('./plots/FitRun%u_%s-excl%u.pdf',obj.RunNr,obj.chi2, obj.exclDataStart);
                     % if strcmp(obj.AnaFlag,'Ring')
                     %     print(fig5,sprintf('./plots/%s_%.0feVrange_%.0f-nRings_%.png',obj.ModelObj.TD,round(obj.ModelObj.qU(obj.exclDataStart)-obj.ModelObj.Q_i),obj.nRings),'-dpng','-r450');
                     % else
                     %print(fig5,[d.savename,'.',saveplot],['-d',saveplot],'-r450');
                     export_fig(fig5,[savename,'.',saveplot],'-q101','-m3');
                     % end
                 else
                     %save_name = sprintf('./plots/FitRun%u_%s-excl%u.%s',obj.RunNr,obj.chi2, obj.exclDataStart,saveplot);
                     export_fig(fig5,[savename,'.pdf']);
                     %publish_figurePDF(fig5,[savename,'.pdf']);
                 end
             end
        end
        
        function PlotResidualsMultiRing(obj,varargin)
            % Plot Fit Results for MultiRing Fit
            % Residuals as function of qU for each Ring + MTD
            % Designed for 4 Pseudo-Rings
            % Call After Fit, otherwise crash...
            % -------------------------------------------------------------%
            p=inputParser;
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF','pdf','png','eps'}));
            p.addParameter('LabelFlag','data',@(x)ismember(x,{'data','simulation','FinalKNM1'})); % not used
            p.addParameter('DisplayStyle','Default',@(x)ismember(x,{'Default','PRL'})); % tested
            p.addParameter('DisplayMTD','ON',@(x)ismember(x,{'ON','OFF'})); % tested
            p.addParameter('YLimRes','',@(x)isfloat(x)); % limits for norm. residuals - not tested
            p.addParameter('XLims','',@(x)isfloat(x));   % x-axis limits - not tested
            p.addParameter('Colors','RGB',@(x)ismember(x,{'RGB','BW'})); % color or back and white - no effect here...
            p.addParameter('qUDisp','Rel',@(x)ismember(x,{'Rel','Abs'}));
            p.addParameter('MaxBkgRange',100,@(x)isfloat(x)); %(eV) maximum range of bkg points shown
            
            p.parse(varargin{:});
            saveplot       = p.Results.saveplot;
            LabelFlag      = p.Results.LabelFlag;
            YLimRes        = p.Results.YLimRes;
            XLims          = p.Results.XLims;
            Colors         = p.Results.Colors;
            DisplayStyle   = p.Results.DisplayStyle;
            qUDisp         = p.Results.qUDisp;
            DisplayMTD     = p.Results.DisplayMTD;
            MaxBkgRange    = p.Results.MaxBkgRange;
            
            BkgEnd = find(obj.RunData.qU(:,1)>=obj.ModelObj.Q+MaxBkgRange,1);
            
            % PlotStyle
            if strcmp(Colors,'RGB')
                obj.GetPlotColor;
            else
                obj.PlotColor = rgb('Silver');
            end
            
            MarkerSize     = 4;
            LocalFontSize  = 17;
            LocalLineWidth = 2;
            
            if ~strcmp(obj.AnaFlag,'Ring') 
                fprintf(2,'PlotResidualsMultiRing: not available for Uniform Fits');
                return;
            end
            
            if contains(obj.DataSet,'FirstTritium')
                fprintf(2,'PlotResidualsMultiRing: not available for First Tritium Analysis');
                return;
            end
            
            ResStyleArg = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',MarkerSize,'Color',rgb('Black')};
            
            % Chi2 Flag
            for r=1:obj.nRings
                [StatCM, ~] = obj.ComputeCM_StatPNP;
                CMdim = ((r-1)*obj.ModelObj.nqU+1):(r*obj.ModelObj.nqU); %ring dimension in covariance matrix
                StatErr(r,:) = diag(StatCM(CMdim,CMdim));
                switch obj.chi2
                    case {'chi2Stat','chi2P'}
                        PlotErr(r,:) = StatErr(r,:);
                    case {'chi2CM','chi2CMFrac'}
                        PlotErr(r,:) = diag(obj.FitCM(CMdim,CMdim));
                    case 'chi2CMShape'
                        PlotErr(r,:) = diag(obj.FitCMShape(CMdim,CMdim));
                end
            end
                        
            RingPreFix = 'Pseudo-';
            % Spectrum + Fit with Residuals
            fig5 = figure('Renderer','painters');
            if obj.nRings<=4
                set(fig5, 'Units', 'normalized', 'Position', [0.001, 0.001,0.6, 0.7]);
                nCol = 1;
                nRow = obj.nRings;
            else
                set(fig5, 'Units', 'normalized', 'Position', [0.001, 0.001,0.7, 0.7]);
                nCol = 2;
                nRow = obj.nRings/2;
                if obj.nRings==12
                    RingPreFix = '';
                end
            end
            
            if strcmp(qUDisp,'Rel')
                qU = obj.ModelObj.qU(obj.exclDataStart:BkgEnd,:)-obj.ModelObj.Q_i;
                xstr = sprintf('Retarding energy - %.0f (eV)',obj.ModelObj.Q_i);
                textx = -38;
            elseif strcmp(qUDisp,'Abs')
                qU = obj.ModelObj.qU(obj.exclDataStart:BkgEnd,:);
                xstr = sprintf('Retarding energy (eV)');
                myxticks = (round(min(qU),0):20:max(qU));
                textx =min(qU)+0.5;
            end
            
            if strcmp(obj.AnaFlag,'Ring') % show only 1 ring
                ring=1;
                qU = qU(:,ring);
            end
            
            if strcmp(qUDisp,'Rel')
                xlim([min(qU)*1.04 max(max(qU))*1.04]);
            elseif strcmp(qUDisp,'Abs')
                xlim([min(min(qU))*0.9999 max(max(qU))*1.0001]);
            end
            
            if strcmp(qUDisp,'Abs')
                xticks(myxticks);
                ax = gca;
                ax.XAxis.Exponent = 0;
            end
            
            % %            ax1 = gca;
            
            for r=1:obj.nRings
 
                if strcmp(DisplayStyle,'PRL') || strcmp(DisplayMTD,'ON')
                    s(r)= subplot(nRow+1,nCol,r);
                else
                    s(r)= subplot(nRow,nCol,r);
                end
                
                DataStat = [qU,zeros(size(qU)),(sqrt(StatErr(r,obj.exclDataStart:BkgEnd)./PlotErr(r,obj.exclDataStart:BkgEnd)))'];
                
                if strcmp(obj.chi2,'chi2Stat')
                    [lstat , pstat]  = boundedline(DataStat(:,1),DataStat(:,2),DataStat(:,3));
                    lstat.LineStyle= '--'; lstat.Color = rgb('DarkSlateGray'); lstat.LineWidth=LocalLineWidth;
                    pstat.FaceColor = obj.PlotColorLight;
                    
                elseif ~strcmp(obj.chi2,'chi2Stat') %&& numel(hdata)==1
                    DataSys = [qU,zeros(numel(qU),1), ones(numel(qU),1)];
                    [l,p] = boundedline(DataSys(:,1),DataSys(:,2),DataSys(:,3),...
                        '-b*',DataStat(:,1),DataStat(:,2),DataStat(:,3),'--ro');
                    lsys = l(1);  lstat = l(2);
                    psys = p(1);  pstat = p(2);
                    
                    pstat.FaceColor = rgb('PowderBlue');
                    if strcmp(Colors,'RGB')
                        psys.FaceColor =obj.PlotColor; %psys.FaceAlpha=0.3;
                        lstat.Color = rgb('Silver');
                    else
                        pstat.FaceColor = rgb('Black')';%obj.PlotColor;
                        psys.FaceColor = rgb('Black'); psys.FaceAlpha=0.4;
                        lstat.Color = rgb('Black');
                    end
                    lsys.LineStyle= 'none';
                    lstat.LineStyle= '--';  lstat.LineWidth=LocalLineWidth;
                    lstat.Marker = 'none'; lsys.Marker = 'none';
                    leg.FontSize = get(gca,'FontSize')+4;
                end
                
                yres = (obj.RunData.TBDIS(obj.exclDataStart:BkgEnd,r)-obj.ModelObj.TBDIS(obj.exclDataStart:BkgEnd,r))...
                    ./sqrt(PlotErr(r,obj.exclDataStart:BkgEnd)');
                
                hold on;
                pRes = errorbar(qU,yres,...
                    zeros(numel(qU),1),...
                    ResStyleArg{:});
                pRes.CapSize = 0;
                
                pleg = plot(0,0,'Color',rgb('White'));
                
                if ~strcmp(obj.chi2,'chi2Stat')
                    leg = legend([pleg,pstat psys],sprintf('%sRing %0.f',RingPreFix,r),'Stat.',sprintf('Stat. and syst.'),'Location','Northeast'); %hsyst
                elseif strcmp(obj.chi2,'chi2Stat')
                    leg = legend([pleg,pstat],sprintf('%sRing %0.f Stat.',RingPreFix,r),'Location','Northeast'); %hsyst
                end
                legend('boxoff');
                %leg.EdgeColor = 'none';
                % set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.5]));
                
                leg.NumColumns = 3;
                
                hold off;
                % xlabel(sprintf('retarding energy - %.1f (eV)',obj.ModelObj.Q_i),'FontSize',LocalFontSize);
                if strcmp(DisplayStyle,'PRL') || strcmp(DisplayMTD,'ON')
                    % xlabel(sprintf('retarding energy - %.0f (eV)',obj.ModelObj.Q_i));
                    if ismember(DisplayStyle,'PRL')
                        PRLFormat;
                    else
                        PrettyFigureFormat
                    end
                    %pstat.delete; psys.delete;
                    lstat.Color = rgb('DarkGray');
                else
                    if r==4
                        xlabel(xstr,'FontSize',LocalFontSize);
                        leg.FontSize = get(gca,'FontSize')+4;
                    end
                end
                
                %leg.FontSize = 16;%get(gca,'FontSize')-2;
                ylabel(sprintf('Res. (\\sigma)'));
                if strcmp(qUDisp,'Rel')
                    xlim([floor(min(qU))-4 max(qU)*1.04]);
                elseif strcmp(qUDisp,'Abs')
                    xlim([min(qU)*0.9999 max(qU)*1.0001]);
                end
                ymin = 1.1*min(yres);
                ymax = 1.8*max(yres);
                if ymin<=-1 || ymax<=1
                    ylim([-2.5 2.5]);% ylim([ymin,ymax]);%ylim([-2 2])
                else
                    ylim([ymin,ymax]);
                end
                
                if ismember(DisplayStyle,'PRL')
                    PRLFormat;
                    xlim([-40 max(qU)+1]);
                else
                    PrettyFigureFormat; set(gca,'FontSize',LocalFontSize);
                    set(gca,'YMinorTick','off');
                    set(gca,'XMinorTick','off');
                    set(gca,'TickLength',[0.01 0.01]);
                    set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);
                    set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
                end
                
                if ~isempty(YLimRes)
                    ylim([min(YLimRes) max(YLimRes)]);
                end
                
                if strcmp(qUDisp,'Abs')
                    xticks(myxticks);
                    ax = gca;
                    ax.XAxis.Exponent = 0;
                end
                
                if strcmp(qUDisp,'Rel')
                    xlim([min(qU)*1.04 max(qU)*1.04]);
                elseif strcmp(qUDisp,'Abs')
                    xticks(myxticks);
                    ax = gca;
                    ax.XAxis.Exponent = 0;
                    xlim([min(qU)*0.9999 max(qU)*1.0001]);
                end
                
               
                    ax1 = gca;
                    tmp = get(ax1,'YLabel');
                    tmp2 = get(gca,'YLabel');
                    set(get(gca,'YLabel'),'Position',[tmp.Position(1),tmp2.Position(2),tmp2.Position(3)]);
                    ax = gca;
                    mypos = ax.Position;
                    if obj.nRings<=4
                        ax.Position = [mypos(1)+0.05 mypos(2)+0.01 mypos(3:4)];
                    else
                        if mod(r,2) % if uneven (left side)
                            ax.Position = [mypos(1)-0.05 mypos(2)-0.01 mypos(3)+0.08,mypos(4)+0.03]; 
                        else
                            ax.Position = [mypos(1) mypos(2)-0.01 mypos(3)+0.08,mypos(4)+0.03];
                        end
                        
                    end
                    if ~isempty(XLims)
                        xlim([min(XLims),max(XLims)])
                    end
              
                
                % MTD - optional
                if strcmp(DisplayStyle,'PRL')  || strcmp(DisplayMTD,'ON')
                    mylim = ylim;
                    %  text(-57,mean(mylim),'b)','FontSize',get(gca,'FontSize')+4,'FontName',get(gca,'FontName'));
                    %text(textx,max(mylim)*0.65,'b)','FontSize',get(gca,'FontSize')+4,'FontName',get(gca,'FontName'));
                    
                    s(5)= subplot(obj.nRings+1,1,obj.nRings+1);
                    ring=1;
                    b1 = bar(qU,obj.RunData.qUfrac(obj.exclDataStart:BkgEnd,ring).*obj.RunData.TimeSec(ring)./(60*60));
                    b1.FaceColor = obj.PlotColor;
                    b1.EdgeColor = obj.PlotColor;
                    xlabel(xstr);
                    ylabel('Time (h)')
                    ylim([0,max(obj.RunData.qUfrac(obj.exclDataStart:obj.exclDataStop,ring).*obj.RunData.TimeSec(ring))/(60*60)*1.05]);
                    if strcmp(DisplayStyle,'PRL')
                        PRLFormat;
                    else
                        PrettyFigureFormat;
                    end
                    tmp = get(ax1,'YLabel');
                    tmp2 = get(gca,'YLabel');
                    set(get(gca,'YLabel'),'Position',[tmp.Position(1),tmp2.Position(2),tmp2.Position(3)]);
                    b1.BarWidth=0.6;
                    
                    if strcmp(qUDisp,'Abs')
                        xticks(myxticks);
                        ax = gca;
                        ax.XAxis.Exponent = 0;
                    end
                    
                  
                    ax = gca;
                    mypos = ax.Position;
                    ax.Position = [mypos(1)+0.05 mypos(2)+0.01 mypos(3:4)];
                    
                    if ~isempty(XLims)
                        xlim([min(XLims),max(XLims)])
                    end
                else
                    if obj.nRings==4 && r==4
                        xlabel(xstr);
                    elseif nCol==2 && (r==obj.nRings  || r==obj.nRings-1 )
                         xlabel(xstr);
                    end
                end
            end
            
            if numel(s)==4
            linkaxes([s(1),s(2),s(3),s(4)],'x');
            end
            if numel(s)==5
            linkaxes([s(1),s(2),s(3),s(4),s(5)],'x');
            end
            
            if ~strcmp(saveplot,'OFF')
                d = obj.GetPlotDescription;
                if strcmp(Colors,'BW')
                    savename = [d.savename,'_PSRresiduals_BW'];
                else
                    savename = [d.savename,'_PSRresiduals'];
                end
                if strcmp(saveplot,'ON') || strcmp(saveplot,'png')
                    export_fig(fig5,[savename,'.',saveplot],'-q101','-m3');
                    % end
                else
                    export_fig(fig5,[savename,'.pdf']);
                end
            end
        end
        
        
        function PlotFitResultRings(obj,varargin)
            p=inputParser;
            p.addParameter('PlotPar','qUOffset',@(x)ischar(x));
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('linFit','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('PlotMode','Rel',@(x)ismember(x,{'Rel','Abs'}));
            
            p.parse(varargin{:});
            
            PlotPar    = p.Results.PlotPar;
            SavePlot   = p.Results.SavePlot;
            linFitFlag = p.Results.linFit;
            PlotMode   = p.Results.PlotMode;
            
            if ~strcmp(obj.AnaFlag,'Ring')
                fprintf('only for ring analysis \n');
                return;
            end
            
           
        end
        function Description=GetPlotDescription(obj,LabelFlag,ring)
            if ~exist('ring','var')
                ring = 0;
            end
            % Gather Relevant information to properly legend the plots
            % -------------------------------------------------------------%
            str = './plots/';
            MakeDir(str);
            if strcmp(obj.AnaFlag,'Ring')
                Description = struct('dataleg','Data','fitleg','Fit result',...
                    'title','Ring fit','minititle','');
                if contains(obj.fixPar,'fix 1 ')  % neutrino mass fixed
                    nu_leg = sprintf(' {\\itm}_\\nu^2 fixed \n');
                else
                    nu_leg = sprintf(' {\\itm}_\\nu^2 = %.2g \\pm %.2f eV^2 \n',...
                        (obj.ModelObj.mnuSq_i+obj.FitResult.par(1)),obj.FitResult.err(1));
                end
                e0_leg = sprintf(' {\\itE}_0^{fit} = %.2f \\pm  %.2f eV',obj.ModelObj.Q_i+obj.FitResult.par(2),obj.FitResult.err(2));
                Description.fitleg = [nu_leg,e0_leg];
                  ErangeLabel = round(obj.ModelObj.Q_i-obj.ModelObj.qU(obj.exclDataStart));
          
                Description.savename = sprintf('./plots/FitRun%s_%s_%.0feV_%s_%.0fsigma_Ring%s-No%.0f',...
                    obj.ModelObj.TD,obj.chi2,ErangeLabel,obj.AnaFlag,obj.ErrorBarScaling,obj.RingMerge,ring);
           
                if strcmp(obj.chi2,'chi2Stat')
                    Description.chi2 = 'stat. only';
                else
                     Description.chi2 = 'stat. and syst.';
                end
                return
            end
            
            if contains(obj.fixPar,'fix 1 ')  % neutrino mass fixed
                nu_leg = sprintf('{\\itm}_\\nu^2 fixed');
            else                          % neutrino mass is being fitted
                if strcmp(LabelFlag,'FinalKNM1')
                    nu_leg = sprintf('{\\itm}_\\nu^2     =  - %.2g _{-  %.2f}^{+ %.2f} eV^2',...
                        0.98,1.06,0.89);
                else
                    nu_leg = sprintf('{\\itm}_\\nu^2 \t \t = \t %.2g \t\\pm %.2f eV^2',...
                        ((obj.ModelObj.mnuSq_i+obj.FitResult.par(1))),obj.FitResult.err(1));
                end
            end
            
            if strcmp(LabelFlag,'FinalKNM1')
                chi2_leg = sprintf('\\chi^2      = \t %.1f (%.0f dof)', 21.4,23);
                e0_leg = sprintf('{\\itE}_0^{fit}  =  %.2f _{-  %.2f}^{+ %.2f} eV',18573.73,0.06,0.06);
                b_leg = sprintf('{\\itB}        =  %.1f _{-  %.1f}^{+ %.1f} mcps',292.3,0.7,0.7);
            else
                chi2_leg = sprintf('\\chi2 / dof \t \t = \t %.1f/%.0f ',obj.FitResult.chi2min,obj.FitResult.dof);
                e0_leg = sprintf('{\\itE}_0^{fit}  =  %.2f \\pm  %.2f eV',obj.ModelObj.Q_i+obj.FitResult.par(2),obj.FitResult.err(2));
                b_leg = sprintf('{\\itB}        =  %.1f \\pm %.1f mcps',(obj.ModelObj.BKG_RateSec_i+obj.FitResult.par(3))*1e3,obj.FitResult.err(3)*1e3);
            end
                myfitleg1 =  sprintf(' %s \n %s \n %s \n %s',...
                    chi2_leg,... %chi2
                    nu_leg,...%nu-mass
                    e0_leg,...%endpoint
                    b_leg);%,...%background
                
            if ~contains(obj.fixPar,sprintf('fix %.0f ;fix %.0f',3+2*obj.ModelObj.nPixels,3+2*obj.ModelObj.nPixels+1)) % in case DT FSD Probabilities are fitted
                myfitleg2=sprintf('\n DT GS-Prob=%.2f%% \t\\pm %.2f%% \n DT ES-Prob=%.2f%% \t\\pm %.2f%%',...
                    100*(obj.ModelObj.DTNormGS_i + obj.FitResult.par(2*obj.ModelObj.nPixels+3)),100*obj.FitResult.err(2*obj.ModelObj.nPixels+3),...
                    100*(obj.ModelObj.DTNormES_i - obj.FitResult.par(2*obj.ModelObj.nPixels+3)), 100*obj.FitResult.err(2*obj.ModelObj.nPixels+4));
            else
                myfitleg2=sprintf('');
            end
            
            if ~contains(obj.fixPar,sprintf('fix %.0f ;fix %.0f',3+2*obj.ModelObj.nPixels+2,3+2*obj.ModelObj.nPixels+3)) % in case HT FSD Probabilities are fitted
                myfitleg3=sprintf('\n HT GS-Prob=%.2f%% \t\\pm %.2f%% \n HT ES-Prob=%.2f%% \t\\pm %.2f%%',...
                    100*(obj.ModelObj.HTNormGS_i + obj.FitResult.par(4*obj.ModelObj.nPixels+3)),100*obj.FitResult.err(4*obj.ModelObj.nPixels+3),...
                    100*(obj.ModelObj.HTNormES_i - obj.FitResult.par(4*obj.ModelObj.nPixels+3)), 100*obj.FitResult.err(4*obj.ModelObj.nPixels+4));
            else
                myfitleg3=sprintf('');
            end
            
            if ~contains(obj.fixPar,sprintf('fix %.0f ;fix %.0f',3+2*obj.ModelObj.nPixels+4,3+2*obj.ModelObj.nPixels+5)) % in case TT FSD Probabilities are fitted
                myfitleg4=sprintf('\n TT GS-Prob=%.2f%% \t\\pm %.2f%% \n TT ES-Prob=%.2f%% \t\\pm %.2f%%',...
                    100*(obj.ModelObj.TTNormGS_i + obj.FitResult.par(6*obj.ModelObj.nPixels+3)),100*obj.FitResult.err(6*obj.ModelObj.nPixels+3),...
                    100*(obj.ModelObj.TTNormES_i - obj.FitResult.par(6*obj.ModelObj.nPixels+3)), 100*obj.FitResult.err(6*obj.ModelObj.nPixels+4));
            else
                myfitleg4=sprintf('');
            end
            
            myfitleg = [myfitleg1,myfitleg2,myfitleg3,myfitleg4];
            switch obj.DataType
                case 'Real'
                    mydatatype = 'data:';
                case 'Twin'
                     mydatatype = 'twin data:';
                case 'Fake'
                    mydatatype = 'fake MC data:';
                case 'FitriumTwin'
                    mydatatype = 'Fitrium twin data';
                case 'KafitTwin'
                    mydatatype = 'Kafit twin data';
            end
            
            if isa(obj,'MultiRunAnalysis')
                runleg =sprintf(' %.0f runs stacked,',numel(obj.StackedRuns));
            else
                runleg = sprintf(' run %.0f,',obj.RunNr);
            end
            pixleg1 = sprintf(' %.0f pixels',numel(obj.PixList));
            switch obj.AnaFlag
                case 'StackPixel'
                    pixelleg2 = sprintf(' stacked ');
                case 'SinglePixel'
                    pixelleg2 = sprintf(' pixel %.0f',obj.pixel);
                case 'MultiPixel'
                    pixelleg2 = sprintf(' %.0f pixel: %.0f - %.0f',numel(obj.PixList), obj.PixList(1),obj.PixList(end));
                case 'Ring'
                    pixelleg2 = sprintf(' %.0f rings',obj.nRings);
            end
            mydataleg = [mydatatype,runleg,pixleg1,pixelleg2];
            
            if isempty(obj.RunNr)   % MultiRunAnalysis
                mytitle =sprintf('Stacked %s Run %u-%u - Samak Fit (%s)',...
                    LabelFlag,obj.RunList(1), obj.RunList(end),obj.chi2);
                minititle = sprintf('Stacked %s Runs %u-%u (%.0f minutes)',...
                    LabelFlag,obj.RunList(1), obj.RunList(end),obj.ModelObj.TimeSec/60);
            else                       % RunAnalysis
                mytitle =sprintf('Run %u %s- Samak Fit (%s)',...% \n\\rho .d = %.2g mol/cm^2',
                    obj.RunNr,LabelFlag,obj.chi2);%,obj.ModelObj.WGTS_CD_MolPerCm2);
                minititle = sprintf('%s Run %u (%.0f minutes)',LabelFlag,obj.RunNr,obj.ModelObj.TimeSec/60);
            end
            ErangeLabel = round(obj.ModelObj.Q_i-obj.ModelObj.qU(obj.exclDataStart));
            savename = sprintf('./plots/FitRun%s_%s_%.0feV_%s_%.0fsigma',obj.ModelObj.TD,obj.chi2,ErangeLabel,obj.AnaFlag,obj.ErrorBarScaling);
            Description = struct('dataleg',mydataleg,'fitleg',myfitleg,...
                'title',mytitle,'minititle',minititle,...
                'savename',savename);
        end
        function PlotFitErrorMatrix(obj,varargin)
            % Display Fit Error Matrix (of Parameters)
            % -------------------------------------------------------------%

            p=inputParser;
            p.addParameter('saveplot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            saveplot = p.Results.saveplot;
            
            % Spectrum + Fit with Residuals
            fig51 = figure(51);
            set(fig51, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
            corplot(obj.FitResult.errmat);
            
            set(gca,'XTickLabel',{'m^2','E_0','B','N','DTgs','DTes','HTgs','HTes','TTgs','TTes','qUoff'});
            set(gca,'YTickLabel',{'m^2','E_0','B','N','DTgs','DTes','HTgs','HTes','TTgs','TTes','qUoff'});
            PrettyFigureFormat; set(gca,'FontSize',16);
            if strcmp(saveplot,'ON')
                save_name = sprintf('./plots/FitRun%u_chi2-%s-excl%u_ErrMat',obj.RunNr,obj.chi2,obj.exclDataStart);
                export_fig(fig51,[save_name,'.pdf']);
                print(fig51,[save_name,'.png'],'-dpng');
            end
        end     
        function PlotFitm2e0Contours(obj,varargin)
            % Construct data points of ellipses representing contour curves of 
            % (E0,m2) assuming Gaussian distributions 
            % Uses the Fit Results & Covariance matrix of the fitted parameters.
            % -------------------------------------------------------------%

            p=inputParser;
            p.addParameter('saveplot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('ErrorOnModel','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            saveplot      = p.Results.saveplot;
            ErrorOnModel  = p.Results.ErrorOnModel;
            
            % Spectrum + Fit with Residuals
            fig111 = figure(111);
            set(fig111, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
            
            if (obj.FitResult.par(1)==0 || obj.FitResult.par(2)==0)
                fprintf(2,'RunAnalysis:PlotFitm2e0Contours - Unaproriate Inputs - Aborting PlotFitm2e0Contours...\n');
                return;
            end
            
            cov_m2e0 = obj.FitResult.errmat(1:2,1:2);
            switch ErrorOnModel
                case 'OFF'
                    fit_m2e0 = [obj.FitResult.par(1) 0];
                case 'ON'
                    fit_m2e0 = [0,obj.FitResult.par(2)];
            end
           
            
           CLsigma(1)='1';
           r_ellipse = DrawEllipseError(obj,fit_m2e0,cov_m2e0,'CL',CLsigma(1));
           h1 = plot(r_ellipse(:,1) + fit_m2e0(1),r_ellipse(:,2) + fit_m2e0(2),'-','LineWidth',3);
           hold on
           hf1 = fill(r_ellipse(:,1) + fit_m2e0(1),r_ellipse(:,2) + fit_m2e0(2),rgb('DarkRed'),'FaceAlpha',0.2);
           CLsigma(2)='2';
           r_ellipse = DrawEllipseError(obj,fit_m2e0,cov_m2e0,'CL',CLsigma(2));
           h2 = plot(r_ellipse(:,1) + fit_m2e0(1),r_ellipse(:,2) + fit_m2e0(2),'-','LineWidth',3);
           hf2 = fill(r_ellipse(:,1) + fit_m2e0(1),r_ellipse(:,2) + fit_m2e0(2),rgb('IndianRed'),'FaceAlpha',0.2);
           CLsigma(3)='3';
           r_ellipse = DrawEllipseError(obj,fit_m2e0,cov_m2e0,'CL',CLsigma(3));
           h3 = plot(r_ellipse(:,1) + fit_m2e0(1),r_ellipse(:,2) + fit_m2e0(2),'-','LineWidth',3);
           hf3 = fill(r_ellipse(:,1) + fit_m2e0(1),r_ellipse(:,2) + fit_m2e0(2),rgb('Amethyst'),'FaceAlpha',0.2);
           h=[h1,h2,h3];
           
            switch ErrorOnModel
                case 'OFF'
                    d=plot(fit_m2e0(1),fit_m2e0(2),'s','Color','Black','MarkerSize',8,'LineWidth',3);
%                     d=errorbar(fit_m2e0(1),fit_m2e0(2),...
%                         sqrt(obj.FitResult.errmat(2,2))/1,sqrt(obj.FitResult.errmat(2,2))/1,...
%                         sqrt(obj.FitResult.errmat(1,1))/1,sqrt(obj.FitResult.errmat(1,1))/1,...
%                        '+','Color','Black','MarkerSize',8,'LineWidth',3);
                    strleg='Data';
                    strleg2='Fit (Gaussian Approximation) ';
                case 'ON'
%                     d=errorbar(obj.FitResult.par(1),obj.FitResult.par(2),...
%                         sqrt(obj.FitResult.errmat(2,2))/1,sqrt(obj.FitResult.errmat(2,2))/1,...
%                         sqrt(obj.FitResult.errmat(1,1))/1,sqrt(obj.FitResult.errmat(1,1))/1,...
%                         'o','Color','Black','MarkerSize',8,'LineWidth',3);
                    d=plot(obj.FitResult.par(1),0,'s','Color','Black','Marker','+','MarkerSize',24,'LineWidth',5,'MarkerFaceColor','Black');
                    strleg='Best Fit to Data';
                    strleg2='Model Uncertainty ';
            end
            hold off
            h1.Color = rgb('DarkRed');
            leg = legend([d,h1],strleg,[strleg2 CLsigma(1) ' \sigma']);
            if numel(h)==2
                h2.Color = rgb('IndianRed');
                leg = legend([d,h1,h2],strleg,[strleg2 CLsigma(1) '\sigma'],[strleg2 CLsigma(2) ' \sigma']);
            elseif numel(h)==3
                h3.Color = rgb('Salmon');
                leg = legend([d,h1,h2,h3],strleg, [strleg2 CLsigma(1) '\sigma'],[strleg2 CLsigma(2) ' \sigma'],[strleg2 CLsigma(3) ' \sigma']);
            elseif numel(h)==4
                h4.Color = rgb('Amethyst');
                leg = legend([d,h1,h2,h3,h4],strleg, [strleg2 CLsigma(1) '\sigma'],[strleg2 CLsigma(2) ' \sigma'],[strleg2 CLsigma(3), ' \sigma'],[strleg2 CLsigma(4) ' \sigma']);
            end
            
            xlabel('m^2 (eV)','FontSize',24);
            ylabel(sprintf('E_0 - %.2f (eV)',obj.ModelObj.Q_i+obj.FitResult.par(2)),'FontSize',24);
            leg.Location = 'northwest';
            legend('boxoff');
            
            grid on;
            PrettyFigureFormat
            set(gca,'FontSize',24);
            if strcmp(saveplot,'ON')
                save_name = sprintf('./plots/FitE0m2Contours_%s-excl%u_ErrMat',obj.chi2,obj.exclDataStart);
                export_fig(fig111,[save_name,'.pdf']);
                print(fig111,[save_name,'.png'],'-dpng');
            end
        end     
         function PlotFitCovCorMatrices(obj,varargin)
            % Plot Total Covariance Matrix / Correlation MAtrix
            
            p=inputParser;
            p.addParameter('Mode','Shape',@(x)ismember(x,{'Shape','Frac','CM'}));
            p.parse(varargin{:});
            Mode      = p.Results.Mode;            
                        
            % Selection Matrix
            switch Mode
                case 'Shape'
            StudyMat = obj.FitCMShape;
                case 'Frac'
            StudyMat = obj.FitCMFrac;
                case 'CM'
            StudyMat = obj.FitCM;
            end
            
            qUmin = round(obj.ModelObj.qU(obj.exclDataStart));
            qUmax = round(obj.ModelObj.qU(end));
            range = round(obj.ModelObj.qU(obj.exclDataStart)-obj.ModelObj.Q_i);


            myMainTitle=[sprintf('KATRIN - KNM1 Uniform - [%.0f - %0.f] eV ',qUmin,qUmax)];
            maintitle=myMainTitle;
            savefile=sprintf('plots/KNM1_FitCM_%s_%.0feVbelowE0-1.png',Mode,abs(range));
            fig1 = figure('Name','KATRIN - KNM1 Uniform  Covariance Matrix','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
            a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
            a.FontSize=18;a.FontWeight='bold';
            
            subplot(1,2,1)
            covplot(StudyMat(obj.exclDataStart:end,obj.exclDataStart:end));
            c = colorbar('southoutside');            c.FontSize = 12;
            axis square;
            set(gca,'xtick',[1 numel(obj.RunData.qU)-obj.exclDataStart+1]),set(gca,'ytick',[]);
            LqUmin = sprintf('qU_{min} = E_0-%.0fV',abs(obj.ModelObj.qU(obj.exclDataStart)-obj.ModelObj.Q_i));
            if (obj.ModelObj.qU(numel(obj.RunData.qU))-obj.ModelObj.Q_i)>=0
                LqUmax = sprintf('qU_{max} = E_0+%.0fV',obj.ModelObj.qU(numel(obj.RunData.qU))-obj.ModelObj.Q_i);
            else
                LqUmax = sprintf('qU_{max} = E_0-%.0fV',abs(obj.ModelObj.qU(numel(obj.RunData.qU))-obj.ModelObj.Q_i));
            end
            set(gca,'xticklabel',{LqUmin,LqUmax}); set(gca,'yticklabel',[]);
            title(sprintf('Fit Stat+Sys Covariance Matrix (%s)',Mode));
            PrettyFigureFormat;set(gca,'FontSize',12);
            
            PrettyFigureFormat;set(gca,'FontSize',12);
            subplot(1,2,2)
            corplot(StudyMat(obj.exclDataStart:end,obj.exclDataStart:end));
            c = colorbar('southoutside');            c.FontSize = 12;
            axis square;
            set(gca,'xtick',[1 numel(obj.RunData.qU)-obj.exclDataStart+1]),set(gca,'ytick',[]);
            LqUmin = sprintf('qU_{min} = E_0-%.0fV',abs(obj.ModelObj.qU(obj.exclDataStart)-obj.ModelObj.Q_i));
            if (obj.ModelObj.qU(numel(obj.RunData.qU))-obj.ModelObj.Q_i)>=0
                LqUmax = sprintf('qU_{max} = E_0+%.0fV',obj.ModelObj.qU(numel(obj.RunData.qU))-obj.ModelObj.Q_i);
            else
                LqUmax = sprintf('qU_{max} = E_0-%.0fV',abs(obj.ModelObj.qU(numel(obj.RunData.qU))-obj.ModelObj.Q_i));
            end
            set(gca,'xticklabel',{LqUmin,LqUmax}); set(gca,'yticklabel',[]);
            title('Fit Stat+Sys Correlation Matrix');
            PrettyFigureFormat;set(gca,'FontSize',12);
            export_fig(gcf,savefile,'-q101','-m3');
            
         end
         
         
         function [parRD, errRD, chi2RD, WGTS_CD_MolPerCm2_local, CD_bestfit] = RhoDScan(obj,varargin)
             % Perform Column Density Scan
             % cRhoD:
             % minRhoD: minimum column density in scan
             % maxRhoD: maximum column density in scan
             % tRhoD:
             % -------------------------------------------------------------%
             
             p=inputParser;
             p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
             p.addParameter('plotFit','OFF',@(x)ismember(x,{'ON','OFF'}));
             p.addParameter('WGTS_CD_MolPerCm2_V',(0.9:0.01:1.1),@(x)isfloat(x)); %vector of relative column densities
             p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
             p.parse(varargin{:});
             saveplot            = p.Results.saveplot;
             WGTS_CD_MolPerCm2_V = p.Results.WGTS_CD_MolPerCm2_V;
             RecomputeFlag       = p.Results.RecomputeFlag;
             plotFit             = p.Results.plotFit;
             
             if strcmp(saveplot,'ON')
                 plotFit = 'ON';
             end
             
             % save column density scan: labeling
             if contains(obj.DataSet,'Knm1')
                 savedir = [getenv('SamakPath'),'knm1ana/knm1_ColumnDensity/results/'];
             elseif contains(obj.DataSet,'FirstTritium')
                 savedir = [getenv('SamakPath'),'first-tritium-studies/ft_ColumnDensity/results/'];
             else
                 savedir = [getenv('SamakPath'),'ColumnDensity/results/'];
             end
             MakeDir(savedir); %create directory if doesnt exist
             savefile = [savedir,sprintf('RDScan_%s_%s_%s_%.0feV_RD%.1f-%.2f-%.1f.mat',obj.RunData.RunName,obj.DataType,obj.chi2,abs(18575-obj.ModelObj.qU(obj.exclDataStart)),...
                 min(WGTS_CD_MolPerCm2_V),WGTS_CD_MolPerCm2_V(2)-WGTS_CD_MolPerCm2_V(1),max(WGTS_CD_MolPerCm2_V))];
             
             if exist(savefile,'file') && strcmp(RecomputeFlag,'OFF')
                 load(savefile,'parRD', 'errRD', 'chi2RD', 'WGTS_CD_MolPerCm2_local','WGTS_CD_MolPerCm2_i',...
                     'CD_bestfit','WGTS_CD_MolPerCm2_ErrUp','WGTS_CD_MolPerCm2_ErrLow','dofRD');
                 obj.CDmin =  CD_bestfit;
                 obj.rhoDupperUnc= WGTS_CD_MolPerCm2_ErrUp;
                 obj.rhoDlowerUnc= WGTS_CD_MolPerCm2_ErrLow ;
             else
                 % save previous settings
                 RD_fixPar     = obj.fixPar;
                 RD_pullFlag   = obj.pullFlag;
                 RD_normFit_i  = obj.ModelObj.normFit_i;
                 
                 WGTS_CD_MolPerCm2_i     = obj.ModelObj.WGTS_CD_MolPerCm2;
                 WGTS_CD_MolPerCm2_local = WGTS_CD_MolPerCm2_i.*WGTS_CD_MolPerCm2_V;
                 parRD                   = zeros(obj.nPar,numel(WGTS_CD_MolPerCm2_local));
                 errRD                   = zeros(obj.nPar,numel(WGTS_CD_MolPerCm2_local));
                 chi2RD                  = zeros(numel(WGTS_CD_MolPerCm2_local),1);
                 
                 progressbar(sprintf('column density scan %.0feV',abs(18575-obj.ModelObj.qU(obj.exclDataStart))));
                 for rd = 1:numel(WGTS_CD_MolPerCm2_local)
                     progressbar(rd/numel(WGTS_CD_MolPerCm2_local));
                     obj.ModelObj.WGTS_CD_MolPerCm2 = WGTS_CD_MolPerCm2_local(rd);
                     obj.ModelObj.AdjustRF;
                     
                     obj.fixPar             = RD_fixPar;
                     obj.pullFlag           = RD_pullFlag;
                     obj.ModelObj.normFit_i = RD_normFit_i;
                     
                     obj.ModelObj.ComputeTBDDS; obj.ModelObj.ComputeTBDIS;
                     if ~strcmp(obj.chi2,'chi2Stat') %Initialize CovMat with statistics
                         SysEffects            = struct('FPDeff','ON','TASR','ON','FSD','ON','RF_RX','OFF','RF_EL','OFF','RF_BF','OFF','BkgShape','ON','TCoff_RAD','ON','TCoff_OTHER','ON','Stack','ON');
                         BkgCM                 = 'ON';
                         obj.ComputeCM('SysEffects',SysEffects,'BkgCM',BkgCM);
                     end
                     obj.Fit();
                     parRD(:,rd) =  obj.FitResult.par;
                     errRD(:,rd) =  obj.FitResult.err;
                     chi2RD(rd)  =  obj.FitResult.chi2min;
                     
                 end
                 dofRD       =  obj.FitResult.dof;
                 
                 %% Find minimum and asymmetric errors
                 [~,minIndex] = min(chi2RD);
                 
                 WGTS_CD_MolPerCm2_percent = WGTS_CD_MolPerCm2_local./WGTS_CD_MolPerCm2_i*100;
                 dataPointsforPlot = max(1,minIndex-4):min(length(chi2RD),minIndex+4);
                 [polyCoeff,structForErrors] = polyfit(WGTS_CD_MolPerCm2_percent(dataPointsforPlot)',chi2RD(dataPointsforPlot),3);
                 
                 r = roots([3*polyCoeff(1) 2*polyCoeff(2) polyCoeff(3)]);
                 [M,I] = min(abs(r-100));
                 CDminPercentage = r(I);
                 
                 obj.CDmin = CDminPercentage*WGTS_CD_MolPerCm2_i/100;
                 chi2RDmin = polyval(polyCoeff,CDminPercentage);
                 
                 chi2err = chi2RDmin + 1;
                 chi2RDLeftofMin = chi2RD(WGTS_CD_MolPerCm2_percent < CDminPercentage);
                 rhoDLeftofMin = WGTS_CD_MolPerCm2_percent(WGTS_CD_MolPerCm2_percent < CDminPercentage);
                 chi2RDRightofMin = chi2RD(WGTS_CD_MolPerCm2_percent > CDminPercentage);
                 rhoDRightofMin = WGTS_CD_MolPerCm2_percent(WGTS_CD_MolPerCm2_percent > CDminPercentage);
                 try
                     obj.rhoDlowerUnc = interp1(chi2RDLeftofMin,rhoDLeftofMin,chi2err,'spline')*WGTS_CD_MolPerCm2_i/100;
                     obj.rhoDupperUnc = interp1(chi2RDRightofMin,rhoDRightofMin,chi2err,'spline')*WGTS_CD_MolPerCm2_i/100;
                 catch
                     try
                         obj.rhoDupperUnc = interp1(chi2RDRightofMin,rhoDRightofMin,chi2err,'spline')*WGTS_CD_MolPerCm2_i/100;
                         obj.rhoDlowerUnc = 2*obj.CDmin - obj.rhoDupperUnc;
                     catch
                         obj.rhoDlowerUnc = interp1(chi2RDLeftofMin,rhoDLeftofMin,chi2err,'spline')*WGTS_CD_MolPerCm2_i/100;
                         obj.rhoDupperUnc = 2*obj.CDmin - obj.rhoDlowerUnc;
                     end
                 end
                 
                 CD_bestfit = obj.CDmin;
                 WGTS_CD_MolPerCm2_ErrUp  = obj.rhoDupperUnc;
                 WGTS_CD_MolPerCm2_ErrLow =  obj.rhoDlowerUnc;
                 save(savefile,'parRD', 'errRD', 'chi2RD', 'WGTS_CD_MolPerCm2_local','WGTS_CD_MolPerCm2_i',...
                     'CD_bestfit','WGTS_CD_MolPerCm2_ErrUp','WGTS_CD_MolPerCm2_ErrLow','dofRD');
                 
                 fprintf('Saving column density scan to file %s \n',savefile);
             end
             
             if strcmp(plotFit,'OFF')
                 return;
             end
             %% Plot Results
             obj.GetPlotColor;
             rhomin = min(WGTS_CD_MolPerCm2_local./WGTS_CD_MolPerCm2_i*100);
             rhomax = max(WGTS_CD_MolPerCm2_local./WGTS_CD_MolPerCm2_i*100);

             x = WGTS_CD_MolPerCm2_local./WGTS_CD_MolPerCm2_i*100;
             PlotArg =  {'o-', 'Color',obj.PlotColor,'LineWidth',2,'MarkerFaceColor',rgb('DarkCyan'),'MarkerSize',5};
             xstr = 'Column density (%)';
             myXticks = (rhomin:10:rhomax);
             
             fig12 = figure('Renderer','painters');
             set(fig12, 'Units', 'normalized', 'Position', [0.1, 0.1, 1 ,0.6]);
             
             s1 = subplot(2,4,1); % Neutrino mass
              y =parRD(1,:)+obj.ModelObj.mnuSq_i;
             ymin = min(y-errRD(1,:)); ymax = max(y+errRD(1,:));
             errorbar(x, y,errRD(1,:),PlotArg{:});
             hold on; 
             plot(100*ones(10,1),linspace(ymin,ymax,10),'Color',rgb('Silver'))
             hold off;
             xlabel(xstr)
             ylabel(sprintf('m^2(\\nu_e) (eV^2)'))
             xlim([rhomin rhomax]); xticks(myXticks);
             ylim([ymin,ymax]);
             PrettyFigureFormat;
             
             s2=subplot(2,4,2);  %Endpoint
             y =parRD(2,:);
             ymin = min(y-errRD(2,:)); ymax = max(y+errRD(2,:));
             errorbar(x, y,errRD(2,:),PlotArg{:});
             hold on;
             plot(100*ones(10,1),linspace(ymin,ymax,10),'Color',rgb('DimGray'))
             hold off;
             xlabel(xstr);
             ylabel(sprintf('E_{0_{fit}} - %0.1f (eV)',obj.ModelObj.Q_i));
             xlim([rhomin rhomax]); xticks(myXticks);
             ylim([ymin,ymax]);
             PrettyFigureFormat;

            
            s3=subplot(2,4,5); % Background 
            y = (obj.ModelObj.BKG_RateSec_i+parRD(3,:))*1e3;
            ymin = min(y-errRD(3,:)*1e3); ymax = max(y+errRD(3,:)*1e3);
            errorbar(x,y,errRD(3,:)*1e3,PlotArg{:})
            hold on;
            plot(100*ones(10,1),linspace(ymin,ymax,10),'Color',rgb('DimGray'))
            hold off;
            xlabel(xstr);
            ylabel(sprintf('B (mcps)'));
            xlim([rhomin rhomax]); xticks(myXticks);
            ylim([ymin,ymax]);
            PrettyFigureFormat;
            
            s4 =subplot(2,4,6); % Normalization
            y = obj.ModelObj.normFit_i+parRD(4,:)+1;
            errorbar(x,y,errRD(4,:),PlotArg{:});
            hold on;
            ymin = min(y-errRD(4,:)); ymax = max(y+errRD(4,:));
            plot(100*ones(10,1),linspace(ymin,ymax,10),'Color',rgb('DimGray'))
            hold off;
            xlim([rhomin rhomax]); xticks(myXticks);
            ylim([ymin,ymax]);
            xlabel(xstr);
            ylabel('N');
            PrettyFigureFormat;
            
            s5 =subplot(2,4,[3,4,7,8]); % chi2
            [~,minIndex] = min(chi2RD);
            
            WGTS_CD_MolPerCm2_percent = WGTS_CD_MolPerCm2_local./WGTS_CD_MolPerCm2_i*100;
            dataPointsforPlot = max(1,minIndex-4):min(length(chi2RD),minIndex+4);
            [polyCoeff,structForErrors] = polyfit(WGTS_CD_MolPerCm2_percent(dataPointsforPlot)',chi2RD(dataPointsforPlot),3);
            
            r = roots([3*polyCoeff(1) 2*polyCoeff(2) polyCoeff(3)]);
            [M,I] = min(abs(r-100));
            CDminPercentage = r(I);
            
            plotchi2 = plot(x, chi2RD,PlotArg{:},'MarkerSize',8);
            hold on;
            plotchi2fit = plot(WGTS_CD_MolPerCm2_local(dataPointsforPlot)./WGTS_CD_MolPerCm2_i*100,...
                polyval(polyCoeff,WGTS_CD_MolPerCm2_local(dataPointsforPlot)./WGTS_CD_MolPerCm2_i*100),...
                'o-','Color',rgb('GoldenRod'),'LineWidth',1.5,'LineStyle','-',...
                'MarkerFaceColor',rgb('DarkGoldenRod'),'MarkerSize',10);
            hold off
            xlabel(xstr);
            ylabel(['\chi^2 (',num2str(dofRD),' dof)']);
            xlim([rhomin rhomax]);  xticks(rhomin:5:rhomax)
            grid on;
            PrettyFigureFormat;
            set(gca,'FontSize',18);

           leg = legend([plotchi2,plotchi2fit],sprintf('Column density (100%%)  =  %.3g mol/cm^2',WGTS_CD_MolPerCm2_i),...
              sprintf('Column density best fit  = %0.3g ^{+ %0.3g}_{- %0.3g} mol/cm^2',obj.CDmin,WGTS_CD_MolPerCm2_ErrUp-CD_bestfit,...
              -1*(WGTS_CD_MolPerCm2_ErrLow-CD_bestfit)));
           leg.Location = 'north';
           legend('boxoff');
            
            linkaxes([s1,s2,s3,s4,s5], 'x');
            
            % reset
            obj.ModelObj.WGTS_CD_MolPerCm2 = WGTS_CD_MolPerCm2_i;
            obj.ModelObj.AdjustRF;
  
            if strcmp(saveplot,'ON')
                plotdir = strrep(savedir,'results','plots');
                MakeDir(plotdir);
                savename = strrep(strrep(savefile,'results','plots'),'.mat','.pdf');
                export_fig(fig12,savename);
                %print(fig12,[save_name,'.png'],'-dpng','-r400');
            end

        end
        function RhoDScanAllRings(obj)
            % Perform RhoD scan fit for All rings
            % -------------------------------------------------------------%
            
            for ri = obj.nRings
                obj.RhoDScan();
            end
            
        end
        
        
        function [parqU, errqU, chi2qU, dofqU,e1] = qUScan(obj,varargin)
            % Perform qUmin Fit Scan
            % -------------------------------------------------------------%
            p=inputParser;
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('qURange',[90,20],@(x)all(isfloat(x)));
            p.addParameter('CorrMean','OFF',@(x)ismember(x,{'ON','OFF'})); % doesnt work atm
            p.addParameter('HoldOn','OFF',@(x)ismember(x,{'ON','OFF','ON1'})); % plotting multiple FSDs...
            p.addParameter('RelFlag','OFF',@(x)ismember(x,{'ON','OFF'})); % show results to mean 
            p.addParameter('ErrorBarScaling',1,@(x)isfloat(x)); % scale error bars in fit results
            p.addParameter('RefLine','OFF',@(x) isfloat(x) || strcmp(x,'OFF')); % reference line for respect to certain range
            p.addParameter('saveStr','',@(x)ischar(x) || isempty(x));
            p.parse(varargin{:});
            saveplot        = p.Results.saveplot;
            qURange         = p.Results.qURange;
            RecomputeFlag   = p.Results.RecomputeFlag;
            CorrMean        = p.Results.CorrMean;
            HoldOn          = p.Results.HoldOn;
            RelFlag         = p.Results.RelFlag;
            ErrorBarScaling = p.Results.ErrorBarScaling;
            saveStr         = p.Results.saveStr; % additional label for result and plots
            RefLine         = p.Results.RefLine;
            
            if strcmp(obj.DataSet,'Knm1')
                savedir = [getenv('SamakPath'),'knm1ana/knm1_qUScan/results/'];
            elseif contains(obj.DataSet,'FirstTritium')
                savedir = [getenv('SamakPath'),'first-tritium-studies/ft_qUScan/results/'];
                MakeDir(savedir);
            elseif strcmp(obj.DataSet,'Knm2')
                savedir = [getenv('SamakPath'),'knm2ana/knm2_qUScan/results/'];
            else
                savedir = './plots/';
            end
            MakeDir(savedir);
            exclDataStart_v = obj.GetexclDataStart(qURange(1)):obj.GetexclDataStart(qURange(2));
            nFits = numel(exclDataStart_v);
            
            freeParStr = ConvertFixPar('Mode','Reverse','freePar',obj.fixPar);
            savename = [savedir,sprintf('qUScan_%s_%s_%sNP%.3f_FitPar%s_%.0feV-%.0feV_%s%s.mat',...
                obj.RunData.RunName,obj.DataType,obj.chi2,obj.NonPoissonScaleFactor,freeParStr,...
                qURange(1),qURange(2),obj.FSDFlag,saveStr)];
            if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
                load(savename,'parqU', 'errqU', 'chi2qU', 'dofqU');
                fprintf('load qU scan from file %s \n',savename)
            else
                parqU                   = zeros(obj.nPar,nFits);
                errqU                   = zeros(obj.nPar,nFits);
                chi2qU                  = zeros(nFits,1);
                dofqU                   = zeros(nFits,1);
                
                progressbar('qU Scan')
                for i=1:nFits
                    progressbar(i/nFits);
                    
                    obj.exclDataStart = exclDataStart_v(i);
                    savename_tmp = strrep(savename,sprintf('%.0feV-%.0feV',qURange(1),qURange(2)),...
                        sprintf('%.0feV',obj.GetRange));
                    
                    if exist(savename_tmp,'file') && strcmp(RecomputeFlag,'OFF')
                        FitResult_tmp = importdata(savename_tmp);
                        obj.FitResult = FitResult_tmp;
                        fprintf('load fit from file %s \n',savename_tmp)
                    else
                        obj.Fit;
                        FitResult = obj.FitResult;
                        save(savename_tmp,'FitResult');
                    end
                    parqU(:,i) = obj.FitResult.par;
                    errqU(:,i) = obj.FitResult.err;
                    chi2qU(i)  = obj.FitResult.chi2min;
                    dofqU(i)   = obj.FitResult.dof;
                end
                save(savename,'parqU', 'errqU', 'chi2qU', 'dofqU');
            end
            
            
            %% Plot qU Scan
            if ismember(HoldOn,{'ON','ON1'})
                fitPar = 1; % hack for nu mass plot with different models
            else
                fitPar = 1:5;
            end
            
            for i=fitPar
              
                yErr = flip(errqU(i,:));
                if i==1
                    y = flip(parqU(i,:));
                    if strcmp(RelFlag,'ON')
                        y = y-wmean(y,1./yErr.^2);
                        ystr = sprintf('{\\itm}_\\nu^2 - \\langle{\\itm}_\\nu^2\\rangle (eV^{ 2})');
                    else
                        ystr = sprintf('{\\itm}_\\nu^2 (eV^{ 2})');
                    end
                elseif i==2
                    y = flip(parqU(i,:));%obj.ModelObj.Q_i-18574;
                    if strcmp(RelFlag,'ON')
                        y = y-wmean(y,1./yErr.^2);
                        ystr = sprintf('{\\itE}_0^{fit} - \\langleE_0^{fit}\\rangle  (eV)');
                    else
                        ystr = sprintf('{\\itE}_0^{fit} - %.1f (eV)',obj.ModelObj.Q_i);
                    end
                elseif i==3
                    ystr = 'Background (mcps)';
                    y =(flip(parqU(i,:))+obj.ModelObj.BKG_RateSec_i)*1e3;
                    yErr = flip(errqU(i,:))*1e3;
                elseif i==4
                    ystr ='Normalization';
                    y =flip(parqU(i,:))+1;
                elseif i==5
                    ystr ='p-value';
                    y =flip(1-chi2cdf(chi2qU,dofqU));
                    yErr = zeros(numel(y),1);
                end
                
                if strcmp(CorrMean,'ON') && strcmp(obj.DataSet,'Knm1')
                    % get correlation matrix with largest number of samples
                    CovMatdir = [getenv('SamakPath'),'knm1ana/knm1_qUScan/results/'];
                    mydir = arrayfun(@(x) x.name,dir(CovMatdir),'UniformOutput',0);
                    Index = cell2mat(cellfun(@(x) contains(x,'CorrMat_qUScan_KNM1'),mydir,'UniformOutput',0));
                    myfiles = mydir(Index);
                    nSamples = str2double(extractBetween(myfiles,[obj.chi2,'_'],'.mat'));
                    thisfile = myfiles{nSamples==max(nSamples)};
                    d = importdata(thisfile);
                    %correlation matrix:
                    % Corr(1,1) describes point farthest away from endpoint
                    % Corr(end,end) describes point closest to endpoint
                    if i==1
                        CorrMat = d.CorrMatmNuSq;
                    elseif i==2
                        CorrMat = d.CorrMatE0;
                    elseif i==3
                        CorrMat = d.CorrMatB;
                    elseif i==4
                        CorrMat = d.CorrMatN;
                    end
                    CovMat = CorrMat.*flip(yErr).*flip(yErr)';
                    Weight = sum(inv(CovMat))./sum(sum(inv(CovMat)));
                    Mean   = sum(flip(y).*Weight);
                    chi2   = sum((flip(y)-Mean)*(CovMat\(flip(y)-Mean)'));
                    pval =  1-chi2cdf(chi2,numel(y)-1);
                else
                    Mean = wmean(y,1./yErr.^2);
                end
                
                yErr = yErr.*ErrorBarScaling;
                if all(y==0) || all(isnan(y))
                    
                else
                    
                 x =flip(obj.RunData.qU(exclDataStart_v(1):exclDataStart_v(end),1))-18575;%obj.ModelObj.Q_i;
                    if (strcmp(HoldOn,'OFF') || strcmp(HoldOn,'ON1')) && ~contains(obj.DataSet,'FirstTritium')
                        fig12345 = figure('Renderer','painters');
                        set(fig12345, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.6]);
                        ColorArg = {'MarkerFaceColor',obj.PlotColor,'Color',obj.PlotColor};
                    elseif contains(obj.DataSet,'FirstTritium')
                        fig12345 = figure('Renderer','painters');
                        set(fig12345, 'Units', 'centimeters', 'Position', [0.1, 0.1, 8.4, 4.5]);
                        ColorArg = {'MarkerFaceColor',rgb('CadetBlue'),'Color',rgb('DarkCyan')};
                    else
                        hold on;
                         ColorArg = {'MarkerFaceColor',obj.PlotColorLight,'Color',obj.PlotColorLight};
                       % ColorArg = {'MarkerFaceColor',rgb('Orange'),'Color',rgb('Orange')};
                        x = x+0.5;
                    end
                    
                    
                    if i~=5 && strcmp(CorrMean,'ON') && ~strcmp(HoldOn,'ON')
                        pref = plot(linspace(min(x)-3,max(x)+3,numel(x)),...
                            Mean.*ones(numel(x),1),':','Color',obj.PlotColor,'LineWidth',3); %rgb('Silver')
                        hold on
                    elseif ~strcmp(RefLine,'OFF')
                        % reference line with respect to certain range
                        RefIndex = find(abs(x+RefLine)==min(abs(x+RefLine)));  
                        pref = plot(linspace(min(x)-3,max(x)+3,numel(x)),...
                            y(RefIndex).*ones(numel(x),1),':','Color',rgb('Silver'),'LineWidth',3);
                        hold on
                    end
                    
                    if contains(obj.DataSet,'FirstTritium') && i==2
                        plot(linspace(min(x)-10,max(x)+10,numel(x)),zeros(numel(x),1),'-','Color',rgb('Black'),'LineWidth',1);
                        hold on;
                    end
                    e1 = errorbar(x, y,yErr,...
                        '.','LineWidth',2.5,'MarkerSize',25,ColorArg{:});
                    e1.CapSize = 0;
                    if ~strcmp(RefLine,'OFF')
                          plot(x(RefIndex),y(RefIndex),...
                        '.','LineWidth',2.5,'MarkerSize',e1.MarkerSize,...
                        'MarkerEdgeColor',obj.PlotColor,'MarkerFaceColor',obj.PlotColorLight);
                    end
                    xlabel(sprintf('Lower fit boundary below {\\itE}_0 (eV)'));
                    ylabel(ystr);
                    xlim([min(x)-2,max(x)+2]);
                    if i==3 %background
                        ylim([min(y-yErr).*0.999,max(y+yErr)*1.001])
                    elseif i==1
                        %ylim([min(y-yErr),max(y+yErr)])
                    end
                    if contains(obj.DataSet,'FirstTritium')
                        FTpaperFormat;
                        set(gca,'FontSize',9);
                        set(get(gca,'XLabel'),'FontSize',9);
                        set(get(gca,'YLabel'),'FontSize',9);
                        e1.LineWidth = 1;
                        e1.MarkerSize = 3;
                        xlim([min(x)-10,max(x)+10]);
                    else
                        PrettyFigureFormat('FontSize',24);
                    end
                    
                    
                    % legend
                    if ~strcmp(obj.chi2,'chi2Stat')
                        legstr='Fit (Stat + syst)';
                    else
                        legstr ='Fit (Stat only)';
                    end
                    
                    if i~=5 && strcmp(CorrMean,'ON') && ~strcmp(HoldOn,'ON')
                        leg = legend([e1,pref],legstr,'correlated weighted mean');
                    elseif ~strcmp(RefLine,'OFF')
                         leg = legend([e1,pref],legstr,sprintf('%.0f eV range',x(RefIndex)));
                    else
                        leg = legend(e1,legstr);
                    end
                    leg.EdgeColor = rgb('Silver');
                    
                    if i==1 && strcmp(CorrMean,'ON')
                        leg.Location = 'southwest';
                    else
                        leg.Location = 'northwest';
                    end
                    
                    if contains(obj.DataSet,'FirstTritium')
                        mypos = leg.Position;
                        leg.Position = [mypos(1)-0.05, mypos(2)+0.02,mypos(3:end)];
                        leg.FontSize = 9;
                    else
                        leg.FontSize = get(gca,'FontSize')+2;
                    end
                    
                    leg.delete;
                    %% save
                    if strcmp(saveplot,'ON')
                        if ErrorBarScaling~=1
                            errStr = sprintf('_%.0fErrScaling',ErrorBarScaling);
                        else
                            errStr = '';
                        end
                        plotdir = strrep(savedir,'results','plots');
                        savename_plot = strrep(strrep(strrep(savename,'results','plots'),'.mat','.png'),...
                            obj.RunData.RunName,[obj.RunData.RunName,'_',num2str(i),errStr]);
                        MakeDir(plotdir);
                        %print(gcf,savename_plot,'-dpng','-r100');
                        %publish_figurePDF(gcf,strrep(savename_plot,'.png','.pdf'));
                        export_fig(gcf,strrep(savename_plot,'.png','.pdf'));
                    end
                    if strcmp(HoldOn,'OFF')
                        close
                    end
                end
            end
            
            
        end
        function [parqU, errqU, chi2qU, dofqU] = qUScanFSD_TT(obj,varargin)
            % Perform qUmin Fit Scan
            % Plot the FSD Probabilities for TT ground state
            % Only valid for fixpar = '5 6 7 8 11'
            % -------------------------------------------------------------%
            
            p=inputParser;
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('firstPoint',8,@(x)isfloat(x) & x>0);
            p.addParameter('lastPoint',25,@(x)isfloat(x) & x>0);
            p.parse(varargin{:});
            saveplot    = p.Results.saveplot;
            firstPoint  = p.Results.firstPoint;
            lastPoint   = p.Results.lastPoint;
            
            parqU                   = zeros(obj.nPar,lastPoint-firstPoint+1);
            errqU                   = zeros(obj.nPar,lastPoint-firstPoint+1);
            chi2qU                  = zeros(lastPoint-firstPoint+1,1);
            dofqU                   = zeros(lastPoint-firstPoint+1,1);
            
            if ~strcmp(obj.fixPar,'fix 5 ;fix 6 ;fix 7 ;fix 8 ;fix 11 ;')
                fprintf(2,'RunAnalysis:qUScanFSD_TT: Wrong Fit Options\n');
                return;
            end
            
            counter=0;
            for i=firstPoint:1:lastPoint
                counter=counter+1;
                obj.exclDataStart = i;
                obj.Fit;
                parqU(:,counter) = obj.FitResult.par;
                errqU(:,counter) = obj.FitResult.err;
                chi2qU(counter)  = obj.FitResult.chi2min;
                dofqU(counter)   = obj.FitResult.dof;
            end
            
            %% Plot qU Scan
            qUmin = min(obj.RunData.qU(firstPoint:lastPoint,1)-obj.ModelObj.Q_i);
            qUmax = max(obj.RunData.qU(firstPoint:lastPoint,1)-obj.ModelObj.Q_i);
            if ~isempty(obj.RunNr) %RunAnalysis
                Runtitle      = sprintf('Run%u',obj.RunNr);
                Runtitle_save = sprintf('Run%u',obj.RunNr);
            else
                Runtitle      = sprintf('%uRuns %.0f - %.0f',numel(obj.StackedRuns),obj.RunList(1),obj.RunList(end));
                Runtitle_save = sprintf('%uRuns_%.0f_%.0f',numel(obj.StackedRuns),obj.RunList(1),obj.RunList(end));
            end
            
            % Plot E0,m2,N,B
            fig12345 = figure(12345); %Endpoint
            set(fig12345, 'Units', 'normalized', 'Position', [0.1, 0.1, 1.7 ,0.88]);
            subplot(2,4,5);
            errorbar(flip(obj.RunData.qU(firstPoint:lastPoint,1))-obj.ModelObj.Q_i, flip(parqU(2,:))+obj.ModelObj.Q_i,flip(errqU(2,:)),...
                's-', 'Color',rgb('Amethyst'),'LineWidth',1.5);
            xlabel(sprintf('qU_{min} - %.1f (eV)',obj.ModelObj.Q_i));
            ylabel('E_{0_{fit}} (eV)');
            xlim([qUmin qUmax]);
            grid on;
            PrettyFigureFormat;
            set(gca,'FontSize',14);
            
            subplot(2,4,[1 4]); % PgsTT
            errorbar(flip(obj.RunData.qU(firstPoint:lastPoint,1))-obj.ModelObj.Q_i, flip(obj.ModelObj.TTNormGS_i+parqU(9,:)),flip(errqU(9,:)),...
                's-', 'Color',rgb('DarkSlateGray'),'LineWidth',2);
            xlabel(sprintf('qU_{min} - %.1f (eV)',obj.ModelObj.Q_i));
            ylabel('P_{gs}(T_2)')
            xlim([qUmin qUmax]);
            grid on;
            PrettyFigureFormat;
            title(sprintf('%s - Fit Range Scan (%s) Results (Samak)',Runtitle,obj.chi2));
            set(gca,'FontSize',14);
            
            subplot(2,4,7); % Background
            errorbar(flip(obj.RunData.qU(firstPoint:lastPoint,1))-obj.ModelObj.Q_i, flip(obj.ModelObj.BKG_RateSec_i+parqU(3,:))*1e3,flip(errqU(3,:))*1e3,...
                's-', 'Color',rgb('Amethyst'),'LineWidth',2);
            xlabel(sprintf('qU_{min} - %.1f (eV)',obj.ModelObj.Q_i));
            ylabel(sprintf('BKG (m cps)'));
            xlim([qUmin qUmax]);
            grid on;
            PrettyFigureFormat;
            set(gca,'FontSize',14);
            
            subplot(2,4,6); % Neutrino Mass Squared
            errorbar(flip(obj.RunData.qU(firstPoint:lastPoint,1))-obj.ModelObj.Q_i, flip(obj.ModelObj.mnuSq_i+parqU(1,:)),flip(errqU(1,:)),...
                's-', 'Color',rgb('Amethyst'),'LineWidth',1.5);
            xlabel(sprintf('qU_{min} - %.1f (eV)',obj.ModelObj.Q_i));
            ylabel('m^2 (eV)');
            xlim([qUmin qUmax]);
            grid on;
            PrettyFigureFormat;
            set(gca,'FontSize',14);
            
            subplot(2,4,8); % p-value
            plot(flip(obj.RunData.qU(firstPoint:lastPoint,1))-obj.ModelObj.Q_i, flip(1-chi2cdf(chi2qU,dofqU))',...
                's-', 'Color',rgb('CadetBlue'),'LineWidth',2);
            xlabel(sprintf('qU_{min} - %.1f (eV)',obj.ModelObj.Q_i));
            ylabel('p-value');
            xlim([qUmin qUmax]);
            grid on;
            PrettyFigureFormat;
            set(gca,'YAxisLocation','right');
            set(gca,'FontSize',14);
            
            if strcmp(saveplot,'ON')
                if  exist('plots/qU_Scan/')==0 %folder doens't exist
                    mkdir plots/qU_Scan
                end
                save_name = sprintf('./plots/qU_Scan/qUScanFSD_TT%s_%s_%s',Runtitle_save,obj.DataType,obj.chi2);
                export_fig(fig12345,[save_name,'.pdf']);
            end
            
        end
        function r_ellipse = DrawEllipseError(obj,fit_m2e0,covariance,varargin)
            
            p=inputParser;
            p.addParameter('CL','1',@(x)ismember(x,{'1','2','3'}));
            p.addParameter('myColor','Black');
            p.parse(varargin{:});
            CL      = p.Results.CL;
            myColor = p.Results.myColor;
            
            % Calculate the eigenvectors and eigenvalues
            [eigenvec, eigenval] = eig(covariance);
            
            % Get the index of the largest eigenvector
            [largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
            largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);
            
            % Get the largest eigenvalue
            largest_eigenval = max(max(eigenval));
            
            % Get the smallest eigenvector and eigenvalue
            if(largest_eigenvec_ind_c == 1)
                smallest_eigenval = max(eigenval(:,2));
                smallest_eigenvec = eigenvec(:,2);
            else
                smallest_eigenval = max(eigenval(:,1));
                smallest_eigenvec = eigenvec(1,:);
            end
            
            % Calculate the angle between the x-axis and the largest eigenvector
            angle = atan2(largest_eigenvec(2), largest_eigenvec(1));
            
            % This angle is between -pi and pi.
            % Let's shift it such that the angle is between 0 and 2pi
            if(angle < 0)
                angle = angle + 2*pi;
            end
            
            % Get the coordinates of the data mean
            avg = fit_m2e0;
            
            % Get the 95% confidence interval error ellipse
            switch CL
                case '99%'
                    %  chisquare_val = sqrt(9.210);
                case '95%'
                    %  chisquare_val = sqrt(5.991);
                case '90%'
                    % chisquare_val = sqrt(4.605);
                case '1'
                    % chisquare_val = sqrt(2.3);
                    chisquare_val = sqrt(1);
                case '2'
                    % chisquare_val = sqrt(6.18);
                    chisquare_val = sqrt(4);
                case '3'
                    % chisquare_val = sqrt(11.83);
                    chisquare_val = sqrt(9);
            end
            
            theta_grid = linspace(0,2*pi);
            phi = angle;
            X0=avg(1);
            Y0=avg(2);
            a=chisquare_val*sqrt(largest_eigenval);
            b=chisquare_val*sqrt(smallest_eigenval);
            
            % the ellipse in x and y coordinates
            ellipse_x_r  = a*cos( theta_grid );
            ellipse_y_r  = b*sin( theta_grid );
            
            %Define a rotation matrix
            R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
            
            %let's rotate the ellipse to some angle phi
            r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
            
            % Draw the error ellipse
            %h = plot(hax2,r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-','LineWidth',2,'Color',myColor);
            
            % Plot the eigenvectors
            %quiver(X0, Y0, largest_eigenvec(1)*sqrt(largest_eigenval), largest_eigenvec(2)*sqrt(largest_eigenval), '-b', 'LineWidth',1);
            %quiver(X0, Y0, smallest_eigenvec(1)*sqrt(smallest_eigenval), smallest_eigenvec(2)*sqrt(smallest_eigenval), '-b', 'LineWidth',1);
            %hold off;
            
        end
    end
    methods % auxillary methods
        function filename = SetTwinOrFakeFileName(obj)
            if ismember(obj.DataType,{'Twin','FitriumTwin','KafitTwin'})
                
                if strcmp(obj.AnaFlag,'Ring')
                    ringfiles = '';%sprintf('_%.0fRings_%.0f-%.0f',obj.nRings,obj.RingList(1),obj.RingList(end));
                else
                    ringfiles ='';
                end
                
                str_bias = '';
                if isa(obj,'MultiRunAnalysis') &&  strcmp(obj.Twin_SameCDFlag,'ON')
                    str_CDbias = sprintf('_sameCD');
                    str_bias = [str_bias,str_CDbias];
                elseif obj.TwinBias_WGTS_CD_MolPerCm2~=1
                    str_CDbias = sprintf('_%.1f-WGTS_CD_MolPerCm2Bias',obj.TwinBias_WGTS_CD_MolPerCm2);
                    str_bias = [str_bias,str_CDbias];
                end
                
                if isa(obj,'MultiRunAnalysis') &&  strcmp(obj.Twin_SameIsotopFlag,'ON')
                    str_Isotop = sprintf('_IsotopSame');
                    str_bias = [str_bias,str_Isotop];
                else
                    if obj.TwinBias_WGTS_MolFrac_TT~=1
                        str_TTbias = sprintf('_%.1f-WGTS_MolFrac_TTBias',obj.TwinBias_WGTS_MolFrac_TT);
                        str_bias = [str_bias,str_TTbias];
                    end
                    
                    if obj.TwinBias_WGTS_MolFrac_DT~=1
                        str_DTbias = sprintf('_%.1f-WGTS_MolFrac_DTBias',obj.TwinBias_WGTS_MolFrac_DT);
                        str_bias = [str_bias,str_DTbias];
                    end
                    
                    if obj.TwinBias_WGTS_MolFrac_HT~=1
                        str_HTbias = sprintf('_%.1f-WGTS_MolFrac_HTBias',obj.TwinBias_WGTS_MolFrac_HT);
                        str_bias = [str_bias,str_HTbias];
                    end
                end
                if isa(obj,'MultiRunAnalysis') &&  strcmp(obj.Twin_SameqUFlag,'ON')
                    str_qUbias = sprintf('_qUSame');
                    str_bias = [str_bias,str_qUbias];
                elseif any(obj.TwinBias_qU~=0)
                    str_qUbias = sprintf('_%.1feV-qUOffset',obj.TwinBias_qU);
                    str_bias = [str_bias,str_qUbias];
                end
                
                if isa(obj,'MultiRunAnalysis') &&  strcmp(obj.Twin_SameqUfracFlag,'ON')
                    str_qUfracbias = sprintf('_qUfracSame');
                    str_bias = [str_bias,str_qUfracbias];
                end
                
                if ~isempty(obj.TwinBias_Time) %&& ~all(obj.TwinBias_Time==obj.ModelObj.TimeSec)
                    if isempty(diff(obj.TwinBias_Time)) || diff(obj.TwinBias_Time)==0
                        str_Timebias = sprintf('_%.0fs-Time',obj.TwinBias_Time(1)); %when all the same
                    else
                        str_Timebias = sprintf('_%.0fs-Time',obj.TwinBias_Time);
                    end
                    str_bias = [str_bias,str_Timebias];
                    
                end
                
                if any(obj.TwinBias_Bkg~=1)
                    if isempty(diff(obj.TwinBias_Bkg)) || diff(obj.TwinBias_Bkg)==0
                        str_Bkgbias = sprintf('_%.2f-BkgBias',obj.TwinBias_Bkg(1));
                    else
                        str_Bkgbias = sprintf('_%.2f-BkgBias',obj.TwinBias_Bkg);
                    end
                    str_bias = [str_bias, str_Bkgbias];
                end
                
                if obj.TwinBias_mnuSq~=0
                    str_bias = [str_bias,sprintf('_mNuSq%.2feV2',obj.TwinBias_mnuSq)];
                end
                
                if ischar(obj.TwinBias_Q)
                    str_bias = [str_bias,sprintf('_E0%s',obj.TwinBias_Q)];
                elseif numel(obj.TwinBias_Q)==1
                    str_bias = [str_bias,sprintf('_E0%.2feV',obj.TwinBias_Q)];
                elseif numel(obj.TwinBias_Q)==obj.nRuns
                    % TwinBias_Q is an array of values: 1 endpoint for each run
                    str_bias = [str_bias,sprintf('_meanE0%.2feV',mean(obj.TwinBias_Q))];
                else
                    % if nothing of the above: do not specify name further
                end
                 if obj.TwinBias_FSDSigma~=0
                    str_bias = [str_bias,sprintf('_FSDsigma%.3geV',obj.TwinBias_FSDSigma)];
                 end
                
                if ~strcmp(obj.KTFFlag,'WGTSMACE')
                    str_bias = [str_bias,'_',obj.KTFFlag];
                end
                
                if obj.TwinBias_BKG_PtSlope~=0
                    str_bias = [str_bias,sprintf('_BkgPtSlope%.1fmuCpsPerS',obj.TwinBias_BKG_PtSlope*1e6)];
                end
                
                filename = [ringfiles,str_bias];
                switch obj.FitNBFlag
                    case 'OFF'
                        filename = [filename,'FitNBFlagOFF'];
                    case 'NormOnly'
                        filename = [filename,'FitBFlagOFF'];
                end
                
                if strcmp(obj.MosCorrFlag,'ON')
                    filename = [filename,'_MosCorr'];
                end
                obj.TwinFakeLabel =filename;
                
            elseif ismember(obj.DataType,{'Fake'})

                if numel(obj.TwinBias_Q)==1
                    filename = sprintf('_E0%.2feV',obj.TwinBias_Q);  
                else
                    filename = '';
                end
                
                if ~isempty(obj.TwinBias_Time) 
                    str_Timebias = sprintf('_%.0fs-Time',obj.TwinBias_Time); %when all the same
                    filename = [filename,str_Timebias]; 
                end
                
                
                obj.TwinFakeLabel = filename;
            else
                filename = '';
                obj.TwinFakeLabel = filename;
            end
        end
        function SetDefaultInit(obj)
            % set some default values if not specified according to DataSet
            if isempty(obj.DataSet)
                return
            end
            
            % energy loss function
            if isempty(obj.ELossFlag)
                switch obj.DataSet
                    case 'Knm1'
                        switch obj.DataType
                            case 'Real'
                                obj.ELossFlag = 'KatrinT2';
                            case 'Twin'
                                obj.ELossFlag = 'KatrinT2';
                                % obj.ELossFlag = 'KatrinD2';
                            case {'FitriumTwin','KafitTwin'}
                                obj.ELossFlag = 'KatrinD2';
                        end
                    case 'FirstTritium.katrin'
                        obj.ELossFlag = 'Abdurashitov';
                    case 'Knm2'
                        obj.ELossFlag = 'KatrinT2A20';
                end
            end
            
            % systematics budget
            if isempty(obj.SysBudget)
                switch obj.DataSet
                    case 'Knm1'
                        obj.SysBudget = 22;
                    case 'FirstTritium.katrin'
                        obj.SysBudget = 0;
                    case 'Knm2'
                        obj.SysBudget = 40;
                end
            end
            
            % Doppler effect
            if isempty(obj.DopplerEffectFlag)
                switch obj.DataSet
                    case 'Knm1'
                        obj.DopplerEffectFlag = 'OFF'; %always included by 'FSD_Doppler'
                    case 'Knm2'
                        obj.DopplerEffectFlag = 'FSD';
                    case 'FirstTritium.katrin'
                        obj.DopplerEffectFlag = 'OFF';
                end
            end
            
        end
        function [TTFSD,DTFSD,HTFSD] = SetDefaultFSD(obj)
            
            % KNM2 blinding: force blind FSDs for all KNM2 analysis
%             if strcmp(obj.DataSet,'Knm2') && ~strcmp(obj.DataType,'Fake')
%                % obj.FSDFlag =  'BlindingKNM2';
%             end
            
            switch obj.FSDFlag
                case 'Sibille'
                    DTFSD = 'Sibille';
                    HTFSD = 'Sibille';
                    TTFSD = 'Sibille';
                case 'Sibille0p5eV'
                    DTFSD = 'Sibille0p5eV';
                    HTFSD = 'Sibille0p5eV';
                    TTFSD = 'Sibille0p5eV';
                case 'SibilleFull'
                    DTFSD = 'SibilleFull';
                    HTFSD = 'SibilleFull';
                    TTFSD = 'SibilleFull';
                case 'BlindingKNM1'
                    DTFSD = 'BlindingKNM1';
                    HTFSD = 'BlindingKNM1';
                    TTFSD = 'BlindingKNM1';
                case 'SAENZ' % to be used for first tritium
                    DTFSD = 'HTFSD';
                    HTFSD = 'SAENZ';
                    TTFSD = 'SAENZ';
                case 'OFF'
                    DTFSD = 'OFF';
                    HTFSD = 'OFF';
                    TTFSD = 'OFF';
                case 'BlindingKNM2'
                    DTFSD = 'BlindingKNM2';
                    HTFSD = 'BlindingKNM2';
                    TTFSD = 'BlindingKNM2';
                case 'KNM2'
                    DTFSD = 'KNM2';
                    HTFSD = 'KNM2';
                    TTFSD = 'KNM2';
                case 'KNM2_0p5eV'
                    DTFSD = 'KNM2_0p5eV';
                    HTFSD = 'KNM2_0p5eV';
                    TTFSD = 'KNM2_0p5eV';
                case 'KNM2_0p1eV'
                    DTFSD = 'KNM2_0p1eV';
                    HTFSD = 'KNM2_0p1eV';
                    TTFSD = 'KNM2_0p1eV';
            end
        end
        function InitFitPar(obj)
            obj.nPar = 4*obj.ModelObj.nPixels+13; % number of avaibale fit parameter
            
            if strcmp(obj.AnaFlag,'StackPixel') % number of FPD segmentations
                nFPDSeg = 1;
            elseif strcmp(obj.AnaFlag,'Ring')
                nFPDSeg = obj.nRings;
            end
            
            if isempty(obj.fixPar)
                % fixPar not specified in parser
                % -> just fit endpoint, background and normalization
                obj.fixPar = ConvertFixPar('freePar','E0 Bkg Norm','nPar',obj.nPar,'nPixels',nFPDSeg);
            elseif isempty(str2num(obj.fixPar)) && ~contains(obj.fixPar,'fix')
                % fixPar is string of free Parameters (characters) e.g. 'E0, Bkg, Norm'
                % -> use new method to convert into fixPar
                obj.fixPar = ConvertFixPar('freePar',obj.fixPar,'nPar',obj.nPar,'nPixels',nFPDSeg);
            elseif contains(obj.fixPar,'fix')
                % do nothing
            else
                % back compatibility: fixPar is string of numbers e.g. '1 2 3 4'
                % -> do nothing but,conversion '1 2' -> 'fix 1; fix 2;'
                all(~isnan(str2double(strrep(obj.fixPar,' ',''))))
                fixParNum = str2num(obj.fixPar);
                obj.fixPar = sprintf('fix %.0f ;',fixParNum);
            end
            
            if obj.exclDataStop==9999
                obj.exclDataStop = obj.ModelObj.nqU;
            end
            
        end
        function SetNPfactor(obj)
            if isempty(obj.NonPoissonScaleFactor)
                if isa(obj,'MultiRunAnalysis') %this routine only works for MultiRunAnalysis (->because of FitRunList...)
                    [obj.NonPoissonScaleFactor, ~]= GetNPfactor(obj,'SlopeCorr','ON','RecomputeFlag','OFF','SanityPlot','OFF');
                else
                    fprintf('WARNING: NonPoissonScaleFactor is set to 1 \n');
                    obj.NonPoissonScaleFactor = ones([1,numel(obj.RunData.MACE_Ba_T)]);
                end
            elseif numel(obj.NonPoissonScaleFactor)~=numel(obj.RunData.MACE_Ba_T)
                % NP scale factor should have same size as number of pseudorings
                obj.NonPoissonScaleFactor = repmat(mean(obj.NonPoissonScaleFactor),[1,numel(obj.RunData.MACE_Ba_T)]);
            end
        end
        function SetROI(obj)
            % change region of interest (only for KNM2 and later)
            % 2 options: Default and [14,32] keV ROI
            if ~(strcmp(obj.DataSet,'Knm2') && strcmp(obj.DataType,'Real'))
                fprintf('ROI change only available for KNM2 real data \n');
                return;
            end
            
            switch obj.ROIFlag
                case 'Default'
                    obj.RunData.TBDIS      = obj.RunData.TBDIS_Default;
                    obj.RunData.TBDIS_RM   = obj.RunData.TBDIS_RM_Default;
                case '14keV'
                    obj.RunData.TBDIS      = obj.RunData.TBDIS14keV;
                    obj.RunData.TBDIS_RM   = obj.RunData.TBDIS14keV_RM;
            end
            
            if isa(obj,'MultiRunAnalysis')
                switch obj.ROIFlag
                    case 'Default'
                        obj.SingleRunData.TBDIS      = obj.SingleRunData.TBDIS_Default;
                        obj.SingleRunData.TBDIS_RM   = obj.SingleRunData.TBDIS_RM_Default;
                    case '14keV'
                        obj.SingleRunData.TBDIS      = obj.SingleRunData.TBDIS14keV;
                        obj.SingleRunData.TBDIS_RM   = obj.SingleRunData.TBDIS14keV_RM;
                end
            end
        end
        function SetMosCorr(obj)
            % correct qU for long term drift -> input from monitor spectrometer
            % cannot be reset -> if you want to change back to no correction
            % -> re-run data import before: ReadData or Data
            
            if strcmp(obj.MosCorrFlag,'OFF')
                % no correction
                return
            end
            
            if ~(strcmp(obj.DataSet,'Knm2') && strcmp(obj.DataType,'Real'))
                fprintf('MoS qU correction only available for KNM2 real data \n');
                return
            end
            
            try  [qUslope, StartTimeMean] = GetMosDriftKNM2('SanityPlot','OFF'); % eV/day
            catch
                fprintf('MoS qU correction file not found \n');
            end
            
            fprintf('Correct qU with monitor spectrometer correction \n');
            
            qUslope = -qUslope; % qU is defined positiv in code -> flip sign of slope
            
            if isa(obj,'MultiRunAnalysis')
                LiveTimeDays            = days(obj.SingleRunData.StartTimeStamp-StartTimeMean);% live time with respect to whole KNM2
                qUCorr                  = repmat(LiveTimeDays.*qUslope,[size(obj.SingleRunData.qU,1),1,148]);
                obj.SingleRunData.qU    = obj.SingleRunData.qU+qUCorr;
                obj.SingleRunData.qU_RM = obj.SingleRunData.qU_RM+LiveTimeDays.*qUslope;
            else
                LiveTimeDays      = days(obj.RunData.StartTimeStamp-StartTimeMean);% live time with respect to whole KNM2
                qUCorr            = LiveTimeDays.*qUslope;
                obj.RunData.qU    = obj.RunData.qU+qUCorr;
                obj.RunData.qU_RM = obj.RunData.qU_RM+LiveTimeDays.*qUslope;
            end
        end
        function range = GetRange(obj)
            % data range below endpoint (eV)
            range = abs(round(obj.RunData.qU(obj.exclDataStart)-obj.ModelObj.Q_i,0));
        end
        function exclDataStart = GetexclDataStart(obj,range)
            exclDataStart = find((obj.ModelObj.qU)>=obj.ModelObj.Q_i-range,1);
        end
    end
        methods % ring methods
            function SaveOrLoadRingFit(obj,mode)
            savedir = [getenv('SamakPath'),sprintf('tritium-data/fit/%s/MultiRing/',obj.DataSet)];
            savename = [savedir,sprintf('MultiRing%s_%s_%s_fixPar%s.mat',...
                obj.RingMerge,obj.RunData.RunName,obj.chi2,strrep(strrep(obj.fixPar,'fix ',''),' ;',''))];
            switch mode
                case 'save'
                    MakeDir(savedir);
                    FitResult = obj.FitResult;
                    save(savename,'FitResult');
                case 'load'
                    if exist(savename,'file')
                        d = importdata(savename);
                        obj.FitResult = d.FitResult;
                        obj.ModelObj.ComputeTBDDS(...
                            'mSq_bias',d.FitResult.par(1),...
                            'E0_bias',d.FitResult.par(2),...
                            'B_bias',d.FitResult.par(3:obj.nRings+2),...
                            'N_bias',d.FitResult.par(obj.nRings+3:(2*obj.nRings+2)),...
                            'DTGS_bias',d.FitResult.par(2*obj.nRings+3),...
                            'DTES_bias',d.FitResult.par(2*obj.nRings+4),...
                            'HTGS_bias',d.FitResult.par(2*obj.nRings+5),...
                            'HTES_bias',d.FitResult.par(2*obj.nRings+6),...
                            'TTGS_bias',d.FitResult.par(2*obj.nRings+7),...
                            'TTES_bias',d.FitResult.par(2*obj.nRings+8),...
                            'qUOffset_bias',d.FitResult.par((2*obj.nRings+9):(3*obj.nRings+8)),...
                            'BSlope_bias',d.FitResult.par(3*obj.nRings+9));
                        obj.ModelObj.ComputeTBDIS;
                    else
                        fprintf(2,'Multi-Ring fit result not avaiable in data bank \n');
                    end
            end
            end
        
        function PlotFitMultiRing(obj,varargin)
            p=inputParser;
            p.addParameter('PlotPar','qU',@(x)ismember(x,{'qU','Norm','Bkg','mTSq'})); %display parameter
            p.addParameter('linFitFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('savePlot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Blind','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            PlotPar = p.Results.PlotPar;
            linFitFlag  = p.Results.linFitFlag;
            savePlot    = p.Results.savePlot;
            Blind       = p.Results.Blind;
            
            switch PlotPar
                case 'qU'
                    yErr  = obj.FitResult.err(2*obj.nRings+9:3*obj.nRings+8);
                    y     = obj.FitResult.par(2*obj.nRings+9:3*obj.nRings+8);
                    ystr = sprintf(' \\DeltaqU (eV)');
                    legPos = 'northeast';
                case 'Bkg'
                    nPixRing = cell2mat(cellfun(@(x) numel(x),obj.RingPixList,'UniformOutput',0));
                    y    = 1e3.*(obj.FitResult.par(3:2+obj.nRings)+obj.ModelObj.BKG_RateSec_i)./nPixRing';
                    yErr = 1e3.*obj.FitResult.err(3:2+obj.nRings)./nPixRing';
                    ystr = 'Background / npixels (mcps)';
                    legPos = 'northwest';
                case 'Norm'
                    y    = obj.FitResult.par(3+obj.nRings:3+2*obj.nRings-1)+1;
                    yErr = obj.FitResult.err(3+obj.nRings:3+2*obj.nRings-1);
                    ystr = sprintf('N_{sig}');%'Signal normalization factor';
                    legPos = 'northeast';
                case 'mTSq'
                    if strcmp(Blind,'ON')
                        fprintf('Blinding is ON -> you cannot look at mTSq \n');
                        return
                    end
                    y    = obj.FitResult.par(3*obj.nRings+10:4*obj.nRings+9);
                    yErr = obj.FitResult.err(3*obj.nRings+10:4*obj.nRings+9);
                    ystr = sprintf('Energy smearing \\sigma^2 (eV^2)');
                    legPos = 'northeast';
            end
            
            PlotStyle = {'o','Color',rgb('DodgerBlue'),'MarkerFaceColor',...
                rgb('DodgerBlue'),'LineWidth',2,'MarkerSize',9,...
                'CapSize',0};
            
            fig1 = figure('Units','normalized','Position',[0.1,0.1,0.6,0.5]);
            if ismember(PlotPar,{'qU','mTSq'})
                plot(linspace(-5,obj.nRings+1,10),zeros(10,1),'LineWidth',2,'Color',rgb('Silver'));
                hold on;
            end
            e1 = errorbar(obj.RingList,y,yErr,PlotStyle{:});
            PrettyFigureFormat('FontSize',24);
            ylabel(ystr);
            xlabel('Ring');
            xticks(obj.RingList); set(gca,'XMinorTick','off');
            if strcmp(obj.RingMerge,'Full')
                xticklabels({'1,2,3','4,5,6','7,8,9','10,11,12'})
            elseif strcmp(obj.RingMerge,'Half')
                xticklabels({'1,2,3,4,5,6','7,8,9,10,11,12'})
            elseif strcmp(obj.RingMerge,'Azi')
                xticklabels({'Pole','NE','SE','SW','NW'})
            elseif strcmp(obj.RingMerge,'AziHalfEW')
                xticklabels({'East','West'})
            elseif strcmp(obj.RingMerge,'AziHalfNS')
                xticklabels({'North','South'})
            end
            % set nice limits
            xlim([min(obj.RingList)-0.2,max(obj.RingList)+0.2]);
            ymin = min(y-yErr);
            if ~ismember(PlotPar,{'Norm','Bkg'})
                if ymin<0 && ~strcmp(PlotPar,'qU')
                    ylim([ymin*1.2,max(y+yErr)*1.5]);
                elseif strcmp(PlotPar,'qU')
                    ylim([-0.02,max(y+yErr)*1.1]);
                else
                    ylim([ymin*0.8,max(y+yErr)*1.5]);
                end
            end
            % linear fit (optional)
            if strcmp(linFitFlag,'ON')
                % test start
                if strcmp(PlotPar,'qU') && yErr(1)==0
                    yErr(1) = obj.FitResult.err(2);
                end
                [linFitpar, linFiterr, linFitchi2min,linFitdof] =...
                    linFit(obj.RingList',y',yErr');
                hold on;
                plot(obj.RingList,linFitpar(1).*obj.RingList+linFitpar(2),'-','Color',e1.Color,'LineWidth',e1.LineWidth);
                t = title(sprintf('Linear fit slope: %.1f \\pm %.1f meV/ring @ \\chi2 = %.1g (%.0f dof)',...
                    linFitpar(1)*1e3,linFiterr(1)*1e3,linFitchi2min,linFitdof));
                t.FontWeight = 'normal';
            end
            
            %legend
                d=obj.GetPlotDescription('data');
                  if strcmp(Blind,'OFF')
                      resultsleg = d.fitleg;
                  else
                      resultsleg = '';  
                  end
                leg = legend(e1,[sprintf(' \\chi^2 = %.1f (%.0f dof) \n',...
                    obj.FitResult.chi2min,obj.FitResult.dof),resultsleg]);
                leg.Title.String = sprintf('%s Multi-ring fit (%s)',upper(obj.DataSet),d.chi2);
                leg.Location = legPos;
                leg.FontSize = get(gca,'FontSize');
                leg.Title.FontWeight = 'normal';
                leg.EdgeColor = rgb('Silver');
                leg.Location = 'northwest';
            
                leg.delete;
            if strcmp(savePlot,'ON')
                %grid on
                plotdir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/MultiRingFit/',obj.DataSet)];
                MakeDir(plotdir);
                 freePar = ConvertFixPar('freePar',obj.fixPar,'nPar',obj.nPar,'nPixels',numel(obj.RunData.MACE_Ba_T),'Mode','Reverse');
                plotname = [plotdir,sprintf('MultiRing%s_%s_%s_freePar%s_%s.pdf',...
                    obj.RingMerge,obj.RunData.RunName,obj.chi2,freePar,PlotPar)];
                export_fig(fig1,plotname,'-painters');
                fprintf('Plot saved to %s \n',plotname)
            end
        end
        end
       methods (Access = protected)
%       function cp = copyElement(obj)
%          % Shallow copy object
%          cp = copyElement@matlab.mixin.Copyable(obj);
%          % Get handle from Prop2
%          hobj = obj.ModelObj.mnu4Sq;
%          % Create default object
%          new_hobj = 0;%eval(class(hobj));
%          % Add public property values from orig object
%           Copy.propValues(new_hobj,hobj);
%          % Assign the new object to property
%          cp.ModelObj.mnu4Sq  = 0;%new_hobj;
%       end
   end
end










