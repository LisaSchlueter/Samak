% Analysis class for the KATRIN Experiment - Building Fake Data Set
%-------------------------------------------------------------------------%
%  This class contains:
%  - Initialization of any KATRIN settings
%  - Compute fake spectra
%  - Distribute in Runs, following the KATRIN data structure
%
% RunData = Structure with all information for run summaries
%-------------------------------------------------------------------------%
%  L.Schlueter              T. Lasserre
%  MPP/TUM                  CEA
%-------------------------------------------------------------------------%
%  Last Update:     (07/02/2019)

classdef FakeRunsGenerator < handle
    
    properties (Access=public)
        
        RefKATRINconfig;      % Reference KATRIN Configuration Flag
        RefKATRINconfigFile;  % Reference KATRIN Configuration File
        FakeRunType;          % Fake1, Fake2, Fake3 to create Runs from scratch or KNM1 to mimick SC data
        ExternalRingCut       % Cut 2 external Rings (label RunXXXex2.mat)
        Debug;                % Printed in blue
        
        % Run Info
        StartRunNumber;       % First Run Number - Then increment by 1 unit
        NumberOfRuns;         % Number of Runs to Generate
        RunTimeSec;           % Elapsed time for single run
        StartRunTime;         % Date / Time of first run
        StartRunTime_LuT;     % Date / Time of first run Look Up Table
        
        % original TD
        TD;
        TDfile;
        % original HV
        qU;
        nqU;
        qUfrac;
        % HV Setting Fluct        
        qUDeltaMinus;
        qUDeltaPlus;
        % Fake TD Table
        qU_LuT;               % qU Look Up Table for all Runs
        
        % KATRIN config
        opt_calc;
        opt_katrin;
        opt_wgts;
        opt_mace;
        opt_wgtsmace;
        opt_fpd;
        opt_bkg;
        opt_fsd;
        opt_doppler;
        opt_integration;
        opt_corr;
            
        % Slow Control
        SCLinearDrift;       % Add linear Drift of Sow Control Parameters
        WGTS_CD_MolPerCm2_i;
        WGTS_CD_MolPerCm2_b;
        WGTS_CD_MolPerCm2_s;
        WGTS_MolFrac_DT_i;
        WGTS_MolFrac_DT_b;
        WGTS_MolFrac_DT_s;
        WGTS_MolFrac_TT_i;
        WGTS_MolFrac_TT_b;
        WGTS_MolFrac_TT_s;
        WGTS_MolFrac_HT_i;
        WGTS_MolFrac_HT_b;
        WGTS_MolFrac_HT_s;
        BKG_RateAllFPDSec;
        
        % Paths / Flags
        FakeRunPath;
        FakeRunSave;
        TDPath;
        TDname_LuT;
        
        % Fake TBD Object
        Asimov;          % 'ON' --> no statistical fluctuations
        SimFakeObj_i;    % Initialization TBD Object
        SimFakeObj_LuT=[];  % TBD Objects Look Up Table
        
        % SC Look Up Tables - Dimension= nqU x n(subruns)
        WGTS_CD_MolPerCm2_SubRun_LuT;
        WGTS_MolFrac_DT_SubRun_LuT;
        WGTS_MolFrac_TT_SubRun_LuT;
        WGTS_MolFrac_HT_SubRun_LuT;
        
        WGTS_CD_MolPerCm2_LuT;
        WGTS_MolFrac_DT_LuT;
        WGTS_MolFrac_TT_LuT;
        WGTS_MolFrac_HT_LuT;
        
        % FPD
        FPD_SegOffPixList;  % List of Pixel (KNM1 by default)
            
    end
    
    methods % Begin Constructor
        
        function obj = FakeRunsGenerator(varargin)
            
            fprintf(2,'-------------------Start MultiRunAnalysis Constructor----------------- \n')
            %--------------------------------- Parser Start----------------------------------------%
            p = inputParser;
            
            % Generic
            p.addParameter('RefKATRINconfig','KITNuScan',@(x)ismember(x,{'FT','KNM1','KNM5y','RefFile','KITNuScan'}));
            p.addParameter('RefKATRINconfigFile','myrefconfigfile.m');
            p.addParameter('Debug','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('FakeRunType','FakeKITNuScan');
            p.addParameter('ExternalRingCut',2,@(x)(isfloat(x) && x>=0)); 
            p.addParameter('Asimov','OFF',@(x)ismember(x,{'ON','OFF'}));

            % Run Wise
            p.addParameter('StartRunNumber',1,@(x)(isfloat(x) && x>0));
            p.addParameter('NumberOfRuns',5,@(x)(isfloat(x) && x>0));
            p.addParameter('RunTimeSec',7200,@(x)(isfloat(x) && x>0));
            p.addParameter('StartRunTime',datetime([2019 03 05 14 00 00]));
            
            % sub-Runs Retarding Voltages
%            p.addParameter('TD','MTDcreator_E018575.0_40eV_B115_Ba7.0_RedF0.7_NuMF0.80_BkgF0.20_B364');
            p.addParameter('TD','KITNuScanex2');
            p.addParameter('qUDeltaMinus',-0.05,@(x)(isfloat(x) && x>0));
            p.addParameter('qUDeltaPlus',+0.05,@(x)(isfloat(x) && x>0));
            
            % Drift of SC Parameters
            p.addParameter('SCLinearDrift','OFF',@(x)ismember(x,{'ON','OFF'}));
            
            % Column Density
            p.addParameter('WGTS_CD_MolPerCm2_i',1e+17,@(x)(isfloat(x) && x>0));
            p.addParameter('WGTS_CD_MolPerCm2_b',0,@(x)(isfloat(x) && x>0));
            p.addParameter('WGTS_CD_MolPerCm2_s',0*1e+17*1e-3,@(x)(isfloat(x) && x>0));
            
            % DT
            p.addParameter('WGTS_MolFrac_DT_i',0.011,@(x)(isfloat(x) && x>0));
            p.addParameter('WGTS_MolFrac_DT_b',0.,@(x)(isfloat(x) && x>0));
            p.addParameter('WGTS_MolFrac_DT_s',0.011*1e-3,@(x)(isfloat(x) && x>0));
            % TT
            p.addParameter('WGTS_MolFrac_TT_i',0.963,@(x)(isfloat(x) && x>0));
            p.addParameter('WGTS_MolFrac_TT_b',0.,@(x)(isfloat(x) && x>0));
            p.addParameter('WGTS_MolFrac_TT_s',0.963*1e-3,@(x)(isfloat(x) && x>0));
            % HT
            p.addParameter('WGTS_MolFrac_HT_i',0.026,@(x)(isfloat(x) && x>0));
            p.addParameter('WGTS_MolFrac_HT_b',0.,@(x)(isfloat(x) && x>0));
            p.addParameter('WGTS_MolFrac_HT_s',0.026*1e-3,@(x)(isfloat(x) && x>0));
            
            % Background
            p.addParameter('BKG_RateAllFPDSec',0.388,@(x)(isfloat(x) && x>0)); %Whole FPD
           
            % Produce FakeRun File
            p.addParameter('FakeRunPath','./FakeTmp');
            p.addParameter('FakeRunSave','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('TDPath','../../simulation/katrinsetup/TD_DataBank/');

            % FPD
            p.addParameter('FPD_SegOffPixList',[]);


            p.parse(varargin{:});
            obj.RefKATRINconfig             = p.Results.RefKATRINconfig;
            obj.RefKATRINconfigFile         = p.Results.RefKATRINconfigFile;
            obj.FakeRunType                 = p.Results.FakeRunType;
            obj.ExternalRingCut             = p.Results.ExternalRingCut;
            obj.Asimov                      = p.Results.Asimov;
            obj.Debug                       = p.Results.Debug;
            obj.StartRunNumber              = p.Results.StartRunNumber;
            obj.NumberOfRuns                = p.Results.NumberOfRuns;
            obj.RunTimeSec                  = p.Results.RunTimeSec;
            obj.StartRunTime                = p.Results.StartRunTime;
            obj.TD                          = p.Results.TD;
            obj.qUDeltaMinus                = p.Results.qUDeltaMinus;
            obj.qUDeltaPlus                 = p.Results.qUDeltaPlus;
            obj.SCLinearDrift               = p.Results.SCLinearDrift;
            obj.WGTS_CD_MolPerCm2_i         = p.Results.WGTS_CD_MolPerCm2_i;
            obj.WGTS_CD_MolPerCm2_b         = p.Results.WGTS_CD_MolPerCm2_b;
            obj.WGTS_CD_MolPerCm2_s         = p.Results.WGTS_CD_MolPerCm2_s;
            obj.WGTS_MolFrac_DT_i           = p.Results.WGTS_MolFrac_DT_i;
            obj.WGTS_MolFrac_DT_b           = p.Results.WGTS_MolFrac_DT_b;
            obj.WGTS_MolFrac_DT_s           = p.Results.WGTS_MolFrac_DT_s;
            obj.WGTS_MolFrac_TT_i           = p.Results.WGTS_MolFrac_TT_i;
            obj.WGTS_MolFrac_TT_b           = p.Results.WGTS_MolFrac_TT_b;
            obj.WGTS_MolFrac_TT_s           = p.Results.WGTS_MolFrac_TT_s;
            obj.WGTS_MolFrac_HT_i           = p.Results.WGTS_MolFrac_HT_i;
            obj.WGTS_MolFrac_HT_b           = p.Results.WGTS_MolFrac_HT_b;
            obj.WGTS_MolFrac_HT_s           = p.Results.WGTS_MolFrac_HT_s;
            obj.BKG_RateAllFPDSec           = p.Results.BKG_RateAllFPDSec;
            obj.FakeRunPath                 = p.Results.FakeRunPath;
            obj.FakeRunSave                 = p.Results.FakeRunSave;
            obj.TDPath                      = p.Results.TDPath;
            obj.FPD_SegOffPixList           = p.Results.FPD_SegOffPixList;
            %------------------------------- Parser End----------------------------------------%
            
            % Load Original TD
            LoadOriginalTD(obj);
            
            % Load Original KATRIN configuration
            LoadKATRINconfig(obj);
            
            % Build TD LookUp Table
            DrawFakeRunsTD(obj);
            SaveFakeRunsTD(obj);
                        
            % Associate Time
            SetRunTimeStart(obj);
            
            % Build SC LookUp Table
            DrawFakeRunsSC(obj);
            
            % Build Fake Runs & Save
            %DrawFakeRunsTBDIS_Uniform(obj);
            
            % Save
            %SaveFakeRuns(obj);
            
            fprintf(2,'-------------------  End MultiRunAnalysis Constructor----------------- \n')
        end
        
    end     % End Constructor
    
    methods %
        
        function LoadOriginalTD(obj)
            % FakeRunsGenerator: Load desired TD
            obj.TDfile            = load([obj.TD '.mat']);
            obj.nqU               = numel(obj.TDfile.qU);
            obj.qU                = obj.TDfile.qU + obj.qUDeltaMinus  + (obj.qUDeltaPlus-obj.qUDeltaMinus).*rand(obj.nqU,1);
            obj.qUfrac            = obj.TDfile.qUfrac;
        end
        
        function DrawFakeRunsTD(obj)
            % qU Look Up Table for all Runs
            obj.qU_LuT   = repmat(obj.qU + obj.qUDeltaMinus,1,obj.NumberOfRuns) ...
                + repmat((obj.qUDeltaPlus-obj.qUDeltaMinus),1,obj.NumberOfRuns).*rand(obj.nqU,obj.NumberOfRuns);
        end
        
        function SaveFakeRunsTD(obj)
            % Save FakeRun TD's in DataBank
            for RunNr=1:1:obj.NumberOfRuns
                TDStruct.RunTime            = obj.RunTimeSec;
                TDStruct.qU                 = obj.qU_LuT(:,RunNr);
                TDStruct.qUfrac             = obj.qUfrac;
                TDStruct.TD                 = [obj.FakeRunType 'Run' num2str(RunNr)];
                save([obj.TDPath obj.FakeRunType 'Run' num2str(RunNr) '.mat'],'-struct','TDStruct');
                obj.TDname_LuT{RunNr}       = TDStruct.TD;
                switch obj.Debug
                    case 'ON'
                    cprintf('blue','Save %s%sRun%s.mat \n',obj.TDPath,obj.FakeRunType,num2str(RunNr));
                end
            end
        end
        
        function CleanFakeRunsTD(obj)
            % Remove FakeRun TD's from DataBank
            command = ['rm ' obj.TDPath obj.FakeRunType 'Run*.mat'];
            system(command);
            for RunNr=1:1:obj.NumberOfRuns
                switch obj.Debug
                    case 'ON'
                    cprintf('blue','Remove %s%sRun%s.mat \n',obj.TDPath,obj.FakeRunType,num2str(RunNr));
                end
            end
        end
        
        
        function LoadKATRINconfig(obj)
            % FakeRunsGenerator: Load KATRING desired settings
            switch obj.RefKATRINconfig
                case {'KNM1','KNM5y','KITNuScan'}
                    nTeBinningFactor = 50;
                    ISCS             = 'Theory'; 
                    Mode             = 'Read';
                    WGTS_B_T         = 2.52;
                    MACE_Bmax_T      = 4.2;
                    MACE_Ba_T        = 6e-4;
                    mnuSq_i          = 0;
                    Q_i              = 18574;
                    DTFSD            = 'HTFSD';
                    HTFSD            = 'SAENZ';
                    TTFSD            = 'SAENZ';
                    ELossFlag        = 'Abdurashitov';
                    DopplerEffectFlag= 'OFF';
                    FPD_MeanEff      = 0.95;
                    FPD_ROIlow       = 14;
                    FPD_ROIEff       = 'OFF';
                    FPD_PileUpEff    = 'OFF';
                    BKG_Flag         = 'ON';
                    BKG_Type         = 'FLAT';
                    FPD_Segmentation = 'OFF';
                    KTFFlag          = 'WGTSMACE';
                    recomputeRF      = 'OFF';
                    UseParallelRF    = 'ON';
                    FPD_Pixel        = 1;
                    FPD_Ring         = 1;
                    RadiativeFlag    = 'ON';
                    NIS              = 11;
                    
                    if isempty(obj.FPD_SegOffPixList)
                        % Pixel List
                        PixelList=1:148;
                        % FPD
                        PixelexclusionFPD =[97, 98, 110, 111, 121, 122]+1;
                        PixelList=PixelList(~ismember(PixelList,PixelexclusionFPD));
                        % FBM
                        PixelexclusionFBM =[100, 112, 123]+1;
                        PixelList=PixelList(~ismember(PixelList,PixelexclusionFBM));
                        % Alignement
                        PixelexclusionAli=[103:105, 113:117, 124:130, 135, 136:142, 147]+1;
                        obj.FPD_SegOffPixList=PixelList(~ismember(PixelList,PixelexclusionAli));
                    end
            end
            
            obj.opt_calc = {...
                'nTeBinningFactor',nTeBinningFactor};
            
            obj.opt_katrin = {...
                'Mode',Mode,...
                'TD',obj.TD,...
                'TimeSec',obj.RunTimeSec,...
                'mnuSq_i',mnuSq_i,...
                'Q_i',Q_i};
            
            obj.opt_wgts = {...
                'WGTS_MolFrac_TT',obj.WGTS_MolFrac_TT_i,...
                'WGTS_MolFrac_TT_SubRun',repmat(obj.WGTS_MolFrac_TT_i,1,obj.nqU),...
                'WGTS_MolFrac_DT',obj.WGTS_MolFrac_DT_i,...
                'WGTS_MolFrac_DT_SubRun',repmat(obj.WGTS_MolFrac_DT_i,1,obj.nqU),...
                'WGTS_MolFrac_HT',obj.WGTS_MolFrac_HT_i,...
                'WGTS_MolFrac_HT_SubRun',repmat(obj.WGTS_MolFrac_HT_i,1,obj.nqU),...
                'WGTS_MolFracRelErr_TT',0,...
                'WGTS_MolFracRelErr_DT',0,...
                'WGTS_MolFracRelErr_HT',0,...
                'WGTS_DTHTr',1,...
                'WGTS_FTR_cm',4.1,...
                'WGTS_CD_MolPerCm2',obj.WGTS_CD_MolPerCm2_i,...
                'WGTS_CD_MolPerCm2_SubRun',repmat(obj.WGTS_CD_MolPerCm2_i,1,obj.nqU),...
                'WGTS_B_T',WGTS_B_T,...
                'WGTS_Temp',30,...
                'ISCS',ISCS,...
                'NIS',NIS,...
                'ELossFlag',ELossFlag,...
                'recomputeRF',recomputeRF,...
                'UseParallelRF',UseParallelRF};
            
            obj.opt_mace = {...
                'MACE_Bmax_T',MACE_Bmax_T,...
                'MACE_Ba_T',MACE_Ba_T};
            
            obj.opt_wgtsmace = {...
                'KTFFlag',KTFFlag};
            
            obj.opt_fpd = {...
                'FPD_Segmentation',FPD_Segmentation,...
                'FPD_SegOffPixList',obj.FPD_SegOffPixList,...
                'FPD_Pixel',FPD_Pixel,...
                'FPD_Ring',FPD_Ring,...
                'FPD_MeanEff',FPD_MeanEff,...
                'FPD_ROIEff',FPD_ROIEff,...
                'FPD_ROIlow',FPD_ROIlow,...
                'FPD_PileUpEff',FPD_PileUpEff};
            
            obj.opt_bkg = {...
                'BKG_Flag',BKG_Flag,...
                'BKG_Type',BKG_Type,...
                'BKG_RateAllFPDSec',obj.BKG_RateAllFPDSec};
            
            obj.opt_fsd= {...
                'TTFSD',TTFSD,...
                'DTFSD',DTFSD,...
                'HTFSD',HTFSD};
            
            obj.opt_doppler = {'DopplerEffectFlag',DopplerEffectFlag};
            
            obj.opt_integration = {'IStype','SIMPFAST'};
            
            obj.opt_corr= {...
                'RadiativeFlag',RadiativeFlag};
            
            % Tritium spectrum definition
            obj.SimFakeObj_i  = TBD(...
                obj.opt_corr{:},...
                obj.opt_calc{:},...
                obj.opt_katrin{:},...
                obj.opt_wgts{:},...
                obj.opt_mace{:},...
                obj.opt_fpd{:},...
                obj.opt_wgtsmace{:},...
                obj.opt_fsd{:},...
                obj.opt_bkg{:},...
                obj.opt_doppler{:},...
                obj.opt_integration{:});
            
            obj.SimFakeObj_i.ComputeTBDDS; obj.SimFakeObj_i.ComputeTBDIS;

        end
        
        function DrawFakeRunsSC(obj)
            % Draw SC values for all Fake Runs
            switch obj.SCLinearDrift
                case 'OFF'
            obj.WGTS_CD_MolPerCm2_SubRun_LuT = obj.GenSC_GaussianRandomBias(obj.WGTS_CD_MolPerCm2_i,obj.WGTS_CD_MolPerCm2_b,obj.WGTS_CD_MolPerCm2_s,obj.nqU,obj.NumberOfRuns);
            obj.WGTS_MolFrac_DT_SubRun_LuT   = obj.GenSC_GaussianRandomBias(obj.WGTS_MolFrac_DT_i,obj.WGTS_MolFrac_DT_b,obj.WGTS_MolFrac_DT_s,obj.nqU,obj.NumberOfRuns);
            obj.WGTS_MolFrac_TT_SubRun_LuT   = obj.GenSC_GaussianRandomBias(obj.WGTS_MolFrac_TT_i,obj.WGTS_MolFrac_TT_b,obj.WGTS_MolFrac_TT_s,obj.nqU,obj.NumberOfRuns);
            obj.WGTS_MolFrac_HT_SubRun_LuT   = obj.GenSC_GaussianRandomBias(obj.WGTS_MolFrac_HT_i,obj.WGTS_MolFrac_HT_b,obj.WGTS_MolFrac_HT_s,obj.nqU,obj.NumberOfRuns);
                case 'ON'
            obj.WGTS_CD_MolPerCm2_SubRun_LuT = obj.GenSC_GaussianRandomBiasDrift(obj.WGTS_CD_MolPerCm2_i,obj.WGTS_CD_MolPerCm2_b,obj.WGTS_CD_MolPerCm2_s,obj.nqU,obj.NumberOfRuns,-1e-05*obj.WGTS_CD_MolPerCm2_i);
            obj.WGTS_MolFrac_DT_SubRun_LuT   = obj.GenSC_GaussianRandomBiasDrift(obj.WGTS_MolFrac_DT_i,obj.WGTS_MolFrac_DT_b,obj.WGTS_MolFrac_DT_s,obj.nqU,obj.NumberOfRuns,0);
            obj.WGTS_MolFrac_TT_SubRun_LuT   = obj.GenSC_GaussianRandomBiasDrift(obj.WGTS_MolFrac_TT_i,obj.WGTS_MolFrac_TT_b,obj.WGTS_MolFrac_TT_s,obj.nqU,obj.NumberOfRuns,0);
            obj.WGTS_MolFrac_HT_SubRun_LuT   = obj.GenSC_GaussianRandomBiasDrift(obj.WGTS_MolFrac_HT_i,obj.WGTS_MolFrac_HT_b,obj.WGTS_MolFrac_HT_s,obj.nqU,obj.NumberOfRuns,0);
            end
                        
            % Average over sub-run values
            obj.WGTS_CD_MolPerCm2_LuT         = mean(obj.WGTS_CD_MolPerCm2_SubRun_LuT,1); 
            obj.WGTS_MolFrac_DT_LuT           = mean(obj.WGTS_MolFrac_DT_SubRun_LuT,1);
            obj.WGTS_MolFrac_TT_LuT           = mean(obj.WGTS_MolFrac_TT_SubRun_LuT,1);
            obj.WGTS_MolFrac_HT_LuT           = mean(obj.WGTS_MolFrac_HT_SubRun_LuT,1);
        end
        
        function DrawFakeRunsTBDIS_Uniform(obj)
            % Draw Random TBDIS Spectra
            % Assuming Uniform FPD (no pixel-wise generation)
            % Renormalise to number of rings excluded
            % (assume uniform background...)
            
            for Run=1:1:obj.NumberOfRuns
                switch obj.Debug
                   case 'ON'
                       cprintf('blue','--- Build %s Data Set - Run %0.0f \n',obj.FakeRunType,Run);
                       cprintf('blue','Use TD: %s \n',obj.TDname_LuT{Run});
                end
               
                tmpTBD = TBD(...
                    obj.opt_corr{:},...
                    obj.opt_calc{:},...
                    obj.opt_katrin{:},...
                    obj.opt_wgts{:},...
                    obj.opt_mace{:},...
                    obj.opt_fpd{:},...
                    obj.opt_wgtsmace{:},...
                    obj.opt_fsd{:},...
                    obj.opt_bkg{:},...
                    obj.opt_doppler{:},...
                    obj.opt_integration{:},...
                    'TD',obj.TDname_LuT{Run},...
                    'TimeSec',obj.RunTimeSec ,...
                    'WGTS_CD_MolPerCm2',obj.WGTS_CD_MolPerCm2_LuT(Run),...
                    'WGTS_MolFrac_TT',obj.WGTS_MolFrac_TT_LuT(Run),...
                    'WGTS_MolFrac_HT',obj.WGTS_MolFrac_HT_LuT(Run),...
                    'WGTS_MolFrac_DT',obj.WGTS_MolFrac_DT_LuT(Run),...
                    'BKG_RateAllFPDSec',obj.BKG_RateAllFPDSec);
                
                tmpTBD.ComputeTBDDS;
                tmpTBD.ComputeTBDIS;
                
                % Apply External ring Cut
                tmpTBD.TBDIS  = tmpTBD.TBDIS .* (4+12*(12-obj.ExternalRingCut))/148;
                
                % Uncertainty Stat
                s = (tmpTBD.TBDIS - obj.RunTimeSec .* obj.qUfrac .* obj.BKG_RateAllFPDSec .* (4+12*(12-obj.ExternalRingCut))/148);
                b = obj.RunTimeSec .* obj.qUfrac .* obj.BKG_RateAllFPDSec .* (4+12*(12-obj.ExternalRingCut))/148;
                %tmpTBD.TBDISE = sqrt(s.^2+b.^2);
                tmpTBD.TBDISE = sqrt(tmpTBD.TBDIS);

                % Asimov Data set ?
                switch obj.Asimov
                    case 'OFF'
                        tmpTBD.TBDIS = tmpTBD.TBDIS + tmpTBD.TBDISE.*randn(obj.nqU,1);
                end
                                                
                obj.SimFakeObj_LuT{Run} = tmpTBD;
                
                switch obj.Debug
                    case 'ON'
                        cprintf('blue','DrawFakeRunsTBDIS_Uniform: TBDIS Run %0.f Computed\n',Run);
                end
                
            end
        end
        
        function SetRunTimeStart(obj)
            % Set Date/Time for each Fake runs
            obj.StartRunTime_LuT = obj.StartRunTime;
            for RunNr=2:1:obj.NumberOfRuns
                obj.StartRunTime_LuT(RunNr) = obj.StartRunTime_LuT(RunNr-1) + seconds(obj.RunTimeSec);
            end
        end
            
        function SaveFakeRuns(obj)
            % Save Fake Runs Files in dedicated folder
            % Can be then used with Ru-n MultiRun- Analyses Classes
            if ~exist(obj.FakeRunPath, 'dir')
                mkdir(obj.FakeRunPath)
            end
            
            for RunNr=1:1:obj.NumberOfRuns
                % Build Fake Data Structure
                TmpRunStruct.PixelList                 = obj.FPD_SegOffPixList;
                TmpRunStruct.TD                        = obj.TDname_LuT{RunNr};
                TmpRunStruct.qU                        = obj.qU_LuT(:,RunNr);
                TmpRunStruct.qUfrac                    = obj.qUfrac;
                TmpRunStruct.RunNr                     = RunNr;
                TmpRunStruct.TimeSec                   = obj.RunTimeSec;
                TmpRunStruct.RunTimeStart              = obj.StartRunTime_LuT(RunNr);
                TmpRunStruct.WGTS_CD_MolPerCm2_SubRun  = obj.WGTS_CD_MolPerCm2_SubRun_LuT(:,RunNr);
                TmpRunStruct.WGTS_CD_MolPerCm2         = obj.WGTS_CD_MolPerCm2_LuT(RunNr);
                TmpRunStruct.WGTS_MolFrac_DT_SubRun    = obj.WGTS_MolFrac_DT_SubRun_LuT(:,RunNr);
                TmpRunStruct.WGTS_MolFrac_DT           = obj.WGTS_MolFrac_DT_LuT(RunNr);
                TmpRunStruct.WGTS_MolFrac_TT_SubRun    = obj.WGTS_MolFrac_TT_SubRun_LuT(:,RunNr);
                TmpRunStruct.WGTS_MolFrac_TT           = obj.WGTS_MolFrac_TT_LuT(RunNr);
                TmpRunStruct.WGTS_MolFrac_HT_SubRun    = obj.WGTS_MolFrac_HT_SubRun_LuT(:,RunNr);
                TmpRunStruct.WGTS_MolFrac_HT           = obj.WGTS_MolFrac_HT_LuT(RunNr);
                TmpRunStruct.TBDIS                     = obj.SimFakeObj_LuT{RunNr}.TBDIS;
                TmpRunStruct.TBDISE                    = sqrt(TmpRunStruct.TBDIS);
                SimFakeObj_i=obj.SimFakeObj_i;
                if obj.ExternalRingCut>0
                    save([obj.FakeRunPath '/' obj.FakeRunType 'Run' num2str(RunNr) 'ex' num2str(obj.ExternalRingCut) '.mat'],'-struct','TmpRunStruct');
                    save([obj.FakeRunPath '/' obj.FakeRunType 'ex' num2str(obj.ExternalRingCut) 'RunObject.mat'],'SimFakeObj_i');
                else
                    save([obj.FakeRunPath '/' obj.FakeRunType 'Run' num2str(RunNr) '.mat'],'-struct','TmpRunStruct');
                    save([obj.FakeRunPath '/' obj.FakeRunType 'RunObject.mat'],'SimFakeObj_i');
                end
            end
        end
        
        function CleanFakeRuns(obj)
            % Remove FakeRun from Storage
            if ~exist(obj.FakeRunPath, 'dir')
                return;
            end
            
            command = ['rm ' obj.FakeRunPath '/' obj.FakeRunType 'Run*.mat'];
            system(command);
            
            for RunNr=1:1:obj.NumberOfRuns
                switch obj.Debug
                    case 'ON'
                        if obj.ExternalRingCut>0
                            cprintf('blue','Remove %s/%sRun%sex%0.0f.mat \n',obj.FakeRunPath,obj.FakeRunType,num2str(RunNr),obj.ExternalRingCut);
                        else
                            cprintf('blue','Remove %s/%sRun%s.mat \n',obj.FakeRunPath,obj.FakeRunType,num2str(RunNr));
                        end
                end
            end
        end
        
        function RandomValues = GenSC_GaussianRandomBiasDrift(obj,m,b,s,nqu,nruns,bdrift)
            % Generate Gaussian Distributed Random SC Value
            % Include a liner Drift / Hours
            % m = expected value
            % b = bias
            % s = 1sigma uncertainty
            % For nqu x nruns
            % bdrift = linear drift of b per hour
            
            bdrift_LuT=[]; bdrift_LuT(1)=b;
            for RunNr=2:1:nruns
                bdrift_LuT(RunNr) =  bdrift_LuT(RunNr-1) ...
                    + bdrift * hours(obj.StartRunTime_LuT(RunNr)-obj.StartRunTime_LuT(RunNr-1));
            end
            RandomValues = squeeze(repmat(m,1,nqu,nruns)) + squeeze(repmat(bdrift_LuT,nqu,1)) + squeeze(randn(1,nqu,nruns)).*s;
        end
        
    end
    
    methods (Static)
        function RandomValues = GenSC_GaussianRandomBias(m,b,s,nqu,nruns)
            % Generate Gaussian Distributed Random SC Value
            % m = expected value
            % b = bias
            % s = 1sigma uncertainty
            % For nqu x nruns
            RandomValues = squeeze(repmat(m + b,1,nqu,nruns) + (randn(1,nqu,nruns)).*s);
        end
    end
end
