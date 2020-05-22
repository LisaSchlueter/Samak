% ------------------------------------------------------------------- % 
    %   Class to Create & Handle Covariance Matrices
    %  
    % - Computation: Options to switch ON/OFF certain Sys. Effects
    % - Computation of Fractional CM out of CM (and vice versa)
    % - Reading CM from File
    % - Plotting CM
    % - Sanity Checks: Rank (& more to come)
    % 
    %   L. Schlueter & Thierry Lasserre
    %   TUM / MPP      CEA
    %   04/2018
    % ------------------------------------------------------------------- % 
    
classdef CovarianceMatrix < handle
    
    properties (Access=public)
        
    CovMat;                   % Covariance Matrix  
    CovMatFrac;               % Fractional CovMat
    CovMatFracShape;          % Shape + Mixed CM (Decomposed)
    CovMatFracNorm;           % Norm CM          (Decomposed)
    CovMatFile;               % path/filename, where CM is stored
    
    MultiCovMat;              % Cell where all CM used are stored
    MultiCovMatFrac;
    MultiCovMatFracShape;         
    MultiCovMatFracNorm;
    
    nTrials;                  % Number of Trials
    StudyObject;              % TBD or Kr Object
    AltObject;                % Alternative TBD Object -> Cell with n Objects  
    SysEffect;                % Struct with Flags (ON/OFF) for all Sys. Effects
    RecomputeFlag;            % Recompute LookupTables    
    nRuns;                    % Number of Runs (only available for column density atm)
    RunFlag;                  % Run Analysis ON/OFF    
    
    % Uncertainties for RF Parameters
    MACE_Bmax_T_RelErr;       % B-Fields
    MACE_Ba_T_RelErr;
    WGTS_B_T_RelErr;
    WGTS_CD_MolPerCm2_RelErr; % Column density
    ISXsection_RelErr;        % Inelastic scattering(IS) cross section  
   
    % Uncertainties for TASR Parameters (subRun activity fluctuation)
    WGTS_TASR_RelErr;         % Relative Error on subRun Tritium activity
    
    % Uncertainties on FSDs
    FSDNorm_RelErr;           % Normalization Error on ground state (sum of GS+ES probability is still always 100%)
    FSDShapeGS_RelErr;        % Bin-to-bin uncorrelated uncertainty of FSD GROUND  STATE
    FSDShapeES_RelErr;        % Bin-to-bin uncorrelated uncertainty of FSD EXCITED STATE
    DTFlag;                   % Isotopologues to be considered
    HTFlag;
    TTFlag;
    
    % Uncertainty FPD
    FPDeff_RelErr;
    
    % Uncertainties on Background Model 1 - Shape - 
    % Decreasing slope with respect to flat background
    BM1_RatePerSec;           % Background Rate Per Second
    BM1_RedFactor;            % Decreasing fraction, with respect to background at endpoint
    BM1_RefSlopeqU;           % Reference potential to fix the background slope of alternative model 
    NonPoissonScaleFactor;    % Non poissonian background fluctuation factor - horizontal Array of Ring-wise
    
    % Uncertainties Stacking 
    Stack_qUErr;
    Stack_qUfracRelErr;
    nStack;
    
    % plasma uncertainties
    RW_SigmaErr;             % absolute uncertainty of plasma/rear wall potential over measurement time (eV)
    RW_MultiPosErr;
    
    % Plasma shifts
    E0Offsets;
    E0OffsetsErr;
    
    % longitidinal
    MACE_VarErr;
    is_EOffsetErr;
    
    %debugging option
    SanityPlots;             % various ON/OFF
    
    TDlabel % information which on binning (usually RunNr)
    
    % Shape-only Covariance Matrix Computation
    ShapeCMnormIndex;
    
    end
    
    methods
        %Constructor
        function obj = CovarianceMatrix(varargin)
         cprintf('blue','-------------------START: CovarianceMatrix Constructor ---------------------- \n')
         p = inputParser; 
         p.addParameter('CovMatFile','' ,@(x)ischar(x)); 
         p.addParameter('nTrials',1000,@(x)isfloat(x) && x>-1);
         p.addParameter('RunFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
         p.addParameter('nRuns',1, @(x)isfloat(x) && x>0); %if RunFlag== OFF ---> nRuns=1 automatically        
         p.addParameter('StudyObject','' , @(x) isa(x,'TBD') || isa(x,'Kr'));     
         p.addParameter('AltObject','' );
         p.addParameter('MACE_Ba_T_RelErr',1e-2,@(x)isfloat(x) && x>=0);
         p.addParameter('MACE_Bmax_T_RelErr',1e-2,@(x)isfloat(x) && x>=0);    
         p.addParameter('WGTS_B_T_RelErr',1e-2,@(x)isfloat(x) && x>=0);
         p.addParameter('WGTS_CD_MolPerCm2_RelErr',0.05,@(x)isfloat(x) && x>=0);
         p.addParameter('ISXsection_RelErr',1e-2,@(x)isfloat(x) && x>=0); 
         p.addParameter('RecomputeFlag','ON',@(x)ismember(x,{'ON','OFF'}));
         p.addParameter('SysEffect','',@(x)isstruct(x));        
         p.addParameter('WGTS_TASR_RelErr',1e-2,@(x)all(isfloat(x))); % can also be vector with subrun values
         p.addParameter('FSDNorm_RelErr',1e-2,@(x)isfloat(x) && x>=0);
         p.addParameter('FSDShapeGS_RelErr',4e-2,@(x)isfloat(x) && x>=0);
         p.addParameter('FSDShapeES_RelErr',21e-2,@(x)isfloat(x) && x>=0);
         p.addParameter('Stack_qUErr',[],@(x)isfloat(x));
         p.addParameter('Stack_qUfracRelErr',[],@(x)isfloat(x));
         p.addParameter('nStack',100,@(x)isfloat(x));
         p.addParameter('DTFlag','ON',@(x)ismember(x,{'ON','OFF'}));
         p.addParameter('TTFlag','ON',@(x)ismember(x,{'ON','OFF'}));
         p.addParameter('HTFlag','ON',@(x)ismember(x,{'ON','OFF'}));
         p.addParameter('BM1_RatePerSec',0.33,@(x)isfloat(x) && x>=0);
         p.addParameter('BM1_RedFactor',0.2,@(x)isfloat(x) && x>=0);
         p.addParameter('BM1_RefSlopeqU',16575,@(x)isfloat(x) && x>=0);
         p.addParameter('NonPoissonScaleFactor',[1.064 1.064],@(x)isfloat(x) && all(x>0));
         p.addParameter('SanityPlots','OFF',@(x)ismember(x,{'ON','OFF'}));
         p.addParameter('ShapeCMnormIndex',2,@(x)isfloat(x) && x>=0);
         p.addParameter('FPDeff_RelErr',0.001,@(x)isfloat(x) && x>=0);
         p.addParameter('RW_SigmaErr',0.1,@(x)isfloat(x) && x>=0);
         p.addParameter('RW_MultiPosErr',0.1,@(x)isfloat(x));
         p.addParameter('MACE_VarErr',0.02,@(x)isfloat(x));
         p.addParameter('is_EOffsetErr',0.02,@(x)isfloat(x));
         
         p.parse(varargin{:});          
         obj.nTrials                  = p.Results.nTrials;
         obj.nRuns                    = p.Results.nRuns;
         obj.RunFlag                  = p.Results.RunFlag;
         obj.SysEffect                = p.Results.SysEffect;
         obj.StudyObject              = p.Results.StudyObject;
         obj.AltObject                = p.Results.AltObject;
         obj.MACE_Ba_T_RelErr         = p.Results.MACE_Ba_T_RelErr;
         obj.MACE_Bmax_T_RelErr       = p.Results.MACE_Bmax_T_RelErr;
         obj.WGTS_B_T_RelErr          = p.Results.WGTS_B_T_RelErr;
         obj.WGTS_CD_MolPerCm2_RelErr = p.Results.WGTS_CD_MolPerCm2_RelErr;
         obj.ISXsection_RelErr        = p.Results.ISXsection_RelErr;                
         obj.CovMatFile               = p.Results.CovMatFile;
         obj.RecomputeFlag            = p.Results.RecomputeFlag;
         obj.WGTS_TASR_RelErr         = p.Results.WGTS_TASR_RelErr;
         obj.FSDNorm_RelErr           = p.Results.FSDNorm_RelErr;   
         obj.FSDShapeGS_RelErr        = p.Results.FSDShapeGS_RelErr;
         obj.FSDShapeES_RelErr        = p.Results.FSDShapeES_RelErr;
         obj.Stack_qUErr              = p.Results.Stack_qUErr;
         obj.Stack_qUfracRelErr       = p.Results.Stack_qUfracRelErr; 
         obj.FPDeff_RelErr            = p.Results.FPDeff_RelErr;
         obj.nStack                   = p.Results.nStack;
         obj.DTFlag                   = p.Results.DTFlag;
         obj.TTFlag                   = p.Results.TTFlag;
         obj.HTFlag                   = p.Results.HTFlag;
         obj.BM1_RatePerSec           = p.Results.BM1_RatePerSec;
         obj.BM1_RedFactor            = p.Results.BM1_RedFactor;
         obj.BM1_RefSlopeqU           = p.Results.BM1_RefSlopeqU;
         obj.NonPoissonScaleFactor    = p.Results.NonPoissonScaleFactor;
         obj.SanityPlots              = p.Results.SanityPlots;
         obj.ShapeCMnormIndex         = p.Results.ShapeCMnormIndex;
         obj.RW_SigmaErr              = p.Results.RW_SigmaErr;
         obj.RW_MultiPosErr           = p.Results.RW_MultiPosErr;
         obj.MACE_VarErr            = p.Results.MACE_VarErr;
         obj.is_EOffsetErr            = p.Results.is_EOffsetErr;
         
         if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
             dim = obj.StudyObject.nqU*obj.StudyObject.nRings;
         else
             dim = obj.StudyObject.nqU;
         end
                    
         obj.CovMat                   = zeros(dim);
         obj.CovMatFrac               = zeros(dim);
         obj.CovMatFracShape          = zeros(dim);
         obj.CovMatFracNorm           = zeros(dim);
         
         GetSamakPath;
         
         SysEffectList = {... % Cell Sys Effects in this Class --> to be updated constantly
             'RF_EL',...      % Response Function(RF) EnergyLoss
             'RF_BF',...      % R F B-Fields
             'RF_RX',...      % RF Column Density, Cross Section
             'TASR',...       % Tritium Activity Stack Runs
             'FSD',...        % FSD Normalization GS/ES
             'DOPoff',...     % Doppler (if switched off)         
             'TCoff_RAD',...  % Radiative corrections (if switched off) 
             'TCoff_OTHER',...% Other corrections (if switched off)
             'BM1S',...       % Background Alternative Model 1 - Shape Only
             'HybStack',...   % Hybrid Stacking (Very First Tritium)
             'Stack',...      % Stacking all qU
             'BkgShape',...   % background shape (slope)
             'FPDeff',...     % detector efficiency from subrun to subrun
             'RW',...        % rear wall /plasma potential over time
             'LongPlasma',...  % longitudinal plasma uncertainty (e-loss shift and variance)
             'PlasmaOffsets'};  % uncertainty on rel. endpoints for rear wall period combi
          
         % Init  SysEffect structure: Default Option: all OFF  
         if isempty(obj.SysEffect)         
             obj.SysEffect = struct(SysEffectList{1},'OFF'); % make SysEffect a structure
         end  
        for i=1:numel(SysEffectList) %if field doesnt exist
         if ~isfield(obj.SysEffect,SysEffectList{i}) 
             obj.SysEffect.(SysEffectList{i}) = 'OFF'; %create and set to OFF
         end 
        end
         
          %shortcut for theoretical corrections: when TC on -> both on
          if isfield(obj.SysEffect,'TC')
              if strcmp(obj.SysEffect.TC,'ON')
                  obj.SysEffect.TCoff_OTHER = 'ON';
                  obj.SysEffect.TCoff_RAD = 'ON';
                  obj.SysEffect = rmfield(obj.SysEffect,'TC');
              end
          end
        
        %shortcut for response function: when RF on -> all RF_ on
          if isfield(obj.SysEffect,'RF')
              if strcmp(obj.SysEffect.RF,'ON')
                  obj.SysEffect.RF_EL = 'ON';
                  obj.SysEffect.RF_BF = 'ON';
                  obj.SysEffect.RF_RX = 'ON';
                  obj.SysEffect = rmfield(obj.SysEffect,'RF');
              end
          end
          
          %shortcut for all: when all on -> set everything to ON
          if isfield(obj.SysEffect,'all')
              if strcmp(obj.SysEffect.all,'ON')
                  obj.SysEffect.RF_EL = 'ON';
                  obj.SysEffect.RF_BF = 'ON';
                  obj.SysEffect.RF_RX = 'ON';
                  obj.SysEffect.TCoff_OTHER = 'ON';
                  obj.SysEffect.TCoff_RAD = 'ON';
                  obj.SysEffect.TASR ='ON';
                  obj.SysEffect.FSD = 'ON';
                  obj.SysEffect.RW = 'ON';
                  obj.SysEffect.LongPlasma = 'ON';
                  obj.SysEffect = rmfield(obj.SysEffect,'all');
              end
          end
          
        % Init MultiCovMat --> Useful when more than 1 CM is used
        obj.MultiCovMat = struct('CM_STAT',zeros(dim),...
            'CM_RF', zeros(dim),... %all RF ON
            'CM_TCoff',zeros(dim));    % all TCoff ON
        obj.MultiCovMatFrac = struct('CM_STAT',zeros(dim),...
            'CM_RF', zeros(dim),...
            'CM_TCoff',zeros(dim));
        obj.MultiCovMatFracShape = struct('CM_STAT',zeros(dim),...
            'CM_RF', zeros(dim),...
            'CM_TCoff',zeros(obj.StudyObject.nqU));
        obj.MultiCovMatFracNorm = struct('CM_STAT',zeros(dim),...
            'CM_RF', zeros(dim),...
            'CM_TCoff',zeros(dim));
                               
          for i=1:numel(SysEffectList)
            obj.MultiCovMat.(strcat('CM_',SysEffectList{i}))          =  zeros(dim);
            obj.MultiCovMatFrac.(strcat('CM_',SysEffectList{i}))      =  zeros(dim);
            obj.MultiCovMatFracShape.(strcat('CM_',SysEffectList{i})) =  zeros(dim);
            obj.MultiCovMatFracNorm.(strcat('CM_',SysEffectList{i}))  =  zeros(dim);   
          end
         
         if strcmp(obj.RunFlag,'OFF')
             obj.nRuns = 1;
         end
         
         cprintf('blue','-------------------END:   CovarianceMatrix Constructor ---------------------- \n')
        end %end constructor
    end
    methods
        function ReadCMFile(obj,varargin)
            % Read CM from file and load all properties
            p = inputParser;
            p.addParameter('filename',obj.CovMatFile,@(x)ischar(x));
            p.addParameter('SysEffect','FSD',@(x)ischar(x));
            
            p.parse(varargin{:});
            filename  = p.Results.filename;
            SysEffect = p.Results.SysEffect;
            
            %Open CM File
            fprintf(2,'CovarianceMatrix:ReadCMFile: Loading CovMatrix File: %s \n',filename)
            %disp(path);
            
            try CM =  importdata(filename);
            catch
                fprintf(2,'--------- ERROR: CovMatrix FILE DOES NOT EXIST/ WRONG FILENAME ---------- \n')
                return
            end
            
            % Always get most current CovMat
            obj.CovMat          = CM.obj.CovMat;
            obj.CovMatFrac      = CM.obj.CovMatFrac;
            obj.CovMatFracShape = CM.obj.CovMatFracShape;
            obj.CovMatFracNorm  = CM.obj.CovMatFracNorm;
            
            % Get corresponding MultiCovMat
            fieldname = ['CM_',SysEffect];
            if strcmp(fieldname,'CM_Combi') % Combi CM
                obj.MultiCovMat = CM.obj.MultiCovMat;
                obj.MultiCovMatFrac = CM.obj.MultiCovMatFrac;
                obj.MultiCovMatFracShape = CM.obj.MultiCovMatFracShape;
                obj.MultiCovMatFracNorm = CM.obj.MultiCovMatFracNorm;
            elseif strcmp(fieldname,'CM_other') % for RF, when 2 effects are switched on
                % dont read multi CM, only obj.CovMat
            else % all others
                obj.MultiCovMat.(fieldname) = CM.obj.MultiCovMat.(fieldname);
                obj.MultiCovMatFrac.(fieldname) = CM.obj.MultiCovMatFrac.(fieldname);
                obj.MultiCovMatFracShape.(fieldname) = CM.obj.MultiCovMatFracShape.(fieldname);
                obj.MultiCovMatFracNorm.(fieldname) = CM.obj.MultiCovMatFracNorm.(fieldname);
            end
            if any(CM.obj.StudyObject.TimeSec ~= obj.StudyObject.TimeSec)
                obj.ComputeCM_STAT;
                fprintf(2,'Different Runtime used im CM (%.0f/%.0f)! launch ComputeCM_STAT ...',...
                    CM.obj.StudyObject.TimeSec,obj.StudyObject.TimeSec);
                
                if strcmp(obj.SanityPlots,'ON')
                    fig17 = figure(17);
                    set(fig17, 'Units', 'normalized', 'Position', [0.1, 0.1, 1 ,0.6]);
                    subplot(1,2,1);
                    imagesc(obj.CovMat); hold on;
                    colorbar;
                    colormap(parula);
                    title(sprintf('CM from File, Time %.2f h',CM.obj.StudyObject.TimeSec/(60*60)));
                end
                
                obj.ComputeFracCM('Mode','Frac2CM');
                
                if strcmp(obj.SanityPlots,'ON')
                    subplot(1,2,2);
                    imagesc(obj.CovMat); hold on;
                    colorbar;
                    title(sprintf('Reweighted CM, Time %.2f h',obj.StudyObject.TimeSec/(60*60)));
                    hold off;
                end
                fprintf(2,' Done\n')
            end
            
        end
        function DisplayCMInfo(obj,varargin)
           p=inputParser;
           p.addParameter('Output','screen',@(x)ismember(x,{'file','screen'}));
           p.addParameter('filename','',@(x)ischar(x));
           p.parse(varargin{:});
           Output   = p.Results.Output;
           filename = p.Results.filename; 
           
%            if ~exist('./CovMatInfo/','dir')
%                mkdir ./CovMatInfo
%            end
%           CMInfoName = sprintf('./CovMatInfo/CovarianceMatrixInfo_%s_%s.txt',obj.StudyObject.TD,filename);
            
            switch Output
                case 'file'
            fileID = fopen(filename,'w+');
                case 'screen'
                    fileID = 1;
            end
            
            fprintf(fileID,'--------------------------------------\n');
            fprintf(fileID,'       Covariance Matrix Info \n');
            fprintf(fileID,'--------------------------------------\n');
            SysEffectNames = fieldnames(obj.SysEffect);
            fprintf(fileID,'Systematic Effects \n');
            fprintf(fileID,'--------------------------------------\n');
            for i=1:numel(SysEffectNames)
                fprintf(fileID,'%s= %s \n',SysEffectNames{i},obj.SysEffect.(SysEffectNames{i}));
            end
            fprintf(fileID,'--------------------------------------\n');
            fprintf(fileID,'Uncertainties (relative) \n');
            fprintf(fileID,'WGTS_B_T_RelErr          = %.3f\n', obj.WGTS_B_T_RelErr);
            fprintf(fileID,'MACE_Ba_T_RelErr         = %.3f\n', obj.MACE_Ba_T_RelErr);
            fprintf(fileID,'MACE_Bmax_T_RelErr       = %.3f\n', obj.MACE_Bmax_T_RelErr);
            fprintf(fileID,'WGTS_CD_MolPerCm2_RelErr = %.3f\n', obj.WGTS_CD_MolPerCm2_RelErr);
            fprintf(fileID,'ISXsection_RelErr        = %.3f\n', obj.ISXsection_RelErr);
            fprintf(fileID,'WGTS_TASR_RelErr         = %.3f\n', mean(obj.WGTS_TASR_RelErr));
            fprintf(fileID,'FSDNorm_RelErr           = %.3f\n', obj.FSDNorm_RelErr);
            fprintf(fileID,'FSDShapeGS_RelErr        = %.3f\n', obj.FSDShapeGS_RelErr);
            fprintf(fileID,'FSDShapeES_RelErr        = %.3f\n', obj.FSDShapeES_RelErr);
            fprintf(fileID,'--------------------------------------\n');
            fprintf(fileID,'Background Systematics \n');
            fprintf(fileID,'Mean Non-Poissoninan Factor   = %.3f\n', mean(obj.NonPoissonScaleFactor));
            fprintf(fileID,'--------------------------------------\n');

         switch Output
             case 'file'
                 %obj.StudyObject.DisplayTDBInfo('Output','file','filename',filename);
                 fclose(fileID);
         end
        end
        function PlotCM(obj,varargin) 
            p = inputParser;
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('savename','',@(x)ischar(x)); % additional label
            p.addParameter('savedir','',@(x)ischar(x));  % save directory
            p.addParameter('ConvergenceTest','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('PlotEffect','',@(x)ischar(x));     % only in title
            p.addParameter('Mode','Frac',@(x)ismember(x,{'CM','Frac','Shape'}));
            p.addParameter('qUWindowIndexMin',90,@(x)isfloat(x)); % start display cm below endpoint (eV)
            p.addParameter('qUWindowIndexMax',10,@(x)isfloat(x)); % stop display cm below endpoint   (eV)
            p.addParameter('CovMatInput','',@(x)isfloat(x)); % plot covmt from outside (optinal)
            
            p.parse(varargin{:});
            saveplot        = p.Results.saveplot;
            savename        = p.Results.savename;
            ConvergenceTest = p.Results.ConvergenceTest; % Convergence Tests on work for a single Covariance Matrix
            PlotEffect      = p.Results.PlotEffect;      % For Labeling only: When manyeffecs are ON, but only 1 is used to compute CovMat
            qUWindowIndexMin= p.Results.qUWindowIndexMin;
            qUWindowIndexMax= p.Results.qUWindowIndexMax;
            Mode            = p.Results.Mode;
            savedir         = p.Results.savedir;
            CovMatInput     = p.Results.CovMatInput;
            
            % convert into an index
            qUWindowIndexMin = find((18574-obj.StudyObject.qU)<qUWindowIndexMin,1);
            qUWindowIndexMax = find((18574-obj.StudyObject.qU)<qUWindowIndexMax,1);
            
            if strcmp(ConvergenceTest,'ON')
                [Trials, Convergence] = obj.ConvergenceTest('Criterium','Cauchy');
            end
            obj.StudyObject.ComputeTBDDS;
            obj.StudyObject.ComputeTBDIS;
            
            %% choose  Covariance Matrix
            % spectrum without background
            TBDIS_NoBKG = obj.nRuns.*obj.StudyObject.TBDIS...
                -(obj.StudyObject.BKG_RateSec.*obj.StudyObject.qUfrac.*obj.StudyObject.TimeSec);
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                TBDIS_NoBKG = reshape(TBDIS_NoBKG,[obj.StudyObject.nqU*obj.StudyObject.nRings,1]);
            end
            
            switch Mode
                case 'CM'
                    CovMatDisp = TBDIS_NoBKG.*obj.CovMatFrac.*TBDIS_NoBKG';
                    CovMatSpec = CovMatDisp; % for sigma sys/ sigma stat plot
                    strTitle   = 'Covariance';
                case 'Frac'
                    CovMatDisp = obj.CovMatFrac; 
                    CovMatSpec = TBDIS_NoBKG.*obj.CovMatFrac.*TBDIS_NoBKG';
                    strTitle   = 'Fractional covariance';
                case 'Shape'
                    CovMatDisp = obj.CovMatFracShape;
                    CovMatSpec = TBDIS_NoBKG(1:end).*obj.CovMatFracShape.*TBDIS_NoBKG(1:end)';
                    strTitle   = 'Fractional covariance';
            end
            
            if ~isempty(CovMatInput)
                CovMatDisp = CovMatInput;
            end
            %% Plot Covariance Matrix
            fig100 = figure('Units', 'normalized', 'Position', [0.1, 0.1, 1 ,0.75]);
            FontSize = 24;
            % --------------------- cov mat ------------------------
            subplot(2,2,[1,3]);
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                nrings = numel(obj.StudyObject.MACE_Ba_T); % always gives correct number of pseudo rings
                nqU_used = size(obj.StudyObject.qU(qUWindowIndexMin:qUWindowIndexMax,:),1);             % number of subruns, which are NOT excluded
                exclIndex = sort(reshape(repmat(qUWindowIndexMin:qUWindowIndexMax,[nrings,1])+[0:nrings-1]'.*obj.StudyObject.nqU,[nqU_used*nrings,1]));
                CovMatDisp = CovMatDisp(exclIndex,exclIndex); 
            else
                CovMatDisp = CovMatDisp(qUWindowIndexMin:qUWindowIndexMax,qUWindowIndexMin:qUWindowIndexMax);
            end
            imagesc(CovMatDisp);
            colormap(parula);
            pbaspect([1 1 1])
            % axis
            xlabel(sprintf('Retarding energy - 18574 (eV)'));
             qUmin = sprintf('%.0f',obj.StudyObject.qU(qUWindowIndexMin)-obj.StudyObject.Q_i);
             qUmax = sprintf('%.0f',obj.StudyObject.qU(qUWindowIndexMax)-obj.StudyObject.Q_i);
             
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                %x-axis
                maxtmp = qUWindowIndexMax-qUWindowIndexMin;
                set(gca,'xtick',[1 maxtmp+1, size(CovMatDisp,1)-maxtmp size(CovMatDisp,1)]);%, 2*maxtmp+3,3*maxtmp+3]);
                myYticks = (1:2:nrings*2).*(1+maxtmp/2)';
                set(gca,'xticklabel',{qUmin,qUmax,qUmin,qUmax}); set(gca,'yticklabel',[]);
                %y-axis
                set(gca,'ytick',myYticks)
                set(gca,'yticklabel',string(1:nrings))
                ylabel('Ring')
                ax = gca;
                ax.YAxisLocation = 'right';
            else
                set(gca,'xtick',[1 qUWindowIndexMax-qUWindowIndexMin+1]); set(gca,'ytick',[])
              
                set(gca,'xticklabel',{qUmin,qUmax}); set(gca,'yticklabel',[]);
            end
            % colorbar
            c = colorbar('westoutside');
            c.Label.String = sprintf(strTitle);
            PrettyFigureFormat('FontSize',FontSize);
            % adjust font sizes of colorbar
            ax1 = gca;
            c.FontSize = ax1.FontSize;
            c.Label.FontSize = ax1.XLabel.FontSize;
            
            %x-axis label position (shift a little bit up)
            xpos = ax1.XLabel.Position;
            ax1.XLabel.Position = [xpos(1), xpos(2)-0.2, xpos(3)];
            %% Plot Spectrum
             if strcmp(ConvergenceTest,'ON')
                subplot(2,2,2)
             else
                 subplot(2,2,[2,4])
             end
             UncorrVarianceWrtStat = diag(CovMatSpec./TBDIS_NoBKG);
             if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                 qU = reshape(obj.StudyObject.qU,[obj.StudyObject.nqU*nrings,1]);
                 nqU = obj.StudyObject.nqU;
                 ringlabel = cell(nrings,1);
                 for i=1:nrings
                     plot(qU((i-1)*nqU+1:i*nqU)...
                         -obj.StudyObject.Q_i,...
                         UncorrVarianceWrtStat((i-1)*nqU+1:i*nqU),...
                         'LineWidth',3);
                     hold on;
                     ringlabel{i} = sprintf('Ring %.0f',i);
                 end
             else
                 qU = obj.StudyObject.qU;
                 plot(qU-obj.StudyObject.Q_i,UncorrVarianceWrtStat,...
                     'Color',rgb('GoldenRod'),'LineWidth',3);
             end
             
             
             
            TimeHour = obj.StudyObject.TimeSec(1)/(60*60)*obj.nRuns; % time in hours
            if TimeHour <= 24
                plotTime = TimeHour;  % in hours
                TimeUnit = 'hours';
            elseif TimeHour > 24 && TimeHour <= 24*365.25
                plotTime = TimeHour/24; % in days
                TimeUnit = 'days';
            elseif  TimeHour > 24*365.25
                plotTime = TimeHour/(24*365.25); % in years
                TimeUnit = 'years';
            end
            time_leg = sprintf('%.0f %s',plotTime,TimeUnit);
            xlabel(sprintf('Retarding energy - 18574 (eV)'));
            ylabel('\sigma_{syst. (uncorr.)}^2/\sigma_{stat.}^2');
            ylim([0 1.1*max(UncorrVarianceWrtStat)]);
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                leg =legend(ringlabel,'Location','northeast');
                leg.Title.String = time_leg;
                leg.EdgeColor = rgb('Silver');
            else
                leg =legend(time_leg,'Location','northeast');
                leg.EdgeColor = rgb('Silver');
            end
            PrettyFigureFormat('FontSize',FontSize);
            leg.FontSize = ax1.FontSize-1;
            xlim([min(qU-obj.StudyObject.Q_i) 0]);
            hold off;
            ax2 = gca;
            mypos = ax2.Position;
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                ax2.Position = [mypos(1)-0.01 mypos(2)+0.003 mypos(3)-0.02 mypos(4)-0.09];
                
            else
                ax2.Position = [mypos(1)-0.05 mypos(2)+0.003 mypos(3)-0.02 mypos(4)-0.09];
            end
            
            %% plot convergence text
            if strcmp(ConvergenceTest,'ON')
                subplot(2,2,4);
                plot(Trials,Convergence,'Color',rgb('DodgerBlue'),'LineWidth',4);
                xlabel('Samples');
                ylabel('Matrix norm');
                leg = legend(sprintf('Convergence test'),'Location','southeast');
                leg.EdgeColor = rgb('Silver');
                PrettyFigureFormat('FontSize',FontSize);
                leg.FontSize = ax2.FontSize;
                ylim([min(Convergence) max(Convergence)]);
                ax3 = gca;
                mypos3 = ax3.Position;
                if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                    ax3.Position = [mypos3(1)-0.01 mypos3(2)+0.089 mypos3(3)-0.02 mypos3(4)-0.09];
                else
                    ax3.Position = [mypos3(1)-0.05 mypos3(2)+0.089 mypos3(3)-0.02 mypos3(4)-0.09];
                end
                xlim([min(Trials), max(Trials)]);
            else
                ax2 = gca;
                mypos = ax2.Position;
               ax2.Position = [mypos(1)+0.02 mypos(2)+0.083 mypos(3)-0 mypos(4)-0.083];
            end

        %% save
            if strcmp(saveplot,'ON')
                %label
                effects_logic = structfun(@(x)strcmp(x,'ON'),obj.SysEffect);
                fields = fieldnames(obj.SysEffect);
                mytitle = '';
                for i=1:numel(fields)
                    %more readable labels:
                    if contains(fields{i},'_')
                        fields{i}=strrep(fields{i},'_','-');
                    end
                    if effects_logic(i)==1 && sum(effects_logic)>=2
                        mytitle = strcat(mytitle,' + ',fields(i));
                    elseif effects_logic(i)==1 && sum(effects_logic)==1
                        mytitle = fields(i);
                    end
                end
                if contains(mytitle, ' +RF_EL +RF-BF +RF_RX')
                    mytitle{:} = strrep(mytitle{:},' +RF-EL +RF-BF +RF-RX','Response Function');
                end
                if contains(mytitle, '+TCoff_RAD +TCoff_OTHER')
                    mytitle{:} = strrep(mytitle{:},'+TCoff_RAD +TCoff_OTHER','+ Theoretical Correction ON/OFF');
                end
                
                if ~isempty(PlotEffect)
                    mytitle = {PlotEffect};
                end
                
                TD = obj.StudyObject.TD;
                
                if isempty(savedir)
                    savedir = './plots/';
                end
                
                MakeDir(savedir);
                if contains(mytitle,' ') || contains(mytitle,'+') || contains(mytitle,'/')
                    mytitle{:} = strrep(mytitle{:},' ','');
                    mytitle{:} = strrep(mytitle{:},'+','_');
                    mytitle{:} = strrep(mytitle{:},'/','_');
                end
                fig100_str = [savedir,sprintf('CM_%s_%s%s.pdf',mytitle{:},TD,savename)];
                if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                    fig100_str = strrep(fig100_str,'.pdf',sprintf('_%.0frings.pdf',nrings));
                end
                export_fig(fig100,fig100_str);
                fprintf('Save plot to %s \n',fig100_str)
            end
        end
        function PlotCorr(obj,varargin)
            p = inputParser;
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('qUWindowIndexMin',90,@(x)isfloat(x)); % start display cm below endpoint (eV)
            p.addParameter('qUWindowIndexMax',10,@(x)isfloat(x)); % stop display cm below endpoint   (eV)
            p.addParameter('CovMatInput','',@(x)isfloat(x) || isempty(x));
            p.addParameter('savename','',@(x)ischar(x)); % additional label
            p.addParameter('savedir','',@(x)ischar(x));  % save directory
            
            p.parse(varargin{:});   
            saveplot         = p.Results.saveplot;
            qUWindowIndexMin = p.Results.qUWindowIndexMin;
            qUWindowIndexMax = p.Results.qUWindowIndexMax;
            CovMatInput      = p.Results.CovMatInput;
            savedir         = p.Results.savedir;
            savename        = p.Results.savename;
            
            % convert into an index
            qUWindowIndexMin = find((18574-obj.StudyObject.qU)<qUWindowIndexMin,1);
            qUWindowIndexMax = find((18574-obj.StudyObject.qU)<qUWindowIndexMax,1);
            
            % plot correlation matrix (of fractional)
            fig1= figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.5 ,0.5]);
            FontSize = 24;
            % --------------------- cov mat ------------------------
            if isempty(CovMatInput)
                CM = obj.CovMatFracShape;
            else
                CM = CovMatInput;
            end
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                nrings = numel(obj.StudyObject.MACE_Ba_T); % always gives correct number of pseudo rings
                nqU_used = size(obj.StudyObject.qU(qUWindowIndexMin:qUWindowIndexMax,:),1);             % number of subruns, which are NOT excluded
                exclIndex = sort(reshape(repmat(qUWindowIndexMin:qUWindowIndexMax,[nrings,1])+[0:nrings-1]'.*obj.StudyObject.nqU,[nqU_used*nrings,1]));
                CovMatDisp =CM(exclIndex,exclIndex);
            else
                CovMatDisp = CM(qUWindowIndexMin:qUWindowIndexMax,qUWindowIndexMin:qUWindowIndexMax);
            end
            
            [~, ~, c ] = corplot(CovMatDisp);
            colormap(parula);
            c.Label.String = sprintf('Correlation coefficient');
            pbaspect([1 1 1])
            % axis
            ax1 = gca;
            xlabel(sprintf('Retarding energy - 18574 (eV)'));
            qUmin = sprintf('%.0f',obj.StudyObject.qU(qUWindowIndexMin)-obj.StudyObject.Q_i);
            qUmax = sprintf('%.0f',obj.StudyObject.qU(qUWindowIndexMax)-obj.StudyObject.Q_i);
            ax1.XAxisLocation = 'bottom';
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                %x-axis
                maxtmp = qUWindowIndexMax-qUWindowIndexMin;
                set(gca,'xtick',[1 maxtmp+1, size(CovMatDisp,1)-maxtmp size(CovMatDisp,1)]);%, 2*maxtmp+3,3*maxtmp+3]);
                myYticks = (1:2:nrings*2).*(1+maxtmp/2)';
                set(gca,'xticklabel',{qUmin,qUmax,qUmin,qUmax}); set(gca,'yticklabel',[]);
                %y-axis
                set(gca,'ytick',myYticks)
                set(gca,'yticklabel',string(1:nrings))
                ylabel('Ring')
                ax = gca;
                ax.YAxisLocation = 'right';
            else
                set(gca,'xtick',[1 qUWindowIndexMax-qUWindowIndexMin+1]); set(gca,'ytick',[])
                set(gca,'xticklabel',{qUmin,qUmax}); set(gca,'yticklabel',[]);
            end
            
            PrettyFigureFormat('FontSize',FontSize);
            % adjust font sizes of colorbar
            c.FontSize = ax1.FontSize;
            c.Label.FontSize = ax1.XLabel.FontSize;
            
            if strcmp(saveplot,'ON')
                %label
                if isempty(savedir)
                    savedir = './plots/';
                end
                MakeDir(savedir);
                
                fig1_str = [savedir,sprintf('CorrMat%s.pdf',savename)];
                if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                    fig1_str = strrep(fig1_str,'.pdf',sprintf('_%.0frings.pdf',nrings));
                end
                export_fig(fig1,fig1_str);
                fprintf('Save plot to %s \n',fig1_str)
            end
            
        end
        function PlotStack(obj,varargin)
            p = inputParser;
            p.addParameter('plotStat','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            plotStat = p.Results.plotStat; %Option to plot stat error
            saveplot = p.Results.saveplot;
            
            effects_logic = structfun(@(x)strcmp(x,'ON'),obj.SysEffect);
            fields = fieldnames(obj.SysEffect);
            
            %exception for TC
            findElement = @(x) find(strcmp(fields,x));
            TCElements = sort([findElement('TCoff_RAD'),findElement('TCoff_OTHER')],'descend');
            if sum(effects_logic(TCElements)) == 2 %if both ON
                fields(TCElements(1)) =[];
                effects_logic(TCElements(1))=[];
                fields{TCElements(2)}= 'TCoff';
            end
            
            %exception for RF
            RFElements = sort([findElement('RF_BF'),findElement('RF_EL'),findElement('RF_RX')],'descend');
            if sum(effects_logic(RFElements)) == 3 %if all ON
                for i=1:2
                    fields(RFElements(i))=[];
                    effects_logic(RFElements(i))=[];
                end
                fields{RFElements(3)} = 'RF';
            end
            
            
            nEffects = numel(find(effects_logic==1)); % Number of Effect 'ON'
            stackdata = zeros(obj.StudyObject.nqU,numel(fields));
            %select CM
            for i=numel(fields):-1:1
                cm_fields{i} = strcat('CM_',fields{i});
                if effects_logic(i)==1 %if Effect is ON
                % stackdata(:,i)= sqrt(diag(obj.MultiCovMat.(cm_fields{i})))./obj.StudyObject.TBDIS; 
                stackdata(:,i)= diag(obj.MultiCovMat.(cm_fields{i}))./(obj.StudyObject.TBDIS);
                else
                 stackdata(:,i)=[];
                end               
            end 
    
           fstack=  figure(999); %stack plot
           set(fstack, 'Units', 'normalized', 'Position', [0.9, 0.9, 1, 1]);
            b1 = bar(obj.StudyObject.qU-obj.StudyObject.Q_i, [stackdata],'stack');         
            %b1 = bar(obj.StudyObject.Q_i-obj.StudyObject.qU, [stackdata],'stack');         
            if strcmp(plotStat,'ON')
                 hold on;
                %bar(obj.StudyObject.qU-obj.StudyObject.Q_i, (obj.StudyObject.TBDISE.^2)./obj.StudyObject.TBDIS, 'FaceAlpha', 0);           
                 lstat = line(obj.StudyObject.qU-obj.StudyObject.Q_i,ones(1,obj.StudyObject.nqU),...
                     'LineStyle','--','Color','black','LineWidth',1.5);               
            end
            PrettyFigureFormat;
            pbaspect([2 1 1])
            mycolors = {rgb('Teal'); rgb('DarkOliveGreen'); rgb('YellowGreen'); rgb('PaleGreen'); rgb('YellowGreen'); rgb('LightSlateGray');};
            set(b1,{'FaceColor'},mycolors(1:nEffects));
            xlabel('retarding potential - 18575 (V)', 'FontSize',20);
            ylabel('\sigma_{sys}^2 / \sigma_{stat}^2', 'FontSize',20);
            set(gca,'YMinorTick','off');
            set(gca,'XMinorTick','off');
            set(gca,'TickLength',[0.015 0.015]);
            xlim([1.05*min(obj.StudyObject.qU-obj.StudyObject.Q_i) 1.1*max(obj.StudyObject.qU-obj.StudyObject.Q_i)]);
            legend_str = cell(nEffects,1);
            counter = 1; %counter
            for i = 1:numel(fields)                
                if effects_logic(i)==1
                    %more readable labels: 
                    if contains(fields{i},'_')
                        fields{i}=strrep(fields{i},'_','-');                       
                    end
                    if contains(fields{i},'TASR')
                        fields{i}='Tritium Activity Fluctuation';
                    elseif contains(fields{i},'TCoff')
                        fields{i}='Theoretical Corrections';
                    elseif contains(fields{i},'RF')
                        fields{i}='Response Function (abs. column density)';
                    elseif contains(fields{i},'FSD')
                        fields{i}='Final State Distribution';
                    end
                    legend_str{counter} = fields{i};
                    counter = counter +1;
                end
                if contains(obj.StudyObject.TD,'RunStack')
                    RunsInTD = (strfind(obj.StudyObject.TD,'_'))+1;
                    firstRun = obj.StudyObject.TD(RunsInTD(1):RunsInTD(2)-2);
                    lastRun  = obj.StudyObject.TD(RunsInTD(end):RunsInTD(end)+2);
                    titleTD = sprintf('Stacked Runs %s - %s',firstRun,lastRun);
                else
                    titleTD = obj.StudyObject.TD;
                end
                title(sprintf('First Tritium Uncorrelated Uncertainties \n %s  (%.0f seconds)',titleTD,obj.StudyObject.TimeSec));
            end
            if strcmp(plotStat,'ON')
                legend_str{counter} = sprintf('Statistics');
            end
            legend(legend_str{:},'Location','northeast');
            set(gca,'FontSize',20);
            hold off;
            if strcmp(saveplot,'ON')
                publish_figure('./plots/eps/StackPlot.eps',fstack);
                publish_figurePDF(fstack,'./plots/pdf/StackPlot.pdf');
            end
        end        
        function ComputeFracCM(obj,varargin)
            % 2 Options: Mode = CM2Frac   &   Frac2CM
            % CM2Frac:   Computes and Saves Fractional Covariance Matrix
            % Frac2CM:   Computes Covariance Matrix out of Fractional CovMat
            p = inputParser;
            p.addParameter('Mode','CM2Frac', @(x)ismember(x,{'CM2Frac', 'Frac2CM'}));            
            p.parse(varargin{:});
            Mode     = p.Results.Mode;
            myfile   = importdata(obj.CovMatFile);
                    
            switch Mode
                case 'CM2Frac' %Compute fractional CovMat out of CovMat
                    cprintf('blue','--------------------------------------------------------------------------\n')
                    cprintf('blue','CovarianceMatrix: Compute Fractional CovMat out of CM \n')                                     
                    
                    %Compute TBDIS with init settings, without background!
                    if strcmp(obj.SysEffect.BM1S,'OFF')
                        myfile.obj.StudyObject.ComputeTBDDS;
                        myfile.obj.StudyObject.ComputeTBDIS;
                        TBDIS = (obj.nRuns.*myfile.obj.StudyObject.TBDIS)-...%integral spectrum without background
                            myfile.obj.StudyObject.BKG_RateSec.*myfile.obj.StudyObject.qUfrac.*myfile.obj.StudyObject.TimeSec;
                      %  myfile.obj.StudyObject.BKG_RateSec_i = BKG_prev; %reset to previous value                        
                    else %Background CM
                        TBDIS = ones(obj.StudyObject.nqU,1).*myfile.Bfit;
                    end
                    
                    %Compute Fraction Covariance Matrix
                    BkgIndex = myfile.obj.StudyObject.qU>myfile.obj.StudyObject.Q_i;
                    TBDIS(BkgIndex) = 0;       %set background to zero
                    TBDIS(TBDIS<=1e-03) = 0 ;  %to avoid huge number for background region
                    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                       TBDIS = reshape(TBDIS,[obj.StudyObject.nqU*obj.StudyObject.nRings,1]);    
                    end
                    obj.CovMatFrac = bsxfun(@rdivide,obj.CovMat,TBDIS);  %divides 1st row of CovMat by TBDIS(1), etc...
                    obj.CovMatFrac = bsxfun(@rdivide,obj.CovMatFrac,TBDIS'); %divides columnwise
                    obj.CovMatFrac(isinf(obj.CovMatFrac)) = 0; %for background region
                    obj.CovMatFrac(isnan(obj.CovMatFrac)) = 0; %for background region
                
                case 'Frac2CM' %Compute CovMat out of CovMatFrac
                    try obj.CovMatFrac = myfile.obj.CovMatFrac;
                    catch
                        fprintf(2,'Fractional CM does not exist! -----> Computing CMFrac!\n')
                        obj.CovMatFrac = obj.ComputeFracCM('Mode','CM2Frac');
                    end
                    
                    if strcmp(obj.SysEffect.BM1S,'OFF')
                        % get original integral spectrum
                        myfile.obj.StudyObject.ComputeTBDDS;
                        myfile.obj.StudyObject.ComputeTBDIS;
                        TBDIS = (obj.nRuns.*myfile.obj.StudyObject.TBDIS)-...%integral spectrum without background
                            myfile.obj.StudyObject.BKG_RateSec.*myfile.obj.StudyObject.qUfrac.*myfile.obj.StudyObject.TimeSec;
                        BkgIndex = myfile.obj.StudyObject.qU>myfile.obj.StudyObject.Q_i;
                        TBDIS(BkgIndex) = 0;       %set background to zero
                        TBDIS(TBDIS<=1e-03) = 0 ;  %to avoid huge number for background region
                        
                        if strcmp(obj.StudyObject.FPD_Segmentation,'RING') % reshape for multiring
                            TBDIS = reshape(TBDIS,[obj.StudyObject.nqU*obj.StudyObject.nRings,1]);
                        end
                    else %Background CM
                        TBDIS =  ones(obj.StudyObject.nqU,1).*obj.nRuns.*(obj.BM1_RatePerSec*obj.StudyObject.qUfrac.*obj.StudyObject.TimeSec);
                    end
                    obj.CovMat = TBDIS.*obj.CovMatFrac.*TBDIS';
                    
                    fields = fieldnames(obj.MultiCovMat);
                    for f=1:numel(fields)
                        obj.MultiCovMat.(fields{f}) = TBDIS.*obj.MultiCovMatFrac.(fields{f}).*TBDIS';
                    end
            end
            
        end   
        
        function Convert2NSPD(obj,varargin)
            % Convert to the nearest (in Frobenius norm) Symmetric Positive Definite matrix to CovMat or CovMatFrac
            p = inputParser;
            p.addParameter('Option','Frac',@(x)ismember(x,{'Frac','CM'}));
            p.parse(varargin{:});
            Option = p.Results.Option; 
            
            %Options: Decompose Fractional or normal CovMat
            if strcmp(Option,'Frac')
                obj.CovMatFrac = nearestSPD(obj.CovMatFrac);
            elseif strcmp(Option,'CM')
                obj.CovMat = nearestSPD(obj.CovMat);
            end
        end
        
        function [CMShapeOnly,CMFracShapeOnly] = DecomposeCM(obj,varargin)
             cprintf('blue','CovarianceMatrix: Decompose CM - Extract Shape Only CM \n');
            % Decompose Covariance Matrix to extract Shape Only component
            % Renormalization to the actual Fit Range
            p = inputParser;
            p.addParameter('CovMatFrac','',@(x)isfloat(x));             % input: fractional covariance matrix 
            p.addParameter('plots','OFF',@(x)ismember(x,{'ON','OFF'})); % some sanity plots
            p.addParameter('exclDataStart',1,@(x)isfloat(x));           % defines range for normalization
            p.addParameter('BkgCM','OFF',@(x)ismember(x,{'ON','OFF'})); % if covmat is background covmat
            
            p.parse(varargin{:});
            plots               = p.Results.plots;
            CovMatFrac_local    = p.Results.CovMatFrac;
            exclDataStart       = p.Results.exclDataStart; 
            BkgCM               = p.Results.BkgCM; 

            % range for normalization: given by exclDataStart
            range = obj.StudyObject.qU(exclDataStart);
            qUindex = obj.StudyObject.qU(:,1)>=range;
            nRings = size(obj.StudyObject.qU,2);
            
            if strcmp(BkgCM,'OFF') 
                % Regular tritium siganl covariance marix
                % Model TBDIS: Subtract background -> signal only!
                TBDIS_NoBKG = obj.StudyObject.TBDIS...
                    -(obj.StudyObject.BKG_RateSec.*obj.StudyObject.qUfrac.*obj.StudyObject.TimeSec);
            elseif strcmp(BkgCM,'ON')
                % Background covariance matrix
                % Use background spectrum 
                TBDIS_NoBKG = obj.StudyObject.BKG_RateSec.*obj.StudyObject.qUfrac.*obj.StudyObject.TimeSec;
            end
            
            %sum of Asiomov spectrum (ringwise)
            TBDIS_SumExpected = sum(TBDIS_NoBKG(qUindex,:));
            
            % Model TBDIS: reshape if FPD segmented
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                TBDIS_NoBKG = reshape(TBDIS_NoBKG,[obj.StudyObject.nqU.*nRings,1]);
            end
            
            %Covariance Matrix: Prepare 
            if isempty(CovMatFrac_local)
                CovMatFrac_local = obj.CovMatFrac;
            end
            CM = TBDIS_NoBKG.*CovMatFrac_local.*TBDIS_NoBKG';
            
            %Draw from multivariate distribution
            TBDIS_Sample       = mvnrnd(TBDIS_NoBKG,CM,obj.nTrials*5);
            
             % Model TBDIS: reshape back if FPD segmented
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                TBDIS_Sample = reshape(TBDIS_Sample,[size(TBDIS_Sample,1),obj.StudyObject.nqU,nRings]); 
            end
            
            % Normliazed to average spectrum (range defined withe exclDataStart and qUindex)
            TBDIS_SumSample    = squeeze(sum(TBDIS_Sample(:,qUindex,:),2));           % sum of MC spectra (sample and ringwise)
            TBDIS_SampleNorm   = permute(permute(TBDIS_Sample,[1,3,2]).*(TBDIS_SumExpected./TBDIS_SumSample),[1,3,2]);  % normalized sample spectra
            
             % Model TBDIS: reshape if FPD segmented
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                TBDIS_SampleNorm = reshape(TBDIS_SampleNorm,[size(TBDIS_Sample,1),obj.StudyObject.nqU*nRings]);
            end
            CMShapeOnly        = cov(TBDIS_SampleNorm);
           
            % compute fractional shape-only
             CMFracShapeOnly = CMShapeOnly./TBDIS_NoBKG./TBDIS_NoBKG';
             CMFracShapeOnly(isnan(CMFracShapeOnly))=0;
             CMFracShapeOnly(isinf(CMFracShapeOnly))=0;
            
             obj.CovMatFracShape = CMFracShapeOnly;
             
            %% covariance matrices plot
            if strcmp(plots,'ON')
                figure(9999);
                subplot(1,2,1);
                corplot(CMShapeOnly);
                colorbar;
                title('CovMatShape');
                PrettyFigureFormat
                pbaspect([1 1 1])
                
                subplot(1,2,2);
                imagesc(CMFracShapeOnly);
                colorbar;colormap(parula);
                title('CovMatFracShape');
                PrettyFigureFormat
                pbaspect([1 1 1])
            end
        end
        function DecomposeCMold(obj,varargin)
            % OLD method to decompose covariance matrix: not used anymore
            p = inputParser;
            p.addParameter('Option','Frac',@(x)ismember(x,{'Frac','CM'}));
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('plots','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('CombiCM','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            saveplot      = p.Results.saveplot;
            plots         = p.Results.plots;
            Option        = p.Results.Option;
            CombiCM       = p.Results.CombiC        
            cprintf('blue','CovarianceMatrix: OLD METHOD Decompose CM - Extract Shape Only CM \n');
            
            %Options: Decompose Fractional or normal CovMat
            if strcmp(Option,'Frac')
                CM = obj.CovMatFrac;
                t_str = sprintf(' FracCovMat');
            elseif strcmp(Option,'CM')
                CM = obj.CovMat;
                t_str = sprintf(' CovMat');
            end
            %init
            Term1 = zeros(length(CM)); Term2 = zeros(length(CM)); Term3 = zeros(length(CM));
            
            cm_file = importdata(obj.CovMatFile);
            if strcmp(CombiCM,'ON')
                obj.StudyObject.ComputeTBDDS; obj.StudyObject.ComputeTBDIS;
                cm_file.TBDIS_V = repmat(obj.StudyObject.TBDIS,1,obj.nTrials);
            end
            if strcmp(obj.SysEffect.BM1S,'OFF')
                Nsum = sum(cm_file.TBDIS_V,2); %Spectrum (summed over all Trials)
                Nsum = sum(Nsum,3); % in case of nRuns>1
            else
                TBDIS_BKG = cm_file.Bfit.*ones(obj.StudyObject.nqU,1);
                Nsum = sum(TBDIS_BKG);
                Nsum = repmat(Nsum,length(CM),1);
            end
            Ntotal = sum(Nsum); % total number of counts in all bins and all Trials
            Msum1 = sum(CM,1);
            Msum2 = sum(CM,2);
            Msumsum = sum(Msum2);
            
            for i=1:length(CM)
                for j=1:length(CM)
                    Term1(i,j) = Nsum(j)/Ntotal*Msum2(i);
                    Term2(i,j) = Nsum(i)/Ntotal*Msum1(j);
                    Term3(i,j) = Nsum(i)*Nsum(j)/Ntotal^2 * Msumsum;
                end
            end
            CMshape = CM - Term1 - Term2 + Term3;
            CMmixed = Term1 + Term2 -2*Term3;
            CMnorm  = Term3;
            
            if strcmp(Option,'Frac')
                obj.CovMatFracShape = CMshape + CMmixed;
                obj.CovMatFracNorm = CMnorm;
            end
            
            %% covariance matrices plot
            switch plots
                case 'ON'
                    fig9999 = figure(9999);
                    subplot(2,2,1);
                    imagesc(CM-(CMshape+CMmixed+CMnorm));
                    colormap(flipud(gray));
                    colorbar;
                    title([t_str,'-\Sigma CM']);
                    PrettyFigureFormat
                    pbaspect([1 1 1])
                    
                    subplot(2,2,2);
                    imagesc(CMshape);
                    colorbar;colormap(flipud(gray));
                    title(strcat('Shape ', t_str));
                    PrettyFigureFormat
                    pbaspect([1 1 1])
                    
                    subplot(2,2,3);
                    imagesc(CMmixed);
                    colorbar;colormap(flipud(gray));
                    title(strcat('Mixed ', t_str));
                    PrettyFigureFormat
                    pbaspect([1 1 1])
                    
                    subplot(2,2,4);
                    imagesc(CMnorm);
                    colorbar; colormap(flipud(gray));
                    title(strcat('Norm ', t_str));
                    PrettyFigureFormat
                    pbaspect([1 1 1])
                    
                    %% correlation plot
                    fig9991 = figure(9991);
                    nqU = obj.StudyObject.nqU;
                    TickStep = ceil(nqU/5);
                    
                    subplot(2,2,1);
                    corplot(CMshape+CMmixed+CMnorm);
                    %ALTERNATIVE PLOT by Thierry
                    % [~, CorMall] = Cov2Corr(nearestSPD(CMshape+CMmixed+CMnorm));
                    %imagesc(CorMall);
                    set(gca,'XTick',1:TickStep:(nqU),'XTickLabel',1:TickStep:nqU);
                    set(gca,'YTick',1:TickStep:(nqU),'YTickLabel',1:TickStep:nqU);
                    title(sprintf(strcat('Correlation\n', t_str)));
                    PrettyFigureFormat
                    pbaspect([1 1 1])
                    
                    subplot(2,2,2);
                    corplot(CMshape);
                    % [~, CorMshape] = Cov2Corr(nearestSPD(CMshape));
                    % imagesc(CorMshape);
                    set(gca,'XTick',1:TickStep:(nqU),'XTickLabel',1:TickStep:nqU);
                    set(gca,'YTick',1:TickStep:(nqU),'YTickLabel',1:TickStep:nqU);
                    title(sprintf(strcat('Correlation\nShape ', t_str)));
                    PrettyFigureFormat
                    pbaspect([1 1 1])
                    
                    subplot(2,2,3);
                    corplot(CMmixed);
                    % [~, CorMmixed] = Cov2Corr(nearestSPD(CMmixed));
                    % imagesc(CorMmixed);
                    set(gca,'XTick',1:TickStep:(nqU),'XTickLabel',1:TickStep:nqU);
                    set(gca,'YTick',1:TickStep:(nqU),'YTickLabel',1:TickStep:nqU);
                    title(sprintf(strcat('Correlation\nMixed ', t_str)));
                    PrettyFigureFormat
                    pbaspect([1 1 1])
                    
                    subplot(2,2,4);
                    corplot(CMnorm);
                    % [~, CorMnorm] = Cov2Corr(nearestSPD(CMnorm));
                    %imagesc(CorMnorm);
                    set(gca,'XTick',1:TickStep:(nqU),'XTickLabel',1:TickStep:nqU);
                    set(gca,'YTick',1:TickStep:(nqU),'YTickLabel',1:TickStep:nqU);
                    title(sprintf(strcat('Correlation\nNorm ', t_str)));
                    PrettyFigureFormat
                    pbaspect([1 1 1])
                    
                    % Plot shape uncertainty estimator
                    fig9992 = figure(9992);
                    boundedline(obj.StudyObject.qU-obj.StudyObject.Q, ones(obj.StudyObject.nqU,1),...
                        sqrt(diag((obj.nRuns.*obj.StudyObject.TBDIS).*(CMmixed).*(obj.nRuns.*obj.StudyObject.TBDIS)'))./(obj.nRuns.*obj.StudyObject.TBDIS));
                    hold on;
                    
                    errorbar(obj.StudyObject.qU-obj.StudyObject.Q, ones(obj.StudyObject.nqU,1),...
                        sqrt(obj.nRuns.*obj.StudyObject.TBDIS)./(obj.nRuns*obj.StudyObject.TBDIS),'x','Color',orange,'MarkerSize',5);
                    xlabel('qU (eV)');
                    ylabel('norm. Counts');
                    leg =legend('sys. err (uncorr)',sprintf('stat. err (%.0f h)',obj.StudyObject.TimeSec/(60*60)),'Location','southwest');
                    leg.FontSize = 10;
                    title('Integral Spectrum');
                    PrettyFigureFormat;
                    xlim([min(obj.StudyObject.qU-obj.StudyObject.Q) max(obj.StudyObject.qU-obj.StudyObject.Q)]);
                    hold off;
                    %save plot
                    if strcmp(saveplot,'ON')
                        effects_logic = structfun(@(x)strcmp(x,'ON'),obj.SysEffect);
                        fields = fieldnames(obj.SysEffect);
                        mytitle = '';
                        for i=1:numel(fields)
                            if effects_logic(i)==1
                                mytitle = strcat(mytitle,' + ',fields(i));
                            end
                        end
                        if strcmp(mytitle, ' +RF_EL +RF_BF +RF_RX')
                            mytitle{:} = 'Response Function';
                        end
                        fig9999_str = sprintf('./plots/DecomposeCM_%s_%s_%gmolPercm2_%.2ferr.eps',mytitle{:},obj.StudyObject.TD,obj.StudyObject.WGTS_CD_MolPerCm2,obj.WGTS_CD_MolPerCm2_RelErr);
                        publish_figure(fig9999,fig9999_str);
                        fig9991_str = sprintf('./plots/DecomposeCorrM_%s_%s_%gmolPercm2_%.2ferr.eps',mytitle{:},obj.StudyObject.TD,obj.StudyObject.WGTS_CD_MolPerCm2,obj.WGTS_CD_MolPerCm2_RelErr);
                        publish_figure(fig9991,fig9991_str);
                    end
            end
        end
        function TBDISmultisim = GenerateTBDISmultisims(obj,varargin)
            % Returns an nqU-by-nMultisims matrix TBDISmultisim of random vectors chosen
            % from the multivariate normal distribution with mean TBDIS, and covariance CovMat.
            p = inputParser;
            p.addParameter('nMultisims',1000,@(x)isfloat(x) && x>0);
            p.parse(varargin{:});
            nMultisims = p.Results.nMultisims;
            TBDISmultisim = mvnrnd(obj.StudyObject.TBDIS,obj.CovMat,nMultisims);
        end       
        function [fscatn,fscatnE]  =  ComputeELossLookupTable(obj,varargin)
            %Computating of empirical energy loss function and convolutions
            p= inputParser;
            p.addParameter('Debug', 'OFF', @(x)ismember(x,{'ON','OFF'}));
            p.addParameter('E','',@(x)isfloat(x));  % energy for eloss binning
            
            p.parse(varargin{:});
            
            Debug   = p.Results.Debug;
            E       = p.Results.E;
            
            if isempty(E) % default binning if non specified
                maxE = 500;% 9288;
                minE=-maxE; NbinE = (maxE-minE)/0.2;
                E = linspace(minE,maxE,NbinE);
            end
            Estep = E(2) - E(1);
            
            dir_eloss = sprintf([getenv('SamakPath'),'inputs/CovMat/RF/LookupTables/ELoss/']);
            MakeDir(dir_eloss);
            file_eloss = [dir_eloss,sprintf('ELossFunctions_%s_%geV-%geV-qU_%.3geV-Binning_%.0fNIS_%.0fTrials.mat',...
                obj.StudyObject.ELossFlag,min(E),max(E),Estep,obj.StudyObject.NIS,obj.nTrials)];
            
            if exist(file_eloss,'file') && strcmp(obj.RecomputeFlag,'OFF')
                tmp = importdata(file_eloss);
            else 
                tmp.E = 0;
            end
            
            if exist(file_eloss,'file') && all(E==tmp.E) && strcmp(obj.RecomputeFlag,'OFF')
                fprintf(2,'loading energy loss function samples from file \n');
                fscatn = tmp.fscatn;
                fscatnE = tmp.fscatnE;  
            else
            % Expectation Values and uncorrelated errors for ELoss Parameterization
            switch obj.StudyObject.ELossFlag
                case {'Aseev','Abdurashitov'}
                    ELossPar_expected = [obj.StudyObject.is_A1,obj.StudyObject.is_A2,obj.StudyObject.is_w1,...
                        obj.StudyObject.is_w2,obj.StudyObject.is_eps1,obj.StudyObject.is_eps2];
                    ELossErr_expected = [obj.StudyObject.is_A1e,obj.StudyObject.is_A2e, obj.StudyObject.is_w1e,...
                        obj.StudyObject.is_w2e,obj.StudyObject.is_eps1e,obj.StudyObject.is_eps2e];
                case 'KatrinD2'
                    ELossPar_expected = obj.StudyObject.parKatrinD2;
                    ELossErr_expected = obj.StudyObject.errKatrinD2;
                case {'KatrinT2','KatrinT2A20'}
                    ELossPar_expected = obj.StudyObject.parKatrinT2;
                    ELossErr_expected = obj.StudyObject.errKatrinT2;
            end
            
            %% Correlation Matrix of ELoss Parameters
            switch obj.StudyObject.ELossFlag
                case 'Aseev'
                    % correlation unknown -> consider uncertainties uncorrelated
                    ELossCorrMat = ELossErr_expected.*diag(ones(numel(ELossErr_expected),1)).*ELossErr_expected';
                case 'Abdurashitov'
                    % correlation taken from A.Nozik (internal communication)
                    d = importdata([getenv('SamakPath'),'/inputs/ELossFunction/Abdurashitov_ELossFunction_CorrMat.mat']);
                    ELossCorrMat =ELossErr_expected.*d.ELossCorrMat_Average.*ELossErr_expected';
                case 'KatrinD2'
                       % correlation taken from V.Hannen (internal communication): May 19
                      %d=importdata([getenv('SamakPath'),'/inputs/ELossFunction/KATRIND2_ELossFunction_CorrMat_old.mat']);     % old
                      % ELossCorrMat = d.ELoss_CorrMat; % old
                      d = importdata([getenv('SamakPath'),'/inputs/ELossFunction/KATRIND2_ELossFunction_CorrMat_May19.mat']);  % new: May 19
                      ELossCorrMat = ELossErr_expected.*d.ELoss_CorrMat.*ELossErr_expected';
                case {'KatrinT2','KatrinT2A20'}
                       % correlation taken from V.Hannen (internal communication): May 19
                      %d=importdata([getenv('SamakPath'),'/inputs/ELossFunction/KATRIND2_ELossFunction_CorrMat_old.mat']);     % old
                      % ELossCorrMat = d.ELoss_CorrMat; % old
                      ElossFit = 'M19'; % 'M19';
                      switch ElossFit
                          case 'M19'
                              d = importdata([getenv('SamakPath'),'/inputs/ELossFunction/KATRINT2_ELossFunction_CorrMat.mat']);  % new: May 19
                          case 'A20'
                              d = importdata([getenv('SamakPath'),'/inputs/ELossFunction/KATRINT2_ELossFunction_CorrMat_April20.mat']);  % new: April 20
                      end
                      ELossCorrMat = ELossErr_expected.*d.ELoss_CorrMat.*ELossErr_expected';
            end
            
            %% draw samples with multivariate normal distribution
            try
            ELossPar = mvnrnd(ELossPar_expected,ELossCorrMat,obj.nTrials);
            catch
                % WARNING: temporary solution
                fprintf(2,'warning: eloss correlation matrix not positive semi-definite \n');
                ELossPar = mvnrnd(ELossPar_expected,nearestSPD(ELossCorrMat),obj.nTrials);
            end
            fscatn  = cell(obj.StudyObject.NIS, obj.nTrials); % cell of function handles
            fscatnE = zeros(obj.StudyObject.NIS,numel(E),obj.nTrials);
            
            switch obj.StudyObject.ELossFlag
                case {'Aseev','Abdurashitov'}
                    is_epsc = zeros(obj.nTrials,1);
            end
            
            %% calculate energy loss and convolutions
            progressbar('Compute ELossFunction Lookup Table');
            for i=1:1:obj.nTrials
                progressbar(i/obj.nTrials)
                switch Debug
                    case 'ON'
                        cprintf('blue','CovarianceMatrix: ComputeRFCovMatrix: Trial %.0f/%.0f - %s Eloss - %.5f %.5f %.5f %.5f %.5f %.5f \n',...
                            i,obj.nTrials,obj.StudyObject.ELossFlag,ELossPar(i,1),ELossPar(i,2),ELossPar(i,3),ELossPar(i,4),ELossPar(i,5),ELossPar(i,6));
                end
                
                switch obj.StudyObject.ELossFlag
                    case {'Aseev','Abdurashitov'}
                        % calculate epsilon c
                        is_epsc_tmp = obj.StudyObject.ComputeELossEpsilonC(ELossPar(i,1),ELossPar(i,2),ELossPar(i,3),...
                            ELossPar(i,4),ELossPar(i,5),ELossPar(i,6));
                        while isempty(is_epsc_tmp)
                            ELossPar(i,:) = mvnrnd(ELossPar_expected,ELossCorrMat,1);
                            is_epsc_tmp = obj.StudyObject.ComputeELossEpsilonC(ELossPar(i,1),ELossPar(i,2),ELossPar(i,3),...
                                ELossPar(i,4),ELossPar(i,5),ELossPar(i,6));
                        end
                        is_epsc(i) = is_epsc_tmp;
                        
                        % calculate convolutions
                        [fscatn(:,i),fscatnE(:,:,i)] =  obj.StudyObject.ComputeELossFunction('is_A1',ELossPar(i,1),'is_A2',ELossPar(i,2),...
                            'is_w1',ELossPar(i,3),'is_w2',ELossPar(i,4),'is_eps1',ELossPar(i,5),'is_eps2',ELossPar(i,6),...
                            'is_epsc',is_epsc(i),'LoadOrSaveEloss','OFF','E',E);
                    case 'KatrinD2'
                        [fscatn(:,i),fscatnE(:,:,i)] =  obj.StudyObject.ComputeELossFunction('parKatrinD2',ELossPar(i,:),'LoadOrSaveEloss','OFF','E',E);
                    case {'KatrinT2','KatrinT2A20'}
                        [fscatn(:,i),fscatnE(:,:,i)] =  obj.StudyObject.ComputeELossFunction('parKatrinT2',ELossPar(i,:),'LoadOrSaveEloss','OFF','E',E);
                end
            end
            
            %save lookup table
            save(file_eloss,'fscatn','fscatnE','E','-mat');
            
            end
        end
        %         function out   =  ComputeISProbLookupTable(obj,varargin)
        %             % Compute: Inelastic Scattering Probability Covariance Maxtrix
        %             % Varied Parameter:
        %             %   - MACE_Bmax_T_RelErr       : Relative Uncertainty on Bmax settings
        %             %   - WGTS_B_T_RelErr          : Relative Uncertainty on B-WGTS settings
        %             %   - WGTS_CD_MolPerCm2_RelErr : Relative Uncertainty on columns density
        %             %   - ISXsection_RelErr        : Relative Uncertainty on IS Cross Section
        %
        %             %% Variation of Parameters
        %             % B-Fields
        %             if strcmp(obj.SysEffect.RF_BF,'ON')
        %                 Bs_v   = obj.VaryParameter('Parameter',obj.StudyObject.WGTS_B_T, 'rel_Error', obj.WGTS_B_T_RelErr);
        %                 Bmax_v = obj.VaryParameter('Parameter',obj.StudyObject.MACE_Bmax_T,'rel_Error', obj.MACE_Bmax_T_RelErr);
        %             elseif strcmp(obj.SysEffect.RF_BF,'OFF')
        %                 Bs_v   = repmat(obj.StudyObject.WGTS_B_T,obj.nTrials,1);
        %                 Bmax_v = repmat(obj.StudyObject.MACE_Bmax_T,obj.nTrials,1);
        %             end
        %
        %             %Rho, XSection
        %             WGTS_CD_MolPerCm2_v = zeros(obj.nTrials,obj.nRuns);
        %             if strcmp(obj.SysEffect.RF_RX,'ON')
        %                 WGTS_CD_MolPerCm2_v(:,1) = obj.VaryParameter('Parameter',obj.StudyObject.WGTS_CD_MolPerCm2,'rel_Error', obj.WGTS_CD_MolPerCm2_RelErr);
        %                 if strcmp(obj.RunFlag, 'ON')
        %                     for i=1:1:obj.nTrials
        %                         WGTS_CD_MolPerCm2_v(i,2:end) = obj.VaryParameter('Parameter',WGTS_CD_MolPerCm2_v(i,1),'rel_Error', 0.05,'Trials', obj.nRuns-1 );
        %                     end
        %                 end
        %             elseif strcmp(obj.SysEffect.RF_RX,'OFF')
        %                 WGTS_CD_MolPerCm2_v(:,1) = repmat(obj.StudyObject.WGTS_CD_MolPerCm2,obj.nTrials,1);
        %             end
        %
        %             %% Compute Scattering Probabilities
        %             Pis_m         = zeros(obj.StudyObject.NIS+1,obj.nTrials,obj.nRuns);
        %             progressbar('Compute ISProb Lookup Table');
        %             for i=1:1:obj.nTrials
        %                 r = 1;
        %                 progressbar(i/obj.nTrials);
        %                 Pis_m(:,i) = obj.StudyObject.ComputeISProb(...
        %                     'MACE_Bmax_T',Bmax_v(i),...
        %                     'WGTS_B_T',Bs_v(i),...
        %                     'ISXsection',obj.StudyObject.ISXsection,...
        %                     'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2_v(i,r),...
        %                     'saveFile','OFF');
        %             end
        %
        %             NIS=obj.StudyObject.NIS+1;
        %             %% Save IS Probabilites Lookup Table
        %             PisCM = zeros(NIS,NIS,obj.nRuns);
        %             Pis_mean = zeros(NIS, obj.nRuns);
        %             for r=1:1:obj.nRuns
        %                 PisCM(:,:,r)  = cov(Pis_m(:,:,r)');   % Covariance Matrix for first run
        %                 Pis_mean(:,r) = mean(Pis_m(:,:,r),2); % Averaged Probabilites for first run
        %             end
        %
        %             if strcmp(obj.SysEffect.RF_BF,'ON') && strcmp(obj.SysEffect.RF_RX,'ON')
        %                 isp_dir =   [getenv('SamakPath'),sprintf('inputs/CovMat/RF/LookupTables/ISProb/')];
        %             else
        %                 isp_dir =   [getenv('SamakPath'),sprintf('inputs/CovMat/RF/LookupTables/PartVar/ISProbPartVar/')];
        %             end
        %
        %             isp_common = sprintf('ISProb-LookupTable_%u-Trials_%.0fNIS',obj.nTrials,obj.StudyObject.NIS+1);
        %             str_cd = sprintf('_%.5g-molPercm2',obj.StudyObject.WGTS_CD_MolPerCm2);
        %             str_is = sprintf('_IsX%.5gm2',obj.StudyObject.ISXsection);
        %             if strcmp(obj.SysEffect.RF_RX,'ON')
        %                 str_cd = [str_cd,sprintf('-%.3gerr',obj.WGTS_CD_MolPerCm2_RelErr)];
        %                 str_is = [str_is,sprintf('-%.3gerr',obj.ISXsection_RelErr)];
        %             end
        %             str_rx = [str_cd,str_is];
        %
        %             str_bs   = sprintf('_%.3gBs',obj.StudyObject.WGTS_B_T);
        %             str_bmax = sprintf('_%.3gBmax',obj.StudyObject.MACE_Bmax_T);
        %             if strcmp(obj.SysEffect.RF_BF,'ON')
        %                 str_bs    = [str_bs,sprintf('-%.3gerr',obj.WGTS_B_T_RelErr)];
        %                 str_bmax  = [str_bmax,sprintf('-%.3gerr',obj.MACE_Bmax_T_RelErr)];
        %             end
        %             str_bf = [str_bs,str_bmax];
        %
        %             isp_filename = [isp_common,str_rx,str_bf,'.mat'];
        %
        %             if strcmp(obj.SysEffect.RF_BF,'OFF') && strcmp(obj.SysEffect.RF_RX,'OFF')
        %                 isp_dir = [getenv('SamakPath'),sprintf('inputs/WGTSMACE/WGTS_ISProb/')];
        %                 isp_filename = sprintf('IS_%.5g-molPercm2_%.5g-Xsection_%.0f-NIS_%.3g-Bmax_%.3g-Bs.mat',...
        %                     obj.StudyObject.WGTS_CD_MolPerCm2,obj.StudyObject.ISXsection,(obj.StudyObject.NIS+1),obj.StudyObject.MACE_Bmax_T, obj.StudyObject.WGTS_B_T);
        %             end
        %             isp_file = [isp_dir,isp_filename];
        %
        %             MakeDir(isp_dir);
        %             save(isp_file,'Pis_mean','Pis_m','PisCM','Bs_v','Bmax_v','WGTS_CD_MolPerCm2_v','ISXsection_v');
        %             out = Pis_m;
        %         end
        function ISProb = ComputeISProbLookupTable(obj,varargin)
            p=inputParser;
            p.addParameter('MACE_Bmax_T','',@(x)isfloat(x)); % samples
            p.addParameter('WGTS_B_T','',@(x)isfloat(x)); % samples
            p.addParameter('WGTS_CD_MolPerCm2','',@(x)isfloat(x)); % samples
            p.parse(varargin{:});
            MACE_Bmax_T       = p.Results.MACE_Bmax_T;
            WGTS_B_T          = p.Results.WGTS_B_T;
            WGTS_CD_MolPerCm2 = p.Results.WGTS_CD_MolPerCm2;
            
            % binning for isProb E==te-qu
            IsProbBinStep = 2;
            maxEis = 1000;
            Energy = 18575+(-maxEis:IsProbBinStep:maxEis);
            
            % label
            savedir = [getenv('SamakPath'),'inputs/CovMat/RF/LookupTables/ISProb/'];
            MakeDir(savedir);
            
            str_cd = sprintf('ISProb_LookupTable_CD%gmolPercm2',obj.StudyObject.WGTS_CD_MolPerCm2);
            if strcmp(obj.StudyObject.ISCS,'Edep')
                str_is = '_IsXEdep';
            else
                str_is = sprintf('_IsX%.5gm2',obj.StudyObject.ISXsection(18575));
            end
            
            if strcmp(obj.SysEffect.RF_RX,'ON')
                str_cd = [str_cd,sprintf('-%.2gerr',obj.WGTS_CD_MolPerCm2_RelErr)];
                str_is = [str_is,sprintf('-%.3gerr',obj.ISXsection_RelErr)];
            end
            str_rx = [str_cd,str_is];
            
            str_bs   = sprintf('_Bs%.2fT',obj.StudyObject.WGTS_B_T);
            str_bmax = sprintf('_Bmax%.2fT',obj.StudyObject.MACE_Bmax_T);
            if strcmp(obj.SysEffect.RF_BF,'ON')
                str_bs = [str_bs,sprintf('-%.2gerr',obj.WGTS_B_T_RelErr)];
                str_bmax = [str_bmax,sprintf('-%.2gerr',obj.MACE_Bmax_T_RelErr)];
            end
            str_bf = [str_bs,str_bmax];
            
            str_common = sprintf('_maxE%.1feV_Estep%.1feV_%gTrials.mat',...
                max(Energy),IsProbBinStep,obj.nTrials);
            
            is_filename = [savedir,str_rx,str_bf,str_common];
            
            
            if exist(is_filename,'file') 
                d = importdata(is_filename);
                ISProb = d.is_Pv;  
                fprintf('load is prob. lookup table from file %s \n',is_filename)
            else
                ISXsection = obj.StudyObject.ISXsection(Energy);
                % load Grid
                if obj.StudyObject.WGTS_CD_MolPerCm2<(2*1e17)
                    DataSet = 'Knm1';
                else
                    DataSet = 'Knm2';
                end
                
                [RhoDSigma,Theta,ISProb] = Compute_InitISProbMeshGrid('DataSet',DataSet);
                
                ThetaFun = @(bmax,bs)  asin(sqrt(bs./bmax));
                
                is_Pv = squeeze(zeros(obj.StudyObject.NIS+1,...
                    numel(squeeze(ISXsection)),obj.nTrials));
                [X,Y] = meshgrid(RhoDSigma,Theta);
                
                ThetaMax = repmat(ThetaFun(MACE_Bmax_T,WGTS_B_T),[1,numel(Energy)]);
                RDISX = (ISXsection.*WGTS_CD_MolPerCm2);

                for i=1:obj.StudyObject.NIS+1
                    is_Pv(i,:,:) = reshape(interp2(X,Y,squeeze(ISProb(i,:,:))',...
                        RDISX,ThetaMax,...
                        'spline')',1,size(is_Pv,2),size(is_Pv,3));
                end
                
                save(is_filename,'is_Pv','MACE_Bmax_T','WGTS_B_T','ISXsection','WGTS_CD_MolPerCm2','Energy');
                ISProb = is_Pv;  
            end
        end
        function ResponseFunction  =  ComputeRFLookupTable(obj,varargin)
            % Compute nTrials Response Functions
            % Each Response Function Matrix has size (nqUxElossSave)
            % do not save anymore, too large
            % pass to ComputeCM_RF
            
            % Parser
            p = inputParser;
            p.addParameter('Debug','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('RFBinStep',0.04,@(x)isfloat(x)); % eV 0.04
            p.addParameter('ELossBinStep',0.2,@(x)isfloat(x)); % for ELoss Convolution
            p.addParameter('maxE_el',500,@(x)isfloat(x)); % maximal energy of ELoss Convolution %9288;
            p.parse(varargin{:});
            Debug               = p.Results.Debug;
            RFBinStep           = p.Results.RFBinStep;
            ELossBinStep        = p.Results.ELossBinStep;
            maxE_el             = p.Results.maxE_el;
            % -----------------------parser end ------------------------------------------------------

            nRings = numel(obj.StudyObject.MACE_Ba_T);
            
            %% Load Lookup Table: Energy Loss Functions
            % energy loss bninning
            minE_el=-maxE_el;
            NbinE_el = (maxE_el-minE_el)/ELossBinStep;
            E_el = linspace(minE_el,maxE_el,NbinE_el);
            
            if strcmp(obj.SysEffect.RF_EL, 'OFF')
                % load or recompute default eloss
                [~,ElossFunctions_default] = obj.StudyObject.ComputeELossFunction('E',E_el,'LoadOrSaveEloss',obj.RecomputeFlag);
                ElossFunctions= repmat(ElossFunctions_default,1,1,obj.nTrials);
            elseif strcmp(obj.SysEffect.RF_EL, 'ON')
                % load or compute eloss samples (lookup table)
                [~, ElossFunctions] = obj.ComputeELossLookupTable('E',E_el);
            end
            
            %% Variation of Parameters  
             % B-fields
            if strcmp(obj.SysEffect.RF_BF,'ON')
                Bs_v   = obj.VaryParameter('Parameter',obj.StudyObject.WGTS_B_T, 'rel_Error', obj.WGTS_B_T_RelErr);
                Bmax_v = obj.VaryParameter('Parameter',obj.StudyObject.MACE_Bmax_T,'rel_Error', obj.MACE_Bmax_T_RelErr);
                Ba_v   = obj.VaryParameter('Parameter',obj.StudyObject.MACE_Ba_T,'rel_Error',obj.MACE_Ba_T_RelErr);
            elseif strcmp(obj.SysEffect.RF_BF,'OFF')
                Bs_v   = repmat(obj.StudyObject.WGTS_B_T,obj.nTrials,1);
                Bmax_v = repmat(obj.StudyObject.MACE_Bmax_T,obj.nTrials,1);
                Ba_v   = repmat(obj.StudyObject.MACE_Ba_T,obj.nTrials,1);
            end
            
             % column density 
            if strcmp(obj.SysEffect.RF_RX,'ON')
                WGTS_CD_MolPerCm2_v = obj.VaryParameter('Parameter',obj.StudyObject.WGTS_CD_MolPerCm2,'rel_Error', obj.WGTS_CD_MolPerCm2_RelErr);
            elseif strcmp(obj.SysEffect.RF_RX,'OFF')
                WGTS_CD_MolPerCm2_v = repmat(obj.StudyObject.WGTS_CD_MolPerCm2,obj.nTrials,1);
            end

            %% scattering probabilities for all samples
            ISProb = obj.ComputeISProbLookupTable('MACE_Bmax_T',Bmax_v,...
                                       'WGTS_B_T',Bs_v,...
                                       'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2_v);
   
            % calculate response functions with varied parameters
            RFfun = @ComputeRF;
            ResponseFunction = zeros(obj.StudyObject.nTe,obj.StudyObject.nqU,obj.nTrials,nRings);

           % ParObj = copy(repmat(obj.StudyObject,obj.nTrials,1));
            progressbar('Compute RF lookup table');
            nTrials_local = obj.nTrials;
            for i = 1:nTrials_local
                progressbar(i/nTrials_local);
                for ri = 1:nRings
                    for ii = 1:obj.StudyObject.nqU
                        ResponseFunction(:,ii,i,ri) = RFfun(obj.StudyObject,obj.StudyObject.Te,...
                            obj.StudyObject.qU(ii,ri),'pixel',ri,...
                            'ELossRange',maxE_el,'ELossBinStep',ELossBinStep,...
                            'ElossFunctions',ElossFunctions(:,:,i),...
                            'RFBinStep',RFBinStep,...
                            'MACE_Bmax_T',Bmax_v(i),...
                            'MACE_Ba_T',Ba_v(i,:),...
                            'WGTS_B_T',Bs_v(i),...
                            'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2_v(i),...
                            'ISProb',squeeze(ISProb(:,:,i)));
                    end
                end
            end
            
            % do not save all samples, just some infos
            if strcmp(obj.SysEffect.RF_BF,'ON') && strcmp(obj.SysEffect.RF_RX,'ON') && strcmp(obj.SysEffect.RF_EL,'ON')
                rf_path = [getenv('SamakPath'),sprintf('inputs/CovMat/RF/LookupTables/RF/')];
                rf_filename = sprintf('RFInfo_%u-Trials_%u-NIS_%s_%g-molPercm2_%.2gerr_xsection%gerr_%.2fT-Bmax_%.2fT-Bs_BT%gerr_%.0fG-Ba_%.2gerr_%s_elossBinning-%.geV-%.2geVStep_RFbinStep-%.2f.mat',...
                    obj.nTrials, obj.StudyObject.NIS, obj.StudyObject.TD, obj.StudyObject.WGTS_CD_MolPerCm2,...
                    obj.WGTS_CD_MolPerCm2_RelErr,obj.ISXsection_RelErr,...
                    obj.StudyObject.MACE_Bmax_T,obj.StudyObject.WGTS_B_T,obj.WGTS_B_T_RelErr,...
                    mean(floor(obj.StudyObject.MACE_Ba_T*1e4)),obj.MACE_Ba_T_RelErr,obj.StudyObject.ELossFlag,...
                    maxE_el,ELossBinStep,RFBinStep);
                rf_file = strcat(rf_path,rf_filename);
            else
                rf_path = [getenv('SamakPath'),sprintf('inputs/CovMat/RF/LookupTables/PartVar/RFPartVar/')];
                rf_filename = sprintf('RFInfo_%u-Trials_%u-NIS_BField-%s_RhoXsection-%s_Eloss-%s_%s_%g-molPercm2_%.2gerr_xsection%gerr_%.2fT-Bmax_%.2fT-Bs_BT%gerr_%.0fG-Ba_%.2gerr_elossBinning-%.geV-%.2geVStep_RFbinStep-%.2f.mat',...
                    obj.nTrials, obj.StudyObject.NIS,obj.SysEffect.RF_BF,obj.SysEffect.RF_RX,obj.SysEffect.RF_EL,...
                    obj.StudyObject.TD,obj.StudyObject.WGTS_CD_MolPerCm2,...
                    obj.WGTS_CD_MolPerCm2_RelErr,obj.ISXsection_RelErr,...
                    obj.StudyObject.MACE_Bmax_T,obj.StudyObject.WGTS_B_T,obj.WGTS_B_T_RelErr,...
                    mean(floor(obj.StudyObject.MACE_Ba_T*1e4)),obj.MACE_Ba_T_RelErr,...
                    maxE_el,ELossBinStep,RFBinStep);
                if strcmp(obj.SysEffect.RF_EL,'ON')
                    rf_filename = strrep(rf_filename,'.mat',sprintf('_%s.mat',obj.StudyObject.ELossFlag));
                end
                rf_file = strcat(rf_path,rf_filename);
            end
            
            MakeDir(rf_path);

            RFmean = squeeze(mean(ResponseFunction,3));
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                RFstd  = squeeze(std(permute(ResponseFunction,[3 1 2 4])));
            else
                RFstd  = squeeze(std(permute(ResponseFunction,[3 1 2])));
            end
            qU = obj.StudyObject.qU;
            Te = obj.StudyObject.Te;
            TD = obj.StudyObject.TD;    
            ISXsection = obj.StudyObject.ISXsection;
            
            save(rf_file,'RFmean','RFstd','qU','Te','TD',...
                'Bmax_v','Bs_v','ISXsection',...
                'WGTS_CD_MolPerCm2_v','Ba_v',...
                '-mat');
        end
        function ComputeCM_RF(obj,varargin)
            fprintf('--------------------------------------------------------------------------\n')
            cprintf('blue','CovarianceMatrix:ComputeCM_RF: Compute Response Function Covariance Matrix \n')

            if ismember(obj.StudyObject.TD,{'FTpaper','StackCD100all','StackCD100_3hours'}) || str2double(obj.StudyObject.TD)<51410 || contains(obj.StudyObject.TD,'FTpaper')% first tritium binning!
                DefaultRFBinStep = 0.4;
                if strcmp(obj.SysEffect.RF_EL,'ON')
                    DefaultmaxE_el = 2500;%(obj.StudyObject.qUmax-obj.StudyObject.qUmin)*1.6;
                else
                    %  DefaultmaxE_el = 2500
                    DefaultmaxE_el = 9288;
                end
            else
                DefaultRFBinStep = 0.04;
                DefaultmaxE_el = 500;
            end 
            DefaultELossBinStep = 0.2;
            
            % Inputs
            p = inputParser;
            p.addParameter('StatFluct','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('RFBinStep',DefaultRFBinStep,@(x)isfloat(x));
            p.addParameter('ELossBinStep',DefaultELossBinStep,@(x)isfloat(x)); % for ELoss Convolution
            p.addParameter('maxE_el',DefaultmaxE_el,@(x)isfloat(x)); % maximal energy of ELoss Convolution %9288;
            p.parse(varargin{:});
            StatFluct  = p.Results.StatFluct;
            RFBinStep           = p.Results.RFBinStep;
            ELossBinStep        = p.Results.ELossBinStep;
            maxE_el             = p.Results.maxE_el;
            % -----------------------parser end ------------------
            
            %Labeling
            if strcmp(obj.SysEffect.RF_BF,'ON') && strcmp(obj.SysEffect.RF_RX,'ON') && strcmp(obj.SysEffect.RF_EL,'ON')
                obj.MultiCovMat.CM_RF = obj.CovMat;
                covmat_path = [getenv('SamakPath'),sprintf('inputs/CovMat/RF/CM/')];
            elseif strcmp(obj.SysEffect.RF_BF,'OFF') || strcmp(obj.SysEffect.RF_RX,'OFF') ||  strcmp(obj.SysEffect.RF_EL,'OFF')
                covmat_path = [getenv('SamakPath'),sprintf('inputs/CovMat/RF/CMPartVar/')];
            end
            
            obj.GetTDlabel;
            
            str_cd = sprintf('WGTSMACE_CovMat_%s_CD%gmolPercm2',obj.TDlabel,obj.StudyObject.WGTS_CD_MolPerCm2);
            if strcmp(obj.StudyObject.ISCS,'Edep')
                str_is = '_IsXEdep';
            else
                str_is = sprintf('_IsX%.5gm2',obj.StudyObject.ISXsection(18575));
            end
            if strcmp(obj.SysEffect.RF_RX,'ON')
                str_cd = [str_cd,sprintf('-%.2gerr',obj.WGTS_CD_MolPerCm2_RelErr)];
                str_is = [str_is,sprintf('-%.3gerr',obj.ISXsection_RelErr)];
            end
            str_rx = [str_cd,str_is];
            
            str_bs   = sprintf('_Bs%.2fT',obj.StudyObject.WGTS_B_T);
            str_bmax = sprintf('_Bmax%.2fT',obj.StudyObject.MACE_Bmax_T);
            str_ba   = sprintf('_Ba%.2fT',round(obj.StudyObject.MACE_Ba_T*1e4));
            if strcmp(obj.SysEffect.RF_BF,'ON')
                str_bs = [str_bs,sprintf('-%.2gerr',obj.WGTS_B_T_RelErr)];
                str_bmax = [str_bmax,sprintf('-%.2gerr',obj.MACE_Bmax_T_RelErr)];
                str_ba = [str_ba,sprintf('-%.2gerr',obj.MACE_Ba_T_RelErr)];
            end
            str_bf = [str_bs,str_bmax,str_ba];
            
            str_el = sprintf('_%s',obj.StudyObject.ELossFlag);
            if strcmp(obj.SysEffect.RF_EL,'OFF')
                str_el = [str_el,'-ELossOFF'];
            end
            
            str_common = sprintf('_elossBinning-%.geV-%.2geVStep_RFbinStep-%.2feV_%gTrials.mat',...
                maxE_el,ELossBinStep,RFBinStep,obj.nTrials);
            
            covmat_filename = [str_rx,str_bf,str_el,str_common];
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                % add ring information for ringwise covariance matrices
                covmat_filename = strrep(covmat_filename,'.mat',sprintf('_Ring%s.mat',obj.StudyObject.FPD_RingMerge));
            end
            obj.CovMatFile = [covmat_path,covmat_filename];
            %% Try reading Covariance Matrix from file
            if exist(obj.CovMatFile,'file') && strcmp(obj.RecomputeFlag,'OFF')
                fprintf('CovarianceMatrix:ComputeCM_RF: Loading WGTS CM from File \n')
                if strcmp(obj.SysEffect.RF_BF,'ON') && strcmp(obj.SysEffect.RF_RX,'ON') && strcmp(obj.SysEffect.RF_EL,'ON')
                    obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','RF');
                    obj.MultiCovMat.CM_RF = obj.CovMat; %in case normalization was changed
                elseif strcmp(obj.SysEffect.RF_BF,'ON') && strcmp(obj.SysEffect.RF_RX,'OFF') && strcmp(obj.SysEffect.RF_EL,'OFF')
                    obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','RF_BF');
                    obj.MultiCovMat.CM_RF_BF = obj.CovMat;
                elseif strcmp(obj.SysEffect.RF_BF,'OFF') && strcmp(obj.SysEffect.RF_RX,'ON') && strcmp(obj.SysEffect.RF_EL,'OFF')
                    obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','RF_RX');
                    obj.MultiCovMat.CM_RF_RX = obj.CovMat;
                elseif strcmp(obj.SysEffect.RF_BF,'OFF') && strcmp(obj.SysEffect.RF_RX,'OFF') && strcmp(obj.SysEffect.RF_EL,'ON')
                    obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','RF_EL');
                    obj.MultiCovMat.CM_RF_EL = obj.CovMat;
                else % when 2 are switched on
                    obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','other');
                end
                return
            end
            %% When CovMat files doesnt exist: compute it
            recompRF_i = obj.StudyObject.recomputeRF;
            RFBinStep_i = obj.StudyObject.RFBinStep;
            Sync_i      =  obj.StudyObject.SynchrotronFlag; 
            AngularTFFlag_i = obj.StudyObject.AngularTFFlag;
            obj.StudyObject.RFBinStep = DefaultRFBinStep;
            obj.StudyObject.SynchrotronFlag = 'OFF';
            obj.StudyObject.AngularTFFlag = 'OFF';
            obj.StudyObject.recomputeRF   = 'OFF';
            
            %Init
            nRings = numel(obj.StudyObject.MACE_Ba_T);
            TBDIS_V  = zeros(obj.StudyObject.nqU,obj.nTrials,nRings);  
      
            % Loop Computing TBDIS Matrix
            RF_Samples = obj.ComputeRFLookupTable('maxE_el',maxE_el,'ELossBinStep',ELossBinStep,'RFBinStep',RFBinStep);  %
            progressbar('ComputeTBDISCovarianceMatrix');
            RF_Samples(isnan(RF_Samples))=0; %easy fix, work in progress!
             
            for i=1:1:obj.nTrials
                progressbar(i/obj.nTrials);
                    obj.StudyObject.RF = squeeze(RF_Samples(:,:,i,:));
                    ComputeTBDDS(obj.StudyObject);
                    ComputeTBDIS(obj.StudyObject);
                    switch StatFluct
                        case 'ON'
                            AddStatFluctTBDIS(obj.StudyObject);
                    end
                    TBDIS_V(:,i,:) = obj.StudyObject.TBDIS;
            end % end nTrials
            
            % Compute Averaged Spectrum
            TBDIS_av = squeeze(mean(TBDIS_V,2));
            
            % reshape
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                TBDIS_V = permute(TBDIS_V,[1 3 2]);
                TBDIS_V = reshape(TBDIS_V,[obj.StudyObject.nqU*nRings,obj.nTrials]);
            end
            
             % Compute Covariance Matrix 
            obj.CovMat = cov(TBDIS_V');
            
            % Save
            if  strcmp(obj.SysEffect.RF_BF,'ON') && strcmp(obj.SysEffect.RF_RX,'OFF') && strcmp(obj.SysEffect.RF_EL,'OFF')
                obj.MultiCovMat.CM_RF_BF = obj.CovMat;
            elseif strcmp(obj.SysEffect.RF_BF,'OFF') && strcmp(obj.SysEffect.RF_RX,'ON') && strcmp(obj.SysEffect.RF_EL,'OFF')
                obj.MultiCovMat.CM_RF_RX = obj.CovMat;
            elseif strcmp(obj.SysEffect.RF_BF,'OFF') && strcmp(obj.SysEffect.RF_RX,'OFF') && strcmp(obj.SysEffect.RF_EL,'ON')
                obj.MultiCovMat.CM_RF_EL = obj.CovMat;
            elseif strcmp(obj.SysEffect.RF_BF,'ON') && strcmp(obj.SysEffect.RF_RX,'ON') && strcmp(obj.SysEffect.RF_EL,'ON')
                obj.MultiCovMat.CM_RF = obj.CovMat;
            end
            
            MakeDir(covmat_path);
            save(obj.CovMatFile,'obj','TBDIS_av','TBDIS_V','-mat');
            
            %Compute and save Fractional CM
            obj.ComputeFracCM('Mode','CM2Frac');
            
            % Compute Decomposed CM
            [~, obj.CovMatFracShape] = obj.DecomposeCM('CovMatFrac',obj.CovMatFrac,'exclDataStart',1);
            
            if  strcmp(obj.SysEffect.RF_BF,'ON') && strcmp(obj.SysEffect.RF_RX,'OFF') && strcmp(obj.SysEffect.RF_EL,'OFF')
                obj.MultiCovMatFrac.CM_RF_BF = obj.CovMatFrac;
                obj.MultiCovMatFracShape.CM_RF_BF = obj.CovMatFracShape;
                obj.MultiCovMatFracNorm.CM_RF_BF = obj.CovMatFracNorm;
            elseif strcmp(obj.SysEffect.RF_BF,'OFF') && strcmp(obj.SysEffect.RF_RX,'ON') && strcmp(obj.SysEffect.RF_EL,'OFF')
                obj.MultiCovMatFrac.CM_RF_RX = obj.CovMatFrac;
                obj.MultiCovMatFracShape.CM_RF_RX = obj.CovMatFracShape;
                obj.MultiCovMatFracNorm.CM_RF_RX  =  obj.CovMatFracNorm;
            elseif strcmp(obj.SysEffect.RF_BF,'OFF') && strcmp(obj.SysEffect.RF_RX,'OFF') && strcmp(obj.SysEffect.RF_EL,'ON')
                obj.MultiCovMatFrac.CM_RF_EL = obj.CovMatFrac;
                obj.MultiCovMatFracShape.CM_RF_EL = obj.CovMatFracShape;
                obj.MultiCovMatFracNorm.CM_RF_EL  =  obj.CovMatFracNorm;
            elseif strcmp(obj.SysEffect.RF_BF,'ON') && strcmp(obj.SysEffect.RF_RX,'ON') && strcmp(obj.SysEffect.RF_EL,'ON')
                obj.MultiCovMatFrac.CM_RF = obj.CovMatFrac;
                obj.MultiCovMatFracShape.CM_RF = obj.CovMatFracShape;
                obj.MultiCovMatFracNorm.CM_RF  =  obj.CovMatFracNorm;
            end
            
            save(obj.CovMatFile, 'obj','-mat','-append');
            
            obj.StudyObject.ElossFunctions  = []; % reset e-loss functions
            obj.StudyObject.RFBinStep       =  RFBinStep_i ;
            obj.StudyObject.SynchrotronFlag =  Sync_i ; 
            obj.StudyObject.recomputeRF     = recompRF_i;
            obj.StudyObject.AngularTFFlag   = AngularTFFlag_i;
        end
        function ComputeCM_TASR(obj,varargin)
            fprintf('--------------------------------------------------------------------------\n')
            cprintf('blue','CovarianceMatrix:ComputeCM_TASR: Compute Subrun Activity Fluctuation Covariance Matrix (TASR)           \n')
            % Inputs
            p = inputParser;
            p.parse(varargin{:});
            
            %Labeling
            covmat_path = [getenv('SamakPath'),sprintf('inputs/CovMat/TASR/CM/')];
            covmat_filename = sprintf('TASR_CovMat_%s_%.4gRelErr.mat',...
                obj.StudyObject.TD,max(max(obj.WGTS_TASR_RelErr)));
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                % add ring information for ringwise covariance matrices
                covmat_filename = strrep(covmat_filename,'.mat',sprintf('_Ring%s.mat',obj.StudyObject.FPD_RingMerge));
            end
            obj.CovMatFile = strcat(covmat_path,covmat_filename);
            
            %Check if CM is already compute 
            if exist(obj.CovMatFile,'file')==2 && strcmp(obj.RecomputeFlag,'OFF')
                fprintf(2,'CovarianceMatrix:ComputeCM_TASR: Loading TASR CM from File \n')
                obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','TASR');
                obj.MultiCovMat.CM_TASR = obj.CovMat;
               return
            end
            
            nRings = numel(obj.StudyObject.MACE_Ba_T);

            % Compute Differential / Integral Spectra: only once!
            obj.StudyObject.BKG_RateSec_i = 1e-9;
            ComputeTBDDS(obj.StudyObject,'B_bias',0);
            ComputeTBDIS(obj.StudyObject);
            TBDIS = obj.StudyObject.TBDIS;
            
            % Compute Fractional Covariance Matrix
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                 obj.CovMatFrac      = repmat(obj.WGTS_TASR_RelErr,[nRings,nRings]);
                 obj.MultiCovMatFrac.CM_TASR = obj.CovMatFrac;
            else
                 obj.CovMatFrac = obj.WGTS_TASR_RelErr;
                 obj.MultiCovMatFrac.CM_TASR = obj.CovMatFrac;
            end
            
            % Save
            MakeDir(covmat_path);
            save(obj.CovMatFile,'obj','TBDIS','-mat');
            
            % Get normalized covariance matrix
            obj.ComputeFracCM('Mode','Frac2CM')
            obj.MultiCovMat.CM_TASR = obj.CovMat;
            
            % Compute Shape-only CM: some problems with KNM2: fluctuations too small
            try
            [~, obj.CovMatFracShape] = obj.DecomposeCM('CovMatFrac',obj.CovMatFrac,'exclDataStart',1);
            catch
                fprintf('Decomposition doesnt work - TASR covariance matrix probably too small (Knm2) \n')
                fprintf('Temporary fix: take regular covariance matrix instead - no shape only \n')
                obj.CovMatFracShape = obj.CovMatFrac;
            end
            obj.MultiCovMatFracShape.CM_TASR = obj.CovMatFracShape;
            
            % Save again
            save(obj.CovMatFile,'obj','-mat','-append');
        end
        function ComputeCM_FPDeff(obj,varargin)
            % systematic uncertainty on FPD efficiency
            % bin-to-bin uncorrelated
            fprintf('--------------------------------------------------------------------------\n')
            cprintf('blue','CovarianceMatrix:ComputeCM_FPDeff: Compute FPD efficiency Covariance Matrix (FPDeff) \n')
            
            %Labeling
            covmat_path = [getenv('SamakPath'),sprintf('inputs/CovMat/FPDeff/CM/')];
            covmat_filename = sprintf('FPDeff_CovMat_%s_%.5g-RelErr_%u-Trials.mat',...
                obj.StudyObject.TD,obj.FPDeff_RelErr,obj.nTrials);
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                % add ring information for ringwise covariance matrices
                covmat_filename = strrep(covmat_filename,'.mat',sprintf('_Ring%s.mat',obj.StudyObject.FPD_RingMerge));
            end 
            obj.CovMatFile = strcat(covmat_path,covmat_filename);
            
            %Check if CM is already computed
            if exist(obj.CovMatFile,'file')==2 && strcmp(obj.RecomputeFlag,'OFF')
                fprintf(2,'CovarianceMatrix:ComputeCM_FPDeff: Loading FPDeff CM from File \n')
                obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','FPDeff');
                obj.MultiCovMat.CM_FPDeff = obj.CovMat;
                return
            end
            
            % Compute Differential / Integral Spectra: only once!
            obj.StudyObject.BKG_RateSec_i = 1e-9;
            obj.StudyObject.ComputeTBDDS('B_bias',0);
            obj.StudyObject.ComputeTBDIS;
            TBDIS_V = obj.StudyObject.TBDIS;
            
            % Compute Covariance Matrix
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                dim = [obj.StudyObject.nqU*obj.StudyObject.nRings,1];
                TBDIS_V = reshape(TBDIS_V,dim);
            else
                dim = [obj.StudyObject.nqU,1];
            end
            FPDeff_RelErr = repmat(obj.FPDeff_RelErr,dim);
            obj.CovMat = TBDIS_V.*diag(FPDeff_RelErr.^2).*TBDIS_V';
            obj.MultiCovMat.CM_FPDeff = obj.CovMat;
            
            %Compute Fractional CM
            obj.CovMatFrac = diag(FPDeff_RelErr.^2);
            obj.MultiCovMatFrac.CM_FPDeff = obj.CovMatFrac;
            
            % Compute Shape-only CM
            obj.CovMatFracShape = obj.CovMatFrac;
            obj.MultiCovMatFracShape.CM_FPDeff =  obj.CovMatFrac;
            
            %save
            MakeDir(covmat_path);
            save(obj.CovMatFile,'obj','TBDIS_V','-mat');
        end
        function ComputeCM_BM1S(obj,varargin)
            %
            % Compute the Shape Error uncertainty of backgrounds
            % by taking the difference between a nominal flat and
            % Case A) A background with a slope
            % Case B) A real background measurement
            %
            fprintf('--------------------------------------------------------------------------\n')
            cprintf('blue','CovarianceMatrix:ComputeCM_BM1S: Compute Background Shape Covariance Matrix (BM1S) - Model 1          \n')
            % Inputs
            p = inputParser;
            p.addParameter('plot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Model','D',@(x)ismember(x,{'S','D'}));% Slope or Data
            p.addParameter('qUStartIndex',2); % For Display Only
            p.parse(varargin{:});
            plot          = p.Results.plot;
            Model         = p.Results.Model;
            qUStartIndex  = p.Results.qUStartIndex;
            
            %Labeling
            covmat_path = [getenv('SamakPath'),sprintf('inputs/CovMat/BM1S/CM/')];
            MakeDir(covmat_path);
            covmat_filename = sprintf('BM1S_CovMat_%s_%u-Runs.mat',...
                obj.StudyObject.TD,obj.nRuns);
            obj.CovMatFile = strcat(covmat_path,covmat_filename);
            
            %Check if CM is already computed
            if exist(obj.CovMatFile,'file')==2
%            if exist(obj.CovMatFile,'file')==2 && strcmp(obj.RecomputeFlag,'OFF')
                fprintf(2,'CovarianceMatrix:ComputeCM_BM1S: Loading BM1S CM from File \n')
                obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','BM1S');
                obj.MultiCovMat.CM_BM1S = obj.CovMat;
                return
            end
            
            %Init
            ComputeTBDDS(obj.StudyObject);
            ComputeTBDIS(obj.StudyObject);
            TBDIS_V    = obj.StudyObject.TBDIS;

            % Model in simulation/fit
            Bfit       = obj.nRuns.*(obj.BM1_RatePerSec*obj.StudyObject.qUfrac.*obj.StudyObject.TimeSec);
            
            switch Model
                case 'S'
            % Alternative Model A
            BM1IS_V    = obj.StudyObject.qU./(obj.StudyObject.qUmax-obj.BM1_RefSlopeqU)*obj.BM1_RedFactor.*Bfit + ...
                Bfit*(1-obj.BM1_RedFactor*obj.StudyObject.qUmax/(obj.StudyObject.qUmax-obj.BM1_RefSlopeqU));
                case 'D'
            % Alternative Model B - Background Measurement on May 15 2018 -
            % Data provided by Anna Pollithy 
            qu=[16575 17175 17575 18175 obj.StudyObject.qUmax ];
            b=[3.325091575091577711e-01 3.433150183150172063e-01 3.396520146520150552e-01 3.427655677655671029e-01 3.364635822262940779e-01];
            be=[3.901894190547565075e-03 3.964789038606086023e-03 3.943581135555766955e-03 3.961615090901693219e-03 3.925926497675827354e-03];
            BM1IS_V = interp1(qu,b,obj.StudyObject.qU);
            BM1IS_V = BM1IS_V./BM1IS_V(obj.StudyObject.nqU).*obj.nRuns.*(obj.BM1_RatePerSec*obj.StudyObject.qUfrac.*obj.StudyObject.TimeSec);
            end
            
            % Compute Covariance Matrix
            obj.CovMat = bsxfun(@times,ones(obj.StudyObject.nqU,obj.StudyObject.nqU),(BM1IS_V-Bfit)); 
            obj.CovMat = bsxfun(@times,obj.CovMat,(BM1IS_V-Bfit)'); 
            obj.MultiCovMat.CM_BM1S  = obj.CovMat;
            
            % Plot 
            switch plot
                case 'ON'
                    figure
                    plt1= Plot(obj.StudyObject.qU-obj.StudyObject.Q,Bfit./(obj.StudyObject.TimeSec*obj.StudyObject.qUfrac),...
                        obj.StudyObject.qU-obj.StudyObject.Q,BM1IS_V./(obj.StudyObject.TimeSec*obj.StudyObject.qUfrac));
                    plt1.LineWidth = 2;
                    plt1.LineStyle = {'-','--'};
                    plt1.Markers   = {'s','d'};
                    plt1.MarkerSpacing = [1,1];
                    pltstr         = sprintf('Background Uncertainty Model %s - %s - %g sec',Model,obj.StudyObject.TD,obj.StudyObject.TimeSec);
                    plt1.Title     = pltstr; % plot title
                    plt1.XLabel    = 'E-qU (V)'; % xlabel
                    plt1.YLabel    = 'Rate Per Second'; %ylabel
                    plt1.YScale    = 'lin'; % 'linear' or 'log'
                    plt1.YLim      = [min(BM1IS_V./(obj.StudyObject.TimeSec*obj.StudyObject.qUfrac))*0.95 max(BM1IS_V./(obj.StudyObject.TimeSec*obj.StudyObject.qUfrac))*1.05];
                    plt1.XScale    = 'lin'; % 'linear' or 'log'
                    plt1.FontSize  = 16;
                    plt1.Legend    = {'Nominal (Flat)','Alternative Model'};
                    plt1.LegendLoc = 'SouthEast';
                    plttitle       = sprintf('BackgroundCovarianceMatrixBM1S_Model%s_1.png',Model);
                    plt1.export(plttitle);
                    figure
                    plt2 = Plot(obj.StudyObject.qU-obj.StudyObject.Q,sqrt(diag(obj.MultiCovMat.CM_BM1S)),...
                        obj.StudyObject.qU-obj.StudyObject.Q,sqrt(TBDIS_V));
                    plt2.LineWidth = 2;
                    plt2.LineStyle = {'-','--'};
                    plt2.Markers   = {'s','d'};
                    pltstr         = sprintf('Background Uncertainty Model %s - %s - %g secs',Model,obj.StudyObject.TD,obj.StudyObject.TimeSec);
                    plt2.Title     = pltstr; % plot title
                    plt2.XLabel    = 'E-qU (V)'; % xlabel
                    plt2.YLabel    = 'Count Difference'; %ylabel
                    plt2.YScale    = 'log'; % 'linear' or 'log'
                    plt2.XScale    = 'lin'; % 'linear' or 'log'
                    plt2.Legend    = {'Syst. Error (background)','Stat. Error (sim data)'};
                    %plt2.YLim      =  [min(sqrt(TBDIS_V))/10 max(sqrt(TBDIS_V))];
                    plt2.FontSize  = 16;                    
                    plttitle       = sprintf('BackgroundCovarianceMatrixBM1S_Model%s_2.png',Model);
                    plt2.export(plttitle);
            end

            % Save         
            save(obj.CovMatFile,'obj','Model','BM1IS_V','Bfit','TBDIS_V','-mat');
            
            %Compute and save Fractional CM
            obj.ComputeFracCM('Mode','CM2Frac');
            obj.MultiCovMatFrac.CM_BM1S = obj.CovMatFrac;
            
            %Compute Decomposed CM
            obj.DecomposeCM('Option','Frac');
            obj.MultiCovMatFracShape.CM_BM1S = obj.CovMatFracShape;
            obj.MultiCovMatFracNorm.CM_BM1S = obj.CovMatFracNorm;
            
            save(obj.CovMatFile, 'obj','-append');
        end  
        
        % LISA's VERSION FOR First Tritium 
        
%         function ComputeCM_Background(obj,varargin)
%             % DataDriven way to compute Background
%             % CovarianceMatrix
%             p = inputParser;
%             p.addParameter('SanityPlots','ON');
%             p.parse(varargin{:});
%             SanityPlots = p.Results.SanityPlots;
% 
%             cm_path = sprintf('../../inputs/CovMat/Background/CM/');
%             cm_name = sprintf('BackgroundCM_%s.mat',obj.StudyObject.TD);
%             obj.CovMatFile = [cm_path,cm_name];
%             if exist(obj.CovMatFile,'file')==2
%             
%             end 
%             % Init
%             par = zeros(2,obj.nTrials);
%             err = zeros(2,obj.nTrials);
%             chi2min = zeros(obj.nTrials,1);
%             BkgPlot_Data = zeros(6,obj.nTrials);              % Store for sanity plot
%             Bkg_Fit = zeros(obj.StudyObject.nqU,obj.nTrials); % Store for covmat computation
%             
%             BKGIndex = min(find(obj.StudyObject.qU>=obj.StudyObject.Q_i)); 
%             BKGnqU = numel(obj.StudyObject.qU(obj.StudyObject.qU>=obj.StudyObject.Q_i));
%             
%             % Fit Background Rate (in counts per second)
%             BKG_Asimov = obj.StudyObject.BKG_RateSec;%.*obj.StudyObject.TimeSec.*obj.StudyObject.qUfrac(BKGIndex:end); % in Counts  
%             for i=1:obj.nTrials  
%                 BKG_RateErr = sqrt(BKG_Asimov./(obj.StudyObject.TimeSec.*obj.StudyObject.qUfrac(BKGIndex:end)));
%                 BKG = BKG_Asimov + BKG_RateErr.*randn(BKGnqU,1);
%                 Data = [obj.StudyObject.qU(BKGIndex:end),BKG, BKG_RateErr];
% %               BKG_i = wmean(BKG,sqrt(BKG)); % WRONG WARNING - weight should be 1/sigma2
%                 BKG_i = wmean(BKG,1./BKG);  
%                 Slope_i = 0;
%                 parInit = [BKG_i, Slope_i];
%                 % Minuit Arguments
%                 tmparg = sprintf(['set pri -10;'...
%                     'migrad minos'],'');
%                 % Minuit Input: Init Fit Parameter, Data, Covariance Matrix
%                 Args = {parInit, Data, '-c', tmparg};
%                 [par(:,i), err(:,i),chi2min(i), ~] = fminuit('chi2BKG',Args{:});
%                 BkgPlot_Data(:,i) = Data(:,2); %for plot in cps
% 
%                 Bkg_Fit(:,i) =  ComputeBkgSlope(par(:,i),obj.StudyObject.qU);
%             end
%             TBDIS_V = Bkg_Fit(:,chi2min<3.8);
%             BkgPlot_Data = BkgPlot_Data(:,chi2min<3.8);
%             obj.CovMat = cov(TBDIS_V');
%             obj.CovMatFrac = obj.CovMat./BKG_Asimov^2;
%             
%             save(obj.CovMatFile,'obj','TBDIS_V')
%             
%             if strcmp(SanityPlots,'ON') 
%                 figure(25);
%                 subplot(2,2,1);
%                 nhist(par(2,chi2min<3.8));
%                 xlabel('Slope: cps/eV');
%                 subplot(2,2,2);
%                 nhist(err(2,chi2min<3.8));
%                 xlabel('Slope Err: cps/eV');
%                 subplot(2,2,3);
%                 nhist(par(1,chi2min<3.8));
%                 xlabel('Flat Background Rate: cps');
%                 subplot(2,2,4);
%                 nhist(err(1,chi2min<3.8));
%                 xlabel('Flat Background Rate Err: cps');
%                 
%                 figure(24) % plot Rate
%                 pfit  = plot(obj.StudyObject.qU-18575,1e3*TBDIS_V(:,1:100)); hold on
%                 BkgPlot_DataErr = sqrt(BkgPlot_Data./(obj.StudyObject.TimeSec.*obj.StudyObject.qUfrac(BKGIndex:end)));
%                 for i=1:100
%                 pdata = errorbar(obj.StudyObject.qU(BKGIndex:end)-18575,1e3.*BkgPlot_Data(:,i),BkgPlot_DataErr(:,i),'x');
%                 end
%                 PrettyFigureFormat;
%                 xlabel('retarding potential - 18575 (eV)');
%                 ylabel('background rate(mcps)');
%                 hold off;
%                 
%                 figure(23);
%                 imagesc(obj.CovMatFrac);
%                 colorbar;
%                 PrettyFigureFormat;
%                 title('Fractional Background Covariance Matrix')     
%             end
%
%         end

function ComputeCM_Background(obj,varargin)
    % Data-Driven method to compute Background Shape
    % Uncertainty Covariance Matrix
    % Include Non-Poissonian fluctuations (scaling factor)
    % Lisa Schlueter, Thierry Lasserre, 10/3/2019
    p = inputParser;
    p.addParameter('Display','ON',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('Chi2CutOff',1e9);
    p.addParameter('qUStartIndex',11);             % here, for Display Only
    p.addParameter('MaxSlopeCpsPereV',0.15E-04);   % qU-Slope constraint - Absolute Slope Value
    p.addParameter('BkgRange',-5,@(x)isfloat(x));  % defines, which points are used to constrain slope. In eV with respect to endpoint
    p.addParameter('RingCorrCoeff',0,@(x)isfloat(x));  % ring to ring correlation, can also be a matrix
    p.addParameter('ScalingOpt',2,@(x)isfloat(x));
    p.addParameter('Mode','SlopeFit',@(x)ismember(x,{'SlopeFit','Gauss'}));
    
    p.parse(varargin{:});
    Display          = p.Results.Display;
    Chi2CutOff       = p.Results.Chi2CutOff;
    qUStartIndex     = p.Results.qUStartIndex;
    MaxSlopeCpsPereV = p.Results.MaxSlopeCpsPereV;
    BkgRange         = p.Results.BkgRange;
    RingCorrCoeff    = p.Results.RingCorrCoeff;
    ScalingOpt       = p.Results.ScalingOpt;
    Mode             = p.Results.Mode;
    fprintf(       '--------------------------------------------------------------------------   \n');
    cprintf('blue','CovarianceMatrix:ComputeCM_Background: Compute Background Covariance Matrix  \n');
    
    MaxSlopeCpsPereV_i = MaxSlopeCpsPereV;
    
    % Number of Trials - Hardcoded
    TrialSave  = obj.nTrials; obj.nTrials = 50000; % BASELINE, termine after obj.nTrials
    
    % Covariance Matrix File
    cm_path        = [getenv('SamakPath'),sprintf('inputs/CovMat/Background/CM/')];
    MakeDir(cm_path);
    if ~strcmp(Mode,'SlopeFit')
        extraStr = sprintf('_%s',Mode);
    else
        extraStr = '';
    end  
    cm_name        = sprintf('BackgroundCM_%s_%.0fmcps_Constrain%.3gCpsPerEv_%.0feVBkgRange_%.0fTrials%s.mat',...
        strrep(obj.StudyObject.TD,'_E018573.73eV2',''),...
        sum(obj.StudyObject.BKG_RateSec)*1e3,sum(MaxSlopeCpsPereV),BkgRange,obj.nTrials,extraStr);
    % add ring information for ringwise covariance matrices
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        if numel(MaxSlopeCpsPereV)~=1 % no scaling, ring by ring input
            ScalingOpt = 99;
        end
        cm_name = strrep(cm_name,'.mat',sprintf('_Ring%s_RingCorrCoeff%.2f_%.0f.mat',...
            obj.StudyObject.FPD_RingMerge,norm(RingCorrCoeff),ScalingOpt));
    end
    obj.CovMatFile = [cm_path,cm_name];
    
    if exist(obj.CovMatFile,'file')==2 && strcmp(obj.RecomputeFlag,'OFF')
        fprintf(2,'CovarianceMatrix::ComputeCM_Background:ReadCMFile: Loading CovMatrix File: %s \n',obj.CovMatFile)
        tmpCov = importdata(obj.CovMatFile);
        if round(sum(tmpCov.BKG_Asimov),3)==round(sum(obj.StudyObject.BKG_RateSec),3)
            obj.CovMat          = tmpCov.obj.CovMat;
            obj.CovMatFrac      = tmpCov.obj.CovMatFrac;
            obj.CovMatFracShape = tmpCov.obj.CovMatFracShape;
            obj.MultiCovMat.CM_BkgShape          = tmpCov.obj.CovMat;
            obj.MultiCovMatFrac.CM_BkgShape      = tmpCov.obj.CovMatFrac;
            obj.MultiCovMatFracShape.CM_BkgShape = tmpCov.obj.CovMatFracShape;
            if strcmp(Display,'OFF')
                return;
            end
            load(obj.CovMatFile);
            calclog = 1; % loaded from file
        else
            calclog = 0; % not loaded from file
            fprintf(2,'CovarianceMatrix:ComputeCM_Background: Recalculate BKG CM for new BKG level \n')
        end
    else
        calclog = 0; % not loaded from file
    end
    
    if calclog==0 % calculate cov mat
        % Check/Prepare for Ring Segmentation
        if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
            nRingLoop = obj.StudyObject.nRings;
        else
            nRingLoop = 1;
        end
        
        TBDIS_V           = zeros(nRingLoop,obj.StudyObject.nqU,obj.nTrials);%cell(nRingLoop,1);      % will be later converted into matrix
        nGoodTrials       = zeros(nRingLoop,1); % number of successfull trials
        BKG_Asimov        = obj.StudyObject.BKG_RateSec;
        
        if nRingLoop>1 && numel(MaxSlopeCpsPereV)== 1% for more than 1 ring
            % normalize MaxSlopeCpsPereV to correct statistics
            switch ScalingOpt   
                case 1
                    MaxSlopeCpsPereV = MaxSlopeCpsPereV.*(BKG_Asimov./sum(BKG_Asimov));
                case 2
                    MaxSlopeCpsPereV = MaxSlopeCpsPereV.*sqrt(BKG_Asimov./sum(BKG_Asimov));
                case 3
                    MaxSlopeCpsPereV = repmat(MaxSlopeCpsPereV,nRingLoop,1);
            end
        end
        
        Slopes      = zeros(nRingLoop,obj.nTrials); % background fit slope
        SlopesExcl  = ones(nRingLoop,obj.nTrials); % 1== included, 0 == excluded
        
        % randomize background
        BKGIndex    = find(obj.StudyObject.qU>=(obj.StudyObject.Q_i+BkgRange),1); % Start 5eV below E0
        BKGnqU      = numel(obj.StudyObject.qU(obj.StudyObject.qU>=(obj.StudyObject.Q_i+BkgRange)));
        BKGnqU      = BKGnqU./nRingLoop;
        BKG_RateErr = obj.NonPoissonScaleFactor.* sqrt(BKG_Asimov./(obj.StudyObject.TimeSec.*obj.StudyObject.qUfrac(BKGIndex:end,:)));
        if nRingLoop>1 % for more than 1 ring
            BKG_RateErr   = reshape(BKG_RateErr,BKGnqU.*nRingLoop,1); % shape back to nqU*nRing x 1
            BKG_AsimovqU  = reshape(repmat(BKG_Asimov,BKGnqU,1),BKGnqU.*nRingLoop,1);
            if numel(RingCorrCoeff)==1
                BKGCorrMat    = diag(ones(BKGnqU.*nRingLoop,1))+RingCorrCoeff.*(repmat(diag(ones(BKGnqU,1) ),nRingLoop,nRingLoop)-diag(ones(BKGnqU.*nRingLoop,1)));
            else
               fprintf('more complex ring to ring correlation matrix to be implemented...\n');
               return
            end
            BKGCovMat     = BKG_RateErr.*BKGCorrMat.*BKG_RateErr';
            
            BKG         = mvnrnd(BKG_AsimovqU,BKGCovMat,obj.nTrials); % ranomize with multivariate distribution
            BKG         = permute(reshape(BKG,obj.nTrials,BKGnqU,nRingLoop),[2,1,3]); % shape back to nqU x nTrials x nRing 
            BKG_RateErr = reshape(BKG_RateErr,BKGnqU,nRingLoop);   % shape back to nqU x nRing
            
            if strcmp(Mode,'Gauss')
                if numel(RingCorrCoeff)==1
                    GaussCorrMat = RingCorrCoeff.*ones(nRingLoop)+(1-RingCorrCoeff).*diag(ones(1,nRingLoop));
                else
                    GaussCorrMat = RingCorrCoeff;
                end
                GaussCovMat  = MaxSlopeCpsPereV.*GaussCorrMat.*MaxSlopeCpsPereV';
                Slopes       = mvnrnd(zeros(nRingLoop,1),GaussCovMat,obj.nTrials)';
            end
        else
            BKG           = BKG_Asimov + BKG_RateErr.*randn(BKGnqU,obj.nTrials);
            if strcmp(Mode,'Gauss')
                Slopes = MaxSlopeCpsPereV.*randn(1,obj.nTrials);% slopes
            end
        end
       
        for rl=1:1:nRingLoop  
            par        = zeros(2,obj.nTrials);
            err        = zeros(2,obj.nTrials);
            chi2min    = zeros(obj.nTrials,1);
            CutOff     = ones(1,obj.nTrials);
            
            % Definitions
            BkgPlot_Data = zeros(BKGnqU,obj.nTrials);              % Store for sanity plot
            Data         = zeros(BKGnqU,3,obj.nTrials);
            Bkg_Fit      = zeros(obj.StudyObject.nqU,obj.nTrials); % Store for covmat computation
           
            % Fit Background Rate (in counts per second)
            progressbar(sprintf('Compute Bkg CM ring %.0f out of %.0f',rl,nRingLoop));
            for i=1:obj.nTrials
                progressbar(i/obj.nTrials);
                Data(:,:,i)    = [obj.StudyObject.qU(BKGIndex:end,rl),BKG(:,i,rl), BKG_RateErr(:,rl)];
                BKG_i         = wmean(BKG(:,i,rl),1./BKG(:,i,rl));
                Slope_i       = 0;
                parInit       = [BKG_i+1e-2*rand(1), Slope_i+1e-4*rand(1)];
                % Call Minuit
                tmparg = sprintf(['set pri -10; min ; migrad '],'');
                Args   = {parInit, squeeze(Data(:,:,i)), '-c', tmparg};
                if strcmp(Mode,'SlopeFit')
                    [par(:,i), err(:,i),chi2min(i), ~] = fminuit('chi2BKG',Args{:});
                    BkgPlot_Data(:,i) = squeeze(Data(:,2,i)); %for plot in cps
                    Slopes(rl,i) = par(2,i);
                    
                    % Exclude poor fits with high chisquare & large slopes
                    % excluded by additional measurements - KNM1 official propaganda -
                    if  chi2min(i)>Chi2CutOff
                        Bkg_Fit(:,i)      = nan(numel(obj.StudyObject.qU),1); CutOff(i) = 0;
                        SlopesExcl(rl,i)  = 0;
                    elseif abs(par(2,i))<= MaxSlopeCpsPereV(rl) 
                        par(1,i) = BKG_Asimov(rl); % do not vary offset
                        Bkg_Fit(:,i)      = ComputeBkgSlope(par(:,i),obj.StudyObject.qU(:,rl));
                    else
                        Bkg_Fit(:,i)      = nan(numel(obj.StudyObject.qU(:,rl)),1); CutOff(i) = 0;
                        SlopesExcl(rl,i)  = 0;
                    end
                elseif strcmp(Mode,'Gauss')
                    par(1,i) = BKG_Asimov(rl);%mean(BKG(:,i,rl));        % offset
                    par(2,i) = Slopes(rl,i);% slopes
                    Bkg_Fit(:,i)  = ComputeBkgSlope(par(:,i),obj.StudyObject.qU(:,rl));
                end
            end
            CutOff                    = logical(CutOff); % transform 0 and 1 into logicals
            nGoodTrials(rl) = sum(SlopesExcl(rl,:));
            
            % Build Integral Spectrum Underlying Background Spectrum
            TBDIS_V(rl,:,:)    = Bkg_Fit.*obj.StudyObject.TimeSec(:,rl).*obj.StudyObject.qUfrac(:,rl);
            BkgPlot_Data      = BkgPlot_Data(:,:);
        end
        
        % convert TBDIS_V back into matrix. for rings: use common good trials   
        if nRingLoop>1  
            SlopesExclCommon =  logical(prod(SlopesExcl));
            nCommonTrials = sum(SlopesExclCommon);
            %  nCommonTrials      = min(nGoodTrials);
            TBDIS_V = permute(TBDIS_V(:,:,SlopesExclCommon),[2,1,3]); % nqU x nRing x n Trials
            %TBDIS_V            = cell2mat(cellfun(@(x) x(:,1:nCommonTrials),TBDIS_V,'UniformOutput',0));
            TBDIS_V = reshape(TBDIS_V,[nRingLoop.*obj.StudyObject.nqU,nCommonTrials]);
        else
            TBDIS_V = squeeze(TBDIS_V(:,:,logical(SlopesExcl)));
        end
        % Compute Covariance Matrix
        obj.CovMat         = cov(TBDIS_V');
        
        % Compute Fractional Covariance Matrix
        BKGIS          = obj.StudyObject.BKG_RateSec.*obj.StudyObject.TimeSec.*obj.StudyObject.qUfrac;
        if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
            % Reshape for Ring Case
            BKGIS = reshape(BKGIS,[obj.StudyObject.nqU*obj.StudyObject.nRings,1]);
        end
        obj.CovMatFrac = bsxfun(@rdivide,obj.CovMat,BKGIS);      %divides 1st row of CovMat by TBDIS(1), etc...
        obj.CovMatFrac = bsxfun(@rdivide,obj.CovMatFrac,BKGIS'); %divides columnwise
        obj.CovMatFrac(isinf(obj.CovMatFrac)) = 0; %for background region
        obj.CovMatFrac(isnan(obj.CovMatFrac)) = 0; %for background region
        
        % Compute Shape Only & Frac Shape Only Covariance Matrices
        [~, obj.CovMatFracShape] = obj.DecomposeCM('CovMatFrac',obj.CovMatFrac,'exclDataStart',1,'BkgCM','ON');
        obj.MultiCovMat.CM_BkgShape   = obj.CovMat;
        obj.MultiCovMatFrac.CM_BkgShape  = obj.CovMatFrac;
        obj.MultiCovMatFracShape.CM_BkgShape = obj.CovMatFracShape;
        save(obj.CovMatFile,'obj','TBDIS_V','BKG_Asimov','par','err','chi2min',...
            'Bkg_Fit','BKG_i','BKG','BKGIndex','BkgPlot_Data','Data','BKGnqU','MaxSlopeCpsPereV',...
            'SlopesExcl','Slopes');
        
        if nRingLoop>1
           save(obj.CovMatFile,'BKGCorrMat','BKGCovMat','-append')
        end
    end

    % Display
    switch Display
        case 'ON'
            savedir = './plots/pdf/';
            savedir2 = './plots/png/';
            MakeDir(savedir);
            MakeDir(savedir2);
            if strcmp(Mode,'SlopeFit')
                % Background (Rate) Spectra
                
                myMainTitle = sprintf('KATRIN - Background Toy Monte Carlo'); maintitle   = myMainTitle;
                %savefile    = [savedir,sprintf('/ComputeCM_Background_Fits_Bkg%.0fmcps_Constrain%.0fMcpsPerKeV.pdf',BKG_i*1e3, MaxSlopeCpsPereV)];
                fig1        = figure('Units','normalized','pos',[0.1 0.1 0.5 0.4]);
                
                Q_ref = 18574;
                pfit = plot(obj.StudyObject.qU(qUStartIndex:end,1)-Q_ref,1e3*Bkg_Fit(qUStartIndex:end,1:1000),...
                    'LineWidth',2,'Color', [rgb('SlateGray') 0.25]);%[1.0000    0.6445      0 0.2]);
                hold on;
                sg = scatter(reshape(squeeze(Data(:,1,:))-Q_ref,[BKGnqU*obj.nTrials,1]),...
                    1e3.*reshape(squeeze(Data(:,2,:)),[BKGnqU*obj.nTrials,1]),'o','filled');
                if size(Data,3)>=5000
                    sg.MarkerFaceAlpha = 0.01;
                elseif size(Data,3)>=1000
                    sg.MarkerFaceAlpha = 0.04;
                else
                    sg.MarkerFaceAlpha = 0.2;
                end
                sg.MarkerFaceColor = rgb('DodgerBlue');
                
                tmp = scatter(1,1,'o','filled'); % just for display in legend
                tmp.MarkerFaceColor=sg.MarkerFaceColor;
                tmp.MarkerFaceAlpha = 0.5;
                
                hold off
                set(gca,'FontSize',20);
                xlabel(sprintf('Retarding Potential qU - %.0f (eV)',Q_ref));
                ylabel('Background (mcps)');
                hold off;
                PrettyFigureFormat('FontSize',24);
                if strcmp(Mode,'SlopeFit')
                    legStr = ' Linear fits';
                else
                    legStr = ' Randomized slopes';
                end
                
                if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                    leg = legend([tmp, pfit(1)],sprintf(' MC background spectra (ring %.0f)',obj.StudyObject.nRings),...
                        legStr,'Location','northeast');
                else
                    leg = legend([tmp, pfit(1)],' MC background spectra',legStr,'Location','northeast');
                end
                
                leg.EdgeColor = rgb('Silver');
                leg.FontSize = get(gca,'FontSize');
                xlim([obj.StudyObject.qU(qUStartIndex,1)-Q_ref,obj.StudyObject.qU(end,1)-Q_ref+5]);
                ylim([min(min(Data(:,2,:)))*1e3-0.005,max(max(Data(:,2,:)))*1e3+0.1])
                % sgtitle(maintitle,'FontSize',22);
                savefile    =  [savedir,strrep(cm_name,'.mat','_ToyMC.png')];
                print(fig1,savefile,'-dpng','-r500');
                % export_fig(fig1,savefile);
                fprintf('Save plot to %s \n',savefile);
            end
          
           if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
               nRingLoop = obj.StudyObject.nRings;
           else
               nRingLoop = 1;
           end
           if nRingLoop >1
               fig2        = figure('Units','normalized','pos',[0.1 0.1 0.7 0.7]);
               s= cell(nRingLoop,1);
               for i=1:nRingLoop
                   s{i} = subplot(floor(nRingLoop/2),ceil(nRingLoop/2),i);
                   h1 = histogram(Slopes(i,:).*1e6,'FaceColor',rgb('LightGray'));
                   if strcmp(Mode,'SlopeFit')
                       hold on;
                       SlopesExclCommon =  logical(prod(SlopesExcl));
                       SlopeOK = Slopes(i,SlopesExclCommon);
                       h22 = histogram(Slopes(i,abs(Slopes(i,:))<MaxSlopeCpsPereV(i)).*1e06,'FaceColor',rgb('PowderBlue'),...
                           'FaceAlpha',1,'BinWidth',h1.BinWidth);
                       h2 = histogram(SlopeOK.*1e06,'FaceColor',rgb('DodgerBlue'),...
                           'FaceAlpha',1,'BinWidth',h1.BinWidth);
                       leg = legend('MC fit slopes',sprintf('MC fit slopes < %.1f mcps / keV',MaxSlopeCpsPereV(i)*1e6),...
                             'MC fit slopes acctepted');
                   else
                        leg = legend(sprintf('Randomized slopes \\sigma = %.1f mcps / keV',MaxSlopeCpsPereV(i)*1e6));
                   end
                   xlabel('Background slope (mcps / keV)');
                   PrettyFigureFormat('FontSize',18)
                   leg.Location='northwest';
                   leg.EdgeColor = rgb('Silver');
                   leg.Title.String = sprintf('Pseudo-Ring %.0f',i);
                   hold off;
                   ylim([0 3200])
                   xlim([-20 20])    
               end
               sgtitle(sprintf('Background slope contrain (all pixels) = %.1f mcps / keV',MaxSlopeCpsPereV_i*1e6));
               linkaxes([s{:}],'x');
               
           else
               fig2        = figure('Units','normalized','pos',[0.1 0.1 0.4 0.4]);
               h1 = histogram(Slopes.*1e6,'FaceColor',rgb('LightGray'));
               
               if strcmp(Mode,'SlopeFit')
                   hold on;
                   SlopeOK = Slopes(logical(SlopesExcl));
                   h2 = histogram(SlopeOK.*1e06,'FaceColor',rgb('DodgerBlue'),...
                       'FaceAlpha',1,'BinWidth',h1.BinWidth);
                   leg = legend('MC fit slopes',sprintf('MC fit slopes < %.1f mcps / keV',MaxSlopeCpsPereV*1e6));
               else
                   h1.FaceColor = rgb('DodgerBlue');
                   leg = legend(sprintf('Randomized slopes \\sigma = %.1f mcps / keV',MaxSlopeCpsPereV*1e6));
               end
               leg.Location='northwest';
               leg.EdgeColor = rgb('Silver');
               xlabel('Background slope (mcps / keV)');
               PrettyFigureFormat('FontSize',22)
               hold off;
               ylim([0 2500])
           end
           savefile    =  [savedir,strrep(cm_name,'.mat','_Hist.pdf')];
           print(fig2,strrep(strrep(savefile,'.pdf','.png'),'plots/','plots/png/'),'-dpng','-r500');
           export_fig(fig2,savefile);
           fprintf('Save plot to %s \n',savefile);
           
           if nRingLoop>1
               xlim([-MaxSlopeCpsPereV_i,MaxSlopeCpsPereV_i]*1e6);
               savefile    =  [savedir,strrep(cm_name,'.mat','_HistZoom.pdf')];
               export_fig(fig2,savefile);
               print(fig2,strrep(strrep(savefile,'.pdf','.png'),'plots/','plots/png/'),'-dpng','-r500');
               fprintf('Save plot to %s \n',savefile);
           end
           %             % Trials Toy MC Features
           %             myMainTitle = sprintf('KATRIN - Background Toy Monte Carlo (Slope Constrain = %.0f mcps/keV)',MaxSlopeCpsPereV*1e6); maintitle   = myMainTitle;
           %             savefile    =  [savedir,sprintf('/ComputeCM_Background_2_Bkg%.0fmcps_Hists_Constrain%.0fMcpsPerKeV.pdf',BKG_i*1e3, MaxSlopeCpsPereV)];
           %             fig2        = figure('Name','KATRIN - Background Toy Monte Carlo','NumberTitle','off','rend','painters','pos',[1 1 1000 1000]);
           %             a = annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
           %             a.FontSize=18;a.FontWeight='bold';
           %
           %             subplot(2,2,2);
           %             h1 = histogram(par(2,CutOff)*1e6,15,'Normalization','probability','FaceColor',rgb('IndianRed'),'LineWidth',2);
           %             xlabel('qU-slope: mcps/keV');ylabel('frequency');
           %             PrettyFigureFormat;set(gca,'FontSize',16);
           %
           %             subplot(2,2,4);
           %             histogram(chi2min(CutOff),15,'Normalization','probability','FaceColor',rgb('IndianRed'),'LineWidth',2);
           %             xlabel('\chi^2');ylabel('frequency');
           %             PrettyFigureFormat;set(gca,'FontSize',16);
           %             subplot(2,2,1);
           %             histogram(par(1,CutOff)*1e3,15,'Normalization','probability','FaceColor',rgb('IndianRed'),'LineWidth',2);
           %             xlabel('Flat Component (mcps)');ylabel('frequency');
           %             PrettyFigureFormat;set(gca,'FontSize',16);
           %             subplot(2,2,3);
           %             histogram(err(1,CutOff)*1e3,15,'Normalization','probability','FaceColor',rgb('IndianRed'),'LineWidth',2);
           %             xlabel('Baseline Component Error (mcps)');ylabel('frequency');
           %             %xlim([0 5]);
           %             set(gca,'YScale','log');
           %             PrettyFigureFormat;set(gca,'FontSize',16);
           %             export_fig(fig2,savefile,'-painters');
           %
           %             % Covariance Matrix
           %             myMainTitle = sprintf('KATRIN - Background Covariance Matrix'); maintitle   = myMainTitle;
           %             savefile    =  [savedir,sprintf('/ComputeCM_Background_Convergence_Bkg%.0fmcps_Constrain%.0fMcpsPerKeV.png',BKG_i*1e3, MaxSlopeCpsPereV)];
           %             fig3        = figure('Name','KATRIN - Background Covariance Matrix','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
           %             a = annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
           %             a.FontSize=18;a.FontWeight='bold';
           %             qULong=obj.StudyObject.qU(qUStartIndex:end);
           %             subplot(2,2,[1,3]);
           %             imagesc(obj.CovMatFrac(qUStartIndex:end,qUStartIndex:end));
           %             PrettyFigureFormat
           %             c = colorbar('northoutside');
           %             c.Label.String = 'fractional covariance (bkg)';
           %             c.FontSize = 18;
           %             colormap(parula);
           %             pbaspect([1 1 1])
           %             set(gca,'xtick',[1 obj.StudyObject.nqU]),set(gca,'ytick',[])
           %             qUmin = sprintf('qU_{min} = E_0 - %.0fV',abs(obj.StudyObject.qU(qUStartIndex)-obj.StudyObject.Q_i));
           %             qUmax = sprintf('qU_{max} = E_0 + %.0fV',obj.StudyObject.qU(end)-obj.StudyObject.Q_i);
           %             set(gca,'xticklabel',{qUmin,qUmax}); set(gca,'yticklabel',[]);
           %             set(gca,'FontSize',16)
           %
           %             subplot(2,2,2);
           %             subruntime  = obj.StudyObject.TimeSec.*obj.StudyObject.qUfrac(qUStartIndex:end);
           %             meanBcounts = mean(Bkg_Fit(qUStartIndex:end,:),2).*subruntime;
%             %                 [SysLine, SysArea] =  boundedline(qULong-obj.StudyObject.Q_i,...
%             %                     1e3.*mean(Bkg_Fit(qUStartIndex:end,:),2),...
%             %                     1e3.*std(permute(Bkg_Fit(qUStartIndex:end,:),[2,1])));
%             [SysLine, SysArea] =  boundedline(qULong-obj.StudyObject.Q_i,...
%                 1e3*meanBcounts./subruntime,...
%                 1e3*sqrt(diag(obj.CovMat(qUStartIndex:end,qUStartIndex:end)))./subruntime); % B counts
%             hold on;
%             %                 pStat = errorbar(obj.StudyObject.qU(qUStartIndex:end)-obj.StudyObject.Q_i,...
%             %                     1e3.*mean(Bkg_Fit(qUStartIndex:end,:),2),...
%             %                     1e3.*sqrt(mean(Bkg_Fit(qUStartIndex:end,:),2)./sqrt((obj.StudyObject.TimeSec.*obj.StudyObject.qUfrac(qUStartIndex:end)))),...
%             %                     'o','Color',rgb('GoldenRod'),'MarkerFaceColor',rgb('RoyalBlue'),'MarkerEdgeColor',rgb('RoyalBlue'),'MarkerSize',2);
%             pStat = errorbar(obj.StudyObject.qU(qUStartIndex:end)-obj.StudyObject.Q_i,...
%                 1e3*meanBcounts./subruntime,...
%                 1e3*sqrt(meanBcounts)./subruntime,...
%                 'o','Color',rgb('GoldenRod'),'MarkerFaceColor',rgb('RoyalBlue'),'MarkerEdgeColor',rgb('RoyalBlue'),'MarkerSize',2);
%             SysArea.FaceColor = rgb('CadetBlue');
%             SysLine.LineWidth = 5; SysLine.Color = rgb('RoyalBlue');
%             SysArea.FaceAlpha = 0.5;
%             xlabel(sprintf('Retarding Potential qU - %.01f (eV)',obj.StudyObject.Q_i));
%             ylabel('Background (mcps)');
%             leg =legend([SysLine, SysArea,pStat],'<B>','1 \sigma_{B,sys}','1 \sigma_{B,stat}','Location','northeast');
%             leg.FontSize = 18;
%             leg.NumColumns = 3;
%             legend boxoff
%             xlim([min(qULong) max(qULong)]-obj.StudyObject.Q_i);
%             ylimdown = (1e3*mean(mean(Bkg_Fit(qUStartIndex:end,:),2)))-1e3*min(sqrt(mean(Bkg_Fit(qUStartIndex:end,:),2)./(obj.StudyObject.TimeSec.*obj.StudyObject.qUfrac(qUStartIndex:end))));%1e3*(mean(mean(Bkg_Fit,2))-max(std(permute(Bkg_Fit,[2,1]))));
%             ylimup   = 1e3*mean(mean(Bkg_Fit(qUStartIndex:end,:),2))+(1e3.*max(sqrt(mean(Bkg_Fit(qUStartIndex:end,:),2)./(obj.StudyObject.TimeSec.*obj.StudyObject.qUfrac(qUStartIndex:end)))));%1e3*(mean(mean(Bkg_Fit,2))+max(std(permute(Bkg_Fit,[2,1]))));
%             if ylimup<=1e3*max(max(mean(Bkg_Fit(qUStartIndex:end,:),2)+std(permute(Bkg_Fit(qUStartIndex:end,:),[2,1]))))
%                 ylimup =1e3*max(max(mean(Bkg_Fit(qUStartIndex:end,:),2)+std(permute(Bkg_Fit(qUStartIndex:end,:),[2,1]))));
%                 ylimdown =1e3*min(min(mean(Bkg_Fit(qUStartIndex:end,:),2)-std(permute(Bkg_Fit(qUStartIndex:end,:),[2,1]))));
%             end
%             ylim([0.97*ylimdown, 1.01*ylimup]);
%             PrettyFigureFormat;
%             set(gca,'FontSize',16);
%             
%             subplot(2,2,4);
%             Nmax = max(size(BkgIS_SampleNorm));
%             StepSize = 10;
%             CovMatnTrials = cell(round(Nmax/StepSize),1);     % CovMat for different nTrials
%             CovMatTrace = zeros(round(Nmax/StepSize),1);      % Trace of CovMats
%             for x = StepSize:StepSize:Nmax                    %Compute CovMats für different nTrials
%                 xn = x/StepSize;
%                 CovMatnTrials{xn,1} = cov(BkgIS_SampleNorm(:,1:x)'); %CovMats
%                 CovMatTrace(xn) = norm(CovMatnTrials{xn,1},'fro');
%             end
%             
%             while max(size([StepSize:StepSize:Nmax]))~=max(size(CovMatTrace))
%                 if max(size([StepSize:StepSize:Nmax]))>=max(size(CovMatTrace))
%                     Nmax = Nmax-StepSize;
%                 elseif max(size([StepSize:StepSize:Nmax]))<=max(size(CovMatTrace))
%                     Nmax = Nmax+StepSize;
%                 end
%             end
%             
%             plot([StepSize:StepSize:Nmax]',CovMatTrace,'Color',rgb('GoldenRod'),'LineWidth',4);
%             xlabel('samples');
%             ylabel('|| M ||');
%             PrettyFigureFormat;
%             set(gca,'FontSize',16);
%             xlim([0 Nmax-2*StepSize]);
%             print(fig3,savefile,'-dpng','-r400');
    end
    
    obj.nTrials = TrialSave;
    
end

function ComputeCM_TCoff(obj,varargin)
    fprintf('--------------------------------------------------------------------------\n')
    cprintf('blue','CovarianceMatrix:ComputeCM_TCoff: Compute Theoretical Corrections (TC) Covariance Matrix (case: TC=OFF) \n')
    % Inputs
    p = inputParser;
    p.addParameter('plot','ON',@(x)ismember(x,{'ON','OFF'}));
    p.parse(varargin{:});
    plot  = p.Results.plot;
    
    % Labeling
    covmat_path = [getenv('SamakPath'),sprintf('inputs/CovMat/TC/CM/')];
    if strcmp(obj.SysEffect.TCoff_RAD,'ON') &&   strcmp(obj.SysEffect.TCoff_OTHER,'ON')
        covmat_filename = sprintf('TCoff_CovMat_%s.mat', obj.StudyObject.TD);
    elseif strcmp(obj.SysEffect.TCoff_RAD,'ON') && strcmp(obj.SysEffect.TCoff_OTHER,'OFF')
        covmat_filename = sprintf('TCoffRad_CovMat_%s.mat', obj.StudyObject.TD);
    elseif strcmp(obj.SysEffect.TCoff_RAD,'OFF') &&   strcmp(obj.SysEffect.TCoff_OTHER,'ON')
        covmat_filename = sprintf('TCoffOTHER_CovMat_%s.mat', obj.StudyObject.TD);
    end
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        % add ring information for ringwise covariance matrices
        covmat_filename = strrep(covmat_filename,'.mat',sprintf('_Ring%s.mat',obj.StudyObject.FPD_RingMerge));
    end
    obj.CovMatFile = strcat(covmat_path,covmat_filename);
    
    MakeDir(covmat_path);
    
    % Check if CM is already computed
    if exist(obj.CovMatFile,'file')==2 && strcmp(obj.RecomputeFlag,'OFF')
        fprintf(2,'CovarianceMatrix:ComputeCM_TCoff: Loading TCoff CM from File \n')
        if strcmp(obj.SysEffect.TCoff_RAD,'ON') &&   strcmp(obj.SysEffect.TCoff_OTHER,'ON')
            obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','TCoff');
            obj.MultiCovMat.CM_TCoff = obj.CovMat; % in case normalization changed
        elseif strcmp(obj.SysEffect.TCoff_RAD,'ON') && strcmp(obj.SysEffect.TCoff_OTHER,'OFF')
            obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','TCoff_RAD');
            obj.MultiCovMat.CM_TCoff_RAD = obj.CovMat;
        elseif strcmp(obj.SysEffect.TCoff_RAD,'OFF') && strcmp(obj.SysEffect.TCoff_OTHER,'OFF')
            obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','TCoff_OTHER');
            obj.MultiCovMat.CM_TCoff_OTHER = obj.CovMat;
        end
        return
    end
    
    % Computing TBDIS TC=off
    ComputeTBDDS(obj.StudyObject);
    ComputeTBDIS(obj.StudyObject);
    TBDIS_V    = obj.StudyObject.TBDIS;
    
    % Computing TBDIS TC=on
    %   .* obj.StudyObject.ComputeScreeningCorr()...
    if strcmp(obj.SysEffect.TCoff_RAD,'ON')
        obj.StudyObject.TBDDS = obj.StudyObject.TBDDS...
            .* obj.StudyObject.ComputeRadiativeCorr();
    end
    if strcmp(obj.SysEffect.TCoff_OTHER,'ON')
        obj.StudyObject.TBDDS = obj.StudyObject.TBDDS...
            .* obj.StudyObject.ComputeFiniteExtChargeCorr()...
            .* obj.StudyObject.ComputeEEexchangeCorr(2)...
            .* obj.StudyObject.ComputeRecoilWmVmACorr()...
            .* obj.StudyObject.ComputeWintFiniteSizeCorr()...
            .* obj.StudyObject.ComputeRecoilCoulombCorr();
    end
    ComputeTBDIS(obj.StudyObject);
    TBDISTC_V  = obj.StudyObject.TBDIS;
    
    % Reshape
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        TBDIS_V   = reshape(TBDIS_V,[obj.StudyObject.nqU*obj.StudyObject.nRings,1]);
        TBDISTC_V = reshape(TBDISTC_V,[obj.StudyObject.nqU*obj.StudyObject.nRings,1]);
    end
    
    % Compute Covariance Matrix
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        obj.CovMat = bsxfun(@times,ones(obj.StudyObject.nqU*obj.StudyObject.nRings,obj.StudyObject.nqU*obj.StudyObject.nRings),(TBDIS_V-TBDISTC_V));
    else
        obj.CovMat = bsxfun(@times,ones(obj.StudyObject.nqU,obj.StudyObject.nqU),(TBDIS_V-TBDISTC_V));
    end
    obj.CovMat = bsxfun(@times,obj.CovMat,(TBDIS_V-TBDISTC_V)');
    
    if strcmp(obj.SysEffect.TCoff_RAD,'ON') && strcmp(obj.SysEffect.TCoff_OTHER,'ON')
        obj.MultiCovMat.CM_TCoff  = obj.CovMat;
    elseif strcmp(obj.SysEffect.TCoff_RAD,'ON') && strcmp(obj.SysEffect.TCoff_OTHER,'OFF')
        obj.MultiCovMat.CM_TCoff_RAD = obj.CovMat;
    elseif strcmp(obj.SysEffect.TCoff_RAD,'OFF') && strcmp(obj.SysEffect.TCoff_OTHER,'ON')
        obj.MultiCovMat.CM_TCoff_OTHER = obj.CovMat;
    end
    
    % Plot
    if strcmp(obj.SanityPlots,'ON')
        figure
        if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
            plt1= stairs(TBDIS_V);
            hold on
            plt1= stairs(TBDISTC_V);
            hold off
        else
            plt1= Plot(obj.StudyObject.qU-obj.StudyObject.Q,TBDIS_V,...
                obj.StudyObject.qU-obj.StudyObject.Q,TBDISTC_V);
        end
        plt1.LineWidth = 2;
        plt1.LineStyle = {'-','--'};
        plt1.Markers   = {'s','d'};
        plt1.MarkerSpacing = [1,1];
        pltstr         = sprintf('Theoretical corrections - %s - %g sec',obj.StudyObject.TD,obj.StudyObject.TimeSec);
        plt1.Title     = pltstr; % plot title
        plt1.XLabel    = 'E-qU (V)'; % xlabel
        plt1.YLabel    = 'Rate'; %ylabel
        plt1.YScale    = 'lin'; % 'linear' or 'log'
        %plt1.YLim      = [min(BM1IS_V./(obj.StudyObject.TimeSec*obj.StudyObject.qUfrac))*0.9 max(BM1IS_V./(obj.StudyObject.TimeSec*obj.StudyObject.qUfrac))*1.1];
        plt1.XScale    = 'lin'; % 'linear' or 'log'
        plt1.FontSize  = 16;
        plt1.Legend    = {'TC OFF','TC ON'};
        plt1.LegendLoc = 'SouthEast';
        plttitle       = sprintf('TheoreticalCorrections_1.png');
        plt1.export(plttitle);
        figure
        plt2 = Plot(obj.StudyObject.qU-obj.StudyObject.Q,sqrt(diag(obj.MultiCovMat.CM_TCoff)),...
            obj.StudyObject.qU-obj.StudyObject.Q,sqrt(TBDIS_V));
        plt2.LineWidth = 2;
        plt2.LineStyle = {'-','--'};
        plt2.Markers   = {'s','d'};
        pltstr         = sprintf('Theoretical corrections - %s - %g secs',obj.StudyObject.TD,obj.StudyObject.TimeSec);
        plt2.Title     = pltstr; % plot title
        plt2.XLabel    = 'E-qU (V)'; % xlabel
        plt2.YLabel    = 'Count Difference'; %ylabel
        plt2.YScale    = 'lin'; % 'linear' or 'log'
        plt2.XScale    = 'lin'; % 'linear' or 'log'
        plt2.Legend    = {'Syst. Error (background)','Stat. Error (sim data)'};
        %plt2.YLim      =  [min(sqrt(TBDIS_V))/10 max(sqrt(TBDIS_V))];
        plt2.FontSize  = 16;
        plttitle       = sprintf('TheoreticalCorrections_2.png');
        plt2.export(plttitle);
    end
    
    % Save
    save(obj.CovMatFile,'obj','TBDIS_V','TBDISTC_V','-mat');
    
    %Compute and save Fractional CM
    obj.ComputeFracCM('Mode','CM2Frac');
    
    CM = {obj.CovMat, obj.CovMatFrac};
    for i=1:numel(CM)
        try mvnrnd(obj.StudyObject.TBDIS,CM{i},1);
        catch
            fprintf(2,'Warning: FitCM not positive semi definite! \n Conversion with nearestSPD \n');
            CM{i} = nearestSPD(CM{i});
        end
    end
    
    %Compute Decomposed CM
    obj.DecomposeCM('CovMatFrac',obj.CovMatFrac,'exclDataStart',1);
    if strcmp(obj.SysEffect.TCoff_RAD,'ON') && strcmp(obj.SysEffect.TCoff_OTHER,'ON')
        obj.MultiCovMatFrac.CM_TCoff  = obj.CovMatFrac;
        obj.MultiCovMatFracShape.CM_TCoff = obj.CovMatFracShape;
        obj.MultiCovMatFracNorm.CM_TCoff = obj.CovMatFracNorm;
    elseif strcmp(obj.SysEffect.TCoff_RAD,'ON') && strcmp(obj.SysEffect.TCoff_OTHER,'OFF')
        obj.MultiCovMatFrac.CM_TCoff_RAD = obj.CovMatFrac;
        obj.MultiCovMatFracShape.CM_TCoff_RAD = obj.CovMatFracShape;
        obj.MultiCovMatFracNorm.CM_TCoff_RAD = obj.CovMatFracNorm;
    elseif strcmp(obj.SysEffect.TCoff_RAD,'OFF') && strcmp(obj.SysEffect.TCoff_OTHER,'ON')
        obj.MultiCovMatFrac.CM_TCoff_OTHER = obj.CovMatFrac;
        obj.MultiCovMatFracShape.CM_TCoff_OTHER = obj.CovMatFracShape;
        obj.MultiCovMatFracNorm.CM_TCoff_OTHER = obj.CovMatFracNorm;
    end
    save(obj.CovMatFile, 'obj','-append');
end
function ComputeCM_DOPoff(obj,varargin)
    fprintf('--------------------------------------------------------------------------\n')
    cprintf('blue','CovarianceMatrix:ComputeCM_DOPoff: Compute Doppler (DOP) Covariance Matrix (case: DOP=OFF) \n')
    % Inputs
    p = inputParser;
    p.addParameter('plot','ON',@(x)ismember(x,{'ON','OFF'}));
    p.parse(varargin{:});
    plot  = p.Results.plot;
    
    %Labeling
    covmat_path = [getenv('SamakPath'),sprintf('/inputs/CovMat/DOP/CM/')];
    MakeDir(covmat_path);
    covmat_filename = sprintf('DOPoff_CovMat_%s_%u-Runs.mat',...
        obj.StudyObject.TD,obj.nRuns);
    obj.CovMatFile = strcat(covmat_path,covmat_filename);
    
    %Check if CM is already computed
    if exist(obj.CovMatFile,'file')==2 && strcmp(obj.RecomputeFlag,'OFF')
        fprintf(2,'CovarianceMatrix:ComputeCM_DOPoff: Loading DOPoff CM from File \n')
        obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','DOPoff');
        obj.MultiCovMat.CM_DOPoff = obj.CovMat;
        return
    end
    
    % Computing TBDIS DOP=off
    ComputeTBDDS(obj.StudyObject);
    ComputeTBDIS(obj.StudyObject); TBDIS_V    = obj.StudyObject.TBDIS;
    
    % Computing TBDIS DOP=on
    obj.AltObject.NormFactorTBDDS = obj.StudyObject.NormFactorTBDDS; % Normalization TO BE FIXED
    ComputeTBDDS(obj.AltObject);
    ComputeTBDIS(obj.AltObject); TBDISDOP_V  = obj.AltObject.TBDIS;
    
    % Compute Covariance Matrix
    obj.CovMat = bsxfun(@times,ones(obj.StudyObject.nqU,obj.StudyObject.nqU),-(TBDIS_V-TBDISDOP_V));
    obj.CovMat = bsxfun(@times,obj.CovMat,-(TBDIS_V-TBDISDOP_V)');
    obj.MultiCovMat.CM_DOPoff  = obj.CovMat;
    
    % Plot
    switch plot
        case 'ON'
            figure
            plt1= Plot(obj.StudyObject.qU-obj.StudyObject.Q,TBDIS_V./TBDIS_V ,...
                obj.StudyObject.qU-obj.StudyObject.Q,TBDISDOP_V./TBDIS_V);
            plt1.LineWidth = 2;
            plt1.LineStyle = {'-','--'};
            plt1.Markers   = {'s','d'};
            plt1.MarkerSpacing = [1,1];
            pltstr         = sprintf('Doppler Effect - %s - %g sec',obj.StudyObject.TD,obj.StudyObject.TimeSec);
            plt1.Title     = pltstr; % plot title
            plt1.XLabel    = 'E-qU (V)'; % xlabel
            plt1.YLabel    = 'Rate'; %ylabel
            plt1.YScale    = 'lin'; % 'linear' or 'log'
            %plt1.YLim      = [min(BM1IS_V./(obj.StudyObject.TimeSec*obj.StudyObject.qUfrac))*0.9 max(BM1IS_V./(obj.StudyObject.TimeSec*obj.StudyObject.qUfrac))*1.1];
            plt1.XScale    = 'lin'; % 'linear' or 'log'
            plt1.FontSize  = 16;
            plt1.Legend    = {'Doppler OFF','Doppler ON'};
            plt1.LegendLoc = 'SouthEast';
            plttitle       = sprintf('Doppler_1.png');
            plt1.export(plttitle);
            figure
            plt2 = Plot(obj.StudyObject.qU-obj.StudyObject.Q,sqrt(diag(obj.MultiCovMat.CM_DOPoff)),...
                obj.StudyObject.qU-obj.StudyObject.Q,sqrt(TBDIS_V));
            plt2.LineWidth = 2;
            plt2.LineStyle = {'-','--'};
            plt2.Markers   = {'s','d'};
            pltstr         = sprintf('Doppler Effect - %s - %g secs',obj.StudyObject.TD,obj.StudyObject.TimeSec);
            plt2.Title     = pltstr; % plot title
            plt2.XLabel    = 'E-qU (V)'; % xlabel
            plt2.YLabel    = 'Count Difference'; %ylabel
            plt2.YScale    = 'lin'; % 'linear' or 'log'
            plt2.XScale    = 'lin'; % 'linear' or 'log'
            plt2.Legend    = {'Syst. Error (background)','Stat. Error (sim data)'};
            %plt2.YLim      =  [min(sqrt(TBDIS_V))/10 max(sqrt(TBDIS_V))];
            plt2.FontSize  = 16;
            plttitle       = sprintf('Doppler_1.png');
            plt2.export(plttitle);
    end
    
    % Save
    save(obj.CovMatFile,'obj','TBDIS_V','TBDISDOP_V','-mat');
    
    % Compute Fractional CM
    obj.ComputeFracCM('Mode','CM2Frac');
    obj.MultiCovMatFrac.CM_DOPoff = obj.CovMatFrac;
    
    % Compute Decomposed CM
    obj.DecomposeCM('Option','Frac');
    obj.MultiCovMatFracShape.CM_DOPoff = obj.CovMatFracShape;
    obj.MultiCovMatFracNorm.CM_DOPoff = obj.CovMatFracNorm;
    
    % Save again
    save(obj.CovMatFile, 'obj','-append');
end
function ComputeCM_FSD(obj,varargin)
    % Covariance Matrix for FSD
    fprintf('--------------------------------------------------------------------------\n')
    cprintf('blue','CovarianceMatrix:ComputeCM_FSD: Compute FSD Covariance Matrix  \n')
    % Inputs
    p = inputParser;
    p.addParameter('StatFluct','OFF',@(x)ismember(x,{'ON','OFF'}));
    
    p.parse(varargin{:});
    StatFluct              = p.Results.StatFluct;
    
    % labeling
    IsoName = '';
    if strcmp(obj.DTFlag,'ON')
        IsoName = [IsoName,'DT-'];
    end
    if strcmp(obj.HTFlag,'ON')
        IsoName = [IsoName,'HT-'];
    end
    if strcmp(obj.TTFlag,'ON')
        IsoName = [IsoName,'TT-'];
    end
    obj.GetTDlabel;
    covmat_path =[getenv('SamakPath'),sprintf('inputs/CovMat/FSD/CM/')];
    MakeDir(covmat_path);
    covmat_filename = sprintf('FSD_%s_%sCovMat_%uTrials_%s_%.2gNormErr_%.2fGS_%.2fES_ShapeErr.mat',...
        obj.StudyObject.TTFSD,IsoName,obj.nTrials, obj.TDlabel,obj.FSDNorm_RelErr, obj.FSDShapeGS_RelErr, obj.FSDShapeES_RelErr);
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        % add ring information for ringwise covariance matrices
        covmat_filename = strrep(covmat_filename,'.mat',sprintf('_Ring%s.mat',obj.StudyObject.FPD_RingMerge));
    end
    obj.CovMatFile = strcat(covmat_path,covmat_filename);
    
    % load if already computed
    if exist(obj.CovMatFile,'file')==2 && strcmp(obj.RecomputeFlag,'OFF')
        fprintf(2,'CovarianceMatrix:ComputeCM_FSD: Loading FSD CM from File \n')
        obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','FSD');
        obj.MultiCovMat.CM_FSD = obj.CovMat; %in case normalization was changed
        return
    end
    
    % FSD from Theory (Reference)
    obj.StudyObject.LoadFSD;
    obj.StudyObject.ComputeTBDDS;
    DT_P_theory = [obj.StudyObject.DTexP_G.*obj.StudyObject.DTNormGS_i,...
        obj.StudyObject.DTexP_E.*obj.StudyObject.DTNormES_i]';
    
    nRings = numel(obj.StudyObject.MACE_Ba_T);
    % Init Variation
    NormGS = zeros(obj.nTrials,nRings); %Normalization Factors Ground + Excited State
    NormES = zeros(obj.nTrials,nRings);
    NormBiasDT = zeros(obj.nTrials,nRings);
    NormBiasHT = zeros(obj.nTrials,nRings);
    NormBiasTT = zeros(obj.nTrials,nRings);
    DT_P_norm = zeros(nRings,size(obj.StudyObject.DTexP,2),obj.nTrials);%Final State Probabilities after norm
    TT_P_norm = zeros(nRings,size(obj.StudyObject.TTexP,2),obj.nTrials);
    HT_P_norm = zeros(nRings,size(obj.StudyObject.HTexP,2),obj.nTrials);
    DT_P = repmat(obj.StudyObject.DTexP,1,1,obj.nTrials);  %Final State Probabilities before norm
    TT_P = repmat(obj.StudyObject.TTexP,1,1,obj.nTrials);
    HT_P = repmat(obj.StudyObject.HTexP,1,1,obj.nTrials);
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        TBDIS_V = zeros(obj.StudyObject.nqU,obj.StudyObject.nRings,obj.nTrials);
    else
        TBDIS_V = zeros(obj.StudyObject.nqU,obj.nTrials);
    end
    
    % Variation
    if strcmp(obj.DTFlag,'ON')
        % Vary Norm
        NormBiasDT = randn(obj.nTrials,1).*obj.StudyObject.DTNormGS_i'.*obj.FSDNorm_RelErr;
        % Bin-to-Bin uncorrelated variation
        DT_P_rand =permute(repmat(randn(size(obj.StudyObject.DTexP,2),obj.nTrials),1,1,nRings),[3,1,2]);
        DT_P(:,1:obj.StudyObject.DTGSTh,:)     = DT_P(:,1:obj.StudyObject.DTGSTh,:).* (1+(DT_P_rand(:,1:obj.StudyObject.DTGSTh,:).*obj.FSDShapeGS_RelErr));
        DT_P(:,obj.StudyObject.DTGSTh+1:end,:) = DT_P(:,obj.StudyObject.DTGSTh+1:end,:).*(1+(DT_P_rand(:,obj.StudyObject.DTGSTh+1:end,:).*obj.FSDShapeES_RelErr)); 
        DT_P(DT_P<0)=0;
    end
    if strcmp(obj.HTFlag,'ON')
        NormBiasHT = randn(obj.nTrials,1).*obj.StudyObject.HTNormGS_i'*obj.FSDNorm_RelErr;
        HT_P_rand = permute(repmat(randn(size(obj.StudyObject.HTexP,2),obj.nTrials),1,1,nRings),[3,1,2]);
        HT_P(:,1:obj.StudyObject.HTGSTh,:)  = HT_P(:,1:obj.StudyObject.HTGSTh,:).* (1+(HT_P_rand(:,1:obj.StudyObject.HTGSTh,:).*obj.FSDShapeGS_RelErr));
        HT_P(:,obj.StudyObject.HTGSTh+1:end,:) = HT_P(:,obj.StudyObject.HTGSTh+1:end,:).*(1+(HT_P_rand(:,obj.StudyObject.HTGSTh+1:end,:).*obj.FSDShapeES_RelErr));
        HT_P(HT_P<0)=0;
    end
    if strcmp(obj.TTFlag,'ON')
        NormBiasTT = randn(obj.nTrials,1).*obj.StudyObject.TTNormGS_i'*obj.FSDNorm_RelErr;
        TT_P_rand = permute(repmat(randn(size(obj.StudyObject.TTexP,2),obj.nTrials),1,1,nRings),[3,1,2]);
        TT_P(:,1:obj.StudyObject.TTGSTh,:)  = TT_P(:,1:obj.StudyObject.TTGSTh,:).* (1+(TT_P_rand(:,1:obj.StudyObject.TTGSTh,:).*obj.FSDShapeGS_RelErr));
        TT_P(:,obj.StudyObject.TTGSTh+1:end,:) = TT_P(:,obj.StudyObject.TTGSTh+1:end,:).*(1+(TT_P_rand(:,obj.StudyObject.TTGSTh+1:end,:).*obj.FSDShapeES_RelErr)); 
        TT_P(TT_P<0)=0;
    end
    
    %Start trials
    progressbar('Computing FSD CovMat')
    for i=1:1:obj.nTrials
        progressbar(i/obj.nTrials)
        obj.StudyObject.DTexP = DT_P(:,:,i); %bin-to-bin uncorrelated variation
        obj.StudyObject.TTexP = TT_P(:,:,i);
        obj.StudyObject.HTexP = HT_P(:,:,i);
        obj.StudyObject.ComputeTBDDS(...      %Inside ComputeTBDDS: SetFitBias
            'DTGS_bias',NormBiasDT(i,:)','DTES_bias', -NormBiasDT(i,:)',...
            'TTGS_bias',NormBiasTT(i,:)','TTES_bias', -NormBiasTT(i,:)',...
            'HTGS_bias',NormBiasHT(i,:)','HTES_bias', -NormBiasHT(i,:)');
        obj.StudyObject.ComputeTBDIS;
        switch StatFluct
            case 'ON'
                ComputeTBDIS(obj.StudyObject);
                AddStatFluctTBDIS(obj.StudyObject);
        end
        if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
            TBDIS_V(:,:,i) = obj.StudyObject.TBDIS;
        else
            TBDIS_V(:,i) = obj.StudyObject.TBDIS;
        end
        NormGS(i,:) = obj.StudyObject.DTNormGS;
        NormES(i,:) = obj.StudyObject.DTNormES;
        DT_P_norm(:,:,i) = [obj.StudyObject.DTexP_G.*NormGS(i,:)',obj.StudyObject.DTexP_E.*NormES(i,:)'];
        TT_P_norm(:,:,i) = [obj.StudyObject.TTexP_G.*NormGS(i,:)',obj.StudyObject.TTexP_E.*NormES(i,:)'];
        HT_P_norm(:,:,i) = [obj.StudyObject.HTexP_G.*NormGS(i,:)',obj.StudyObject.HTexP_E.*NormES(i,:)'];
    end %end Trials
    
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        TBDIS_av = mean(TBDIS_V,3);
    else
        TBDIS_av = mean(TBDIS_V,2);
    end
    
    % Reshape
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        % make sure nTrials is always last dimension before reshape!
        TBDIS_V = reshape(TBDIS_V,[obj.StudyObject.nqU*obj.StudyObject.nRings,obj.nTrials]);
    end
    
    % Compute Covariance Matrix
    obj.CovMat = cov(TBDIS_V');
    obj.MultiCovMat.CM_FSD = obj.CovMat;
    
    % Save
    save(obj.CovMatFile,'obj','TBDIS_av','TBDIS_V',...
        'DT_P_norm','TT_P_norm','HT_P_norm','-mat');
    
    % Compute Fractional CM
    obj.ComputeFracCM('Mode','CM2Frac');
    obj.MultiCovMatFrac.CM_FSD = obj.CovMatFrac;
    
    % Compute Decomposed CM
    [~, obj.CovMatFracShape] = obj.DecomposeCM('CovMatFrac',obj.CovMatFrac,'exclDataStart',1);
    obj.MultiCovMatFracShape.CM_FSD = obj.CovMatFracShape;
    obj.MultiCovMatFracNorm.CM_FSD  = obj.CovMatFracNorm;
    
    % Save again
    save(obj.CovMatFile, 'obj','-append');
    
    %% some plots for sanity check
    
    if strcmp(obj.SanityPlots,'ON')
        fig55 = figure(55);
        set(fig55, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.8, 1.2]);
        p1 = plot(obj.StudyObject.DTexE, DT_P_theory*100,'k','LineWidth',3);
        hold on;
        p2 = plot(obj.StudyObject.DTexE, DT_P_norm*100,'--','Linewidth',1.5);
        set(gca,'XScale','log')
        hold off;
        xlabel('Energy (eV)');
        ylabel('Probability (%)');
        title(sprintf('DT Final State Distribution'));
        sample_leg = sprintf('%u Sample Distributions with \nnormalization uncertainty %u%%  \nbin-to-bin uncorrelated uncertainty GS:%u%% + ES: %u %% ',...
            obj.nTrials, obj.FSDNorm_RelErr*100, obj.FSDShapeGS_RelErr*100, obj.FSDShapeGS_RelErr);
        theo_leg = sprintf('Theory Distribution (%s)',obj.StudyObject.DTFSD);
        legend([p1 p2(4)],theo_leg, sample_leg);
        PrettyFigureFormat;
        xlim([min(obj.StudyObject.DTexE) max(obj.StudyObject.DTexE)]);
        ylim([min(DT_P_theory*100) max(DT_P_theory*100)]);
        set(gca,'FontSize',18)
        
        %                 % weighted mean and weighted variance
        %                 FSDmean = zeros(obj.nTrials,1); %weighted mean energy   (eV)
        %                 FSDvar  = zeros(obj.nTrials,1);%weighted mean variance (eV^2)
        %                 for i=1:obj.nTrials
        %                     FSDmean(i) = wmean(obj.StudyObject.DTexE,DT_P_norm(:,i)');
        %                     FSDvar(i) = var(obj.StudyObject.DTexE,DT_P_norm(:,i)');
        %                 end
        %                 fig66 = figure(66);
        %                 [t1 n1 x1] = nhist(FSDmean,'color',rgb('CadetBlue'));
        %                 mean_leg1 = sprintf('%.2f eV \n',mean(x1));
        %                 mean_leg2 = sprintf('%.2f eV²',var(x1));
        %                 legend(['<\mu_E>  = ',mean_leg1,'var(\mu_E)=',mean_leg2]);
        %                 xlabel('mean energy \mu_E (eV)'); ylabel('');
        %                 title('DT Final State Distribution - Mean FSD Energy Sampling');
        %                 set(gca,'FontSize',18);
        %                 PrettyFigureFormat;
        %
        %                 fig67 = figure(67);
        %                 [t2 n2 x2] = nhist(FSDvar,'color',rgb('CadetBlue'));
        %                 var_leg1 = sprintf('%.2f eV² \n',mean(x2));
        %                 var_leg2 = sprintf('%.2f eV²',var(x2));
        %                 theo_leg = sprintf('%.2f eV²',var(obj.StudyObject.DTexE,DT_P_theory'));
        %                 legend(['\sigma_{theo} = ',theo_leg],['<\sigma²>  = ',var_leg1,'var(\sigma²)= ',var_leg2]);
        %                 xlabel('variance \sigma^2 (eV²)'); ylabel('');
        %                 title('DT Final State Distribution - FSD Variance Sampling');
        %                 set(gca,'FontSize',18);
        %                 PrettyFigureFormat;
        %                 %end
    end
end
function ComputeCM_HybridStacking(obj,varargin)
    %Stacking CM for Very First Tritium (8 Runs)
    p = inputParser;
    p.addParameter('RecomputeMTD','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('StackqUBinStart',obj.StudyObject.nqU-13,@(x)isfloat(x) && x>0);
    p.addParameter('StackqUBinStop',obj.StudyObject.nqU,@(x)isfloat(x) && x>0);
    p.addParameter('RunList',[40603 40604 40610], @(x)isfloat(x));
    
    p.parse(varargin{:});
    
    RecomputeMTD    = p.Results.RecomputeMTD;
    StackqUBinStart = p.Results.StackqUBinStart;
    StackqUBinStop  = p.Results.StackqUBinStop;
    RunList         = p.Results.RunList;
    
    % Labeling
    covmat_path =[getenv('SamakPath'),sprintf('/inputs/CovMat/HybridStacking/CM/')];
    MakeDir(covmat_path);
    covmat_name = sprintf('CMSTackingRuns%s_%uTrials',num2str(RunList),obj.nTrials);
    obj.CovMatFile = [covmat_path,covmat_name,'.mat'];
    
    %Check if CM is already computed
    if exist(obj.CovMatFile,'file')==2 && strcmp(obj.RecomputeFlag,'OFF')
        fprintf(2,'CovarianceMatrix:ComputeCM_HybridStacking: Loading Hybrid Stacking CM from File \n')
        obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','HybStack');
        obj.MultiCovMat.CM_HybStack = obj.CovMat;
        return
    end
    
    %% Create nTrial MTDs
    if strcmp(RecomputeMTD,'ON')
        VFTRuns = numel(RunList);
        qURuns = zeros(VFTRuns,obj.StudyObject.nqU); %qU from the 8 runs
        for i=1:VFTRuns
            VFT_TD = importdata(['Run',num2str(RunList(i)),'.mat']);
            qURuns(i,:) = VFT_TD.qU;
        end
        qUmin = min(qURuns)';
        qUmax = max(qURuns)';
        offset_l = qUmin - obj.StudyObject.qU;
        offset_r = qUmax - obj.StudyObject.qU;
        FT_SimMTD_multiruns('nRuns',obj.nTrials,'TDi',obj.StudyObject.TD,...%Create MTDs
            'qUoffset_l', offset_l, 'qUoffset_r',offset_r,...
            'ClearTDruns','OFF','SanityPlots','ON');
    end
    %%  Compute Spectra
    TBDIS_V  = zeros(obj.StudyObject.nqU,obj.nTrials);
    
    obj.AltObject = cell(obj.nTrials,1);
    progressbar('Compute Hybrid Stacking CM')
    for i=1:1:obj.nTrials
        progressbar(i/obj.nTrials)
        obj.AltObject{i} = ref_run40263('TD',['Run40263_run',num2str(i)]);
        obj.AltObject{i}.ComputeTBDDS;
        obj.AltObject{i}.ComputeTBDIS;
        TBDIS_V(:,i) = obj.AltObject{i}.TBDIS;
        figure(2);
        plot(obj.AltObject{1}.qU, TBDIS_V(:,i)-obj.AltObject{1}.TBDIS);
        hold on
    end
    
    CM_tmp= cov((TBDIS_V)');
    obj.CovMat = zeros(obj.StudyObject.nqU);
    obj.CovMat(StackqUBinStart:StackqUBinStop,StackqUBinStart:StackqUBinStop) = CM_tmp(StackqUBinStart:StackqUBinStop,StackqUBinStart:StackqUBinStop);
    obj.MultiCovMat.CM_HybStack = obj.CovMat;
    figure(3)
    imagesc(obj.CovMat);
    colorbar;
    
    %save
    save(obj.CovMatFile,'obj','TBDIS_V','-mat');
    
    % Compute Decomposed CM
    obj.DecomposeCM('Option','Frac');
    obj.MultiCovMatFracShape.CM_HybStack = obj.CovMatFracShape;
    obj.MultiCovMatFracNorm.CM_HybStack = obj.CovMatFracNorm;
    
    % Save again
    obj.MultiCovMatFrac.CM_HybStack = obj.CovMatFrac;
    save(obj.CovMatFile, 'obj','-append');
end
function ComputeCM_Stacking(obj,varargin)
    % Stacking Covariance Matrix: HV fluctuations only
    cprintf('blue','CovarianceMatrix:ComputeCM_Stacking: Compute Stacking Covariance Matrix  \n')
    
    %----------------------------------------- --------parser start
    p = inputParser;
    p.addParameter('SanityPlot','OFF');
    p.parse(varargin{:});
    SanityPlot      = p.Results.SanityPlot;
    %----------------------------------------- --------parser end
    
    % Labeling
    obj.GetTDlabel;
    covmat_path =[getenv('SamakPath'),sprintf('inputs/CovMat/Stacking/CM/')];
    covmat_name = sprintf('CMStack_%s_%.0fRuns_%.0fmV-qUErr_%.3f-qUfracRelErr_%.0fTrials',...
        obj.TDlabel,obj.nStack,mean(mean(obj.Stack_qUErr))*1e3,mean(mean(obj.Stack_qUfracRelErr)),obj.nTrials);
    
    if contains(covmat_name,' ')
        covmat_name=strrep(covmat_name,'  ','_');
    end
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        % add ring information for ringwise covariance matrices
        covmat_name = [covmat_name,sprintf('_Ring%s',obj.StudyObject.FPD_RingMerge)];
    end
    obj.CovMatFile = [covmat_path,covmat_name,'.mat'];
    
    %Check if CM is already computed
    if exist(obj.CovMatFile,'file')==2 && strcmp(obj.RecomputeFlag,'OFF')
        fprintf(2,'CovarianceMatrix:ComputeCM_Stacking: Loading Stacking Covariance Matrix from File \n')
        obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','Stack');
        obj.MultiCovMat.CM_Stack = obj.CovMat;
        return
    end
    
    %% calculate covariance matrix
    obj.StudyObject.ComputeTBDDS;
    obj.StudyObject.ComputeTBDIS;
    
    % vary qU and qU-frac
    nRings = numel(obj.StudyObject.MACE_Ba_T);
    qUSamples     = obj.StudyObject.qU+squeeze(randn(obj.StudyObject.nqU,nRings,obj.nTrials,obj.nStack)).*obj.Stack_qUErr;
    %qUfracSamples = obj.StudyObject.qUfrac.*(1+squeeze(randn(obj.StudyObject.nqU,nRings,obj.nTrials,obj.nStack)).*obj.Stack_qUfracRelErr);
    %qUfracSamples = qUfracSamples./sum(qUfracSamples);
    
    % mean integral spectrum (rate!) without background
    meanTime = obj.StudyObject.qUfrac.*obj.StudyObject.TimeSec;
    MeanRate = obj.StudyObject.TBDIS./meanTime-obj.StudyObject.BKG_RateSec;
    
    % compute integral spectra
    TBDIS_SingleRun  = squeeze(zeros(obj.StudyObject.nqU,nRings,obj.nTrials,obj.nStack));
    progressbar('compute stacking CM')
    for i=1:obj.nTrials
        %interpolate rate not counts!
        progressbar(i/obj.nTrials)
        for s=1:obj.nStack
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                for r=1:nRings
                    TBDIS_SingleRun(:,r,i,s) = interp1(obj.StudyObject.qU(:,r),MeanRate(:,r),...
                        qUSamples(:,r,i,s),'spline','extrap');
                end
            else
                TBDIS_SingleRun(:,i,s) = interp1(obj.StudyObject.qU,MeanRate,...
                    qUSamples(:,i,s),'spline','extrap');
            end
        end
    end
    
    TBDIS_SingleRun = TBDIS_SingleRun.*meanTime; %convert back to counts
    
    % normalize to mean time
    TBDIS_SingleRun = TBDIS_SingleRun./obj.nStack;%.*(obj.StudyObject.qUfrac./qUfracSamples);
    
    % stack spectra (runwise)
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        TBDIS_V = sum(TBDIS_SingleRun,4);
        TBDIS_V = reshape(TBDIS_V,[obj.StudyObject.nqU*nRings,obj.nTrials]);
    else
        TBDIS_V = sum(TBDIS_SingleRun,3);
    end
    
    if contains(obj.StudyObject.TD,'KNM1')
        qUrange = 2:obj.StudyObject.nqU;
        TBDIS_V(1,:) = 0;
    else
        qUrange = 1:obj.StudyObject.nqU;
    end
    
    % calculate sample covariance matrix
    obj.CovMat = cov(TBDIS_V');
    obj.MultiCovMat.CM_Stack = obj.CovMat;
    
    %save
    MakeDir(covmat_path);
    save(obj.CovMatFile,'obj','TBDIS_V','-mat');
    
    %Compute Fractional CM
    obj.ComputeFracCM('Mode','CM2Frac');
    obj.MultiCovMatFrac.CM_Stack = obj.CovMatFrac;
    
    % Compute Decomposed CM
    [~, obj.CovMatFracShape] = obj.DecomposeCM('CovMatFrac',obj.CovMatFrac,'exclDataStart',1);
    obj.MultiCovMatFracShape.CM_Stack = obj.CovMatFracShape;
    
    % Save again
    save(obj.CovMatFile, 'obj','-append');
    
    %% sanity Plots
    switch SanityPlot
        case 'ON'
            
            f1 = figure(32);
            set(f1, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.7, 0.7]);
            [l,a] = boundedline(obj.StudyObject.qU(qUrange)-18575,mean(TBDIS_V(qUrange,:),2),1000.*std(TBDIS_V(qUrange,:)')');
            l.LineWidth =3; legend([l,a],'mean','1 \sigma band x 1000'); legend boxoff
            l.Color = rgb('DarkRed');
            a.FaceColor=rgb('IndianRed');a.FaceAlpha=0.5;
            PrettyFigureFormat;
            xlabel('retarding potential - 18575 (V)');
            ylabel('stacked counts');
            ylim([0 max(mean(TBDIS_V(qUrange,:),2)+100.*std(TBDIS_V(qUrange,:)')')]);
            xlim([obj.StudyObject.qU(qUrange(1))-18575,obj.StudyObject.qU(qUrange(end))-18575]);
            set(gca,'FontSize',22);
            publish_figurePDF(f1,sprintf('StackingCM_%s_spectrum.pdf',obj.StudyObject.TD));
            
            
            f2 = figure(33);
            subrun = 14;
            set(f2, 'Units', 'normalized', 'Position', [0.9, 0.9, 1, 0.6]);
            subplot(1,2,1);
            h= histogram(qUSamples(subrun,:,:)-18575); h.FaceColor=rgb('IndianRed');
            PrettyFigureFormat;
            xlabel('qU (eV) samples')
            leg = legend(sprintf('\\sigma = %.0f mV',1e3*std(reshape(qUSamples(subrun,:,:),obj.nTrials.*obj.nStack,1))));
            legend boxoff;
            
            subplot(1,2,2);
            h=histogram(100*qUfracSamples(subrun,:,:));  h.FaceColor=rgb('IndianRed');
            PrettyFigureFormat;
            xlabel('qUfrac (%) samples')
            leg = legend(sprintf('\\sigma = %.0f %%',1e2.*obj.Stack_qUfracRelErr(subrun)));
            legend boxoff;
            
            a=annotation('textbox', [0.35 0.88 1 0.12], ...
                'String', sprintf('HV and time fluctuations within subrun %.0f',subrun), ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'left');
            a.FontSize=20;a.FontWeight='bold';
            
            print(f2,sprintf('StackingCM_%s_qU.png',obj.StudyObject.TD),'-dpng','-r450');
            
    end
end
function ComputeCM_STAT(obj,varargin)
    fprintf('--------------------------------------------------------------------------\n')
    cprintf('blue','CovarianceMatrix:ComputeCM_STAT: Compute Statistics Covariance Matrix  \n')
    obj.StudyObject.ComputeTBDDS; obj.StudyObject.ComputeTBDIS;
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        TBDIS = reshape(obj.StudyObject.TBDIS,[obj.StudyObject.nqU*obj.StudyObject.nRings,1]);
    else
        TBDIS = obj.StudyObject.TBDIS;
    end
    obj.MultiCovMat.CM_STAT = diag(TBDIS);
    obj.MultiCovMatFrac.CM_STAT = (obj.MultiCovMat.CM_STAT./TBDIS)./TBDIS';
end

function ComputeCM_LongPlasma(obj,varargin)
    % longitudinal plasma uncertainty
    % effectively described by e-loss shift and variance
    p=inputParser;
    p.addParameter('CorrCoeff',0,@(x)isfloat(x));
    p.addParameter('NegSigma','Troitsk',@(x)ismember(x,{'Abs','Troitsk'}));
    p.addParameter('SanityPlot','OFF');
    p.parse(varargin{:});
    CorrCoeff  = p.Results.CorrCoeff;
    NegSigma   = p.Results.NegSigma;    % how to deal with negative sigmas
    SanityPlot = p.Results.SanityPlot;

    switch NegSigma
        case 'Troitsk'
            NegSigmaStr = '_Troitsk'; % troitsk formula
        case 'Abs'
            NegSigmaStr = '';         % absolute value of sigma
    end
    nRings = numel(obj.StudyObject.MACE_Ba_T);
    
    % initial longi plasma parameters
    is_EOffset_i = obj.StudyObject.is_EOffset;
    MACE_Sigma_i = mean(obj.StudyObject.MACE_Sigma);
    
    obj.GetTDlabel;
    covmat_path =[getenv('SamakPath'),sprintf('inputs/CovMat/LongPlasma/CM/')];
    MakeDir(covmat_path);
    covmat_filename = sprintf('LongPlasma_%s_%.0fTrials_%.0fmeVEloss_%.0fmeVElossErr_%.0fmeVSigma_%.0fmeVSigmaErr_%.0fCorrCoeff%s.mat',...
        obj.TDlabel,obj.nTrials,...
        is_EOffset_i*1e3,obj.is_EOffsetErr*1e3,...
        MACE_Sigma_i*1e3,sqrt(obj.MACE_VarErr).*1e3,...
        CorrCoeff*100,NegSigmaStr);
    
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        % add ring information for ringwise covariance matrices
        covmat_filename = strrep(covmat_filename,'.mat',sprintf('_Ring%s.mat',obj.StudyObject.FPD_RingMerge));
    end
    obj.CovMatFile = strcat(covmat_path,covmat_filename);
    
    loadSuccess = 0;
    % load if already computed
    if exist(obj.CovMatFile,'file')==2 && strcmp(obj.RecomputeFlag,'OFF')
        fprintf(2,'CovarianceMatrix:ComputeCM_LongPlasma: Loading LongPlasma CM from File \n')
        obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','LongPlasma');
        obj.MultiCovMat.CM_LongPlasma = obj.CovMat; %in case normalization was changed
        loadSuccess = 1;
        if strcmp(SanityPlot,'OFF')
            return
        end
    end
    
    if loadSuccess==0
        
        recompRF_i = obj.StudyObject.recomputeRF;
        RFBinStep_i = obj.StudyObject.RFBinStep;
        Sync_i      =  obj.StudyObject.SynchrotronFlag;
        AngularTFFlag_i = obj.StudyObject.AngularTFFlag;
        obj.StudyObject.RFBinStep = 0.04;
        obj.StudyObject.SynchrotronFlag = 'OFF';
        obj.StudyObject.AngularTFFlag = 'OFF';
        obj.StudyObject.recomputeRF   = 'OFF';
        
        % longi plasma uncertainty
        PlasmaErr =  [obj.MACE_VarErr,obj.is_EOffsetErr];
        PlasmaCovMat      = PlasmaErr.*[1,CorrCoeff;CorrCoeff,1].*PlasmaErr';
        
        % randomize longi plasma parameters
        PlasmaPar_expected = [MACE_Sigma_i.^2,is_EOffset_i];
        PlasmaPar =  mvnrnd(PlasmaPar_expected,PlasmaCovMat,obj.nTrials);
        
        MACE_Var_v = PlasmaPar(:,1);                 % negative broadenings dealt with later
        MACE_Var_v(sqrt(abs(MACE_Var_v))<1e-3) = 0;  % too small to resolve anyway -> saves from time
        is_EOffset_v  = PlasmaPar(:,2);
        
        % e-loss shift:
        ELossRange = 500;
        ELossBinStep =0.2;
        maxE = ELossRange;
        minE=-maxE; NbinE = (maxE-minE)/ELossBinStep;
        E = linspace(minE,maxE,NbinE);
        [ElossFunc,~] = obj.StudyObject.ComputeELossFunction('E',E); % load if already exists, otherwise compute
        
        if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
            TBDIS_V = zeros(obj.StudyObject.nqU,obj.nTrials,obj.StudyObject.nRings);
        else
            TBDIS_V = zeros(obj.StudyObject.nqU,obj.nTrials);
        end
        
        % calculate response functions with varied parameters
        RFfun = @ComputeRF;
        ResponseFunction = zeros(obj.StudyObject.nTe,obj.StudyObject.nqU,obj.nTrials,nRings);
        ElossFunc_v       = zeros(obj.StudyObject.NIS,numel(E),obj.nTrials);
        obj.StudyObject.ElossFunctions = []; % reset

        progressbar('Compute longitudinal plasma cov mat')
        for i=1:obj.nTrials
            progressbar(i/obj.nTrials);
            % eloss shift
            ElossFunc_v(:,:,i) =  cell2mat(cellfun(@(x) x(E-is_EOffset_v(i)), ElossFunc,'UniformOutput',0));
            for ri = 1:nRings
                for ii = 1:obj.StudyObject.nqU
                    ResponseFunction(:,ii,i,ri) = RFfun(obj.StudyObject,obj.StudyObject.Te,obj.StudyObject.qU(ii,ri),'pixel',ri,...
                        'ELossRange',ELossRange,'ELossBinStep',ELossBinStep,...
                        'ElossFunctions',ElossFunc_v(:,:,i),...
                        'RFBinStep',0.04);
                end
            end
            obj.StudyObject.RF = squeeze(ResponseFunction(:,:,i,:));
            
            % variance
            obj.StudyObject.LoadFSD('Sigma',sqrt(abs(MACE_Var_v(i)))); % take absolute sigma
            obj.StudyObject.ComputeTBDDS;
            obj.StudyObject.ComputeTBDIS;
            TBDIS_V(:,i,:) = obj.StudyObject.TBDIS;
            
            if strcmp(NegSigma,'Troitsk') % Use Troitsk formula
                % TBDIS(sigma<0) = 2*(sigma=0)-(abs(sigma))
                if MACE_Var_v(i)<0
                    obj.StudyObject.LoadFSD; % sigma = 0
                    obj.StudyObject.ComputeTBDDS;
                    obj.StudyObject.ComputeTBDIS;
                    TBDIS_V(:,i,:) = 2.*obj.StudyObject.TBDIS - squeeze(TBDIS_V(:,i,:));
                end
            end
        end
        
        % Reshape
        if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
            % make sure nTrials is always last dimension before reshape!
            TBDIS_V = permute(TBDIS_V,[1 3 2]);
            TBDIS_V = reshape(TBDIS_V,[obj.StudyObject.nqU*obj.StudyObject.nRings,obj.nTrials]);
        end
        
        % Compute Covariance Matrix
        obj.CovMat = cov(TBDIS_V');
        obj.MultiCovMat.CM_LongPlasma = obj.CovMat;
        
        TBDIS_av = mean(TBDIS_V);
        ResponseFunction_mean = mean(ResponseFunction,3);
        ResponseFunction_std = std(ResponseFunction,0,3);
        ElossFunc_mean       = mean(ElossFunc_v,3);
        ElossFunc_std        = std(ElossFunc_v,0,3);
        
        % Save
        save(obj.CovMatFile,'obj','TBDIS_av','TBDIS_V',...
            'MACE_Var_v','is_EOffset_v',...
            'ResponseFunction_mean','ResponseFunction_std',...
            'ElossFunc_mean','ElossFunc_std','-mat');
        
        % Compute Fractional CM
        obj.ComputeFracCM('Mode','CM2Frac');
        obj.MultiCovMatFrac.CM_LongPlasma = obj.CovMatFrac;
        
        % Compute Decomposed CM
        [~, obj.CovMatFracShape] = obj.DecomposeCM('CovMatFrac',obj.CovMatFrac,'exclDataStart',1);
        obj.MultiCovMatFracShape.CM_LongPlasma = obj.CovMatFracShape;
        obj.MultiCovMatFracNorm.CM_LongPlasma  = obj.CovMatFracNorm;
        
        % Save again
        save(obj.CovMatFile, 'obj','-append');
        obj.StudyObject.ElossFunctions  = []; % reset e-loss functions
        obj.StudyObject.RFBinStep       =  RFBinStep_i ;
        obj.StudyObject.SynchrotronFlag =  Sync_i ;
        obj.StudyObject.recomputeRF     = recompRF_i;
        obj.StudyObject.AngularTFFlag   = AngularTFFlag_i;
    end
    
    if strcmp(SanityPlot,'ON')
        
        if loadSuccess==1
            d = importdata(obj.CovMatFile);
            MACE_Var_v = d.MACE_Var_v;
            is_EOffset_v = d.is_EOffset_v;
        end
        fLP = figure('Units','normalized','Position',[0.1,0.1,0.75,0.5]);
        subplot(1,2,1);
        h1 = histogram(MACE_Var_v);
        PrettyFigureFormat('FontSize',24);
        xlabel(sprintf('Broadening \\sigma^2 (eV^2)'));
        ylabel('Occurrence');
        h1.FaceAlpha = 1;
        h1.FaceColor = rgb('PowderBlue');
        
        subplot(1,2,2);
        h2 = histogram(is_EOffset_v);
        PrettyFigureFormat('FontSize',24);
        xlabel(sprintf('Energy-loss shift \\Delta\\epsilon (eV)'));
        ylabel('Occurrence');
        h2.FaceAlpha = 1;
        h2.FaceColor = rgb('PowderBlue');
        
        plotdir = './plots/';
        MakeDir(plotdir);
        plotname = [plotdir,strrep(covmat_filename,'.mat','.pdf')];
        export_fig(plotname);
        fprintf('save plot to %s \n',plotname)
    end
end
function ComputeCM_PlasmaOffsets(obj,varargin)
    % longitudinal plasma uncertainty
    % effectively described by e-loss shift and variance
    p=inputParser;
    p.addParameter('NegSigma','Troitsk',@(x)ismember(x,{'Abs','Troitsk'}));
    p.parse(varargin{:});
    NegSigma = p.Results.NegSigma;    % how to deal with negative sigmas
    
    switch NegSigma
        case 'Troitsk'
            NegSigmaStr = '_Troitsk'; % troitsk formula
        case 'Abs'
            NegSigmaStr = '';         % absolute value of sigma
    end
    
    nRings = numel(obj.StudyObject.MACE_Ba_T);
   
    obj.GetTDlabel;
    covmat_path =[getenv('SamakPath'),sprintf('inputs/CovMat/LongPlasma/CM/')];
    MakeDir(covmat_path);
    covmat_filename = sprintf('PlasmaOffsets_%s_%.0fTrials_%.2feVmeanE0_%.0fmeVmeanE0Err%s.mat',...
        obj.TDlabel,obj.nTrials,...
        mean(obj.E0Offsets),...   % average endpoint
        mean(obj.E0OffsetsErr)*1e3,...
        NegSigmaStr); % average endpoint uncertainty
    
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        % add ring information for ringwise covariance matrices
        covmat_filename = strrep(covmat_filename,'.mat',sprintf('_Ring%s.mat',obj.StudyObject.FPD_RingMerge));
    end
    obj.CovMatFile = strcat(covmat_path,covmat_filename);
    
    % load if already computed
    if exist(obj.CovMatFile,'file')==2 && strcmp(obj.RecomputeFlag,'OFF')
        fprintf(2,'CovarianceMatrix:ComputeCM_PlasmaOffsets: Loading PlasmaOffsets CM from File \n')
        obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','PlasmaOffsets');
        obj.MultiCovMat.CM_PlasmaOffsets = obj.CovMat; %in case normalization was changed
        return
    end
    
    % randomize endpoints 
    E0Offset_v = (obj.E0Offsets+obj.E0OffsetsErr.*randn(obj.nTrials,1))';
    FSDSigma_v = std(E0Offset_v);           
    FSDSigma_v(abs(FSDSigma_v)<1e-3) = 0;  % too small to resolve anyway
    
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        TBDIS_V = zeros(obj.StudyObject.nqU,obj.nTrials,obj.StudyObject.nRings);
    else
        TBDIS_V = zeros(obj.StudyObject.nqU,obj.nTrials);
    end
    
    progressbar('Compute Plasma Offsets cov mat')
    for i=1:obj.nTrials
        progressbar(i/obj.nTrials);

        % variance
        obj.StudyObject.LoadFSD('Sigma',abs(FSDSigma_v(i))); % take absolute sigma
        obj.StudyObject.ComputeTBDDS;
        obj.StudyObject.ComputeTBDIS;
        TBDIS_V(:,i,:) = obj.StudyObject.TBDIS;
        
        if strcmp(NegSigma,'Troitsk') % Use Troitsk formula
            % TBDIS(sigma<0) = 2*(sigma=0)-(abs(sigma))
            if FSDSigma_v(i)<0
                obj.StudyObject.LoadFSD; % sigma = 0
                obj.StudyObject.ComputeTBDDS;
                obj.StudyObject.ComputeTBDIS;
                TBDIS_V(:,i,:) = 2.*obj.StudyObject.TBDIS - squeeze(TBDIS_V(:,i,:));
            end
        end 
    end
    
    % Reshape
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        % make sure nTrials is always last dimension before reshape!
        TBDIS_V = permute(TBDIS_V,[1 3 2]);
        TBDIS_V = reshape(TBDIS_V,[obj.StudyObject.nqU*obj.StudyObject.nRings,obj.nTrials]);
    end
    
    % Compute Covariance Matrix
    obj.CovMat = cov(TBDIS_V');
    obj.MultiCovMat.CM_PlasmaOffsets = obj.CovMat;
    
    TBDIS_av = mean(TBDIS_V);

    % Save
    save(obj.CovMatFile,'obj','TBDIS_av','TBDIS_V','FSDSigma_v','-mat');
    
    % Compute Fractional CM
    obj.ComputeFracCM('Mode','CM2Frac');
    obj.MultiCovMatFrac.CM_PlasmaOffsets = obj.CovMatFrac;
    
    % Compute Decomposed CM
    [~, obj.CovMatFracShape] = obj.DecomposeCM('CovMatFrac',obj.CovMatFrac,'exclDataStart',1);
    obj.MultiCovMatFracShape.CM_PlasmaOffsets = obj.CovMatFracShape;
    obj.MultiCovMatFracNorm.CM_PlasmaOffsets  = obj.CovMatFracNorm;
    
    % Save again
    save(obj.CovMatFile, 'obj','-append');

end
function ComputeCM_RW(obj,varargin)
    % Goal:  Rear wall potential systematics
    %
    % Method:
    % FSDs are replaced with superposition of gaussians or rectangles or 3 Gaussians
    % sigma       == Width of gaussians or width of rectangle
    % sigmaErr    == Absolute uncertainty on sigma
    %
    % Lisa,  December 2019
    %
    % ---------------------------------parser start  ----------------------------------------------
    p=inputParser;
    p.addParameter('Sigma',obj.StudyObject.FSD_Sigma,@(x)isfloat(x));
    p.addParameter('SanityPlot','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('MultiPos',obj.StudyObject.FSD_MultiPos,@(x) isfloat(x) || isempty(x));     %3 gaussians instead of using 1 gaussian per energy (for 3 RW settings)
    p.addParameter('MultiWeights',obj.StudyObject.FSD_MultiWeights,@(x) isfloat(x) || isempty(x)); %3 gaussians instead of using 1 gaussian per energy
    p.addParameter('Dist',obj.StudyObject.FSD_Dist,@(x)ismember(x,{'Gauss','Rect'}));
    p.addParameter('BinningFactor','',@(x) isfloat(x) || isempty(x)); % multiplies number of bins in rebinning
    
    p.parse(varargin{:});
    
    Sigma              = p.Results.Sigma;
    SanityPlot         = p.Results.SanityPlot;
    MultiPos           = p.Results.MultiPos;
    MultiWeights       = p.Results.MultiWeights;
    BinningFactor      = p.Results.BinningFactor;
    Dist               = p.Results.Dist;
    
    % ---------------------------------parser end  ----------------------------------------------
    % label
    obj.GetTDlabel;
    covmat_path =[getenv('SamakPath'),sprintf('inputs/CovMat/RearWall/CM/')];
    MakeDir(covmat_path);
    cm_filename = sprintf('RearWall_%s_%s_%.0fmeVSigma_%.0fmeVSigmaErr_%.0fTrials.mat',...
        obj.TDlabel,Dist,Sigma*1e3,obj.RW_SigmaErr*1e3,obj.nTrials);
    
    if BinningFactor>1 % add binning factor
        cm_filename = strrep(cm_filename,'.mat',sprintf('_%.0fBinFactor.mat',BinningFactor));
    end
    
    if ~isempty(MultiPos) % add multi positions
        cm_filename = strrep(cm_filename,'.mat',sprintf('_%.0fmeV-%.0fmeV-%.0fmeVMultiPos_%.0fmeVPosErr.mat',...
            MultiPos(1)*1e3,MultiPos(2)*1e3,MultiPos(3)*1e3,obj.RW_MultiPosErr*1e3));
    end
    
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        % add ring information for ringwise covariance matrices
        cm_filename = strrep(cm_filename,'.mat',sprintf('_Ring%s.mat',obj.StudyObject.FPD_RingMerge));
    end
    obj.CovMatFile = strcat(covmat_path,cm_filename);
    
    if exist(obj.CovMatFile,'file')
        fprintf(2,'CovarianceMatrix:ComputeCM_RearWall: Loading RearWall CM from File \n')
        obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','RW');
        obj.MultiCovMat.CM_RW = obj.CovMat; %in case normalization was changed
        return
    end
    
    % Common FSD argument for all samples
    FSDarg = {'Dist',Dist,...
        'MultiWeights',MultiWeights,...
        'BinningFactor',BinningFactor,...
        'SanityPlot',SanityPlot};
    
    % Vary parameter according to uncertainty: sigma
    Sigma_V = Sigma + obj.RW_SigmaErr.*randn(obj.nTrials,1);
    Sigma_V(Sigma_V<0.001) = 0.001; % no negative broadening possible
    
    if ~isempty(MultiPos) % add multi positions uncertainty
        MultiPos_V = MultiPos+obj.RW_MultiPosErr.*randn(obj.nTrials,1);
        Sigma_V = 0.001.*ones(obj.nTrials,1); % small for multi-gauss
    else
        MultiPos_V = '';
    end
    
    TBDIS_V = zeros(obj.StudyObject.nqU,obj.nTrials);
    
    % Calculate integral spectra with varied FSDs
    progressbar('Compute CM RW');
    for i=1:obj.nTrials
        progressbar(i/obj.nTrials)
        if ~isempty(MultiPos) % add multi positions
            obj.StudyObject.LoadFSD(FSDarg{:},'MultiPos',MultiPos_V(i,:),'Sigma',Sigma_V(i));
        else
            obj.StudyObject.LoadFSD(FSDarg{:},'MultiPos',MultiPos,'Sigma',Sigma_V(i));
        end
        
        obj.StudyObject.ComputeTBDDS;
        obj.StudyObject.ComputeTBDIS;
        TBDIS_V(:,i) = obj.StudyObject.TBDIS;
    end
    
    % Reshape
    
    if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
        % make sure nTrials is always last dimension before reshape!
        TBDIS_V = reshape(TBDIS_V,[obj.StudyObject.nqU*obj.StudyObject.nRings,obj.nTrials]);
    end
    
    % Compute Covariance Matrix
    obj.CovMat = cov(TBDIS_V');
    obj.MultiCovMat.CM_RW = obj.CovMat;
    
    % Save
    save(obj.CovMatFile,'obj','TBDIS_V','Sigma_V','MultiPos_V','-mat');
    
    % Compute Fractional CM
    obj.ComputeFracCM('Mode','CM2Frac');
    obj.MultiCovMatFrac.CM_RW = obj.CovMatFrac;
    
    % Compute Decomposed CM
    [~, obj.CovMatFracShape] = obj.DecomposeCM('CovMatFrac',obj.CovMatFrac,'exclDataStart',1);
    obj.MultiCovMatFracShape.CM_RW = obj.CovMatFracShape;
    obj.MultiCovMatFracNorm.CM_RW  = obj.CovMatFracNorm;
    
    % Save again
    save(obj.CovMatFile, 'obj','-append');
    
end
    end
    
    methods %combi cov mat
        function ComputeCM(obj,varargin)
            p=inputParser;
            p.addParameter('PlotSaveCM','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SysBudget_Label','',@(x)ischar(x));
            p.parse(varargin{:});
            PlotSaveCM = p.Results.PlotSaveCM;
            SysBudget_Label = p.Results.SysBudget_Label;
            
            fprintf('--------------------------------------------------------------------------\n')
            cprintf('blue','CovarianceMatrix:ComputeCM: Compute Combi Covariance Matrix  \n')
           %  obj.nTrials = 5000;
            %Labeling
            combi_path= [getenv('SamakPath'),sprintf('/inputs/CovMat/Combi/')];
            effects_logic = structfun(@(x)strcmp(x,'ON'),obj.SysEffect);
            fields = fieldnames(obj.SysEffect);
            mytitle = '';
            for i=1:numel(fields)
                if effects_logic(i)==1
                    mytitle = strcat(mytitle,'-',fields(i));
                end
            end
            IsoName = '';
            if strcmp(obj.DTFlag,'ON')
                IsoName = [IsoName,'DT-'];
            end
            if strcmp(obj.HTFlag,'ON')
                IsoName = [IsoName,'HT-'];
            end
            if strcmp(obj.TTFlag,'ON')
                IsoName = [IsoName,'TT-'];
            end
            obj.GetTDlabel;
            mytitle{:} = strrep(mytitle{:},'-RF_EL-RF_BF-RF_RX',['RF_',num2str(obj.WGTS_CD_MolPerCm2_RelErr),'CDErr']);
            mytitle{:} = strrep(mytitle{:},'FSD',['FSD_',IsoName,num2str(obj.FSDNorm_RelErr),'NormErr',num2str(obj.FSDShapeGS_RelErr),'ShapeGSErr',num2str(obj.FSDShapeES_RelErr),'ShapeESErr']);
            mytitle{:} = strrep(mytitle{:},'TASR',['TASR_',num2str(mean(obj.WGTS_TASR_RelErr)),'Err']);
            combi_name = sprintf('CombiCM_%s_%u-Trials_%s.mat', obj.TDlabel,obj.nTrials,mytitle{:});
            obj.CovMatFile = strcat(combi_path,combi_name);
            
            %             %Import Option
            %              if exist(obj.CovMatFile,'file')==2 && strcmp(obj.RecomputeFlag,'OFF')
            %                  fprintf('Loading Combi CM from File \n')
            %                  obj.ReadCMFile('filename',obj.CovMatFile,'SysEffect','Combi');
            %                  return
            %              end
            
            %Init
            if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                CovMatFracCombi = zeros(obj.StudyObject.nqU*obj.StudyObject.nRings,obj.StudyObject.nqU*obj.StudyObject.nRings);
            else
                CovMatFracCombi = zeros(obj.StudyObject.nqU,obj.StudyObject.nqU);
            end
            
            %Load / Compute CM of SysEffects
            %% Response Function
           % obj.nTrials = 1000;
            if strcmp(obj.SysEffect.RF_EL,'ON') && strcmp(obj.SysEffect.RF_BF,'ON') && strcmp(obj.SysEffect.RF_RX,'ON') % all RF Effects ON
                %all 'ON'
                obj.ComputeCM_RF;
                CovMatFracCombi = obj.MultiCovMatFrac.CM_RF;
                if strcmp(PlotSaveCM,'ON')
                    obj.PlotCM('PlotEffect','Response Function (all)','savePlot','ON','savename',SysBudget_Label);
                end
            elseif strcmp(obj.SysEffect.RF_EL,'ON') || strcmp(obj.SysEffect.RF_BF,'ON') || strcmp(obj.SysEffect.RF_RX,'ON')
                %one or two 'ON'
                nRings = numel(obj.StudyObject.MACE_Ba_T);
                obj.MultiCovMatFrac.CM_RF_EL = zeros(obj.StudyObject.nqU*nRings);
                obj.MultiCovMatFrac.CM_RF_BF = zeros(obj.StudyObject.nqU*nRings);
                obj.MultiCovMatFrac.CM_RF_RX = zeros(obj.StudyObject.nqU*nRings);
                obj.ComputeCM_RF;
                a1 = strcmp(obj.SysEffect.RF_EL,'ON');  a2=strcmp(obj.SysEffect.RF_BF,'ON'); a3 =strcmp(obj.SysEffect.RF_RX,'ON');
                if sum([a1,a2,a3])==1 %if only 1 is ON
                    CovMatFracCombi = obj.MultiCovMatFrac.CM_RF_EL + obj.MultiCovMatFrac.CM_RF_BF + obj.MultiCovMatFrac.CM_RF_RX;
                else  % if two are ON
                    CovMatFracCombi = obj.CovMatFrac;
                end
                if strcmp(PlotSaveCM,'ON')
                    BF = '';EL = '';RX = '';
                    if strcmp(obj.SysEffect.RF_BF,'ON') BF='B-Field'; end
                    if strcmp(obj.SysEffect.RF_EL,'ON') EL='Energy Loss'; end
                    if strcmp(obj.SysEffect.RF_RX,'ON') RX='CD x-section'; end
                    obj.PlotCM('PlotEffect',sprintf('Response Function %s %s %s',BF,EL,RX),'savePlot','ON','savename',SysBudget_Label);
                end
            end
            
          %   obj.nTrials = 5000;
            %% FSD   
            if strcmp(obj.SysEffect.FSD,'ON')
                obj.ComputeCM_FSD;
                CovMatFracCombi = CovMatFracCombi + obj.MultiCovMatFrac.CM_FSD;
                if strcmp(PlotSaveCM,'ON')
                    obj.PlotCM('PlotEffect','Final States','savePlot','ON','savename',SysBudget_Label);
                end
            end
            %% TASR
            
            if strcmp(obj.SysEffect.TASR,'ON')
                obj.ComputeCM_TASR;
                CovMatFracCombi = CovMatFracCombi + obj.MultiCovMatFrac.CM_TASR;
                if strcmp(PlotSaveCM,'ON')
                    obj.PlotCM('PlotEffect','Tritium Activity Fluctuations','savePlot','ON','savename',SysBudget_Label);
                end
            end
            %% Background: should be used seperately
            if strcmp(obj.SysEffect.BM1S,'ON')
                obj.ComputeCM_BM1S;
                CovMatFracCombi = CovMatFracCombi + obj.MultiCovMatFrac.CM_BM1S;
                if strcmp(PlotSaveCM,'ON')
                    obj.PlotCM('PlotEffect','Background','savePlot','ON','savename',SysBudget_Label);
                end
            end
            %% Doppler Effect ON/OFF
            if strcmp(obj.SysEffect.DOPoff,'ON')
                obj.ComputeCM_DOPoff;
                CovMatFracCombi = CovMatFracCombi + obj.MultiCovMatFrac.CM_DOPoff;
                if strcmp(PlotSaveCM,'ON')
                    obj.PlotCM('PlotEffect','Doppler Effect','savePlot','ON','Convergence','OFF','savename',SysBudget_Label);
                end
            end
            %% Theoretical Corrections (TC)
            if strcmp(obj.SysEffect.TCoff_RAD,'ON') &&  strcmp(obj.SysEffect.TCoff_RAD,'ON')
                %both ON
                obj.ComputeCM_TCoff;
                CovMatFracCombi = CovMatFracCombi + obj.MultiCovMatFrac.CM_TCoff;
                if strcmp(PlotSaveCM,'ON')
                    obj.PlotCM('PlotEffect','Theoretical Corrections','savePlot','ON','Convergence','OFF','savename',SysBudget_Label);
                end
            elseif strcmp(obj.SysEffect.TCoff_RAD,'ON') ||  strcmp(obj.SysEffect.TCoff_OTHER,'ON')
                %only one ON
                if strcmp(obj.StudyObject.FPD_Segmentation,'RING')
                    dim = obj.StudyObject.nqU*obj.StudyObject.nRings;
                else
                    dim = obj.StudyObject.nqU;
                end
                obj.MultiCovMatFrac.CM_TCoff_RAD = zeros(dim);
                obj.MultiCovMatFrac.CM_TCoff_OTHER = zeros(dim);
                obj.ComputeCM_TCoff;
                CovMatFracCombi = CovMatFracCombi + obj.MultiCovMatFrac.CM_TCoff_RAD + obj.MultiCovMatFrac.CM_TCoff_OTHER;
                if strcmp(PlotSaveCM,'ON')
                    obj.PlotCM('PlotEffect','Theoretical Corrections','savePlot','ON','Convergence','OFF','savename',SysBudget_Label);
                end
            end
            %% Stacking CM
            if strcmp(obj.SysEffect.Stack,'ON')
                obj.ComputeCM_Stacking;
                CovMatFracCombi = CovMatFracCombi + obj.MultiCovMatFrac.CM_Stack;
                if strcmp(PlotSaveCM,'ON')
                    obj.PlotCM('PlotEffect','Stacking','savePlot','ON','savename',SysBudget_Label);
                end
            end
            %% FPD efficiency
            if strcmp(obj.SysEffect.FPDeff,'ON')
                obj.ComputeCM_FPDeff;
                CovMatFracCombi = CovMatFracCombi + obj.MultiCovMatFrac.CM_FPDeff;
                if strcmp(PlotSaveCM,'ON')
                    obj.PlotCM('PlotEffect','FPDeff','savePlot','ON','savename',SysBudget_Label,...
                        'ConvergenceTest','OFF');
                end
            end
            
            %% Plasma time evolution
            if strcmp(obj.SysEffect.RW,'ON')
                obj.ComputeCM_RW;
                CovMatFracCombi = CovMatFracCombi + obj.MultiCovMatFrac.CM_RW;
                if strcmp(PlotSaveCM,'ON')
                    obj.PlotCM('PlotEffect','RW','savePlot','ON','savename',SysBudget_Label,...
                        'ConvergenceTest','OFF');
                end
            end
            
            %% Plasma longitudinal
            if strcmp(obj.SysEffect.LongPlasma,'ON')
                nTrials_i = obj.nTrials;
                if numel(obj.StudyObject.MACE_Ba_T)>1 && obj.nTrials>1000% multiring
                    obj.nTrials = 1000;
                end
                obj.ComputeCM_LongPlasma;
                CovMatFracCombi = CovMatFracCombi + obj.MultiCovMatFrac.CM_LongPlasma;
                if strcmp(PlotSaveCM,'ON')
                    obj.PlotCM('PlotEffect','LongPlasma','savePlot','ON','savename',SysBudget_Label,...
                        'ConvergenceTest','OFF');
                end
                
                obj.nTrials = nTrials_i;
            end
            
            %%Combined Fractional CM
            obj.CovMatFrac = CovMatFracCombi;
            
            %Labeling (one more time, because overwritten in the meantime)
            combi_path= [getenv('SamakPath'),sprintf('/inputs/CovMat/Combi/')];
            effects_logic = structfun(@(x)strcmp(x,'ON'),obj.SysEffect);
            fields = fieldnames(obj.SysEffect);
            mytitle = '';
            for i=1:numel(fields)
                if effects_logic(i)==1
                    mytitle = strcat(mytitle,'-',fields(i));
                end
            end
            mytitle{:} = strrep(mytitle{:},'-RF_EL-RF_BF-RF_RX',['RF_',num2str(obj.WGTS_CD_MolPerCm2_RelErr),'CDErr']);
            mytitle{:} = strrep(mytitle{:},'FSD',['FSD_',num2str(obj.FSDNorm_RelErr),'NormErr',num2str(obj.FSDShapeGS_RelErr),'ShapeGSErr',num2str(obj.FSDShapeES_RelErr),'ShapeESErr']);
            mytitle{:} = strrep(mytitle{:},'TASR',['TASR_',num2str(mean(obj.WGTS_TASR_RelErr)),'Err']);
            %             combi_name = sprintf('CombiCM_%s_%u-Trials_%s.mat', obj.StudyObject.TD,obj.nTrials,mytitle{:});
            %             obj.CovMatFile = strcat(combi_path,combi_name);
            %             %if strcmp(PlotSaveCM,'ON')
            %              %   obj.PlotCM('PlotEffect','Combined','savePlot','ON','Convergence','OFF');
            %             %end
            %             save(obj.CovMatFile, 'obj','CovMatFracCombi');
            %
            %             %Compute normal CovMat
            %             obj.ComputeFracCM('Mode','Frac2CM');
            %             save(obj.CovMatFile, 'obj','-append');
            %
            %             %Compute Decomposition
            %             [~, obj.CovMatFracShape] = obj.DecomposeCM('CovMatFrac',obj.CovMatFrac);
            %             save(obj.CovMatFile, 'obj','-append');
        end
    end 
    methods % auxillary methods
        function out = VaryParameter(obj,varargin)
            p = inputParser;
            p.addParameter('Parameter', obj.StudyObject.HVRipplesP2PV, @(x)isfloat(x));
            p.addParameter('rel_Error',0.1,@(x)isfloat(x));
            p.addParameter('Trials', obj.nTrials, @(x)isfloat(x));
            p.parse(varargin{:});
            
            Parameter = p.Results.Parameter;
            rel_Error = p.Results.rel_Error;
            Trials    = p.Results.Trials;
            
            out = Parameter.*(1+randn(Trials,1).*rel_Error);
        end
        %Sanity Checks
        function out = ComputeRank(obj)
            %Computes Rank of Covariance Matrix
            %Rank has to be <= (nTrials -1)
            %Otherswise CM is singular (cannot be inverted)
            
            CovMatRank = rank(obj.CovMat);
            if CovMatRank <= (obj.nTrials-1)
                fprintf('--------------------------------------------------------------------------\n')
                cprintf('blue','CovarianceMatrix: CovMatRank\n');
                fprintf('Rank of CovMat: %u  <= %u (Trials-1) \n', CovMatRank,obj.nTrials-1)
                fprintf('Santity Check successful! \n');
            else
                fprintf('--------------------------------------------------------------------------\n')
                fprintf('Rank of CovMat: %u  > %u (Trials-1) \n', CovMatRank,obj.nTrials-1)
                fprintf(2,'WARNING: CovMat singular! \n');
            end
            out = CovMatRank;
        end
        function [Trials, Convergence] = ConvergenceTest(obj,varargin)
            % Approximate Convergence Test Study of Sample Covariance Matrices
            % Convergence Test after Cauchy Convergence Criterium
            % For all natural numbers k it holds, that if....
            % ||CovMat(nTrials+k)-CovMat(nTrials)||< eps    ; ||CovMat|| = sqrt(trace(CovMat))
            % ...for all arbitrary eps
            
            p = inputParser;
            p.addParameter('filename',obj.CovMatFile,@(x)ischar(x));
            p.addParameter('Nmax',obj.nTrials,@(x)isfloat(x) && x>0);
            p.addParameter('StepSize',10,@(x)isfloat(x) && x>0);
            p.addParameter('Criterium','Cauchy', @(x)ismember(x,{'Cauchy','FracVar','Determinant'}));
            p.addParameter('PlotFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            
            p.parse(varargin{:});
            filename  = p.Results.filename;
            Nmax      = p.Results.Nmax;
            Criterium = p.Results.Criterium;
            StepSize  = p.Results.StepSize;
            PlotFlag  = p.Results.PlotFlag;
            
            myfile = importdata(filename);
            if isfield(myfile,'TBDIS_V')
                TBDIS_V = sum(myfile.TBDIS_V,3); % in case of nRuns >1
            else
                message = 'Problem: Covariance Matrix File doesnt contain sample spectra (TBDIS_V) \nReminder: Convergence only works for single systematic effects - not for combi CovMat! \nSuggestion: To Plot CovMat use option "Convergence"-->"OFF"';
                error('MyFunction:fileNotFound',message);
                
            end
            
            if Nmax>size(TBDIS_V,2) % for background CM not all Trials are sucessfull
                Nmax = size(TBDIS_V,2);
            end
            
            %Init
            nStep = floor(Nmax/StepSize);
            CovMatnTrials = cell(nStep,1);     % CovMat for different nTrials
            CovMatFracnTrials = cell(nStep,1); % Fractional Covmats
            CovMatTrace = zeros(nStep,1);      % Trace of CovMats
            sumfrac = zeros(nStep,1);          % Sum of fractional CovMat Elements
            detCM   = zeros(nStep,1);          % determinant of CovMat
            TBDIS_NoBKG = (obj.nRuns.*myfile.obj.StudyObject.TBDIS)-...%integral spectrum without background
                myfile.obj.StudyObject.BKG_RateSec.*myfile.obj.StudyObject.qUfrac.*myfile.obj.StudyObject.TimeSec;
            
            for x = StepSize:StepSize:Nmax %Compute CovMats für different nTrials
                xn = x/StepSize;
                
                CovMatnTrials{xn,1} = cov(TBDIS_V(:,1:x)'); %CovMats
                CovMatTrace(xn) = norm(CovMatnTrials{xn,1},'fro'); %Frobenius Norm
                
                if strcmp(Criterium, 'FracVar')  %Compute Fractional CovMats for different convergence test
                    CMFrac = CovMatnTrials{xn,1}./TBDIS_NoBKG./TBDIS_NoBKG';
                    CMFrac(isnan(CMFrac)) = 0;
                    CMFrac(isinf(abs(CMFrac))) = 0;
                    CovMatFracnTrials{xn,:} = CMFrac;
                    %TBDIS_Sample = mvnrnd(TBDIS_NoBKG,CovMatnTrials{xn,1},obj.nTrials*5);
                    %TBDIS_SumExpected = sum(TBDIS_NoBKG);
                    %TBDIS_SumSample = sum(TBDIS_Sample,2);
                    % TBDIS_SampleNorm = TBDIS_Sample.*(TBDIS_SumExpected./TBDIS_SumSample);
                    %CovMatShape = cov(TBDIS_SampleNorm);
                    %CovMatFracShape = CovMatShape./TBDIS_NoBKG./TBDIS_NoBKG';
                    %CovMatFracShape(isnan(CovMatFracShape))=0;
                    %CovMatFracShape(isinf(CovMatFracShape))=0;
                    %CovMatFracnTrials{xn,:} = CovMatFracShape;
                    sumfrac(xn) = sum(sum(CovMatFracnTrials{xn,:}));
                elseif strcmp(Criterium, 'Determinant')
                    detCM(xn) = det(CovMatnTrials{xn,1});
                end
            end
            
            %plots
            Trials = StepSize:StepSize:Nmax;
            switch Criterium
                case 'Cauchy'
                    Convergence = CovMatTrace;
                case 'FracVar'
                    Convergence = sumfrac;
                case 'Determinant'
                    Convergence = detCM;
            end
            if strcmp(PlotFlag,'ON')
                if strcmp(Criterium,'Cauchy')
                    f123 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
                    plot(StepSize:StepSize:Nmax,CovMatTrace,'-','LineWidth',3,'Color',rgb('DodgerBlue'));
                    xlabel('Trials','FontSize',15);
                    ylabel('|| CM ||', 'FontSize',15);
                    leg = legend(sprintf('Convergence Test (Cauchy Criterium)')); %\n Response Function CM for First Tritium
                    legend boxoff;
                    PrettyFigureFormat('FontSize',24);
                    export_fig(f123,'./plots/ConvergenceTestCauchy.pdf','-painters');
                end
                
                if strcmp(Criterium,'FracVar')
                    figure(11)
                    plot(StepSize:StepSize:Nmax,sumfrac,'x');
                    xlabel('Trials');
                    ylabel('$\mathbf{\Sigma}$ Fractional CM','Interpreter','latex', 'FontSize',15);
                    title('Convergence Test of Sample Covariance Matrix: Fractional Variance Criterium')
                end
                %save
                % save('results/CovMatnTrialStudy.mat',...
                %'obj', 'myfile.TBDIS_V'','CovMatnTrialsdia', CovMatnTrials,'-mat');
            end
        end
        function CombineTrials(obj,varargin)
            nTotal = 1000;
            set1 = importdata('WGTSMACE_CovMat_1500Trials_FT2_Sec_86400_1-Runs.mat');
            %set2 = importdata('WGTSMACE_CovMat_1000Trials_FT1_Sec_86400.mat');
            %set3 = importdata('WGTSMACE_CovMat_1500Trials_FT1_Sec_86400.mat');
            
            TBDIS_V = set1.TBDIS_V(:,1:1000)*10;%zeros(set1.obj.StudyObject.nqU, nTotal);
            %TBDIS_V(:,1:500)     = set1.TBDIS_V;
            %TBDIS_V(:,501:1500)  = set2.TBDIS_V;
            %TBDIS_V(:,1501:3000) = set3.TBDIS_V;
            %Object = set1.obj.StudyObject;
            %TBDIS_V = TBDIS_V(:,randperm(length(TBDIS_V)));
            %save('../../inputs/CovMat/RF/CM/WGTSMACE_CovMat_3000Trials_FT1_Sec_86400.mat','TBDIS_V', 'Object', '-mat');
            save('WGTSMACE_CovMat_1000Trials_FT2_Sec_86400_1-Runs.mat', 'TBDIS_V','-mat');
        end
        
        function GetTDlabel(obj)
             obj.TDlabel = strrep(obj.StudyObject.TD,'_','');
            
            if contains(obj.StudyObject.TD,'mNuSq')
                obj.TDlabel = extractBefore(obj.TDlabel,'mNuSq');
            end
             
            if contains(obj.TDlabel,'E0')
                obj.TDlabel = extractBefore(obj.TDlabel,'E0');
            end
            
            if ~isnan(str2double(obj.StudyObject.TD)) % if fake run
                 obj.TDlabel = ['FakeRun',obj.TDlabel];
            end
            
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
            dsdW     = @(y,w) -1.4135 ./ (1.200943945173595e+03+1) .* ( y + 1./(1.200943945173595e+03-w) ) + 1.4135 .* (y.^2 + 1./(1.200943945173595e+03-w).^2) + 3.545431573764978.*y.* dfdwH2(y);
            f = (e <=  0) .* 0 + ...
                (e <  obj.Ei) .*   (amp1 .* gaussian(e , pos1, sig1)  + amp2 .* gaussian(e, pos2, sig2) + amp3 .* gaussian(e , pos3, sig3)) + ...
                (e >= obj.Ei) .*  ((amp1 .* gaussian(obj.Ei, pos1, sig1)  + amp2 .* gaussian(obj.Ei, pos2, sig2) + amp3 .* gaussian(obj.Ei, pos3, sig3))  .* dsdW(y(e),w(e)) ./ dsdW(1,0));
            f(isnan(f)) = 0;
        end
    end
    
    
    
end
