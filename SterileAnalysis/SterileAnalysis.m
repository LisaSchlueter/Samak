classdef SterileAnalysis < handle
    % class to collect and organise sterile analysis routines
    %%
    properties (Access=public)
        RunAnaObj; % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        nGridSteps;   
        SmartGrid;
        RecomputeFlag;
        SysEffect; 
        range; % in eV
        LoadGridArg;
        
        % randomized MC
        RandMC;       % 1) if 'OFF' -> Asimov twins. 2) if @(x)isfloat(x) && numel(x)==1 -> randomize twins 
        RandMC_TBDIS; % external (e.g. randomized) spectrum for grid search -> use together with RandMC isfloat(x) for label
        Twin_mNu4Sq; % for randomized MCs with sterile-nu hypothesis
        Twin_sin2T4; % for randomized MCs with sterile-nu hypothesis
        
        % grid parameters
        mNu4Sq;
        sin2T4;   
        chi2;
        mNuSq;
        E0;
        
        % contour
        mNu4Sq_contour;
        sin2T4_contour;
        
        % best fit
        chi2_ref;
        mNu4Sq_bf;
        sin2T4_bf; 
        chi2_bf;
        mNuSq_bf;
        E0_bf;
        
        % null hypothesis
        chi2_Null; 
        DeltaChi2;
        ConfLevel;
        dof;
        
        
        % plot
        PlotColors;
        PlotLines;
        InterpMode; %lin or spline
        NullHypothesis;
    end
    
    methods % constructor
        function obj = SterileAnalysis(varargin)
            p = inputParser;
            p.addParameter('RunAnaObj','', @(x)  isa(x,'RunAnalysis') || isa(x,'MultiRunAnalysis'));
            p.addParameter('nGridSteps',50,@(x)isfloat(x)); % time on server (KSN1): 25 points ~ 7-10 minutes, 50 points ~ 40 minutes on csltr server
            p.addParameter('SmartGrid','OFF',@(x)ismember(x,{'ON','OFF'})); % work in progress
            p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SysEffect','all',@(x)ischar(x)); % if chi2CMShape: all or only 1
            p.addParameter('RandMC','OFF',@(x)ischar(x) || isfloat(x)); % randomize twins if RandMC is float
            p.addParameter('RandMC_TBDIS','',@(x)isfloat(x) || isempty(x));
            p.addParameter('range',65,@(x)isfloat(x));
            p.addParameter('ConfLevel',95,@(x)isfloat(x));
            p.addParameter('InterpMode','spline',@(x)ismember(x,{'lin','spline'}));
            p.addParameter('Twin_mNu4Sq',0,@(x)isfloat(x)); % for randomized MCs with sterile-nu hypothesis
            p.addParameter('Twin_sin2T4',0,@(x)isfloat(x)); % for randomized MCs with sterile-nu hypothesis
            p.addParameter('LoadGridArg','',@(x)iscell(x) ||isempty(x));
             p.addParameter('NullHypothesis','OFF',@(x)ismember(x,{'ON','OFF'})); % use NH instead of global minimum
            p.parse(varargin{:});
            
            obj.RunAnaObj = p.Results.RunAnaObj;
            obj.nGridSteps    = p.Results.nGridSteps;
            obj.SmartGrid     = p.Results.SmartGrid;
            obj.RecomputeFlag = p.Results.RecomputeFlag;
            obj.SysEffect     = p.Results.SysEffect;
            obj.RandMC        = p.Results.RandMC;
            obj.RandMC_TBDIS  = p.Results.RandMC_TBDIS;
            obj.range         = p.Results.range;
            obj.ConfLevel     = p.Results.ConfLevel;
            obj.InterpMode    = p.Results.InterpMode;
            obj.Twin_sin2T4   = p.Results.Twin_sin2T4;
            obj.Twin_mNu4Sq   = p.Results.Twin_mNu4Sq;
            obj.LoadGridArg   = p.Results.LoadGridArg;
            obj.NullHypothesis = p.Results.NullHypothesis;
            
            if isempty(obj.RunAnaObj)
                fprintf(2,'RunAnaObj has to be specified! \n');
                return
            end
            GetSamakPath; %sets current samak path as enviromental variable
            
            obj.SetNPfactor; % -> stat or syst
            obj.RunAnaObj.exclDataStart = obj.RunAnaObj.GetexclDataStart(obj.range);
            obj.InitPlotArg; % some plotting defaults
        end
    end
    
    methods
        function GridSearch(obj,varargin)
            p = inputParser;
            p.addParameter('Negsin2T4','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('NegmNu4Sq','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Extsin2T4','OFF',@(x)ismember(x,{'ON','OFF'})); %extended sin2T2 (up to 1)
            p.addParameter('ExtmNu4Sq','OFF',@(x)ismember(x,{'ON','OFF','0.01'})); %extended m4Sq (from 0.1)
            p.addParameter('mNu4SqTestGrid','OFF',@(x) strcmp(x,'OFF') || isfloat(x)); % different grid..
            p.addParameter('FixmNuSq',0,@(x)isfloat(x)); % if light nu-mass fixed (eV^2)
             p.parse(varargin{:});
            Negsin2T4     = p.Results.Negsin2T4;
            NegmNu4Sq     = p.Results.NegmNu4Sq;
            Extsin2T4     = p.Results.Extsin2T4;
            ExtmNu4Sq     = p.Results.ExtmNu4Sq;
            mNu4SqTestGrid = p.Results.mNu4SqTestGrid;
            FixmNuSq      = p.Results.FixmNuSq;

            savefile = obj.GridFilename('Negsin2T4',Negsin2T4,'NegmNu4Sq',NegmNu4Sq,...
                                        'Extsin2T4',Extsin2T4,'ExtmNu4Sq',ExtmNu4Sq,...
                                        'FixmNuSq',FixmNuSq,'mNu4SqTestGrid',mNu4SqTestGrid);
            
            if exist(savefile,'file') && strcmp(obj.RecomputeFlag,'OFF')
                fprintf('Grid file already exists: %s \nIf you want to load grid into memory call: obj.LoadGridFile() \n',savefile)
            else
                tStart = cputime;
                
                % set range
                obj.RunAnaObj.exclDataStart = obj.RunAnaObj.GetexclDataStart(obj.range);
                
                % get covariance matrix
                if strcmp(obj.RunAnaObj.chi2,'chi2CMShape') && strcmp(obj.SysEffect,'Bkg')
                    obj.RunAnaObj.ComputeCM('SysEffect',struct('FSD','OFF'),'BkgCM','ON','BkgPtCM','OFF');
                elseif strcmp(obj.RunAnaObj.chi2,'chi2CMShape') && strcmp(obj.SysEffect,'BkgPT') 
                    obj.RunAnaObj.ComputeCM('SysEffect',struct('FSD','OFF'),'BkgCM','OFF','BkgPtCM','ON');
                elseif strcmp(obj.RunAnaObj.chi2,'chi2CMShape') && ~strcmp(obj.SysEffect,'all')
                    obj.RunAnaObj.ComputeCM('SysEffect',struct(obj.SysEffect,'ON'),'BkgCM','OFF','BkgPtCM','OFF');
                elseif strcmp(obj.RunAnaObj.chi2,'chi2CMShape') && strcmp(obj.SysEffect,'all')
                    obj.RunAnaObj.ComputeCM('BkgCM','ON','BkgPtCM','ON');
                end
                
                %% ranomized mc data
                if strcmp(obj.RunAnaObj.DataType,'Twin') && isfloat(obj.RandMC) && numel(obj.RandMC_TBDIS)==obj.RunAnaObj.ModelObj.nqU
                   TBDIS_mc = obj.RandMC_TBDIS;
                   obj.RunAnaObj.RunData.TBDIS = TBDIS_mc;
                   obj.RunAnaObj.RunData.TBDISE = sqrt(TBDIS_mc);
                elseif strcmp(obj.RunAnaObj.DataType,'Twin') && (isfloat(obj.RandMC)  || (obj.Twin_mNu4Sq~=0 || obj.Twin_sin2T4~=0))
                    % 1) twin and isfloat(obj.RandMC) -> randomized data
                    % 2) twin and (obj.Twin_mNu4Sq~=0 || obj.Twin_sin2T4~=0) -> asimov twin with sterile nu-signal
                    obj.RunAnaObj.InitModelObj_Norm_BKG('RecomputeFlag','ON');
                    if obj.Twin_mNu4Sq~=0 || obj.Twin_sin2T4~=0
                        obj.RunAnaObj.ModelObj.BKG_RateSec_i = obj.RunAnaObj.ModelObj.BKG_RateSec;
                        obj.RunAnaObj.ModelObj.normFit_i = obj.RunAnaObj.ModelObj.normFit;
                        obj.RunAnaObj.ModelObj.SetFitBiasSterile(obj.Twin_mNu4Sq,obj.Twin_sin2T4);
                        obj.RunAnaObj.ModelObj.ComputeTBDDS;
                        obj.RunAnaObj.ModelObj.ComputeTBDIS;
                        TBDIS_i = obj.RunAnaObj.ModelObj.TBDIS';
                    else
%                         obj.RunAnaObj.ModelObj.ComputeTBDDS;
%                         obj.RunAnaObj.ModelObj.ComputeTBDIS;
                        TBDIS_i = obj.RunAnaObj.ModelObj.TBDIS';
                    end
                    
                    if isfloat(obj.RandMC)
                        % change to randomized MC data
                        TBDIS_mc = zeros(obj.RunAnaObj.ModelObj.nqU,1);
                        TBDIS_mc(obj.RunAnaObj.exclDataStart:end) = ...
                            mvnrnd(TBDIS_i(obj.RunAnaObj.exclDataStart:end),...
                            obj.RunAnaObj.FitCMShape(obj.RunAnaObj.exclDataStart:end,obj.RunAnaObj.exclDataStart:end),1)';
                        %TBDIS_mc = mvnrnd(TBDIS_i,obj.RunAnaObj.FitCMShape,1)';
                    else
                        TBDIS_mc = TBDIS_i';
                    end
                    obj.RunAnaObj.RunData.TBDIS = TBDIS_mc;
                    obj.RunAnaObj.RunData.TBDISE = sqrt(TBDIS_mc);
                end
                
                %% if light nu-mass is fixed in fit, but model nu-mass shall not be something else than 0 eV^2
                if FixmNuSq~=0 && contains(obj.RunAnaObj.fixPar,'fix 1 ;')
                    obj.RunAnaObj.ModelObj.mnuSq_i = FixmNuSq;
                end
                
                %% null hypothesis : no steriles
                obj.RunAnaObj.Fit;
                FitResults_Null = obj.RunAnaObj.FitResult;
                
                if strcmp(obj.RunAnaObj.DataType,'Twin') && (isfloat(obj.RandMC)  || (obj.Twin_mNu4Sq~=0 || obj.Twin_sin2T4~=0)) 
                    % re-set sterile parameters in model
                    if obj.Twin_mNu4Sq~=0 || obj.Twin_sin2T4~=0
                        obj.RunAnaObj.SimulateStackRuns;
                    end
                end
                %% define grid
                if strcmp(Extsin2T4,'ON')
                    sin2T4Max = 1;
                else
                    sin2T4Max = 0.5;
                end
                
                if strcmp(ExtmNu4Sq,'ON')
                    mnu4Sq_ex = [0.1;0.35;0.7];%;logspace(0,log10((obj.range)^2),obj.nGridSteps-3)'];
                    nGridSteps_i = obj.nGridSteps;
                    obj.nGridSteps = nGridSteps_i-3;
                elseif strcmp(ExtmNu4Sq,'0.01')
                    mnu4Sq_ex = [0.01;0.05;0.1;0.35;0.7];%;logspace(0,log10((obj.range)^2),obj.nGridSteps-3)'];
                    nGridSteps_i = obj.nGridSteps;
                    obj.nGridSteps = nGridSteps_i-5;
                end
                
                if mNu4SqTestGrid==1
                    mnu4SqSmall  = logspace(0,log10((obj.range-11)^2),obj.nGridSteps-5)';
                    mnu4SqLarge  = logspace(log10((obj.range-10)^2),log10((obj.range)^2),5)';
                    mnu4Sq       = sort([mnu4SqSmall;mnu4SqLarge]);
                elseif mNu4SqTestGrid==2
                    mnu4SqSmall  = logspace(0,log10((obj.range-11)^2),obj.nGridSteps-7)';
                    mnu4SqLarge  = logspace(log10((obj.range-10)^2),log10((obj.range)^2),7)';
                    mnu4Sq       = sort([mnu4SqSmall;mnu4SqLarge]);
                elseif mNu4SqTestGrid==3
                    mnu4SqSmall  = logspace(0,log10((obj.range-11)^2),obj.nGridSteps-10)';
                    mnu4SqLarge  = logspace(log10((obj.range-10)^2),log10((obj.range)^2),10)';
                    mnu4Sq       = sort([mnu4SqSmall;mnu4SqLarge]);
                elseif mNu4SqTestGrid==4
                    mnu4SqSmall  = logspace(0,log10((obj.range-11)^2),obj.nGridSteps-5)';
                    mnu4SqLarge  = logspace(log10((obj.range-10)^2),log10((obj.range)^2),5)';
                    mnu4Sq       = sort([mnu4SqSmall;mnu4SqLarge]);
                elseif mNu4SqTestGrid==5 || mNu4SqTestGrid==5.5
                    mnu4SqSmall  = logspace(0,log10((obj.range-11)^2),obj.nGridSteps-15)';
                    mnu4SqLarge  = logspace(log10((obj.range-10)^2),log10((obj.range)^2),15)';
                    mnu4Sq       = sort([mnu4SqSmall;mnu4SqLarge]);
                else
                    mnu4Sq      = logspace(0,log10((obj.range)^2),obj.nGridSteps)';
                end
                
                if ismember(ExtmNu4Sq,{'ON','0.01'})
                    mnu4Sq = [mnu4Sq_ex;mnu4Sq];%;logspace(0,log10((obj.range)^2),obj.nGridSteps-3)'];
                    obj.nGridSteps = nGridSteps_i;
                end
                
                if mNu4SqTestGrid==5.5
                    sin2T4      = logspace(log10(0.5),log10(1),obj.nGridSteps);
                elseif strcmp(Extsin2T4,'ON') && obj.nGridSteps<50
                    % add more points are large mixing
                    sin2T4      = sort([logspace(-3,log10(sin2T4Max),obj.nGridSteps-7),0.5,0.6,0.7,0.85,0.9,0.95,0.98]);
                else
                    sin2T4      = logspace(-3,log10(sin2T4Max),obj.nGridSteps);
                end
                
                mnu4Sq      = repmat(mnu4Sq,1,obj.nGridSteps);
                sin2T4      = repmat(sin2T4,obj.nGridSteps,1);
                
                %% make copy of model for parallel computing
                D = copy(repmat(obj.RunAnaObj,obj.nGridSteps^2,1));
                D = reshape(D,numel(D),1);
                %% scan over msq4-sin2t4 grid
                chi2Grid       = zeros(obj.nGridSteps^2,1);
                FitResultsGrid = cell(obj.nGridSteps^2,1);
                mnu4Sq_Grid    = reshape(mnu4Sq',obj.nGridSteps^2,1);
                sin2T4_Grid    = reshape(sin2T4',obj.nGridSteps^2,1);
                
                if strcmp(NegmNu4Sq,'ON')
                    mnu4Sq_Grid = -mnu4Sq_Grid;
                    mnu4Sq = - mnu4Sq;
                end
                
                if strcmp(Negsin2T4,'ON')
                    sin2T4_Grid = -sin2T4_Grid;
                    sin2T4 = -sin2T4;
                end
                
                fixPartmp = obj.RunAnaObj.fixPar;
                DataTypetmp = obj.RunAnaObj.DataType;
                
                parfor i= 1:(obj.nGridSteps^2)
                    if strcmp(DataTypetmp,'Fake')
                        D(i).SimulateRun;
                    else
                        D(i).SimulateStackRuns;
                    end
                    if FixmNuSq~=0 && contains(fixPartmp,'fix 1 ;')
                        % if light nu-mass is fixed and shall not be fixed to 0 eV^2
                        D(i).ModelObj.mnuSq_i = FixmNuSq;
                    end
                    D(i).ModelObj.SetFitBiasSterile(mnu4Sq_Grid(i),sin2T4_Grid(i));
                    D(i).Fit
                    chi2Grid(i) = D(i).FitResult.chi2min;
                    FitResultsGrid{i} = D(i).FitResult;
                end
                
                chi2       = reshape(chi2Grid,obj.nGridSteps,obj.nGridSteps);
                FitResults = reshape(FitResultsGrid,obj.nGridSteps,obj.nGridSteps);
                
                if min(min(chi2))<obj.RunAnaObj.FitResult.chi2min
                    chi2_ref = min(min(chi2));
                else
                    chi2_ref = obj.RunAnaObj.FitResult.chi2min;
                end
                
                tCpuHour = (cputime-tStart)/60; % cpu time in hours
                %% save
                save(savefile,'chi2_ref',...%'FitResults_ref'
                    'chi2','mnu4Sq','sin2T4','FitResults','FitResults_Null','tCpuHour');
                fprintf('save file to %s \n',savefile);
                
                if strcmp(obj.RunAnaObj.DataType,'Twin') && isfloat(obj.RandMC)
                    save(savefile,'TBDIS_mc','-append');
                end
            end
        end
    end
    
    methods
        function Interp1Grid(obj,varargin)
            % mesh grid interpolation
            % Finer binning, nicer appearance!
            p = inputParser;
            p.addParameter('nInter',1e3,@(x)isfloat(x));
            p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'})); % force recompute
            p.addParameter('Maxm4Sq','',@(x)isfloat(x) || isempty(x));
            p.addParameter('Minm4Sq',min(min(obj.mNu4Sq)),@(x)isfloat(x) || isempty(x));
           
            p.parse(varargin{:});
            nInter        = p.Results.nInter;
            Maxm4Sq       = p.Results.Maxm4Sq;
            Minm4Sq       = p.Results.Minm4Sq;
            
            RecomputeFlag_loc = p.Results.RecomputeFlag;
            
            if size(obj.mNu4Sq,1)>=nInter && strcmp(RecomputeFlag_loc,'OFF')
                fprintf('Interp1 stoped - mNuSq size already large than interpolation \n')
                return
            end
            
            if strcmp(obj.InterpMode,'Mix')
                % test new method:
                % spline interpolation up to Maxm4Sq
                % lin interpolation up to range^2
                if isempty(Maxm4Sq)
                    Maxm4Sq = (obj.range-4)^2;
                end
                
                nInter1 = 700;
                nInter2 = 300;
                
                % first part: spline
                LogIdx = obj.mNu4Sq(:,1)<Maxm4Sq;
                [X1,Y1] = meshgrid(obj.mNu4Sq(LogIdx,1),obj.sin2T4(1,:));
                Z1 = obj.chi2(:,LogIdx);
                mNu4tmp = logspace(log10(min(obj.mNu4Sq(LogIdx,1))),log10(max(obj.mNu4Sq(LogIdx,1))),nInter1);
                mNu4Sq_1 = repmat(mNu4tmp,nInter,1);
                sin2T4_1 = repmat(logspace(log10(min(min(obj.sin2T4))),log10(max(max(obj.sin2T4))),nInter),nInter1,1)'; 
                chi2_1  = reshape(interp2(X1,Y1,Z1,mNu4Sq_1,sin2T4_1,'spline'),nInter,nInter1);
                chi2_1(chi2_1<0) = NaN;
                chi2_ref_1 = min(min(chi2_1));
                mNuSq_1 = reshape(interp2(X1,Y1,obj.mNuSq(:,LogIdx),mNu4Sq_1,sin2T4_1,'spline'),nInter,nInter1);
                E0_1 = reshape(interp2(X1,Y1,obj.E0(:,LogIdx),mNu4Sq_1,sin2T4_1,'spline'),nInter,nInter1);
               
                
                % second part: lin
                LogIdx = obj.mNu4Sq(:,1)>=Maxm4Sq;
                [X2,Y2] = meshgrid(obj.mNu4Sq(LogIdx,1),obj.sin2T4(1,:));
                Z2 = obj.chi2(:,LogIdx);
                mNu4tmp = logspace(log10(min(obj.mNu4Sq(LogIdx,1))),log10(obj.range.^2),nInter2);%max(obj.mNu4Sq(LogIdx,1))
                mNu4Sq_2 = repmat(mNu4tmp,nInter,1);
                sin2T4_2 = repmat(logspace(log10(min(min(obj.sin2T4))),log10(max(max(obj.sin2T4))),nInter),nInter2,1)'; 
                chi2_2  = reshape(interp2(X2,Y2,Z2,mNu4Sq_2,sin2T4_2,'lin'),nInter,nInter2);
                chi2_2(chi2_2<0) = NaN;
                chi2_ref_2 = min(min(chi2_2));
                mNuSq_2 = reshape(interp2(X2,Y2,obj.mNuSq(:,LogIdx),mNu4Sq_2,sin2T4_2,'lin'),nInter,nInter2);
                E0_2 = reshape(interp2(X2,Y2,obj.E0(:,LogIdx),mNu4Sq_2,sin2T4_2,'lin'),nInter,nInter2);
              
                
                % merge
                obj.mNu4Sq   = [mNu4Sq_1,mNu4Sq_2];
                obj.sin2T4   = [sin2T4_1,sin2T4_2];
                obj.chi2     = [chi2_1,chi2_2];
                obj.chi2_ref = min([chi2_ref_1,chi2_ref_2]);
                obj.mNuSq    = [mNuSq_1,mNuSq_2];    % best fit results (nuisance parameter)
                obj.E0       = [E0_1,E0_2];          % best fit results (nuisance parameter)
     
            else
                %% define maximum m4:
                if isempty(Maxm4Sq)
                    Maxm4Sq = (obj.range)^2;
                end
                
                [X,Y] = meshgrid(obj.mNu4Sq(:,1),obj.sin2T4(1,:));
                
                % define fine interpolation grid
                % negative quadrants make this complicated
                if  obj.mNu4Sq(1,1)<0
                    Minm4Sq = min(min(-obj.mNu4Sq));
                    mNu4tmp = logspace(log10(Minm4Sq),log10(Maxm4Sq),nInter);
                    mNu4tmp = - mNu4tmp;
                else
                    mNu4tmp = logspace(log10(Minm4Sq),log10(Maxm4Sq),nInter);
                end
                obj.mNu4Sq = repmat(mNu4tmp,nInter,1);
                if obj.sin2T4(1,1)<0
                    obj.sin2T4 = repmat(logspace(log10(min(min(-obj.sin2T4))),log10(max(max(-obj.sin2T4))),nInter),nInter,1)';
                    obj.sin2T4 = - obj.sin2T4;
                else
                    obj.sin2T4 = repmat(logspace(log10(min(min(obj.sin2T4))),log10(max(max(obj.sin2T4))),nInter),nInter,1)';
                end
                obj.chi2   = reshape(interp2(X,Y,obj.chi2,obj.mNu4Sq,obj.sin2T4,obj.InterpMode),nInter,nInter);
                
                % best fit results (nuisance parameter)
                obj.mNuSq  = reshape(interp2(X,Y,obj.mNuSq,obj.mNu4Sq,obj.sin2T4,obj.InterpMode),nInter,nInter);
                obj.E0  = reshape(interp2(X,Y,obj.E0,obj.mNu4Sq,obj.sin2T4,obj.InterpMode),nInter,nInter);
                
                obj.chi2(obj.chi2<0) = NaN;
                obj.chi2_ref = min(min(obj.chi2));
            end
            
            fprintf('Grid interpolation sucessfull \n');
        end
        function FindBestFit(obj,varargin)
            % best fit parameters of m4 and sin2t4
            % based on interpolated grid
            p = inputParser;
            p.addParameter('Mode','Def',@(x)ismember(x,{'Def','Imp'}));
            p.parse(varargin{:});
            Mode = p.Results.Mode;
           switch Mode 
               case 'Def'
            if size(obj.mNu4Sq,1)<1e3
                obj.Interp1Grid;
            end
               case 'Imp'
                   % improved best fit finding
                   % WARNING: still in testing phase
                   % do interpolation only in vicinity of best fit
                   % to be called after Mode "Def" -> needs first estimateof best fit
                   
                   % location of old best fit
                   [row, col]    = find(obj.chi2 == min(obj.chi2(:)));
                   rowIdx_min = row-2;
                   rowIdx_max = row+2;
                   colIdx_min = col-2;
                   colIdx_max = col+2;

                   if row > size(obj.mNu4Sq,1)-2
                       % upper edge Ue4^2
                       rowIdx_max = size(obj.mNu4Sq,1);
                   elseif row<=2
                       % lower edge Ue4^2
                      rowIdx_min = 1;
                   end
                   
                   if col <=2
                      colIdx_min = 1; 
                   elseif col > size(obj.mNu4Sq,1)-2
                        colIdx_max = size(obj.mNu4Sq,1);   
                   end
                   
                   mNu4Sq_inter_min = obj.mNu4Sq(row,colIdx_min);%obj.mNu4Sq(col,row);
                   mNu4Sq_inter_max = obj.mNu4Sq(row,colIdx_max);
                   sin2T4_inter_min = obj.sin2T4(rowIdx_min,col);
                   sin2T4_inter_max = obj.sin2T4(rowIdx_max,col);
                   
                   obj.LoadGridFile(obj.LoadGridArg{:});
                   nInter = 1e3;
                   [X,Y] = meshgrid(obj.mNu4Sq(:,1),obj.sin2T4(1,:));
                   %mNu4tmp = logspace(log10(obj.mNu4Sq_bf-1),log10(obj.mNu4Sq_bf+1),nInter);
                   
                   mNu4tmp = linspace(mNu4Sq_inter_min,mNu4Sq_inter_max,nInter);
                   obj.mNu4Sq = repmat(mNu4tmp,nInter,1);
                   % obj.sin2T4 = repmat(logspace(log10(obj.sin2T4_bf-0.01),log10(obj.sin2T4_bf+0.01),nInter),nInter,1)';
                   obj.sin2T4 = repmat(linspace(sin2T4_inter_min,sin2T4_inter_max,nInter),nInter,1)';
                   obj.chi2   = reshape(interp2(X,Y,obj.chi2,obj.mNu4Sq,obj.sin2T4,'spline'),nInter,nInter);
                   obj.mNuSq  = reshape(interp2(X,Y,obj.mNuSq,obj.mNu4Sq,obj.sin2T4,'spline'),nInter,nInter);
                   obj.E0  = reshape(interp2(X,Y,obj.E0,obj.mNu4Sq,obj.sin2T4,'spline'),nInter,nInter);
           end
           
           [row, col]    = find(obj.chi2 == min(obj.chi2(:)));
            
            obj.mNu4Sq_bf = obj.mNu4Sq(row,col);%obj.mNu4Sq(col,row);
            obj.sin2T4_bf = obj.sin2T4(row,col);
            
            obj.chi2_ref = min(min(obj.chi2));
            obj.chi2_bf   = obj.chi2_ref;
            obj.mNuSq_bf  = obj.mNuSq(row,col);
            obj.E0_bf     = obj.E0(row,col);
            
            if strcmp(Mode,'Imp')
                % reload grid and do normal interpolation
%                 obj.LoadGridFile;
%                 obj.Interp1Grid;
            end
        end
        function [DeltaChi2, SignificanceBF] = CompareBestFitNull(obj,varargin)
            if isempty(obj.chi2_bf)
                obj.FindBestFit;
            end
            
            %  best fit
            fprintf('Best fit sinTsq = %.3f and m4sq = %.1f eV^2 \n',obj.sin2T4_bf,obj.mNu4Sq_bf);
            fprintf('Best fit:        chi2 = %.3f (%.0f dof) -> p-value = %.2f\n',obj.chi2_bf,obj.dof,1-chi2cdf(obj.chi2_bf,obj.dof));
            
            % null 
            fprintf('Null hypothesis: chi2 = %.3f (%.0f dof) -> p-value = %.2f\n',obj.chi2_Null,obj.dof+2,1-chi2cdf(obj.chi2_Null,obj.dof));
          %  x = linspace(10,99,1e2);
          %  y = GetDeltaChi2(x,2);
            DeltaChi2 = obj.chi2_Null-obj.chi2_bf;
            SignificanceBF = chi2cdf(DeltaChi2,2);%interp1(y,x, DeltaChi2,'spline');
            
            fprintf('Delta chi2 = %.2f -> %.1f%% C.L. significance \n',obj.chi2_Null-obj.chi2_bf,100.*SignificanceBF);
          
        end
        function [DeltamNu41Sq,sin2T4Sq] = Convert2Osci(obj,varargin)
            p = inputParser;
            p.addParameter('m4',obj.mNu4Sq,@(x)isfloat(x));
            p.addParameter('sinT4',obj.sin2T4,@(x)isfloat(x));
            p.parse(varargin{:});
            m4    = p.Results.m4;
            sinT4 = p.Results.sinT4;
            % convert KATRIN parameters into oscillation experiment parameter space
            % (sin(t4)^2,m4^2) --> (sin(2t4)^2,Delta(m41)^2)
            DeltamNu41Sq = m4 - obj.RunAnaObj.ModelObj.mnuSq;
            sin2T4Sq     = 4*sinT4.*(1-sinT4);
        end
        function [sin2T4_Stat, sin2T4_Sys, sin2T4_Tot, mNu4SqCommon, StatDomFraction] = StatOverSys(obj,varargin)
            p = inputParser;
            p.addParameter('Ranges',[95:-5:65],@(x)isfloat(x));
            p.parse(varargin{:});
            Ranges   = p.Results.Ranges;
            chi2_i = obj.RunAnaObj.chi2;
            
            DeltaChi2_1Par = GetDeltaChi2(0.6827,1);         
            mNu4SqCommon = cell(numel(Ranges),1);
            sin2T4_Sys = cell(numel(Ranges),1);
            sin2T4_Stat = cell(numel(Ranges),1);
            sin2T4_Tot = cell(numel(Ranges),1);
            
            if strcmp(obj.RunAnaObj.DataSet,'Knm1')
            obj.nGridSteps = 25;
            obj.InterpMode  = 'lin';
            end
            NP_prev = obj.RunAnaObj.NonPoissonScaleFactor;
            
            for i=1:numel(Ranges)
                %  progressbar(i/numel(Ranges));
                obj.range = Ranges(i);
                
                obj.RunAnaObj.chi2 = 'chi2Stat';
                obj.RunAnaObj.NonPoissonScaleFactor = 1;
                obj.LoadGridFile('CheckSmallerN','OFF','CheckLargerN','OFF',obj.LoadGridArg{:});
                obj.Interp1Grid('RecomputeFlag','ON','nInter',1e3);      
                [M,c]= contour(obj.sin2T4,obj.mNu4Sq,obj.chi2-obj.chi2_ref,...
                    [DeltaChi2_1Par DeltaChi2_1Par]);
                sinStattmp = M(1,2:end);
                mNu4Stattmp = M(2,2:end);
                c.delete;

                
                obj.RunAnaObj.chi2 = 'chi2CMShape';
                obj.SetNPfactor;
                obj.LoadGridFile('CheckSmallerN','OFF','CheckLargerN','OFF',obj.LoadGridArg{:});
                
                obj.Interp1Grid('RecomputeFlag','ON','nInter',1e3);
                [M,c]= contour(obj.sin2T4,obj.mNu4Sq,obj.chi2-obj.chi2_ref,...
                    [DeltaChi2_1Par DeltaChi2_1Par]);
                
                sinCMtmp = M(1,2:end);
                mNu4CMtmp = M(2,2:end);
                c.delete;
                
                % find common
                [LogicCM,LogicStat] = ismember(mNu4CMtmp,mNu4Stattmp);
                mNu4SqCommontmp     = mNu4CMtmp(logical(LogicCM)); %equivalent to mNuStattmp(LogicStat,i); 
                sin2T4_Stattmp      = sinStattmp(logical(LogicStat));
                sin2T4_totTmp       = sinCMtmp(logical(LogicCM));
                
                mNu4SqCommon{i}     = mNu4SqCommontmp(sin2T4_Stattmp<sin2T4_totTmp);
                sin2T4_Stat{i}      = sinStattmp(sin2T4_Stattmp<sin2T4_totTmp);
                sin2T4_Tot{i}      = sin2T4_totTmp(sin2T4_Stattmp<sin2T4_totTmp);   
                sin2T4_Sys{i}        = sqrt(sin2T4_Tot{i}.^2-sin2T4_Stat{i}.^2);
            end
            
             obj.RunAnaObj.NonPoissonScaleFactor =  NP_prev;
            %% get ratio syst % stat
            StatDomFraction = zeros(numel(Ranges),1); % fraction of m4 that are stat. dominated;
            mnu4SqCommonLin   = cell(numel(Ranges),1);
            
            for i=1:numel(Ranges)
                mnu4SqCommonLin{i}   = linspace(min(mNu4SqCommon{i}),max(mNu4SqCommon{i}),1e5);
                sin2T4Syst_tmp = interp1(mNu4SqCommon{i},sin2T4_Sys{i},mnu4SqCommonLin{i},'lin');
                sin2T4Stat_tmp = interp1(mNu4SqCommon{i},sin2T4_Stat{i},mnu4SqCommonLin{i},'lin');
                
                Ratio = sin2T4Syst_tmp.^2./sin2T4Stat_tmp.^2;
                StatDomFraction(i) = sum(Ratio<=1)./numel(Ratio); 
            end
            
             obj.RunAnaObj.chi2 = chi2_i;
        end  
        function  [Ratio,StatDomFraction,mNu4Sq_inter,sin2t4_Stat_inter,sin2t4_Tot_inter,sin2t4_Sys_inter] = StatOverSysKsn2(obj,varargin)
             p = inputParser;
            p.addParameter('RasterScan','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            RasterScan = p.Results.RasterScan;
            chi2_i = obj.RunAnaObj.chi2;

            NP_prev = obj.RunAnaObj.NonPoissonScaleFactor;
            
            obj.RunAnaObj.chi2 = 'chi2Stat';
            obj.RunAnaObj.NonPoissonScaleFactor = 1;
            obj.LoadGridFile(obj.LoadGridArg{:});
            obj.Interp1Grid('RecomputeFlag','ON','nInter',1e3,'Maxm4Sq',38^2);
            obj.ContourPlot('RasterScan',RasterScan); close;
            mNu4Sq_Stat = obj.mNu4Sq_contour;
            sin2t4_Stat = obj.sin2T4_contour;
            
            if ~strcmp(obj.SysEffect,'NP')
                obj.RunAnaObj.chi2 = 'chi2CMShape';
            end
            obj.RunAnaObj.NonPoissonScaleFactor = NP_prev;
            obj.LoadGridFile(obj.LoadGridArg{:});
            obj.Interp1Grid('RecomputeFlag','ON','nInter',1e3,'Maxm4Sq',38^2);
            obj.ContourPlot('RasterScan',RasterScan); close;
            mNu4Sq_Tot = obj.mNu4Sq_contour;
            sin2t4_Tot = obj.sin2T4_contour;
            
            %% interpolate at common mNu4Sq
            mNu4Sq_inter = linspace(min([mNu4Sq_Stat,mNu4Sq_Stat]),min([max(mNu4Sq_Stat),max(mNu4Sq_Stat)]),1e4);
            sin2t4_Stat_inter = interp1(mNu4Sq_Stat,sin2t4_Stat,mNu4Sq_inter,'spline');
            sin2t4_Tot_inter  = interp1(mNu4Sq_Tot,sin2t4_Tot,mNu4Sq_inter,'spline');

            %% get ratio syst % stat
            sin2t4_Sys_inter = sqrt(sin2t4_Tot_inter.^2-sin2t4_Stat_inter.^2);
            Ratio = sin2t4_Sys_inter.^2./sin2t4_Tot_inter.^2;
            
            StatDomFraction = sum(Ratio<=1)./numel(Ratio); 
        
            obj.RunAnaObj.NonPoissonScaleFactor =  NP_prev;
            obj.RunAnaObj.chi2 = chi2_i;
        end
    end
    
    %% Plotting
    methods
        function [pHandle,sin2T4_min,PlotHandleBf] = ContourPlot(obj,varargin)
            p = inputParser;
            p.addParameter('BestFit','OFF',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('ReCalcBF','ON',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('CL',obj.ConfLevel,@(x)isfloat(x));   
            p.addParameter('HoldOn','OFF',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('Color',obj.PlotColors{1},@(x) isfloat(x));
            p.addParameter('LineStyle',obj.PlotLines{1},@(x) ischar(x));
            p.addParameter('PlotSplines','OFF',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','png'}));
            p.addParameter('ExtraStr','',@(x)ischar(x)); % for plot label    
            p.addParameter('RasterScan','OFF',@(x) ismember(x,{'ON','OFF'})); % 1 par instead of 2 par
            p.addParameter('NullHypothesis',obj.NullHypothesis,@(x)ismember(x,{'ON','OFF'}));
            
            p.parse(varargin{:});  
            CL          = p.Results.CL;      % also works with vector
            HoldOn      = p.Results.HoldOn;
            myColor     = p.Results.Color;
            myLineStyle = p.Results.LineStyle;
            BestFit     = p.Results.BestFit;
            PlotSplines = p.Results.PlotSplines;
            SavePlot    = p.Results.SavePlot;
            ExtraStr    = p.Results.ExtraStr;
            RasterScan  = p.Results.RasterScan;
            ReCalcBF    = p.Results.ReCalcBF;
            NullHypothesis_local = p.Results.NullHypothesis;
            
            if strcmp(HoldOn,'ON')
                hold on;
            elseif strcmp(HoldOn,'OFF')
                GetFigure;
            end
            
            if strcmp(NullHypothesis_local,'ON')
                chi2_ref_local = obj.chi2_Null;
            else
                chi2_ref_local = obj.chi2_ref;
            end
            if strcmp(RasterScan,'ON')
                obj.DeltaChi2 = GetDeltaChi2(CL,1);
                if strcmp(obj.RunAnaObj.DataType,'Real')
                    chi2_ref_local = min(obj.chi2); % minimum for each m4^2
                end
            else
                obj.DeltaChi2 = GetDeltaChi2(CL,2);
            end
            
             % contour plot
            PlotArg = {'LineWidth',2.5,'LineStyle',myLineStyle};
            if numel(CL)==1
                PlotArg = [PlotArg,{'LineColor',myColor}];
            end
            
            [M,pHandle]= contour(obj.sin2T4,obj.mNu4Sq,obj.chi2-chi2_ref_local,...
                [obj.DeltaChi2 obj.DeltaChi2],...
                PlotArg{:});
            
            sin2T4_min = min(M(1,:));
            ExclIndex = find(M(1,:)==obj.DeltaChi2(1));
            nContour = numel(ExclIndex);
            if nContour>1
                PlotSplines = 'OFF';
                fprintf('set plot splines to off \n')
            end
            
            if strcmp(PlotSplines,'ON')
                sin2t4_tmp = M(1,2:end);  % remove information about contour
                mnu4sq_tmp = M(2,2:end);
                
                mnu4Sq_contour = logspace(log10(min(mnu4sq_tmp)),log10(max(mnu4sq_tmp)),1e4);
                obj.mNu4Sq_contour = sort([mnu4Sq_contour,logspace(log10(8e2),log10(2e3),1e3)]);
                
                obj.sin2T4_contour = interp1(mnu4sq_tmp,sin2t4_tmp,obj.mNu4Sq_contour,'spline');
                
                pHandle.delete;
                PlotArg = {'LineWidth',2.5,'LineStyle',myLineStyle};
                if numel(CL)==1
                    PlotArg = [PlotArg,{'Color',myColor}];
                end
                pHandle = plot(obj.sin2T4_contour,obj.mNu4Sq_contour,PlotArg{:});
            else
                obj.sin2T4_contour = M(1,2:end);  % remove information about contour
                obj.mNu4Sq_contour = M(2,2:end); 
            end
            
            if numel(CL)>1
                cm = colormap('cool');
                bf_color = cm(1,:);
            else
               bf_color = pHandle.LineColor;
            end
            
            % best fit
            if strcmp(BestFit,'ON')
                if strcmp(ReCalcBF,'ON')
                    obj.FindBestFit;
                    % obj.FindBestFit('Mode','Imp');
                end
                hold on;
                PlotHandleBf= plot(obj.sin2T4_bf,obj.mNu4Sq_bf,'x','MarkerSize',9,'Color',bf_color,'LineWidth',pHandle.LineWidth);
            else 
                PlotHandleBf = 0;
            end
            set(gca,'YScale','log');
            set(gca,'XScale','log');
            xlabel(sprintf('|{\\itU}_{e4}|^2'));
            ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
            PrettyFigureFormat('FontSize',22);
            title(sprintf('%s',obj.GetPlotTitle),'FontWeight','normal','FontSize',get(gca,'FontSize'));
            if any(CL<1)
                CL = CL*1e2;
            end
            legStr = sprintf(' %.3g%% C.L. -',CL);  
            legend(legStr(1:end-1),'EdgeColor',rgb('Silver'),'Location','southwest');
            
              if ~strcmp(SavePlot,'OFF')
                  plotname = sprintf('%s_Contour_%.2gCL%s',obj.DefPlotName,obj.ConfLevel,ExtraStr);       
                 if strcmp(SavePlot,'ON')
                     plotname = [plotname,'.pdf'];
                     export_fig(gcf,plotname);
                 elseif strcmp(SavePlot,'png')
                     plotname = [plotname,'.png'];
                     print(gcf,plotname,'-dpng','-r450');
                 end
                 fprintf('save plot to %s \n',plotname);
             end
             
        end
        function [legHandle,legStr] = ContourPlotOsci(obj,varargin)
            % contour plot in osicllation parameter space
            p = inputParser;
            p.addParameter('BestFit','OFF',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('CL',obj.ConfLevel,@(x)isfloat(x));   
            p.addParameter('HoldOn','OFF',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('Color',obj.PlotColors{1},@(x) isfloat(x));
            p.addParameter('LineStyle',obj.PlotLines{1},@(x) ischar(x));
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','png'}));
            p.addParameter('RAA','ON',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('Mainz','ON',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('Troitsk','ON',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('Neutrino4','ON',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('Prospect','ON',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('DANSS','ON',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('DayaBay','ON',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('DoubleChooz','ON',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('Stereo','ON',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('Style','Reg',@(x) ismember(x,{'Reg','PRL'}));
            p.addParameter('FinalSensitivity','OFF',@(x) ismember(x,{'ON','OFF'}));
            
            p.parse(varargin{:});  
            CL          = p.Results.CL;      % also works with vector
            HoldOn      = p.Results.HoldOn;
            myColor     = p.Results.Color;
            myLineStyle = p.Results.LineStyle;
            SavePlot    = p.Results.SavePlot;
            RAA         = p.Results.RAA;
            Mainz       = p.Results.Mainz;
            Troitsk     = p.Results.Troitsk;
            Neutrino4   = p.Results.Neutrino4;
            Prospect    = p.Results.Prospect;
            DANSS       = p.Results.DANSS;
            DayaBay     = p.Results.DayaBay;
            DoubleChooz = p.Results.DoubleChooz;
            Stereo      = p.Results.Stereo;
            Style       = p.Results.Style;
            FinalSensitivity = p.Results.FinalSensitivity;
            
            [DeltamNu41Sq,sin2T4Sq] = obj.Convert2Osci;
            
            if strcmp(HoldOn,'ON')
                hold on;
            elseif strcmp(HoldOn,'OFF')
               pHandle =  figure('Units','normalized','Position',[0.1,0.1,0.382,0.68]);%0.618]);
            end
            
            obj.DeltaChi2 = GetDeltaChi2(CL,2);
            legHandle = cell(0,0);
            legStr = '';
            savedirOther = [getenv('SamakPath'),'SterileAnalysis/GridSearchFiles/Knm1/Others/'];
            %% Mainz
            if strcmp(Mainz,'ON')
                filenameMainz = sprintf('%scoord_Mainz_95CL.mat',savedirOther);
                dMainz = importdata(filenameMainz);
                pMainz = plot(dMainz.SinSquare2Theta_X,dMainz.DmSquare41_Y,'-.','LineWidth',1.5,'Color',rgb('Salmon'));
                legHandle{numel(legHandle)+1} = pMainz;
                legStr = [legStr,{sprintf('Mainz 95%% C.L.')}];
                hold on;
            end
            %% Troitsk
            if strcmp(Troitsk,'ON')
                filenameTroitsk = sprintf('%scoord_Troitsk_95CL.mat',savedirOther);
                dTroitsk = importdata(filenameTroitsk);
                pTroitsk = plot(dTroitsk.SinSquare2Theta_X,dTroitsk.DmSquare41_Y,'--','LineWidth',1.5,'Color',rgb('Black'));
                legHandle{numel(legHandle)+1} = pTroitsk;
                legStr = [legStr,{sprintf('Troitsk 95%% C.L.')}];
                hold on;
            end
            %% Prospect
            if strcmp(Prospect,'ON')
               % outdated filenameProspect = sprintf('%scoord_Prospect_95CL.mat',savedirOther);
                filenameProspect = sprintf('%scoord_Prospect2020_95CL.mat',savedirOther);
                dProspect = importdata(filenameProspect);
                pProspect = plot(dProspect.SinSquare2Theta_X,dProspect.DmSquare41_Y,'-','LineWidth',1,'Color',rgb('PowderBlue'));
                legHandle{numel(legHandle)+1} = pProspect;
                legStr = [legStr,{sprintf('Prospect 95%% C.L.')}];
                hold on;
            end
            %% DANSS
            if strcmp(DANSS,'ON')
                filenameDANSS = sprintf('%scoord_DANSS_95CL.mat',savedirOther);
                dDANSS = importdata(filenameDANSS);
                pDANSS = plot(dDANSS.SinSquare2Theta_X,dDANSS.DmSquare41_Y,'-','LineWidth',1,'Color',rgb('LightGray'));
                legHandle{numel(legHandle)+1} = pDANSS;
                legStr = [legStr,{sprintf('DANSS 95%% C.L.')}];
                hold on;
            end
            %% DayaBay
            if strcmp(DayaBay,'ON')
                filenameDayaBay = sprintf('%scoord_DayaBay1230_90CL.mat',savedirOther);
                dDayaBay = importdata(filenameDayaBay);
                pDayaBay = plot(dDayaBay.SinSquare2Theta_X,dDayaBay.DmSquare41_Y,'-','LineWidth',1,'Color',rgb('IndianRed'));
                legHandle{numel(legHandle)+1} = pDayaBay;
                legStr = [legStr,{sprintf('Daya Bay 90%% C.L.')}];
                hold on;
            end
            %% DoubleChooz
            if strcmp(DoubleChooz,'ON')
                filenameDoubleChooz = sprintf('%scoord_DoubleChooz5y_95CL.mat',savedirOther);
                dDoubleChooz = importdata(filenameDoubleChooz);
                pDoubleChooz = plot(dDoubleChooz.SinSquare2Theta_X,dDoubleChooz.DmSquare41_Y,'-','LineWidth',1,'Color',rgb('CadetBlue'));
                legHandle{numel(legHandle)+1} = pDoubleChooz;
                legStr = [legStr,{sprintf('Double Chooz 95%% C.L.')}];
                hold on;
            end
            %% Stereo
            if strcmp(Stereo,'ON')
            % outdated    filenameStereo = sprintf('%scoord_Stereo_95CL.mat',savedirOther); % old stereo 2019
                filenameStereo = sprintf('%scoord_STEREOprd102_95CL.mat',savedirOther); % new stereo 2020
                dStereo = importdata(filenameStereo);
                pStereo = plot(dStereo.SinSquare2Theta_X,dStereo.DmSquare41_Y,'-','LineWidth',1,'Color',rgb('Orange'));

                legHandle{numel(legHandle)+1} = pStereo;
                legStr = [legStr,{sprintf('STEREO 95%% C.L.')}];
                hold on;
            end
            %% RAA
            if strcmp(RAA,'ON')
                hold on;
                filenameRAA1 = sprintf('%scoord_RAA_95_A.mat',savedirOther);
                filenameRAA2 = sprintf('%scoord_RAA_95_B.mat',savedirOther);
                dRAA1 = importdata(filenameRAA1,'file');
                dRAA2 = importdata(filenameRAA2,'file');  
                pRAA = plot([dRAA1.sith4_X(1),dRAA1.sith4_X,dRAA1.sith4_X(end)],...
                         [1e5,dRAA1.m4_Y,1e5],'-','LineWidth',2,'Color',rgb('ForestGreen'));
                plot(dRAA2.sith4_X,dRAA2.m4_Y,'-','LineWidth',pRAA.LineWidth,'Color',pRAA.Color);
                legHandle{numel(legHandle)+1} = pRAA;
                legStr = [legStr,{sprintf('RAA + GA 95%% CL')}];%-PRD 83, 073006 (2011) -
            end
            %% Neutrino 4
            if strcmp(Neutrino4,'ON')
                filenameN4 = sprintf('%scoord_Neutrino4_123sigma.mat',savedirOther);
                dN4 = importdata(filenameN4);
                pN4 = plot(dN4.SinSquare2Theta_X_2sigma,dN4.DmSquare41_Y_2sigma,'-','LineWidth',1.5,'Color',rgb('FireBrick'));
                legHandle{numel(legHandle)+1} = pN4;
                legStr = [legStr,{sprintf('Neutrino-4 2\\sigma')}];
                hold on;
            end
            %% KATRIN
            PlotArg = {'LineWidth',2.5,'LineStyle',myLineStyle};
            if numel(CL)==1
                PlotArg = [PlotArg,{'LineColor',myColor}];
            end
            [~,legHandle{numel(legHandle)+1}]= contour(sin2T4Sq,DeltamNu41Sq,obj.chi2-obj.chi2_ref,...
                [obj.DeltaChi2 obj.DeltaChi2],...
                PlotArg{:});
             legStr = [legStr,{sprintf('KATRIN %.0f%% C.L.',obj.ConfLevel)}];
             %(%s)',obj.ConfLevel,obj.GetPlotTitle('Mode','chi2'))}];
            
             %% final sensitivity
            if strcmp(FinalSensitivity,'ON')
                DataType_i = obj.RunAnaObj.DataType;
                AngTF = obj.RunAnaObj.AngularTFFlag;
                ElossFlag = obj.RunAnaObj.ELossFlag;
                Budget = obj.RunAnaObj.SysBudget;
                FakeInitFile = @ref_KSNX_KATRIN_Final;
                obj.RunAnaObj.DataType = 'Fake';
                obj.RunAnaObj.DataSet = 'Knm2';
                
                obj.RunAnaObj.FakeInitFile = FakeInitFile;
                obj.RunAnaObj.fixPar = 'E0 Norm Bkg';
                obj.RunAnaObj.InitFitPar;
                obj.RunAnaObj.pullFlag = 99;
                obj.RunAnaObj.AngularTFFlag = 'ON';
                obj.RunAnaObj.ELossFlag = 'KatrinT2A20';
                obj.RunAnaObj.SysBudget = 66;
                obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
                obj.Interp1Grid('RecomputeFlag','ON');
                
                PlotArg = {'LineWidth',2,'LineStyle',':','LineColor',rgb('SteelBlue')};
                [DeltamNu41Sq,sin2T4Sq] = obj.Convert2Osci;
                [~,legHandle{numel(legHandle)+1}]= contour(sin2T4Sq,DeltamNu41Sq,obj.chi2-obj.chi2_ref,...
                    [obj.DeltaChi2 obj.DeltaChi2],...
                    PlotArg{:});
               
                obj.RunAnaObj.DataType= DataType_i ;
                obj.RunAnaObj.DataSet = 'Knm1';
                obj.RunAnaObj.AngularTFFlag = AngTF;
                obj.RunAnaObj.ELossFlag = ElossFlag;
                obj.RunAnaObj.SysBudget = Budget;
                legStr = [legStr,'Projected KATRIN final',sprintf('sensitivity %.0f%% C.L.',obj.ConfLevel)];
                
                pNone = plot(NaN,NaN,'Color',rgb('White'));
                legHandle{numel(legHandle)+1} = pNone;
                
            end
            %% leg + appearance
            set(gca,'YScale','log');
            set(gca,'XScale','log');
        %   xlabel(sprintf('{sin}^2(2\\theta_{ee})'));
             xlabel(sprintf('{sin}^2(2\\theta_{14})'));
            ylabel(sprintf('\\Delta{\\itm}_{41}^2 (eV^2)'));    
            leg = legend([legHandle{:}],legStr{:},'EdgeColor','none','Location','northoutside');
            
            if strcmp(Style,'Reg')
                PrettyFigureFormat('FontSize',22);
                leg.FontSize = 12;
            else
                PRLFormat;
                set(gca,'FontSize',20);
                leg.FontSize = 16;%13;
            end
           
          xlim([1.2e-02 1]);
          if obj.range==65
                ylim([0.1 6e3]);       
            elseif obj.range==95
                ylim([0.1 1e4]);
            elseif obj.range==40
                ylim([0.1 2e3]);   
                yticks([1e-01 1e0 1e1 1e2 1e3])
          end

          if numel(legStr)>4
              leg.NumColumns = 2;
          end
          
          if strcmp(FinalSensitivity,'ON')
              xlim([6e-03 1]);
              extraStr = '_FinalSensitivity';
          else
              extraStr = '';
          end
          
        
          if ~strcmp(SavePlot,'OFF')
              if strcmp(SavePlot,'ON')
                  plotname = sprintf('%s_OsciContour_%.2gCL%s.pdf',obj.DefPlotName,obj.ConfLevel,extraStr);
                  export_fig(gcf,plotname);
              elseif strcmp(SavePlot,'png')
                  plotname = sprintf('%s_OsciContour_%.2gCL%s.png',obj.DefPlotName,obj.ConfLevel,extraStr);
                  print(gcf,plotname,'-dpng','-r450');
              end
              fprintf('save plot to %s \n',plotname);
          end
        end
        function GridPlot(obj,varargin)
            p = inputParser;
            p.addParameter('CL',obj.ConfLevel,@(x)isfloat(x));
            p.addParameter('HoldOn','OFF',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('BestFit','OFF',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('Contour','OFF',@(x)ismember(x,{'ON','OFF'}));  
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','png'}));
            p.addParameter('ExtraStr','',@(x)ischar(x));
            p.addParameter('xLim','',@(x)isfloat(x));
            
            p.parse(varargin{:});
            CL       = p.Results.CL;
            HoldOn   = p.Results.HoldOn;
            BestFit  = p.Results.BestFit;
            Contour  = p.Results.Contour;
            SavePlot = p.Results.SavePlot;
            ExtraStr = p.Results.ExtraStr;
            xLim     = p.Results.xLim;
            
            if CL<1
                CL = CL*100;
            end
            obj.DeltaChi2 = GetDeltaChi2(CL,2);
            chi2grid = obj.chi2;
            chi2grid((chi2grid-obj.chi2_ref)>obj.DeltaChi2) =  NaN;%DeltaChi2+chi2_ref;% NaN;
            zlimMax = obj.DeltaChi2;
            
            if strcmp(HoldOn,'ON')
                hold on
            else
            GetFigure;
            end
            surf(obj.sin2T4,obj.mNu4Sq,chi2grid-obj.chi2_ref,'EdgeColor','interp','FaceColor','interp');
            %% best fit
            if strcmp(BestFit,'ON')
                obj.FindBestFit;
                hold on;
                pbf = plot3(obj.sin2T4_bf,obj.mNu4Sq_bf,obj.DeltaChi2,...
                    'x','MarkerSize',9,'Color',rgb('White'),'LineWidth',3);
            end
            
            PrettyFigureFormat('FontSize',22)
            %% contour
            if strcmp(Contour,'ON')
                hold on;
                [~,pContour] = contour(obj.sin2T4,obj.mNu4Sq,obj.chi2-obj.chi2_ref,...
                    [obj.DeltaChi2 obj.DeltaChi2],...
                    'LineWidth',3,'LineStyle','-','Color',rgb('Black'));
            end
               
            if strcmp(Contour,'ON')  && strcmp(BestFit,'ON')
                leg = legend([pContour,pbf],sprintf('%.3g%% C.L.',CL),'Best fit',...
                   'Location','southwest','FontSize',get(gca,'FontSize')); 
               PrettyLegendFormat(leg,'alpha',0.6);
            elseif strcmp(Contour,'ON')
                 leg = legend(pContour,sprintf('%.3g%% C.L.',CL),...
                   'Location','southwest','FontSize',get(gca,'FontSize'));   
               PrettyLegendFormat(leg,'alpha',0.6);
            end
            
            zlim([0 zlimMax])
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            c =colorbar;
            c.Label.String = sprintf('\\Delta\\chi^2');
            c.Label.FontSize = get(gca,'FontSize')+2;
            c.Limits=[0 zlimMax];
            xlabel(sprintf('|{\\itU}_{e4}|^2'));
            ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
            zlabel(sprintf('\\Delta\\chi^2'))
            grid off
            view([0 0 1])
            
            
            ylim([min(min(obj.mNu4Sq)),max(max(obj.mNu4Sq))]);
            if isempty(xLim)
                xlim([min(min(obj.sin2T4)), max(max(obj.sin2T4))]);
            else
                xlim(xLim);
            end
            
            title(sprintf('%s',obj.GetPlotTitle),'FontWeight','normal','FontSize',get(gca,'FontSize'));
            
         if obj.range == 65
             xlim([2e-03,0.5]);
         end
             %% save
           if ~strcmp(SavePlot,'OFF')
               name_i = obj.DefPlotName;
               if strcmp(SavePlot,'ON')
                   plotname = sprintf('%s_GridPlot_%.2gCL%s.pdf',name_i,obj.ConfLevel,ExtraStr);
                   export_fig(gcf,plotname);
               elseif strcmp(SavePlot,'png')
                   plotname = sprintf('%s_GridPlot_%.2gCL%s.png',name_i,obj.ConfLevel,ExtraStr);
                   print(gcf,plotname,'-dpng','-r450');
               end
               fprintf('save plot to %s \n',plotname);
           end
           
           
        end
        function GridPlotFitPar(obj,varargin)
             p = inputParser;
            p.addParameter('CL',obj.ConfLevel,@(x)isfloat(x));
            p.addParameter('HoldOn','OFF',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('BestFit','OFF',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('Contour','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('ContourTxt','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','png'}));
            p.addParameter('FitPar','mNuSq',@(x)ismember(x,{'mNuSq','E0'}));
            p.addParameter('ContourVec','',@(x)isfloat(x) || isempty(x));
            p.parse(varargin{:});
            CL       = p.Results.CL;
            HoldOn   = p.Results.HoldOn;
            BestFit  = p.Results.BestFit;
            Contour  = p.Results.Contour;
            SavePlot = p.Results.SavePlot;
            FitPar   = p.Results.FitPar;     
           ContourVec = p.Results.ContourVec;
           ContourTxt = p.Results.ContourTxt;
           
            if strcmp(HoldOn,'ON')
                hold on
            else
            GetFigure;
            end
      
            if strcmp(FitPar,'mNuSq')
                if isempty(ContourVec)
                    if strcmp(obj.RunAnaObj.DataSet,'Knm1')
                        ContourVec = [-10 -2 -1 -0.5 1 10];
                    else
                        ContourVec = [-10 -2 -1 0 0.5 1 5 10,50];
                    end
                end
                PlotPar = obj.mNuSq;
                Plot_bf = obj.mNuSq_bf;
                zStr = sprintf('{\\itm}_\\nu^2 (eV^2)');
            elseif strcmp(FitPar,'E0')
                PlotPar = obj.E0-obj.E0_bf;
                Plot_bf = obj.E0_bf;
                if isempty(ContourVec)
                    ContourVec = [-0.2 -0.5 0 0.5];%+obj.RunAnaObj.ModelObj.Q_i;
                end
                zStr = sprintf('{\\itE}_0 - \\langle{\\itE}_0^{bf}\\rangle (eV)');
            end
            
            surf(obj.sin2T4,obj.mNu4Sq,PlotPar,'EdgeColor','interp','FaceColor','interp');
            colormap('jet')
            %% best fit
            %% contour
           
            if strcmp(Contour,'ON')
                hold on;
                if strcmp(ContourTxt,'ON')
                    [M,pContour] = contour3(obj.sin2T4,obj.mNu4Sq,PlotPar,...
                        [ContourVec],...
                        'LineWidth',1.5,'LineStyle','-','Color',rgb('Black'),...
                        'ShowText','on');
                else
%                     for i=1:numel(ContourVec)
%                         [M,pContour] = contour3(obj.sin2T4,obj.mNu4Sq,PlotPar,...
%                             [ContourVec(i) ContourVec(i)],...
%                             'LineWidth',1.5,'LineStyle','-','Color',rgb('Black'),...
%                             'ShowText','on');
%                         %InclIdx = find(M(1,:)~=ContourVec(i));% & M(1,:)>0);
%                         % pC = plot3(M(1,InclIdx),M(2,InclIdx),99.*ones(size(InclIdx,2),1),'-','LineWidth',1.5,'LineStyle','-','Color',rgb('Black'));
%                         
%                         
%                        % NewContourIdx = find(M(1,:)==ContourVec(i));
%                         
% %                         for j=1:numel(NewContourIdx)-1
% %                             nPoints = NewContourIdx(j+1)-NewContourIdx(j)-1;
% %                             plot3(M(1,NewContourIdx(j)+1:NewContourIdx(j+1)-1),M(2,NewContourIdx(j)+1:NewContourIdx(j+1)-1),1e5.*ones(nPoints,1),...
% %                                 '-','LineWidth',1.5,'LineStyle','-','Color',rgb('Black'));
% %                         end
% %                         pContour.delete;
%                     end
                    
                    [M,pContour] = contour3(obj.sin2T4,obj.mNu4Sq,PlotPar,...
                        [ContourVec],...
                        'LineWidth',1.5,'LineStyle','-','Color',rgb('Black'),...
                        'ShowText','on');
                  %  t =clabel(M,'FontSize',12,'Rotation',-45,'LineWidth',2);
                end
            end
            
            if strcmp(BestFit,'ON')
                obj.FindBestFit;
                hold on;
                pbf = plot3(obj.sin2T4_bf,obj.mNu4Sq_bf,Plot_bf+1e5,...
                    'x','MarkerSize',12,'Color',rgb('White'),'LineWidth',3);
            end
            
            PrettyFigureFormat('FontSize',22)
            
            
            if strcmp(Contour,'ON')  && strcmp(BestFit,'ON')
                leg = legend([pbf],'Best fit',...
                    'EdgeColor','none','Location','southwest','Color',rgb('White'),...
                    'FontSize',get(gca,'FontSize'));
                set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.3]));
            end
            
            set(gca,'XScale','log')
            set(gca,'YScale','log')     
            c =colorbar;
            c.Label.String = zStr;
            c.Label.FontSize = get(gca,'FontSize')+2;
            xlabel(sprintf('|{\\itU}_{e4}|^2'));
            ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
            zlabel(zStr)
            grid off
            view([0 0 1])
         
            
            ylim([1,max(max(obj.mNu4Sq))]);
            xlim([min(min(obj.sin2T4)), max(max(obj.sin2T4))]);
            
            title(sprintf('%s',obj.GetPlotTitle),'FontWeight','normal','FontSize',get(gca,'FontSize'));
            
         if obj.range == 65
             xlim([2e-03,0.5]);
         end
             %% save
           if ~strcmp(SavePlot,'OFF')
               name_i = obj.DefPlotName;
               if strcmp(ContourTxt,'OFF')
                   name_i = [name_i,'_TxtOff'];
               end
               if strcmp(SavePlot,'ON')
                   plotname = sprintf('%s_GridPlot_%.2gCL.pdf',name_i,obj.ConfLevel);
                   export_fig(gcf,plotname);
               elseif strcmp(SavePlot,'png')
                   plotname = sprintf('%s_GridPlot_%.2gCL.png',name_i,obj.ConfLevel);
                   print(gcf,plotname,'-dpng','-r450');
               end
               fprintf('save plot to %s \n',plotname);
           end
        end
        function IsoPlotFitPar(obj,varargin)
             p = inputParser;
            p.addParameter('CL',obj.ConfLevel,@(x)isfloat(x));
            p.addParameter('HoldOn','OFF',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('BestFit','OFF',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','png'}));
            p.addParameter('FitPar','mNuSq',@(x)ismember(x,{'mNuSq','E0'}));
            p.addParameter('ContourVec','',@(x)isfloat(x) || isempty(x));
            p.parse(varargin{:});
            CL       = p.Results.CL;
            HoldOn   = p.Results.HoldOn;
            SavePlot = p.Results.SavePlot;
            FitPar   = p.Results.FitPar;     
           ContourVec = p.Results.ContourVec;
           BestFit = p.Results.BestFit;
           
            if strcmp(HoldOn,'ON')
                hold on
            else
            GetFigure;
            end
      
            if strcmp(FitPar,'mNuSq')
                if isempty(ContourVec)
                    if strcmp(obj.RunAnaObj.DataSet,'Knm1')
                        ContourVec = [-10 -2 -1 -0.5 1 10];
                    else
                        ContourVec = [-10 -2 -1 0 0.5 1 5 10,50];
                    end
                end
                PlotPar = obj.mNuSq;
                Plot_bf = obj.mNuSq_bf; 
                LegStr = sprintf('{\\itm}_\\nu^2');
            elseif strcmp(FitPar,'E0')
                PlotPar = obj.E0-obj.E0_bf;
                Plot_bf = obj.E0_bf;
                if isempty(ContourVec)
                    ContourVec = [-0.2 -0.5 0 0.5];%+obj.RunAnaObj.ModelObj.Q_i;
                end 
                LegStr = sprintf('{\\itE}_0^{fit}');
            end
            
            %% best fit
            %% contour
            [M,pContour] = contour3(obj.sin2T4,obj.mNu4Sq,PlotPar,...
                [ContourVec],...
                'LineWidth',1.5,'LineStyle','-','Color',rgb('Black'),...
                'ShowText','on');
            clabel(M,pContour,'FontSize',15);
            
            if strcmp(BestFit,'ON')
                obj.FindBestFit;
                hold on;
                pbf = plot3(obj.sin2T4_bf,obj.mNu4Sq_bf,Plot_bf+1e5,...
                    'x','MarkerSize',12,'Color',rgb('White'),'LineWidth',3);
            end
            
            PrettyFigureFormat('FontSize',22)

            if strcmp(BestFit,'ON')
                leg = legend([pContour,pbf],sprintf('Isoline %s best fit',LegStr),'Best fit',...
                    'Location','northwest', 'FontSize',get(gca,'FontSize'));
            else
                leg = legend(pContour,sprintf('Isoline %s best fit',LegStr),...
                    'Location','northwest', 'FontSize',get(gca,'FontSize'));
            end
            PrettyLegendFormat(leg,'alpha',0.95);
            
            set(gca,'XScale','log')
            set(gca,'YScale','log')       
            xlabel(sprintf('|{\\itU}_{e4}|^2'));
            ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
            grid off
            view([0 0 1])
         
            
            ylim([1,max(max(obj.mNu4Sq))]);
            xlim([min(min(obj.sin2T4)), max(max(obj.sin2T4))]);
            
            title(sprintf('%s',obj.GetPlotTitle),'FontWeight','normal','FontSize',get(gca,'FontSize'));
            
         if obj.range == 65
             xlim([2e-03,0.5]);
         end
             %% save
           if ~strcmp(SavePlot,'OFF')
               name_i = obj.DefPlotName;
               if strcmp(ContourTxt,'OFF')
                   name_i = [name_i,'_TxtOff'];
               end
               if strcmp(SavePlot,'ON')
                   plotname = sprintf('%s_IsoPlot_%.2gCL.pdf',name_i,obj.ConfLevel);
                   export_fig(gcf,plotname);
               elseif strcmp(SavePlot,'png')
                   plotname = sprintf('%s_IsoPlot_%.2gCL.png',name_i,obj.ConfLevel);
                   print(gcf,plotname,'-dpng','-r450');
               end
               fprintf('save plot to %s \n',plotname);
           end
        end
        function PlotqUScan(obj,varargin)
         % plot contours with same settings but different fit ranges
         p = inputParser;
         p.addParameter('Ranges',[95:-5:45,41,40],@(x)isfloat(x));
         p.addParameter('SavePlot','ON',@(x)ismember(x,{'ON','OFF','png'}));
         p.addParameter('BestFit','OFF',@(x)ismember(x,{'ON','OFF'}));
         p.parse(varargin{:});
         Ranges   = p.Results.Ranges;
         SavePlot = p.Results.SavePlot;
         BestFit  = p.Results.BestFit;
         
         legStr = cell(numel(Ranges),1);
         pl     = cell(numel(Ranges),1);
         range_i = obj.range;
         
         if numel(Ranges)>3
             Colors = parula(numel(Ranges));
         else
             Colors = cell2mat(obj.PlotColors');
             
         end
         
         for i=1:numel(Ranges)
             progressbar(i/numel(Ranges));
             obj.range = Ranges(i);
             obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
             obj.Interp1Grid('RecomputeFlag','ON');
             PlotArg = {'Color',Colors(i,:),'LineStyle',obj.PlotLines{i},'BestFit',BestFit};
             if i>1
                 pl{i} = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',PlotArg{:});
             else
                 pl{i} = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','OFF',PlotArg{:});
             end
             legStr{i} = sprintf('%.0f eV range',Ranges(i));
         end
         
         leg = legend([pl{:}],legStr{:},'EdgeColor',rgb('Silver'),'Location','southwest');
         %leg.Title.String = 'Lower fit boundary';
         %leg.Title.FontWeight = 'normal';
         ylim([1 1e4])
         if numel(Ranges)>5 &&numel(Ranges)<10
             leg.NumColumns=2;
         elseif numel(Ranges)>=10
              leg.NumColumns=3; 
         end
%          xlim([5e-03,0.4]);
%          ylim([1 3e4]);
          set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.3]));
          
         title(sprintf('%s (%s) %.0f%% C.L.',obj.GetPlotTitle('Mode','data'),obj.GetPlotTitle('Mode','chi2'),obj.ConfLevel),...
             'FontWeight','normal','FontSize',get(gca,'FontSize'));

             if ~strcmp(SavePlot,'OFF')
                 name_i = strrep(obj.DefPlotName,sprintf('_%.0feVrange',Ranges(end)),'');
                 if strcmp(SavePlot,'ON')
                     plotname = sprintf('%s_qUScan_%.2gCL.pdf',name_i,obj.ConfLevel);
                     export_fig(gcf,plotname);
                 elseif strcmp(SavePlot,'png')
                     plotname = sprintf('%s_qUScan_%.2gCL.png',name_i,obj.ConfLevel);
                     print(gcf,plotname,'-dpng','-r450');
                 end
                 fprintf('save plot to %s \n',plotname);
             end
             
             obj.range = range_i;
        end    
        function [pFree,pFix] = PlotmNuSqOverview(obj,varargin)
            p = inputParser;
            p.addParameter('BestFit','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('PullmNuSq','ON',@(x)ismember(x,{'ON','OFF'})); % show also contour with constrained m^2
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','png'}));
            p.addParameter('HoldOn','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            BestFit  = p.Results.BestFit;
            SavePlot = p.Results.SavePlot;
            PullmNuSq = p.Results.PullmNuSq;
            HoldOn   = p.Results.HoldOn;
            
            fixPar_i = obj.RunAnaObj.fixPar;
            pull_i = obj.RunAnaObj.pullFlag;

            %% 1. nuissance nu-mass without pull
            obj.RunAnaObj.fixPar = 'mNu E0 Norm Bkg'; obj.RunAnaObj.InitFitPar;
            obj.RunAnaObj.pullFlag = 99;
            obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
            obj.Interp1Grid('RecomputeFlag','ON');
            pFree = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn',HoldOn,...
                'Color',rgb('ForestGreen'),'LineStyle','-','BestFit',BestFit);

            %%  nuissance nu-mass + pull
            if strcmp(PullmNuSq,'ON')
                obj.RunAnaObj.pullFlag = 12;
                obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
                obj.Interp1Grid('RecomputeFlag','ON');
                pPull = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                    'Color',rgb('Orange'),'LineStyle','-.','BestFit',BestFit);
                
                if strcmp(obj.RunAnaObj.DataSet,'Knm2')
                     obj.RunAnaObj.pullFlag = 26;
                obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
                obj.Interp1Grid('RecomputeFlag','ON');
                pPullK = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                    'Color',rgb('FireBrick'),'LineStyle','--','BestFit',BestFit);
                end
            end
            
            %% fixed nu-mass
            obj.RunAnaObj.fixPar = 'E0 Norm Bkg'; obj.RunAnaObj.InitFitPar;
            obj.RunAnaObj.pullFlag = 99;
            if strcmp(obj.RunAnaObj.DataSet,'Knm2') &&  strcmp(obj.RunAnaObj.DataType,'Real')
                % best fit is at sin2t4=1, but grid looks better for sin2t4 up to 0.5
                % load first sin2t4=1 file, get best fit, then load regular file
                obj.LoadGridFile(obj.LoadGridArg{:},'Extsin2T4','ON','IgnoreKnm2FSDbinning','ON');
                obj.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',34^2);
                obj.FindBestFit;
                chi2_ref_ExtSin2T4 = obj.chi2_ref;
                mNuSq_bf_ExtSin2T4 = obj.mNuSq_bf;
                sin2T4_bf_ExtSin2T4 = obj.sin2T4_bf;
                mNu4Sq_bf_ExtSin2T4 = obj.mNu4Sq_bf;
            end
            
            obj.LoadGridFile(obj.LoadGridArg{:});
            obj.Interp1Grid('RecomputeFlag','ON');
            

            if strcmp(obj.RunAnaObj.DataSet,'Knm2') &&  strcmp(obj.RunAnaObj.DataType,'Real')
                obj.chi2_ref  = chi2_ref_ExtSin2T4;
                obj.mNuSq_bf  = mNuSq_bf_ExtSin2T4 ;
                obj.sin2T4_bf = sin2T4_bf_ExtSin2T4;
                obj.mNu4Sq_bf = mNu4Sq_bf_ExtSin2T4;
            else
                obj.FindBestFit;
            end
            
            [pFix,pfixmin] = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                'Color',rgb('DodgerBlue'),'LineStyle',':','BestFit',BestFit,'ReCalcBF','OFF');
            
            PrettyFigureFormat('FontSize',22);
            if strcmp(PullmNuSq,'ON') && strcmp(obj.RunAnaObj.DataSet,'Knm1')
                leg = legend([pFree,pPull,pFix],...
                    sprintf('Free {\\itm}_\\nu^2 unconstrained'),...
                    sprintf('Free {\\itm}_\\nu^2 with pull term \\sigma({\\itm}_\\nu^2) = 1.94 eV^2'),...
                    sprintf('Fixed {\\itm}_\\nu^2 = 0 eV^2'),...
                    'Location','southwest');
            elseif strcmp(PullmNuSq,'ON') && strcmp(obj.RunAnaObj.DataSet,'Knm2')
                  leg = legend([pFree,pPull,pPullK,pFix],...
                    sprintf('Free {\\itm}_\\nu^2 unconstrained'),...
                    sprintf('Free {\\itm}_\\nu^2 with pull term \\sigma({\\itm}_\\nu^2) = 1.94 eV^2 (Mainz/Troitsk)'),...
                    sprintf('Free {\\itm}_\\nu^2 with pull term \\sigma({\\itm}_\\nu^2) = 1.1 eV^2 (KATRIN KNM-1)'),...
                    sprintf('Fixed {\\itm}_\\nu^2 = 0 eV^2'),...
                    'Location','southwest');
            else
                leg =  legend([pFree,pFix],...
                    sprintf('Free {\\itm}_\\nu^2 without pull term'),...   
                    sprintf('Fixed {\\itm}_\\nu^2 = 0 eV^2'),...
                    'Location','southwest');
            end
            
                PrettyLegendFormat(leg);
                legend boxoff
                obj.RunAnaObj.fixPar = fixPar_i;
                obj.RunAnaObj.pullFlag = pull_i ;
               
                if obj.range==65
                    ylim([1 6e3]);
                    xlim([2e-03 0.5]);
                else
                    xlim([pfixmin*0.5,0.5]);
                end
                %% save
                if ~strcmp(SavePlot,'OFF')
                    name_i = strrep(obj.DefPlotName,'_mNuE0BkgNorm','');
                    if strcmp(SavePlot,'ON')
                        plotname = sprintf('%s_mNuSqOverview_%.2gCL.pdf',name_i,obj.ConfLevel);
                        export_fig(gcf,plotname);
                    elseif strcmp(SavePlot,'png')
                        plotname = sprintf('%s_mNuSqOverview_%.2gCL.png',name_i,obj.ConfLevel);
                        print(gcf,plotname,'-dpng','-r450');
                    end
                    fprintf('save plot to %s \n',plotname);
                end
                
        end
        function TestCoverageImpact(obj,varargin)
            p = inputParser;
            p.addParameter('BestFit','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','png'}));
            p.parse(varargin{:});
            BestFit  = p.Results.BestFit;
            SavePlot = p.Results.SavePlot;
            fixPar_i = obj.RunAnaObj.fixPar;
            pull_i   = obj.RunAnaObj.pullFlag;

            %%  95CL - Wilks Theorem
            obj.RunAnaObj.fixPar = 'E0 Norm Bkg'; obj.RunAnaObj.InitFitPar;
            obj.RunAnaObj.pullFlag = 99;
            obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
            obj.Interp1Grid('RecomputeFlag','ON');
            pFix1 = obj.ContourPlot('CL',95,'HoldOn','OFF',...
                'Color',rgb('Orange'),'LineStyle','-.','BestFit',BestFit);
            
            %% 95CL - Wilks Theorem Corrected wia MC simulation
            obj.RunAnaObj.fixPar = 'E0 Norm Bkg'; obj.RunAnaObj.InitFitPar;
            obj.RunAnaObj.pullFlag = 99;
            obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
            obj.Interp1Grid('RecomputeFlag','ON');
            
            if obj.range==95
                CL_coverage = 96.45;
            elseif obj.range==65   
                CL_coverage = 95.7;
            elseif obj.range==40   
                CL_coverage = 94.82;
            end
            DeltaChi2Crit = GetDeltaChi2(CL_coverage,2);
            pFix2 = obj.ContourPlot('CL',CL_coverage,'HoldOn','ON',...
                'Color',rgb('DodgerBlue'),'LineStyle',':','BestFit',BestFit);
            
            PrettyFigureFormat('FontSize',22);
                legend([pFix1,pFix2],...
                    sprintf('Fixed {\\itm}_\\nu^2 = 0 eV^2 - \\Delta\\chi^2 = 5.99'),...
                    sprintf('Fixed {\\itm}_\\nu^2 = 0 eV^2 - \\Delta\\chi^2 = %.2f',DeltaChi2Crit),...
                    'EdgeColor',rgb('Silver'),'Location','southwest');
                obj.RunAnaObj.fixPar = fixPar_i;
                obj.RunAnaObj.pullFlag = pull_i ;
               
                if obj.range==65
                    ylim([1 6e3]);
                    xlim([2e-03 0.5]);
                elseif obj.range==95
                    ylim([1 1e4]);
                    xlim([1e-03 0.5]);
                elseif obj.range==40
                    ylim([1 2e3]);
                    xlim([1e-02 0.5]);
                end
                %% save
                if ~strcmp(SavePlot,'OFF')
                    name_i = strrep(obj.DefPlotName,'_mNuE0BkgNorm','');
                    if strcmp(SavePlot,'ON')
                        plotname = sprintf('%s_TestCoverageImpact_%.2gCL.pdf',name_i,obj.ConfLevel);
                        export_fig(gcf,plotname);
                    elseif strcmp(SavePlot,'png')
                        plotname = sprintf('%s_TestCoverageImpact_%.2gCL.png',name_i,obj.ConfLevel);
                        print(gcf,plotname,'-dpng','-r450');
                    end
                    fprintf('save plot to %s \n',plotname);
                end
                
        end      
        function PlotTwinData(obj,varargin)
            p=inputParser;
            p.addParameter('BestFit','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            BestFit  = p.Results.BestFit;
            SavePlot = p.Results.SavePlot;
            
            DataType_i = obj.RunAnaObj.DataType;
            
            % twins
            obj.RunAnaObj.DataType = 'Twin';
            obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
            obj.Interp1Grid('RecomputeFlag','ON');
            [pTwin,sin2T4_Twinmin] = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','OFF',...
                'Color',rgb('Orange'),'LineStyle','-.','BestFit','OFF');
            
            % data
            obj.RunAnaObj.DataType = 'Real';
            if strcmp(obj.RunAnaObj.DataSet,'Knm2') && contains(obj.RunAnaObj.fixPar,'fix 1 ;')
                InterpMode_i = obj.InterpMode;
                obj.InterpMode = 'Mix';
                obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:},'Extsin2T4','ON');
                obj.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',36^2);
                obj.InterpMode = InterpMode_i;
            else
                obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
                obj.Interp1Grid('RecomputeFlag','ON');
            end
            
            [pData,sin2T4_Datamin] = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                'Color',rgb('DodgerBlue'),'LineStyle','-','BestFit',BestFit);
            
            titleStr = sprintf('%.0f%% C.L. , %.0feV range (%s)',obj.ConfLevel,obj.range,obj.GetPlotTitle('Mode','chi2'));
            title(titleStr,'FontWeight','normal','FontSize',get(gca,'FontSize'));
            leg = legend([pData,pTwin],sprintf('Data (%s)',obj.GetPlotTitle('Mode','chi2')),...
                sprintf('Twin (%s)',obj.GetPlotTitle('Mode','chi2')),'Location','southwest');
            PrettyLegendFormat(leg);
            xlim([min([sin2T4_Datamin,sin2T4_Twinmin])*0.5,0.5]);
            
            obj.RunAnaObj.DataType = DataType_i;
            
               %% save
                if ~strcmp(SavePlot,'OFF')
                    name_i = strrep(obj.DefPlotName,sprintf('%s_',DataType_i),'');
                    if strcmp(SavePlot,'ON')
                        plotname = sprintf('%s_DataTwin_%.2gCL.pdf',name_i,obj.ConfLevel);
                        export_fig(gcf,plotname);
                    elseif strcmp(SavePlot,'png')
                        plotname = sprintf('%s_DataTwin_%.2gCL.png',name_i,obj.ConfLevel);
                        print(gcf,plotname,'-dpng','-r450');
                    end
                    fprintf('save plot to %s \n',plotname);
                end
                
        end     
        function PlotFitriumSamak(obj,varargin)
            p = inputParser;
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','png'}));
            p.addParameter('PlotStat','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('PlotTot','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('PlotKaFit','OFF',@(x)ismember(x,{'ON','OFF'}));
                   
            p.parse(varargin{:});
            SavePlot = p.Results.SavePlot;
            PlotStat = p.Results.PlotStat;
            PlotTot  = p.Results.PlotTot;
            PlotKaFit = p.Results.PlotKaFit;
             
            chi2_i   = obj.RunAnaObj.chi2;
            if strcmp(obj.RunAnaObj.DataType,'Real')
                BestFit = 'ON';
            else
                BestFit = 'OFF';
            end
        
            LineWidth = 2.5;
            %% load samak
             HoldOn = 'OFF';
            if strcmp(PlotStat,'ON')
                obj.RunAnaObj.chi2 = 'chi2Stat';
                obj.SetNPfactor;
                obj.LoadGridFile(obj.LoadGridArg{:});
                obj.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',40^2);
                pStat = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn',HoldOn,...
                    'Color',rgb('FireBrick'),'LineStyle','-','BestFit',BestFit,...
                    'PlotSplines','OFF');
                hold on;
                HoldOn = 'ON';
            end
            
            if strcmp(PlotTot,'ON')
                obj.RunAnaObj.chi2 = 'chi2CMShape';
                obj.SetNPfactor;
                obj.LoadGridFile(obj.LoadGridArg{:});
                obj.Interp1Grid('RecomputeFlag','ON');
                pSys = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn',HoldOn,...
                    'Color',rgb('FireBrick'),'LineStyle','-','BestFit',BestFit,...
                    'PlotSplines','OFF');
                hold on;
                HoldOn = 'ON';
            end
            %% load fitrium
            savedirF = [getenv('SamakPath'),'SterileAnalysis/GridSearchFiles/',obj.RunAnaObj.DataSet,'/Others/'];
            if strcmp(obj.RunAnaObj.DataSet,'Knm1')
                fstat = sprintf('%scontour_KSN1_Fitrium_%s_%.0feV_stat_95CL_0.txt',savedirF,obj.RunAnaObj.DataType,obj.range);
                fsys = sprintf('%scontour_KSN1_Fitrium_%s_%.0feV_total_95CL_0.txt',savedirF,obj.RunAnaObj.DataType,obj.range);
            else
                if contains(ConvertFixPar('freePar',obj.RunAnaObj.fixPar,'Mode','Reverse'),'mNu')
                    fstat = sprintf('%scontour_KSN2_Fitrium_%s_%.0feV_mNuFree_stat_95CL.dat',savedirF,obj.RunAnaObj.DataType,obj.range);
                    fsys = NaN;             
                elseif strcmp(obj.NullHypothesis,'ON')
                    fstat = sprintf('%scontour_KSN2_Fitrium_%s_%.0feV_stat_95CL_NH.dat',savedirF,obj.RunAnaObj.DataType,obj.range);
                    fsys = sprintf('%scontour_KSN2_Fitrium_%s_%.0feV_total_95CL_NH.dat',savedirF,obj.RunAnaObj.DataType,obj.range);
                 else
                     fstat = sprintf('%scontour_KSN2_Fitrium_%s_%.0feV_stat_95CL.dat',savedirF,obj.RunAnaObj.DataType,obj.range);
                     fsys = sprintf('%scontour_KSN2_Fitrium_%s_%.0feV_total_95CL.dat',savedirF,obj.RunAnaObj.DataType,obj.range);
                end
                 
              
            end
            
           if strcmp(PlotStat,'ON')
               dfStat = importdata(fstat);
               %                  if strcmp(obj.RunAnaObj.DataSet,'Knm2')
               %                      dfStat.data(:,2) = dfStat.data(:,2);
               %                  end %PowderBlue
               if strcmp(obj.RunAnaObj.DataSet,'Knm2') && strcmp(obj.RunAnaObj.DataType,'Twin') && contains(ConvertFixPar('freePar',obj.RunAnaObj.fixPar,'Mode','Reverse'),'mNu')
                   pFStat = plot(dfStat(:,1),dfStat(:,2),'LineStyle','-.','Color',rgb('Orange'),'LineWidth',LineWidth);  
               elseif strcmp(obj.RunAnaObj.DataSet,'Knm2') && strcmp(obj.RunAnaObj.DataType,'Twin')        
                   pFStat = plot(dfStat(:,1),dfStat(:,2).^2,'LineStyle','-.','Color',rgb('Orange'),'LineWidth',LineWidth);
               else
                   pFStat = plot(dfStat(:,1),dfStat(:,2),'LineStyle','-.','Color',rgb('Orange'),'LineWidth',LineWidth);
               end
           end
           
           if strcmp(PlotTot,'ON')
               dfSys = importdata(fsys);
         %      pFSys  = plot(dfSys.data(:,1),dfSys.data(:,2),'LineStyle','-.','Color',rgb('Orange'),'LineWidth',LineWidth);
                pFSys  = plot(dfSys(:,1),dfSys(:,2),'LineStyle','-.','Color',rgb('Orange'),'LineWidth',LineWidth);
             
             %  pFSys = plot(NaN,NaN,'LineStyle','-.','Color',rgb('Orange'),'LineWidth',LineWidth);
               if obj.range==95 && strcmp(obj.RunAnaObj.DataType,'Real')
                   fsys1 = sprintf('%scontour_KSN1_Fitrium_%s_%.0feV_total_95CL_1.txt',savedirF,obj.RunAnaObj.DataType,obj.range);
                   dfSys1 = importdata(fsys1);
                   fsys2 = sprintf('%scontour_KSN1_Fitrium_%s_%.0feV_total_95CL_2.txt',savedirF,obj.RunAnaObj.DataType,obj.range);
                   dfSys2 = importdata(fsys2);
                   plot(dfSys1.data(:,1),dfSys1.data(:,2),'LineStyle','-.','Color',rgb('Orange'),'LineWidth',LineWidth);
                   plot(dfSys2.data(:,1),dfSys2.data(:,2),'LineStyle','-.','Color',rgb('Orange'),'LineWidth',LineWidth);
               end
           end
           
           if strcmp(obj.RunAnaObj.DataType,'Real') && strcmp(BestFit,'ON') && strcmp(obj.RunAnaObj.DataSet,'Knm1')
               if obj.range==65    
                   if strcmp(PlotStat,'ON')
                       pF_bfStat = plot(2.532e-02,7.466e+01,'o','MarkerSize',8,'Color',pFStat.Color,'LineWidth',pFStat.LineWidth);
                   end
                   if strcmp(PlotTot,'ON')
                       pF_bfSys = plot(2.532e-02,7.466e+01,'o','MarkerSize',8,'Color',pFSys.Color,'LineWidth',pFSys.LineWidth);
                   end
               elseif obj.range==95
                   if strcmp(PlotStat,'ON')
                       pF_bfStat = plot(0.015401, 3942.813341,'o','MarkerSize',8,'Color',pFStat.Color,'LineWidth',pFStat.LineWidth);
                   end
                   if strcmp(PlotTot,'ON')
                       pF_bfSys = plot(0.013601, 3942.813341,'o','MarkerSize',8,'Color',pFSys.Color,'LineWidth',pFSys.LineWidth);
                   end
               elseif obj.range==40
                     if strcmp(PlotStat,'ON')
                      % pF_bfStat = plot(3.676e-02, 7.218e+01,'o','MarkerSize',8,'Color',pFStat.Color,'LineWidth',pFStat.LineWidth);
                       pF_bfStat = plot(0.036, 73.3,'o','MarkerSize',8,'Color',pFStat.Color,'LineWidth',pFStat.LineWidth);
                     end
                   if strcmp(PlotTot,'ON')
                       pF_bfSys = plot(0.034, 74.3,'o','MarkerSize',8,'Color',pFSys.Color,'LineWidth',pFSys.LineWidth);
                   end
               end
           end
          
           %% plot kafit
           if strcmp(PlotKaFit ,'ON')
               if contains(ConvertFixPar('freePar',obj.RunAnaObj.fixPar,'Mode','Reverse'),'mNu')
                    kstat = sprintf('%scontour_KSN2_Kafit_%s_%.0feV_mNuFree_stat_95CL.txt',savedirF,obj.RunAnaObj.DataType,obj.range);
            
               else
                   kstat = sprintf('%scontour_KSN2_Kafit_%s_%.0feV_stat_95CL.txt',savedirF,obj.RunAnaObj.DataType,obj.range);
               end
               dkStat = importdata(kstat);
               if strcmp(PlotStat,'ON')
                   if strcmp(obj.RunAnaObj.DataType,'Real')
                       pKStat = plot(dkStat(:,1),dkStat(:,2),'LineStyle',':','Color',rgb('DodgerBlue'),'LineWidth',LineWidth+0.5);
                   else
                       pKStat = plot(dkStat.data(:,1),dkStat.data(:,2),'LineStyle',':','Color',rgb('DodgerBlue'),'LineWidth',LineWidth+0.5);
                       
                   end
               end
               if strcmp(PlotTot,'ON')
                   ksyst = sprintf('%scontour_KSN2_Kafit_%s_%.0feV_statsyst_95CL.txt',savedirF,obj.RunAnaObj.DataType,obj.range);
                   dkSyst = importdata(ksyst);
                   pKSys = plot(dkSyst(:,1),dkSyst(:,2),'LineStyle',':','Color',rgb('DodgerBlue'),'LineWidth',LineWidth+0.5);
               end
           else
               pKStat = 0;
                pKSys = 0;
           end
           
           if strcmp(PlotStat,'ON') && strcmp(PlotTot,'ON')
               if strcmp(PlotKaFit ,'OFF')
               legStr = {'Samak (stat. only)','Fitrium (stat. only)','Samak (stat. and syst.)','Fitrium (stat. and syst.)'};
               legend([pStat,pFStat,pSys,pFSys],legStr,'EdgeColor',rgb('Silver'),'Location','southwest');
               else
               legStr = {'Fitrium (stat. only)','KaFit (stat. only)','Samak (stat. only)','Fitrium (stat. and syst.)','KaFit (stat. and syst.)','Samak (stat. and syst.)'};
               legend([pFStat,pKStat,pStat,pFSys,pKSys,pSys],legStr,'EdgeColor',rgb('Silver'),'Location','southwest');   
               end
               extraStr = '';
               
           elseif strcmp(PlotStat,'ON') || strcmp(PlotTot,'ON')
               if strcmp(PlotTot,'ON')
                   pF = pFSys; pK = pKSys; pS = pSys;
                   extraStr = '_Tot';
                   legTitle = sprintf('Stat. and syst.');
               elseif strcmp(PlotStat,'ON')
                   pF = pFStat; pK = pKStat; pS = pStat;
                    extraStr = '_StatOnly';
                    legTitle = sprintf('Stat. only');
               end
               
               if strcmp(PlotKaFit ,'ON')
                   legStr = {'Fitrium','KaFit','Samak'};
                   legPlt = [pF,pK,pS];   
               else
                   legStr = {'Samak','Fitrium'};
                   legPlt = [pF,pS];     
               end
                  leg = legend(legPlt,legStr,'Location','southwest'); 
                  PrettyLegendFormat(leg)
                  leg.Title.String = legTitle;
                  leg.Title.FontWeight = 'normal';
           end
           
           obj.RunAnaObj.chi2 = chi2_i;
           if obj.range==65
               xlim([4e-03 0.5])
               ylim([1 1e4])
           elseif obj.range==40
               if strcmp(obj.RunAnaObj.DataSet,'Knm1')
                   xlim([1e-02 0.5])
                   ylim([1 3e3])
               else
                   xlim([3e-03 0.5])
                   ylim([1 40^2])
               end   
            elseif obj.range==95
                xlim([3e-03 0.5])
                ylim([1 2e4]) 
            end
            
            %title(sprintf('%s , %.0f eV range , %.0f%% C.L.',obj.GetPlotTitle('Mode','data'),obj.range,obj.ConfLevel),'FontWeight','normal','FontSize',get(gca,'FontSize'));
            title(obj.GetPlotTitle,'FontWeight','normal','FontSize',get(gca,'FontSize'));
            if contains(ConvertFixPar('freePar',obj.RunAnaObj.fixPar,'Mode','Reverse'),'mNu')
                ylim([6 40^2]);
                xlim([8e-03 0.5]);
                
            end
           %% save
           if ~strcmp(SavePlot,'OFF')
               name_i = strrep(obj.DefPlotName,sprintf('_%s',chi2_i),'');
               
               if strcmp(obj.NullHypothesis,'ON')
                   extraStr = [extraStr,'_NH'];
               end
               
               if strcmp(SavePlot,'ON')
                   plotname = sprintf('%s_CompareFitter_%.2gCL%s.pdf',name_i,obj.ConfLevel,extraStr);
                   export_fig(gcf,plotname);
               elseif strcmp(SavePlot,'png')
                   plotname = sprintf('%s_CompareFitter_%.2gCL%s.png',name_i,obj.ConfLevel,extraStr);
                   print(gcf,plotname,'-dpng','-r450');
               end
               fprintf('save plot to %s \n',plotname);
           end
        end   
        function PlotStatandSys(obj,varargin)
            % plot for a given range: stat. only and stat + syst
            p = inputParser;
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','png'}));
            
            p.parse(varargin{:});
            SavePlot = p.Results.SavePlot;
            
            chi2_i   = obj.RunAnaObj.chi2;
        
            if strcmp(obj.RunAnaObj.DataType,'Real')
                BestFit = 'ON';
            else
                BestFit = 'OFF';
            end
        
            LineWidth = 2.5;
            %% load stat and syst
            obj.RunAnaObj.chi2 = 'chi2Stat';
            obj.SetNPfactor;
            obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
            obj.Interp1Grid('RecomputeFlag','ON');
            pStat = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','OFF',...
                'Color',rgb('DodgerBlue'),'LineStyle','-','BestFit',BestFit,'PlotSplines','OFF');
            
            obj.RunAnaObj.chi2 = 'chi2CMShape';
            obj.SetNPfactor;
            obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
            obj.Interp1Grid('RecomputeFlag','ON');
            [pSys,pSysMin] = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                'Color',rgb('Orange'),'LineStyle','-','BestFit',BestFit,'PlotSplines','OFF');

        %% legend
        
           
               legStr = {'Stat. only','All syst. combined'};
               legend([pStat,pSys],legStr,'EdgeColor',rgb('Silver'),'Location','southwest');
               extraStr = '';
           
           obj.RunAnaObj.chi2 = chi2_i;
           
            if strcmp(obj.RunAnaObj.DataSet,'Knm1') && obj.range==65
                xlim([4e-03 0.5])
                ylim([1 1e4])
            elseif strcmp(obj.RunAnaObj.DataSet,'Knm1') && obj.range==40
                xlim([1e-02 0.5])
                ylim([1 3e3])
            elseif strcmp(obj.RunAnaObj.DataSet,'Knm1') && obj.range==95
                xlim([3e-03 0.5])
                ylim([1 2e4]) 
            else
                xlim([pSysMin*0.5,0.5]);
            end
            
            title(sprintf('%s , %.0f eV range , %.0f%% C.L.',obj.GetPlotTitle('Mode','data'),obj.range,obj.ConfLevel),'FontWeight','normal','FontSize',get(gca,'FontSize'));
          
           %% save
           if ~strcmp(SavePlot,'OFF')
               name_i = strrep(obj.DefPlotName,sprintf('_%s',chi2_i),'');
               if strcmp(SavePlot,'ON')
                   plotname = sprintf('%s_StatandSyst_%.2gCL%s.pdf',name_i,obj.ConfLevel,extraStr);
                   export_fig(gcf,plotname);
               elseif strcmp(SavePlot,'png')
                   plotname = sprintf('%s_StatandSyst_%.2gCL%s.png',name_i,obj.ConfLevel,extraStr);
                   print(gcf,plotname,'-dpng','-r450');
               end
               fprintf('save plot to %s \n',plotname);
           end
        end
        function PlotPRL1(obj,varargin)
            % prl plot 1: comparison with mainz & troitsk
            p = inputParser;
            p.addParameter('BestFit','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','png'}));
            p.addParameter('Mainz','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Troitsk','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Style','Reg',@(x)ismember(x,{'Reg','PRL'}));
            p.addParameter('FinalSensitivity','OFF',@(x)ismember(x,{'OFF','ON'}))
            p.addParameter('Sensitivity','OFF',@(x)ismember(x,{'OFF','ON'}))
            p.addParameter('FreemNuSq','OFF',@(x)ismember(x,{'OFF','ON','Pull'})) % without pull
            p.addParameter('AddPull','', @(x)isfloat(x) || isempty(x)); % with additional pull
            
            p.parse(varargin{:});
            BestFit          = p.Results.BestFit;
            SavePlot         = p.Results.SavePlot;
            Troitsk          = p.Results.Troitsk;
            Mainz            = p.Results.Mainz;
            Style            = p.Results.Style;
            FinalSensitivity = p.Results.FinalSensitivity;
            Sensitivity      = p.Results.Sensitivity;
            FreemNuSq        = p.Results.FreemNuSq;
            AddPull          = p.Results.AddPull;
            
            fixPar_i = obj.RunAnaObj.fixPar;
            pull_i = obj.RunAnaObj.pullFlag;
            
            fPRL = figure('Units','normalized','Position',[0.1,0.1,0.4,0.6]);
            legHandle = cell(0,0);
            legStr = '';
            savedirOther = [getenv('SamakPath'),'SterileAnalysis/GridSearchFiles/Knm1/Others/'];
            
            if strcmp(Mainz,'ON')
                filenameMainz = sprintf('%scoord_Mainz_95CL.mat',savedirOther);
                dMainz = importdata(filenameMainz);
                sinTsq = 0.5*(1-sqrt(1-dMainz.SinSquare2Theta_X));
                pMainz = plot(sinTsq,dMainz.DmSquare41_Y,'-.','LineWidth',1.5,'Color',rgb('Salmon')); %Red
                legHandle{numel(legHandle)+1} = pMainz;
                legStr = [legStr,{sprintf('Mainz 95%% C.L. :  {\\itm}_\\nu^2 = 0 eV^2')}];
                hold on;
            end
            %% Troitsk
            if strcmp(Troitsk,'ON')
                filenameTroitsk = sprintf('%scoord_Troitsk_95CL.mat',savedirOther);
                dTroitsk = importdata(filenameTroitsk);
                pTroitsk = plot(dTroitsk.SinSquareTheta_X,dTroitsk.m4Square_Y,'--','LineWidth',1.5,...
                    'Color',rgb('Black')); % Orange
                legHandle{numel(legHandle)+1} = pTroitsk;
                legStr = [legStr,{sprintf('Troitsk 95%% C.L. :  {\\itm}_\\nu^2 = 0 eV^2')}];
                hold on;
            end
            
            %% KSN1 sensitivity
            if strcmp(Sensitivity,'ON') && strcmp(obj.RunAnaObj.DataType,'Real')
                obj.RunAnaObj.DataType = 'Twin';
                obj.RunAnaObj.fixPar = 'E0 Norm Bkg'; obj.RunAnaObj.InitFitPar;
                obj.RunAnaObj.pullFlag = 99;
                obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
                obj.Interp1Grid('RecomputeFlag','ON');
                pFixSensi = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                    'Color',rgb('Silver'),'LineStyle','-','BestFit','OFF');
                pFixSensi.LineWidth = 2;
                legHandle{numel(legHandle)+1} = pFixSensi;
                legStr = [legStr,{sprintf('KATRIN sensitivity %.0f%% C.L. :  {\\itm}_\\nu^2 = 0 eV^2',obj.ConfLevel)}];
                
                obj.RunAnaObj.DataType = 'Real';
            end
            
            %% fixed nu-mass
            obj.RunAnaObj.fixPar = 'E0 Norm Bkg'; obj.RunAnaObj.InitFitPar;
            obj.RunAnaObj.pullFlag = 99;
            obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
            obj.Interp1Grid('RecomputeFlag','ON');
            pFix = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                'Color',rgb('DodgerBlue'),'LineStyle','-','BestFit',BestFit);
            legHandle{numel(legHandle)+1} = pFix;
            legStr = [legStr,{sprintf('KATRIN %.0f%% C.L. :  {\\itm}_\\nu^2 = 0 eV^2',obj.ConfLevel)}];
            
            
            %% free unconstrained nu mass
            if strcmp(FreemNuSq,'ON')
                if strcmp(obj.RunAnaObj.DataSet,'Knm1')
                obj.InterpMode = 'lin';
                end
                obj.RunAnaObj.fixPar = 'mNu E0 Norm Bkg'; obj.RunAnaObj.InitFitPar;
                obj.RunAnaObj.pullFlag = 99;
                obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
                obj.Interp1Grid('RecomputeFlag','ON');
                pfree = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                    'Color',rgb('Navy'),'LineStyle',':','BestFit',BestFit);
                
                legHandle{numel(legHandle)+1} = pfree;
                
                legStr = [legStr,{...
                    sprintf('KATRIN %.0f%% C.L. :  {\\itm}_\\nu^2 free',obj.ConfLevel)}];
                
            end
            
            if strcmp(FreemNuSq,'Pull')
                %%  nuissance nu-mass + pull
                obj.RunAnaObj.fixPar = 'mNu E0 Norm Bkg'; obj.RunAnaObj.InitFitPar;
                obj.RunAnaObj.pullFlag = 12;
                obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
                obj.Interp1Grid('RecomputeFlag','ON');
                pPull = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                    'Color',rgb('Navy'),'LineStyle',':','BestFit',BestFit);
                
                legHandle{numel(legHandle)+1} = pPull;
                
                legStr = [legStr,{sprintf('KATRIN %.0f%% C.L. :  {\\itm}_\\nu^2 free , \\sigma({\\itm}_\\nu^2) = 1.94 eV^2',obj.ConfLevel)}];
            end
            
            %% additional contour with some pull
            if ~isempty(AddPull)
                %%  nuissance nu-mass + pull
                obj.RunAnaObj.fixPar = 'mNu E0 Norm Bkg'; obj.RunAnaObj.InitFitPar;
                obj.RunAnaObj.pullFlag = AddPull;
                obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
                obj.Interp1Grid('RecomputeFlag','ON');
                pAddPull = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                    'Color',rgb('CadetBlue'),'LineStyle','-.','BestFit',BestFit);
                
                legHandle{numel(legHandle)+1} = pAddPull;
                
                switch AddPull
                    case 12
                        pullStr = sprintf('\\sigma({\\itm}_\\nu^{2 }) = 1.94 eV^2');
                    case 15
                        pullStr = sprintf('\\sigma({\\itm}_\\nu^{2 }) = 1 eV^2');
                    case 16
                        pullStr = sprintf('\\sigma({\\itm}_\\nu^{2 }) = 2 eV^2');
                    case 17
                        pullStr = sprintf('\\sigma({\\itm}_\\nu^{2 }) = 3 eV^2');
                    case 18
                        pullStr = sprintf('\\sigma({\\itm}_\\nu^{2 }) = 0.5 eV^2');
                    case 26
                        pullStr = sprintf('\\sigma({\\itm}_\\nu^{2 }) = 1.1 eV^2');
                    otherwise
                        pullStr = '';
                end
                legStr = [legStr,{sprintf('KATRIN %.0f%% C.L. :  {\\itm}_\\nu^2 free , %s',obj.ConfLevel,pullStr)}];
            end
            
            %% final sensitivity
            if strcmp(FinalSensitivity,'ON')
                DataType_i = obj.RunAnaObj.DataType;
                AngTF = obj.RunAnaObj.AngularTFFlag;
                ElossFlag = obj.RunAnaObj.ELossFlag;
                Budget = obj.RunAnaObj.SysBudget;
                FakeInitFile = @ref_KSNX_KATRIN_Final;
                obj.RunAnaObj.DataType = 'Fake';
                obj.RunAnaObj.DataSet = 'Knm2';
                
                obj.RunAnaObj.FakeInitFile = FakeInitFile;
                obj.RunAnaObj.fixPar = 'E0 Norm Bkg';
                obj.RunAnaObj.InitFitPar;
                obj.RunAnaObj.pullFlag = 99;
                obj.RunAnaObj.AngularTFFlag = 'ON';
                obj.RunAnaObj.ELossFlag = 'KatrinT2A20';
                obj.RunAnaObj.SysBudget = 66;
                obj.LoadGridFile('CheckSmallerN','ON',obj.LoadGridArg{:});
                obj.Interp1Grid('RecomputeFlag','ON');
                pFull = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                    'Color',rgb('LightGreen'),'LineStyle','-','BestFit',BestFit);
                
                legHandle{numel(legHandle)+1} = pFull;
                obj.RunAnaObj.DataType= DataType_i ;
                obj.RunAnaObj.DataSet = 'Knm1';
                obj.RunAnaObj.AngularTFFlag = AngTF;
                obj.RunAnaObj.ELossFlag = ElossFlag;
                obj.RunAnaObj.SysBudget = Budget;
                legStr = [legStr,sprintf('KATRIN Final %.0f%% C.L. :  {\\itm}_\\nu^2 = 0 eV^2)',obj.ConfLevel)];
            end
            
            
            
            
            %% appearance + legend
            leg =  legend([legHandle{:}],legStr{:},...
                'EdgeColor',rgb('Silver'),'Location','southwest');
            set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
            
            if strcmp(Style,'Reg')
                PrettyFigureFormat('FontSize',20);
                leg.FontSize = get(gca,'FontSize')-2;
            else
                PRLFormat
                set(gca,'FontSize',20);
                leg.FontSize = 17;
            end
            
            legend boxoff
            if strcmp(obj.RunAnaObj.DataSet,'Knm1')
                if obj.range==65
                    ylim([1 6e3]);
                    xlim([2e-03 0.5]);
                elseif obj.range==95
                    ylim([1 1e4]);
                    xlim([9e-04 0.5]);
                elseif obj.range==40
                    %ylim([1 2e3]);
                    ylim([0.2 2e3]);
                    xlim([6e-03 0.5]);
                end
            else
                ylim([1 1600])
                xlim([4e-03 0.5])
            end
            title('');%remove title
            %grid on
            
            if strcmp(FinalSensitivity,'ON')
                xlim([1e-03 0.5]);
            end
            %% save
            if ~strcmp(SavePlot,'OFF')
                name_i = extractBefore(strrep(obj.DefPlotName,'_mNuE0BkgNorm',''),'_pull');
                if strcmp(FinalSensitivity,'ON')
                    name_i = [name_i,'_FinalSensi'];
                end
                if strcmp(SavePlot,'ON')
                    plotname = sprintf('%s_PRL1_%.2gCL.pdf',name_i,obj.ConfLevel);
                    export_fig(gcf,plotname);
                elseif strcmp(SavePlot,'png')
                    plotname = sprintf('%s_PRL1_%.2gCL.png',name_i,obj.ConfLevel);
                    print(gcf,plotname,'-dpng','-r450');
                end
                fprintf('save plot to %s \n',plotname);
            end
            
            obj.RunAnaObj.fixPar = fixPar_i;
            obj.RunAnaObj.pullFlag = pull_i ;
        end
        function PlotStatOverSys(obj,varargin)
            p = inputParser;
            p.addParameter('Ranges',[95:-5:65],@(x)isfloat(x));
            p.addParameter('SavePlot','ON',@(x)ismember(x,{'ON','OFF','png'}));
            p.parse(varargin{:});
            Ranges   = p.Results.Ranges;
            SavePlot = p.Results.SavePlot;
            PlotStatDom = 'OFF';
            [sin2T4_Stat, sin2T4_Sys, sin2T4_Tot, mmNu4SqCommon, StatDomFraction] = obj.StatOverSys('Ranges',Ranges);
            pl = cell(numel(Ranges),1);
            GetFigure;
            pref = plot(logspace(0,4,10),ones(10,1),'k-','LineWidth',2);
            hold on;
            for i=1:numel(Ranges)
                %     pl{i}= plot(mnu4SqCommonPlot{i},sin2T4StatOnlyPlot{i}.^2,...
                %         LineStyles{i},'LineWidth',2.5,'Color',rgb(Colors{i}));
                %     hold on;
                %     pl{i}= plot(mnu4SqCommonPlot{i},sin2T4SystOnlyPlot{i}.^2,...
                %         ':','LineWidth',2.5,'Color',rgb(Colors{i}));
                pl{i}= plot(mmNu4SqCommon{i},sin2T4_Sys{i}.^2./sin2T4_Stat{i}.^2,...
                    obj.PlotLines{i},'LineWidth',2.5,'Color',obj.PlotColors{i});
                %     set(gca,'YScale','lin');
                set(gca,'XScale','log');
                ylabel(sprintf('\\sigma^2_{syst.}/\\sigma^2_{stat.}(|U_{e4}|^2)'));
                xlabel(sprintf('{\\itm}_4^2 (eV^2)'));
                
                if strcmp(PlotStatDom,'ON')
                    legStr{i+1} = sprintf('%.0f eV range (%.0f%% stat. dominated)',...
                        Ranges(i),100*StatDomFraction(i));
                else
                    legStr{i+1} = sprintf('%.0f eV range',Ranges(i));
                end
            end
            
            if strcmp(obj.RunAnaObj.DataSet,'Knm1')
                ylim([0 6.5])
            end
            PrettyFigureFormat;
            %  leg = legend([pref,pl{:}],legStr{:},'EdgeColor',rgb('Silver'),'Location','northwest');
            
            
        end
        function PlotQuadrant(obj,varargin)
            % plot 4 quadrant grid  plot of physical and nonphysical parameter space
             p = inputParser;
            p.addParameter('SavePlot','ON',@(x)ismember(x,{'ON','OFF','png'}));
            p.parse(varargin{:});
            SavePlot = p.Results.SavePlot;
            PltxTics = [1e-03,1e-02,1e-01,0.5,1];
            PltyTics = [0.1,1,10,1e2,1e3];
            
            CL = 95;
            %obj.GridPlot;
             if CL<1
                CL = CL*100;
            end
            obj.DeltaChi2 = GetDeltaChi2(CL,2);
            zlimMax = obj.DeltaChi2;
                
            mNuSq_bf  = zeros(4,1);
            sin2T4_bf = zeros(4,1);
            chi2_bf   = zeros(4,1);
            pbf = cell(4,1);
            
            GetFigure; 
            % 1. NW
            obj.LoadGridFile(obj.LoadGridArg{:},'Negsin2T4','ON');
            obj.Interp1Grid;
            chi2grid1 = obj.chi2;
            chi2grid1((chi2grid1-obj.chi2_ref)>obj.DeltaChi2) =  NaN;
            
            s1 = subplot(2,2,1); 
            surf(obj.sin2T4,obj.mNu4Sq,chi2grid1-obj.chi2_ref,'EdgeColor','interp','FaceColor','interp');
            PrettyFigureFormat; 
             zlim([0 zlimMax])
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            xlim([-0.5 -1e-03]);
            ylim([0.1,40^2]);
            view([0 0 1])
            grid off
            ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
            ax1 = gca;
            xticks(sort(-PltxTics));
            yticks(PltyTics);
            yticklabels({'10^{-1}','10^0','10^1','10^2','10^3'});
            xticklabels({''});
            
            obj.LoadGridFile(obj.LoadGridArg{:},'Negsin2T4','ON');
            obj.InterpMode = 'spline';
            obj.Interp1Grid('Maxm4Sq',34^2); 
            obj.FindBestFit; 
            obj.InterpMode = 'lin';
            mNuSq_bf(1)  = obj.mNu4Sq_bf;
            sin2T4_bf(1) = obj.sin2T4_bf;
            chi2_bf(1)    = obj.chi2_bf;
            if chi2_bf(1)<obj.chi2_Null
                hold on;
                pbf{1} = plot3(sin2T4_bf(1),mNuSq_bf(1),2,'x','Color',rgb('White'),'LineWidth',2,'MarkerSize',8);
                leg = legend(pbf{1},sprintf('Best fit \\chi^2 = %.1f',chi2_bf(1)),'Location','southeast');
                PrettyLegendFormat(leg);   
            end
            % 2. NE
            obj.LoadGridFile(obj.LoadGridArg{:},'Extsin2T4','ON');
            obj.Interp1Grid;
            chi2grid1 = obj.chi2;
            chi2grid1((chi2grid1-obj.chi2_ref)>obj.DeltaChi2) =  NaN;
                
            s2 = subplot(2,2,2); 
            surf(obj.sin2T4,obj.mNu4Sq,chi2grid1-obj.chi2_ref,'EdgeColor','interp','FaceColor','interp');
             PrettyFigureFormat; 
             xlim([1e-03,1]);
             ylim([0.1,40^2]);
            
            zlim([0 zlimMax])
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            grid off
            view([0 0 1])
            ax2 = gca;
            xticks(PltxTics);
            yticks(PltyTics);
            yticklabels({'',''});
            xticklabels({''});
            
            obj.LoadGridFile(obj.LoadGridArg{:},'Extsin2T4','ON');
            obj.InterpMode = 'spline';
            obj.Interp1Grid('Maxm4Sq',34^2); 
            obj.FindBestFit; 
            obj.InterpMode = 'lin';
            
            mNuSq_bf(2)  = obj.mNu4Sq_bf;
            sin2T4_bf(2) = obj.sin2T4_bf;
            chi2_bf(2)    = obj.chi2_bf;
            if chi2_bf(2)<obj.chi2_Null
                hold on;
                pbf{2} = plot3(sin2T4_bf(2),mNuSq_bf(2),2,'x','Color',rgb('White'),'LineWidth',2,'MarkerSize',8);
                leg = legend(pbf{2},sprintf('Best fit \\chi^2 = %.1f',chi2_bf(2)),'Location','southwest');
                PrettyLegendFormat(leg);
            end
            
            % 3. SW
            obj.LoadGridFile(obj.LoadGridArg{:},'NegmNu4Sq','ON','Negsin2T4','ON');
            obj.Interp1Grid;
            chi2grid1 = obj.chi2;
            chi2grid1((chi2grid1-obj.chi2_ref)>obj.DeltaChi2) =  NaN;
 
            s3 = subplot(2,2,3); 
            surf(obj.sin2T4,obj.mNu4Sq,chi2grid1-obj.chi2_ref,'EdgeColor','interp','FaceColor','interp');
             PrettyFigureFormat; 
         
            zlim([0 zlimMax])
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            grid off
            view([0 0 1])
            ax3 = gca;
            xlim([-0.5,-1e-03]);
            ylim([-40^2,-0.1]);
             xticks(sort(-PltxTics));
             yticks(sort(-PltyTics));
             yticklabels({'-10^3','-10^2','-10^1','-10^0',''})
             xticklabels({'','','-10^{-1}','-10^{-2}',''});
            
             obj.LoadGridFile(obj.LoadGridArg{:},'NegmNu4Sq','ON','Negsin2T4','ON');
            obj.InterpMode = 'spline';
            obj.Interp1Grid('Maxm4Sq',34^2); 
            obj.FindBestFit; 
            obj.InterpMode = 'lin';
            
            mNuSq_bf(3)  = obj.mNu4Sq_bf;
            sin2T4_bf(3) = obj.sin2T4_bf;
            chi2_bf(3)    = obj.chi2_bf;
             if chi2_bf(3)<obj.chi2_Null
            hold on;
            pbf{3} = plot3(sin2T4_bf(3),mNuSq_bf(3),2,'x','Color',rgb('White'),'LineWidth',2,'MarkerSize',8);
             leg = legend(pbf{3},sprintf('Best fit \\chi^2 = %.1f',chi2_bf(3)),'Location','northeast');
                PrettyLegendFormat(leg); 
             end
             
            % 4.SE
            obj.LoadGridFile(obj.LoadGridArg{:},'NegmNu4Sq','ON','Extsin2T4','ON');
            obj.Interp1Grid;
            chi2grid1 = obj.chi2;
            chi2grid1((chi2grid1-obj.chi2_ref)>obj.DeltaChi2) =  NaN;
            s3 = subplot(2,2,4);
            surf(obj.sin2T4,obj.mNu4Sq,chi2grid1-obj.chi2_ref,'EdgeColor','interp','FaceColor','interp');
            PrettyFigureFormat;
            xlim([1e-03,1]);
            ylim([-40^2,-0.1]);        
            zlim([0 zlimMax])
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            grid off
            view([0 0 1])
            ax4 = gca;
            xlabel(sprintf('|{\\itU}_{e4}|^2'));
            xticks(PltxTics);
            yticks(sort(-PltyTics));
            yticklabels({''});
            xticklabels({'10^{-3}','10^{-2}','10^{-1}','','10^0'});
            
            obj.LoadGridFile(obj.LoadGridArg{:},'NegmNu4Sq','ON','Extsin2T4','ON');
            obj.InterpMode = 'spline';
            obj.Interp1Grid('Maxm4Sq',34^2); 
            obj.FindBestFit; 
            obj.InterpMode = 'lin'; 
            mNuSq_bf(4)  = obj.mNu4Sq_bf;
            sin2T4_bf(4) = obj.sin2T4_bf;
            chi2_bf(4)    = obj.chi2_bf;
             if chi2_bf(4)<obj.chi2_Null
            hold on;
            pbf{4} = plot3(sin2T4_bf(4),mNuSq_bf(4),2,'x','Color',rgb('White'),'LineWidth',2,'MarkerSize',8);
             leg = legend(pbf{1},sprintf('Best fit \\chi^2 = %.1f',chi2_bf(4)),'Location','southwest');
                PrettyLegendFormat(leg); 
             end
            %% find global minimum
            MinIdx = find(chi2_bf==min(chi2_bf));
            pbf{MinIdx}.Color = rgb('IndianRed');
            %% move closer together
            ax1.Position(2) = 0.61;
            ax3.Position(2) = 0.245;%0.245;
            ax4.Position(2) = 0.245;
            ax2.Position(2) = 0.61;
            ax4.Position(1) = 0.49;
            ax2.Position(1) = 0.49;
            ax4.Position(3) = ax1.Position(3);
            
            ax1.YLabel.Position(2) = 0.05;
          %  ax1.YLabel.Position(1) = 1.5;
            ax1.YLabel.FontSize = 21;
            ax4.XLabel.Position(1) = 0.9e-03;
            ax4.XLabel.FontSize = 21;

            c =colorbar;
            c.Label.String = sprintf('\\Delta\\chi^2');
            c.Label.FontSize = 21;
            c.Limits=[0 zlimMax]; 
            c.Position(1) = 0.86;
            c.Position(2) = 0.24;
            c.Position(3) = 0.023;
            c.Position(4) = 0.71;
        
              if ~strcmp(SavePlot,'OFF')
                 if strcmp(SavePlot,'ON')
                     plotname = sprintf('%s_Contour_%.2gCL_Quadrant.pdf',obj.DefPlotName,obj.ConfLevel);
                     export_fig(gcf,plotname);
                 elseif strcmp(SavePlot,'png')
                     plotname = sprintf('%s_Contour_%.2gCL_Quadrant.png',obj.DefPlotName,obj.ConfLevel);
                     print(gcf,plotname,'-dpng','-r450');
                 end
                 fprintf('save plot to %s \n',plotname);
              end
             
        end
        function PlotTririumSpectrumImprint(obj,varargin)
            % look at imprint of (mnu4Sq_Grid,sin2T4) in tritium spectrum
            % define chi2 = (H0-H1)^2/sigma(H0)^2
           p = inputParser;
           p.addParameter('mnu4Sq','',@(x)isfloat(x) || isempty(x));
           p.addParameter('sin2T4',0.01,@(x)isfloat(x));
           p.addParameter('Mode','Residuals',@(x)ismember(x,{'Residuals','Ratio'}));
           p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF'}));
           p.parse(varargin{:});
           mNu4Sq_Plot = p.Results.mnu4Sq;
           sin2T4_Plot = p.Results.sin2T4;
           Mode        = p.Results.Mode;
           SavePlot    = p.Results.SavePlot;
           
           qU = obj.RunAnaObj.ModelObj.qU;
           Time = obj.RunAnaObj.ModelObj.qUfrac.*obj.RunAnaObj.ModelObj.TimeSec;
           
           % null hypothesis (no sterile)
           obj.RunAnaObj.ModelObj.SetFitBiasSterile(0,0);
           obj.RunAnaObj.ModelObj.ComputeTBDDS;
           obj.RunAnaObj.ModelObj.ComputeTBDIS;
           H0 = obj.RunAnaObj.ModelObj.TBDIS./Time;
           H0Err = sqrt(obj.RunAnaObj.ModelObj.TBDIS)./Time;
           
           qU_plt = linspace(min(qU),max(qU),1e3)';
           H0 = interp1(qU,H0,qU_plt,'spline');
           H0Err = interp1(qU,H0Err,qU_plt,'spline');
           
          % alt. hypothesis (3+1 sterile)
          if isempty(mNu4Sq_Plot)
              mNu4Sq_Plot = [0,5,10,20,30,40].^2;
              nmNuSq = numel(mNu4Sq_Plot);
              H1 = zeros(numel(qU_plt),nmNuSq);
              H1Err = zeros(numel(qU_plt),nmNuSq);
              for i=1:nmNuSq
                 progressbar(i/nmNuSq)
                  obj.RunAnaObj.ModelObj.SetFitBiasSterile(mNu4Sq_Plot(i),sin2T4_Plot);
                  obj.RunAnaObj.ModelObj.ComputeTBDDS;
                  obj.RunAnaObj.ModelObj.ComputeTBDIS;
                  H1tmp = obj.RunAnaObj.ModelObj.TBDIS./Time;
                  H1Errtmp = sqrt(obj.RunAnaObj.ModelObj.TBDIS)./Time;    
                  H1(:,i) = interp1(qU,H1tmp,qU_plt,'spline');
                  H1Err(:,i) = interp1(qU,H1Errtmp,qU_plt,'spline');
              end
              
          else
              nmNuSq = 1;
              obj.RunAnaObj.ModelObj.SetFitBiasSterile(mNu4Sq_Plot,sin2T4_Plot);
              obj.RunAnaObj.ModelObj.ComputeTBDDS;
              obj.RunAnaObj.ModelObj.ComputeTBDIS;
              H1 = obj.RunAnaObj.ModelObj.TBDIS./Time;
              H1Err = sqrt(obj.RunAnaObj.ModelObj.TBDIS)./Time;
              H1 = interp1(qU,H1,qU_plt,'spline');
              H1Err = interp1(qU,H1Err,qU_plt,'spline');
          end
          
        
          if strcmp(Mode,'Residuals') 
              y = (H1-H0)./H1Err;
              yStr = sprintf('Residuals (\\sigma)');
          else
              y = H1./H0;
              yStr = sprintf('Ratio H_1/H_0');
          end
          GetFigure;
          a = area(qU(1:obj.RunAnaObj.exclDataStart)-obj.RunAnaObj.ModelObj.Q_i-0.5,...
              50.*ones(numel(qU(1:obj.RunAnaObj.exclDataStart)),1),-50,...
              'FaceColor',rgb('LightGray'),'EdgeColor','none');
          hold on;
          Colors = flipud(colormap('jet'));
         %  Colors = flipud(colormap('hsv'));
         p = cell(nmNuSq,1);
          for i=1:nmNuSq
          p{i} = plot(qU_plt-obj.RunAnaObj.ModelObj.Q_i,y(:,i),'-.','LineWidth',2,...
              'Color',Colors(i*floor(256/nmNuSq),:),'LineStyle',obj.PlotLines{i},...
              'MarkerSize',15);
          hold on;
          end
          PrettyFigureFormat('FontSize',20);
          ylabel(yStr);
          xlabel(sprintf('Retarding energy - {\\itE}_0 (eV)'));
          
          legStr = arrayfun(@(x) sprintf('{\\itm}_4 = %.0f eV',x),sqrt(mNu4Sq_Plot),'UniformOutput',false);
          leg = legend([p{:}],legStr);
          leg.Title.String = sprintf('|{\\itU}_{e4}|^2 = %.3g',sin2T4_Plot);
          leg.Title.FontWeight = 'normal';
          PrettyLegendFormat(leg);
          leg.Location = 'southeast';
          xlim([-54,5])
          if min(y)<0
              ylim([min(min(y)).*1.01 1.01.*max(max(y))]);
          else
              ylim([min(min(y)).*0.996 1.004.*max(max(y))]);
          end
          text(-46.8,min(min(ylim))+0.9*(max(max(ylim))-min(min(ylim))),...
              sprintf('Outside of \nROI'),'FontSize',get(gca,'FontSize'),'HorizontalAlignment','center')
          if strcmp(SavePlot,'ON')
              plotdir = [getenv('SamakPath'),sprintf('SterileAnalysis/plots/%s/spectra/',obj.RunAnaObj.DataSet)];
              MakeDir(plotdir)
              plotname = sprintf('%sTririumSpectrumImprint_%s_%s.png',plotdir,obj.RunAnaObj.DataSet,Mode);
              print(plotname,'-dpng','-r300');
              fprintf('save plot to %s \n',plotname);
          end 
        end
        function PlotTririumSpectrum(obj,varargin)
            % look at imprint of (mnu4Sq_Grid,sin2T4) in tritium spectrum
            % define chi2 = (H0-H1)^2/sigma(H0)^2
           p = inputParser;
           p.addParameter('mnu4Sq',20^2,@(x)isfloat(x) || isempty(x));
           p.addParameter('sin2T4',0.1,@(x)isfloat(x));
           p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF'}));
           p.parse(varargin{:});
           mNu4Sq_Plot = p.Results.mnu4Sq;
           sin2T4_Plot = p.Results.sin2T4;
           SavePlot    = p.Results.SavePlot;
           
           qU = obj.RunAnaObj.ModelObj.qU;
           Time = obj.RunAnaObj.ModelObj.qUfrac.*obj.RunAnaObj.ModelObj.TimeSec;
           
           % null hypothesis (no sterile)
           obj.RunAnaObj.ModelObj.SetFitBiasSterile(0,0);
           obj.RunAnaObj.ModelObj.ComputeTBDDS;
           obj.RunAnaObj.ModelObj.ComputeTBDIS;
           H0 = obj.RunAnaObj.ModelObj.TBDIS./Time;

           % only. hypothesis (3+1 sterile)
           obj.RunAnaObj.ModelObj.SetFitBiasSterile(mNu4Sq_Plot,1);%sin2T4_Plot);
           obj.RunAnaObj.ModelObj.ComputeTBDDS;
           obj.RunAnaObj.ModelObj.ComputeTBDIS;
           H1 = obj.RunAnaObj.ModelObj.TBDIS./Time;
   
           BKG_RateSec = obj.RunAnaObj.ModelObj.BKG_RateSec+0.5.*obj.RunAnaObj.ModelObj.BKG_PtSlope.*Time./obj.RunAnaObj.nRuns;
           GetFigure;
           p0 = plot(qU-obj.RunAnaObj.ModelObj.Q,(1-sin2T4_Plot).*(H0-BKG_RateSec)+BKG_RateSec,'LineWidth',2.5,'Color',rgb('DodgerBlue'));
           hold on;
           p1 = plot(qU-obj.RunAnaObj.ModelObj.Q,sin2T4_Plot.*(H1-BKG_RateSec)+BKG_RateSec,'-.','LineWidth',2.5,'Color',rgb('ForestGreen'));
           p2 = plot(qU-obj.RunAnaObj.ModelObj.Q,sin2T4_Plot.*(H1-BKG_RateSec)+BKG_RateSec+...
               (1-sin2T4_Plot).*(H0-BKG_RateSec),':','LineWidth',2.5,'Color',rgb('Orange'));
           set(gca,'YScale','log')
           PrettyFigureFormat('FontSize',22);
           ylabel('Rate');
           xlabel(sprintf('Retarding energy - {\\itE}_0 (eV)'));
           leg = legend([p0,p1,p2],sprintf('Active branch'),...
               sprintf('Sterile branch: {\\itm}_4 = %.3g eV , |{\\itU}_{e4}|^2 = %.3g',sqrt(mNu4Sq_Plot),sin2T4_Plot),...
               'Active + Sterile branch');
           PrettyLegendFormat(leg);
          
           ylim([0.1,1e2])
           xlim([-42,5])
           title('Integral spectrum','FontWeight','normal','FontSize',get(gca,'FontSize'));
           
           if strcmp(SavePlot,'ON')
              plotdir = [getenv('SamakPath'),sprintf('SterileAnalysis/plots/%s/spectra/',obj.RunAnaObj.DataSet)];
              MakeDir(plotdir)
              plotname = sprintf('%sTririumSpectrum_%s_%.3geVm4_%.3gsint4Sq.png',...
                  plotdir,obj.RunAnaObj.DataSet,sqrt(mNu4Sq_Plot),sin2T4_Plot);
              print(plotname,'-dpng','-r300');
              fprintf('save plot to %s \n',plotname);
           end 
          
        end
       
    end
    
    % Data stream: labels, loading, saving
    methods
        function f = LoadGridFile(obj,varargin)
            p = inputParser;
            p.addParameter('CheckLargerN','ON',@(x)ismember(x,{'ON','OFF'})); % if specified ngrids do not exist - look also for larger n grids
            p.addParameter('CheckSmallerN','OFF',@(x)ismember(x,{'ON','OFF'})); % if specified ngrids  +larger do not exist - look also for smaller n grids
            p.addParameter('CheckExtmNu4Sq','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Negsin2T4','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('NegmNu4Sq','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('IgnoreKnm2FSDbinning','OFF',@(x)ismember(x,{'ON','OFF'})); % if OFF -> load file with other binning if available
            p.addParameter('Extsin2T4','OFF',@(x)ismember(x,{'ON','OFF'})); %extended sin2T2 (up to 1)
            p.addParameter('ExtmNu4Sq','OFF',@(x)ismember(x,{'ON','OFF'})); %extended m4Sq (from 0.1)
            p.addParameter('mNu4SqTestGrid','OFF',@(x)strcmp(x,'OFF') || isfloat(x));
            p.addParameter('FixmNuSq',0,@(x)isfloat(x)); % if light nu-mass fixed (eV^2)
            p.parse(varargin{:});
            CheckLargerN  = p.Results.CheckLargerN;
            CheckSmallerN = p.Results.CheckSmallerN;
            CheckExtmNu4Sq= p.Results.CheckExtmNu4Sq;
            Negsin2T4     = p.Results.Negsin2T4;
            NegmNu4Sq     = p.Results.NegmNu4Sq;
            Extsin2T4     = p.Results.Extsin2T4;
            ExtmNu4Sq     = p.Results.ExtmNu4Sq;
            mNu4SqTestGrid = p.Results.mNu4SqTestGrid;
            FixmNuSq      = p.Results.FixmNuSq;
            IgnoreKnm2FSDbinning = p.Results.IgnoreKnm2FSDbinning;
            
            filename = obj.GridFilename('Negsin2T4',Negsin2T4,'NegmNu4Sq',NegmNu4Sq,...
                                        'Extsin2T4',Extsin2T4,'ExtmNu4Sq',ExtmNu4Sq,...
                                        'FixmNuSq',FixmNuSq,'mNu4SqTestGrid',mNu4SqTestGrid);
            
            loadSuccess = 0;
            
            if exist(filename,'file')
                f = importdata(filename);
                fprintf('load grid from file %s \n',filename)
                loadSuccess = 1;
            end
            
            if strcmp(CheckExtmNu4Sq,'ON') && loadSuccess == 0
                if strcmp(ExtmNu4Sq,'OFF')
                    TestFile =  strrep(filename,'.mat','_ExtmNu4Sq.mat');
                    if exist(TestFile,'file')
                     f = importdata(TestFile);
                    fprintf('change local ExtmNu4Sq to ON - load grid from file %s \n',TestFile)
                    loadSuccess = 1;
                    end
                end
            end
            
            if strcmp(IgnoreKnm2FSDbinning,'ON')
                FSD_i = obj.RunAnaObj.FSDFlag;
                obj.RunAnaObj.FSDFlag = 'KNM2';
                TestFile = obj.GridFilename('Negsin2T4',Negsin2T4,'NegmNu4Sq',NegmNu4Sq,...
                    'Extsin2T4',Extsin2T4,'ExtmNu4Sq',ExtmNu4Sq,...
                    'FixmNuSq',FixmNuSq,'mNu4SqTestGrid',mNu4SqTestGrid);
                TestFile2 = strrep(TestFile,'FSDKNM2','FSDKNM2_0p1eV');
                TestFile3 = strrep(TestFile,'FSDKNM2','FSDKNM2_0p5eV');
                if exist(TestFile,'file')
                    f = importdata(TestFile);
                    fprintf('change local FSD to KNM2 - load grid from file %s \n',TestFile)
                    loadSuccess = 1;
                elseif exist(TestFile2,'file')
                    f = importdata(TestFile2);
                    fprintf('change local FSD to KNM2_0p1eV - load grid from file %s \n',TestFile2)
                    loadSuccess = 1;
                elseif exist(TestFile3,'file')
                    f = importdata(TestFile3);
                    fprintf('change local FSD to KNM2_0p5eV - load grid from file %s \n',TestFile3)
                    loadSuccess = 1;
                end
               obj.RunAnaObj.FSDFlag = FSD_i;
            end
            
            if strcmp(CheckLargerN,'ON') && loadSuccess == 0
                Nmax = 50;
                TestnGrid = (obj.nGridSteps+5):5:Nmax;
                TestFiles = arrayfun(@(x) strrep(filename,sprintf('%.0fnGrid',obj.nGridSteps),...
                    sprintf('%.0fnGrid',x)),TestnGrid,'UniformOutput',0);
                FindFile = find(cellfun(@(x) exist(x,'file'),TestFiles),1,'last'); %largest existing file
                
                if ~isempty(FindFile)
                    f = importdata(TestFiles{FindFile});
                    fprintf('change local grid size to %.0f - load grid from file %s \n',TestnGrid(FindFile),TestFiles{FindFile})
                    loadSuccess = 1;
                end
            end
            
            if strcmp(CheckSmallerN,'ON') && loadSuccess == 0
                Nmin = 10;
                TestnGrid = (obj.nGridSteps-5):-5:Nmin;
                TestFiles = arrayfun(@(x) strrep(filename,sprintf('%.0fnGrid',obj.nGridSteps),...
                    sprintf('%.0fnGrid',x)),TestnGrid,'UniformOutput',0);
                FindFile = find(cellfun(@(x) exist(x,'file'),TestFiles),1,'first'); %largest existing file
                
                if ~isempty(FindFile)
                    f = importdata(TestFiles{FindFile});
                    fprintf('change local grid size to %.0f - load grid from file %s \n',TestnGrid(FindFile),TestFiles{FindFile})
                    loadSuccess = 1;
                end
            end
            
            if loadSuccess == 0
                fprintf(2,'Cannot find grid %s \n',filename);
                f = 0;
            else
                obj.mNu4Sq = f.mnu4Sq;
                obj.sin2T4 = f.sin2T4;
                obj.mNuSq  = cell2mat(cellfun(@(x) x.par(1),f.FitResults,'UniformOutput',0));
                obj.E0     = cell2mat(cellfun(@(x) x.par(2),f.FitResults,'UniformOutput',0))+obj.RunAnaObj.ModelObj.Q_i; 
                obj.chi2   = f.chi2;
                if min(min(obj.chi2)) < f.chi2_ref
                    obj.chi2_ref = min(min(obj.chi2));
                else
                    obj.chi2_ref = f.chi2_ref;
                end
                
                if isfield(f,'FitResults_Null')
                    obj.chi2_Null = f.FitResults_Null.chi2min;
                    obj.dof = f.FitResults_Null.dof-2;
                end
            end
            
        end
        function filename = GridFilename(obj,varargin)
            p = inputParser;
            p.addParameter('Negsin2T4','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('NegmNu4Sq','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('Extsin2T4','OFF',@(x)ismember(x,{'ON','OFF'})); %extended sin2T2 (up to 1)
            p.addParameter('ExtmNu4Sq','OFF',@(x)ismember(x,{'ON','OFF','0.01'}));
            p.addParameter('mNu4SqTestGrid','OFF',@(x)strcmp(x,'OFF') || isfloat(x));
            p.addParameter('FixmNuSq',0,@(x)isfloat(x)); % if light nu-mass fixed (eV^2)
            
            p.parse(varargin{:});
            Negsin2T4     = p.Results.Negsin2T4;
            NegmNu4Sq     = p.Results.NegmNu4Sq;
            Extsin2T4     = p.Results.Extsin2T4;
            ExtmNu4Sq     = p.Results.ExtmNu4Sq;
            mNu4SqTestGrid= p.Results.mNu4SqTestGrid;
            FixmNuSq      = p.Results.FixmNuSq;
            
            %% label
            if strcmp(obj.SmartGrid,'ON')
                AddSin2T4 = 0.1;
                extraStr = sprintf('_SmartGrid%.0e',AddSin2T4);
            else
                extraStr = '';
            end
            if strcmp(obj.RunAnaObj.chi2,'chi2CMShape')
                extraStr = [extraStr,sprintf('_Budget%.0f',obj.RunAnaObj.SysBudget)];
                if ~strcmp(obj.SysEffect,'all')
                    extraStr = [extraStr,sprintf('_%s',obj.SysEffect)];
                end
            end

            savedir = sprintf('%sSterileAnalysis/GridSearchFiles/%s/%s/',...
                     getenv('SamakPath'),obj.RunAnaObj.DataSet,obj.RunAnaObj.DataType);

                 if ~strcmp(obj.RunAnaObj.ELossFlag,'KatrinT2')
                     extraStr = [extraStr,sprintf('_%s',obj.RunAnaObj.ELossFlag)];
                 end
                 
                 if ~strcmp(obj.RunAnaObj.AngularTFFlag,'OFF')
                     extraStr = [extraStr,'_AngTF'];
                 end
                 
                 if obj.Twin_sin2T4~=0
                     extraStr = [extraStr,sprintf('_sinT4Sq%.3g',obj.Twin_sin2T4)];
                 end
                 
                 if obj.Twin_mNu4Sq~=0
                     extraStr = [extraStr,sprintf('_mNu4Sq%.3g',obj.Twin_mNu4Sq)];
                 end
                 
                 if isfloat(obj.RandMC) && strcmp(obj.RunAnaObj.DataType,'Twin') && numel(obj.RandMC_TBDIS)==obj.RunAnaObj.ModelObj.nqU
                     extraStr = sprintf('%s_RandMC%.0f_SumCounts%.0f',extraStr,obj.RandMC,sum(obj.RandMC_TBDIS));
                     savedir = strrep(savedir,'Twin/','TwinExternalMC/');
                 elseif isfloat(obj.RandMC) && strcmp(obj.RunAnaObj.DataType,'Twin')
                     extraStr = sprintf('%s_RandMC%.0f',extraStr,obj.RandMC);
                     savedir = strrep(savedir,'Twin/','TwinRandomizedMC/');
                 end
                 
                 if obj.RunAnaObj.pullFlag<99
                     extraStr = sprintf('%s_pull%.0f',extraStr,obj.RunAnaObj.pullFlag);
                 end
                 
                 if strcmp(NegmNu4Sq,'ON') 
                     extraStr = [extraStr,'_NegmNu4Sq'];
                 end
                 
                 if strcmp(Negsin2T4,'ON')
                     extraStr = [extraStr,'_Negsin2T4'];
                 end
                 
                 if strcmp(Extsin2T4,'ON')
                     extraStr = [extraStr,'_Extsin2T4'];
                 end
                 
                  if strcmp(ExtmNu4Sq,'ON') && (strcmp(obj.RunAnaObj.DataType,'Real') || isfloat(obj.RandMC) )
                      extraStr = [extraStr,'_ExtmNu4Sq'];
                  end
                  
                  if strcmp(ExtmNu4Sq,'0.01') && strcmp(obj.RunAnaObj.DataType,'Real')
                      extraStr = [extraStr,'_ExtmNu4Sq0.01'];
                  end
                  
                  if isfloat(mNu4SqTestGrid)
                      extraStr = [extraStr,sprintf('_mNu4SqTestGrid%.2g',mNu4SqTestGrid)];
                  end
                  
                  if FixmNuSq~=0
                      extraStr = [extraStr,sprintf('_FixmNuSq%.2feV2',FixmNuSq)];
                  end
                  
                  if strcmp(obj.RunAnaObj.DataSet,'Knm1') 
                      if  obj.RunAnaObj.NonPoissonScaleFactor~=1 && obj.RunAnaObj.NonPoissonScaleFactor~=1.064
                      extraStr = [extraStr,sprintf('_NP%.3f',obj.RunAnaObj.NonPoissonScaleFactor)];
                      end
                      if obj.RunAnaObj.ModelObj.BKG_PtSlope ~=0
                          % twin (KNM-2): twin or model BKG_PtSlope not default 3e-06
                          extraStr = [extraStr,sprintf('_BkgPtSlope%.3gmuCpsS',...
                              1e6.*obj.RunAnaObj.ModelObj.BKG_PtSlope)];
                      end
                  elseif strcmp(obj.RunAnaObj.DataSet,'Knm2')
                      if strcmp(obj.RunAnaObj.chi2,'chi2Stat') && obj.RunAnaObj.NonPoissonScaleFactor~=1
                          extraStr = [extraStr,sprintf('_NP%.3g',obj.RunAnaObj.NonPoissonScaleFactor)];
                      elseif strcmp(obj.RunAnaObj.chi2,'chi2CMShape') && obj.RunAnaObj.NonPoissonScaleFactor~=1.112
                          extraStr = [extraStr,sprintf('_NP%.3g',obj.RunAnaObj.NonPoissonScaleFactor)];
                      end
                      % if knm2
                      if strcmp(obj.RunAnaObj.DataType,'Twin')
                          % twin
                          if obj.RunAnaObj.ModelObj.BKG_PtSlope==obj.RunAnaObj.TwinBias_BKG_PtSlope && obj.RunAnaObj.ModelObj.BKG_PtSlope ~=3*1e-06
                              % twin (KNM-2): model and twin BKG_PtSlope same, but not default 3e-06
                              extraStr = [extraStr,sprintf('_BkgPtSlope%.3gmuCpsS',1e6.*obj.RunAnaObj.ModelObj.BKG_PtSlope)];
                          elseif obj.RunAnaObj.ModelObj.BKG_PtSlope ~=3*1e-06 || obj.RunAnaObj.TwinBias_BKG_PtSlope ~=3*1e-06
                              % twin (KNM-2): twin or model BKG_PtSlope not default 3e-06
                              extraStr = [extraStr,sprintf('_BkgPtSlope%.3gmuCpsS_TwinPtSlope%.3gmuCpsS',...
                                  1e6.*obj.RunAnaObj.ModelObj.BKG_PtSlope,1e6.*obj.RunAnaObj.TwinBias_BKG_PtSlope)];
                          end
                      elseif obj.RunAnaObj.ModelObj.BKG_PtSlope ~=3*1e-06
                          % data (KNM-2): model BKG_PtSlope not default 3e-06
                          extraStr = [extraStr,sprintf('_BkgPtSlope%.3gmuCpsS',1e6.*obj.RunAnaObj.ModelObj.BKG_PtSlope)];
                      end
                  end
                  
                  MakeDir(savedir);
                  
                 % get runlist-name
                 RunList = extractBefore(obj.RunAnaObj.RunData.RunName,'_E0');
                 if isempty(RunList)
                     RunList = obj.RunAnaObj.RunData.RunName;
                 end
                 
                 freeParStr =  ConvertFixPar('freePar',obj.RunAnaObj.fixPar,'Mode','Reverse');
                 if strcmp(obj.RunAnaObj.DataType,'Fake')
                     filename = sprintf('%sKSNX_GridSearch_%s_%s_%s_%.0feVrange_%s_%.0fnGrid%s.mat',...
                         savedir,func2str(obj.RunAnaObj.FakeInitFile),obj.RunAnaObj.DataType,strrep(freeParStr,' ',''),...
                         obj.range,obj.RunAnaObj.chi2,obj.nGridSteps,extraStr);
                 else
                     filename = sprintf('%sKSN%.0f_GridSearch_%s_%s_%s_%.0feVrange_FSD%s_%s_%.0fnGrid%s.mat',...
                         savedir,str2double(obj.RunAnaObj.DataSet(end)),RunList,obj.RunAnaObj.DataType,strrep(freeParStr,' ',''),...
                         obj.range,obj.RunAnaObj.FSDFlag,obj.RunAnaObj.chi2,obj.nGridSteps,extraStr);
                 end
                 
        end
        function plotname = DefPlotName(obj)
            % generic plot name
            filename = obj.GridFilename;
            filename = extractBetween(filename,sprintf('KSN%.0f_',str2double(obj.RunAnaObj.DataSet(end))),'.mat');
            plotdir = sprintf('%sSterileAnalysis/plots/%s/%s/',...
                     getenv('SamakPath'),obj.RunAnaObj.DataSet,obj.RunAnaObj.DataType);
            MakeDir(plotdir);
            if isfloat(obj.RandMC) && strcmp(obj.RunAnaObj.DataType,'Twin')
                plotdir = strrep(savedir,'Twin/','TwinRandomizedMC/'); 
            end
            
            if strcmp(obj.NullHypothesis,'ON')
                filename = [filename,'_NH'];
            end
            
             plotname = sprintf('%s%s',plotdir,filename{:});
        end
        function titleStr = GetPlotTitle(obj,varargin)
            p=inputParser;
            p.addParameter('Mode','all',@(x)ismember(x,{'all','chi2','data'}));
            p.parse(varargin{:});
            Mode = p.Results.Mode;
            
            % some useful default title
            if strcmp(obj.RunAnaObj.DataType,'Real')
                DataStr = 'Data';
            elseif strcmp(obj.RunAnaObj.DataType,'Twin')
                DataStr = 'Twin';
            elseif strcmp(obj.RunAnaObj.DataType,'Fake')
                DataStr = 'Simulation';
            end
            
            if strcmp(obj.RunAnaObj.chi2,'chi2Stat')
                chi2Str = 'stat. only';
            else
                chi2Str = 'stat. and syst.';
            end
            
            if ~contains(ConvertFixPar('freePar',obj.RunAnaObj.fixPar,'Mode','Reverse'),'mNu')
                fitparStr = sprintf('{\\itm}_\\nu^2 = %.3g eV^2',obj.RunAnaObj.ModelObj.mnuSq_i);
            elseif obj.RunAnaObj.pullFlag == 99
                 fitparStr = sprintf('{\\itm}_\\nu^2 free');
            else
                 fitparStr = sprintf('{\\itm}_\\nu^2 free + pull');
            end
            
            switch Mode
                case 'all'
            titleStr = sprintf('%s , %.0f eV range (%s) , %s',DataStr,obj.range,chi2Str,fitparStr);
                case 'chi2'
                    titleStr = chi2Str;
                case 'data'
                    titleStr = DataStr;
            end
            
            
        end      
        function InitPlotArg(obj)
            obj.PlotColors =  {rgb('DodgerBlue'),rgb('Orange'),rgb('DarkSlateGray'),rgb('FireBrick'),...
                rgb('Magenta'),rgb('LimeGreen'),rgb('CadetBlue'),rgb('Navy'),...
                rgb('ForestGreen'),rgb('PowderBlue'),rgb('Pink'),rgb('DarkOrange'),rgb('Black'),...
                rgb('ForestGreen'),rgb('PowderBlue'),rgb('Pink'),rgb('DarkOrange')};
            obj.PlotLines = {'-','-.',':','--','-','-.',':','--','-','-.',':','--','-','-.',':','--'};
        end
    end
    
    % small auxillary methods
    methods
        function SetNPfactor(obj)
            if ~strcmp(obj.RunAnaObj.chi2,'chi2Stat') && ~strcmp(obj.RunAnaObj.chi2,'chi2Stat+')
                switch obj.RunAnaObj.DataSet
                    case 'Knm1'
                        obj.RunAnaObj.NonPoissonScaleFactor= 1.064;
                    case 'Knm2'
                        obj.RunAnaObj.NonPoissonScaleFactor= 1.112;
                    otherwise
                        obj.RunAnaObj.NonPoissonScaleFactor= 1;
                end
            elseif strcmp(obj.RunAnaObj.chi2,'chi2Stat')
                obj.RunAnaObj.NonPoissonScaleFactor=1;
            elseif strcmp(obj.RunAnaObj.chi2,'chi2Stat+')
                % keep factor as it is
                obj.RunAnaObj.chi2 = 'chi2Stat';
            end 
        end
    end
    
end

