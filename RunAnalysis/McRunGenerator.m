classdef McRunGenerator < handle
    properties(Access=public)
        MCFlag;
        RunObj;
        RealData;
        TwinObj;

        WGTS_CD_MolPerCm2;
        WGTS_CD_MolPerCm2_SubRun;
        WGTS_MolFrac_TT   ;
        WGTS_MolFrac_TT_SubRun  ;
        WGTS_MolFrac_DT   ;
        WGTS_MolFrac_DT_SubRun   ;
        WGTS_MolFrac_HT          ;
        WGTS_MolFrac_HT_SubRun   ;
        qU               ;
        qUfrac;
        TimeSec;
        
        AsimovData;
        
        %fake run properties
        InitFile;
    end
    methods
        function obj = McRunGenerator(varargin)
            p = inputParser;
            p.addParameter('MCFlag','Twin',@(x)ismember(x,{'Twin','Fake'})); %later more features
            p.addParameter('RunObj','',@(x)isa(x,'RunAnalysis'));
            p.addParameter('InitFile','',@(x)isa(x,'function_handle'));
            p.parse(varargin{:});
            
            obj.MCFlag    = p.Results.MCFlag;
            obj.RunObj    = p.Results.RunObj;
            obj.InitFile  = p.Results.InitFile;
            
            if strcmp(obj.MCFlag,'Twin') && isempty(obj.RunObj)
                fprintf(2,'Error: give a RunAnalysis object for MCFlag "Twin" \n');
                return
            elseif strcmp(obj.MCFlag,'Fake') && isempty(obj.InitFile)
                fprintf(2,'Error: give a InitFile for MCFlag "Fake" \n');
                return
            else
                fprintf('start Monte Carlo run generation \n');
            end
        end % end constructor
    end
    
    methods
        function ComputeTwinRun(obj,varargin)
            % label
            obj.RunObj.RunData.matFilePath = [getenv('SamakPath'),'/tritium-data/mat/Twin',obj.RunObj.DataSet,'/'];
            MakeDir(obj.RunObj.RunData.matFilePath);

            % retrieve real data
            chi2_prev = obj.RunObj.chi2;
            obj.RunObj.chi2 = 'chi2Stat';
            
            obj.ReadRealData;
            
            if ~isa(obj.RunObj,'MultiRunAnalysis')
                RunListAll = obj.RunObj.RunNr;
            else
                RunListAll = obj.RunObj.RunList;
            end

            % fit real data
            switch obj.RunObj.FitNBFlag
                case 'ON'
                    [FitBkg, FitN, FitE0] = obj.FitRealData;
                case 'OFF'
                    FitBkg = zeros(numel(RunListAll),1);
                    FitN   = ones(numel(RunListAll),1);
            end
            
            if strcmp(obj.RunObj.TwinBias_Q,'Fit')
                TwinQ = FitE0;
            elseif numel(obj.RunObj.TwinBias_Q)==1
                TwinQ = repmat(obj.RunObj.TwinBias_Q,[numel(RunListAll),1]);
            elseif numel(obj.RunObj.TwinBias_Q)==numel(RunListAll) 
                TwinQ = obj.RunObj.TwinBias_Q;
            else
                fprintf('Invalid TwinBias_Q \n')
                return
            end
            
            % Set Biases (optional)
            obj.SetSameTwinBias; % only for MultiRunAnalysis
            obj.SetTwinBias;
            
            %% simulate twin run
            ref_func = @ref_RunAnalysis;
             
            obj.TwinObj = cell(numel(RunListAll),1);
            
            progressbar('computing twins');
            for i=1:numel(RunListAll)
                progressbar(i/numel(RunListAll));
               savepath = [getenv('SamakPath'),sprintf('tritium-data/mat/Twin%s/',obj.RunObj.DataSet)];
               MakeDir(savepath);
               savename_twin =  [savepath,'Twin',num2str(RunListAll(i)),obj.RunObj.SetTwinOrFakeFileName,'.mat'];
               if exist(savename_twin,'file')
                   continue
               end
               [TTFSD,DTFSD,HTFSD] = obj.RunObj.SetDefaultFSD;
               
               % get real data of single run
               if numel(RunListAll)==1
                   savename_real = [obj.RealData.matFilePath,obj.RealData.matFileName];
               else
                   savename_real = [obj.RealData.matFilePath,obj.RealData.matFileName{i}];
               end
               RunData_real             = load(savename_real);
               RunData_real.matFilePath = obj.RealData.matFilePath;
               RunData_real.RunName      = num2str(RunListAll(i));
               
               if strcmp(obj.RunObj.KTFFlag,'WGTSMACE_NIS1')
                   KTFFlag = 'WGTSMACE';
                   NIS = 1;
               else
                   KTFFlag =obj.RunObj.KTFFlag;
                   NIS = 7;
               end

               TBDarg = {'ISCS',obj.RunObj.ModelObj.ISCS,...
                    'recomputeRF','OFF',...
                    'RadiativeFlag',obj.RunObj.RadiativeFlag,...
                    'ELossFlag',obj.RunObj.ELossFlag,...
                    'SynchrotronFlag',obj.RunObj.SynchrotronFlag,...
                    'AngularTFFlag',obj.RunObj.AngularTFFlag,...
                    'WGTS_CD_MolPerCm2',obj.WGTS_CD_MolPerCm2(i),...
                    'WGTS_CD_MolPerCm2_SubRun',obj.WGTS_CD_MolPerCm2_SubRun(:,i),...
                    'WGTS_MolFrac_TT',obj.WGTS_MolFrac_TT(i),...
                    'WGTS_MolFrac_DT',obj.WGTS_MolFrac_DT(i),...
                    'WGTS_MolFrac_HT',obj.WGTS_MolFrac_HT(i),...
                    'WGTS_MolFrac_TT_SubRun',obj.WGTS_MolFrac_TT_SubRun(:,i),...
                    'WGTS_MolFrac_DT_SubRun',obj.WGTS_MolFrac_DT_SubRun(:,i),...
                    'WGTS_MolFrac_HT_SubRun',obj.WGTS_MolFrac_HT_SubRun(:,i),...
                    'MACE_Ba_T',obj.RealData.MACE_Ba_T(i),...
                    'MACE_Bmax_T',obj.RealData.MACE_Bmax_T(i),...
                    'qU',obj.qU(:,i),...
                    'qUfrac',obj.qUfrac(:,i),...
                    'TimeSec',obj.TimeSec(i),...
                    'BKG_RateAllFPDSec',0,...% added later (pixelwise)
                    'PixList',obj.RunObj.PixList,...
                    'RingList',obj.RunObj.RingList,...
                    'RingMerge',obj.RunObj.RingMerge...
                    'DopplerEffectFlag',obj.RunObj.DopplerEffectFlag,...
                    'DTFSD',DTFSD,...
                    'HTFSD',HTFSD,...
                    'TTFSD',TTFSD,...
                    'Q_i',TwinQ(i),...
                    'KTFFlag',KTFFlag,...
                    'NIS',NIS,...
                    'FSD_Sigma',obj.RunObj.TwinBias_FSDSigma,...
                    'BKG_PtSlope',obj.RunObj.TwinBias_BKG_PtSlope,...
                    'nRuns',1};

                switch obj.RunObj.AnaFlag
                    case 'StackPixel'
                        obj.TwinObj{i} =  ref_func(RunData_real,TBDarg{:},...
                            'FPD_Segmentation','OFF');
                    case 'SinglePixel'
                        obj.TwinObj{i} = ref_func(RunData_real,TBDarg{:},...
                            'FPD_Segmentation','SINGLEPIXEL');
                    case 'MultiPixel'
                        obj.TwinObj{i} = ref_func(RunData_real,TBDarg{:},...
                            'FPD_Segmentation','MULTIPIXEL','nTeBinningFactor',5);
                    case 'Ring'
                        obj.TwinObj{i} =  ref_func(RunData_real,TBDarg{:},...
                            'FPD_Segmentation','RING');
                end
               
                %% Calculate Asimov pixel-wise spectra
                obj.TwinObj{i}.ComputeTBDDS('N_bias',FitN(i),'mSq_bias',obj.RunObj.TwinBias_mnuSq);
                obj.TwinObj{i}.ComputeTBDIS;
                TBDIS_NoBkg   = obj.TwinObj{i}.TBDIS;
                %TBDIS = cell2mat(cellfun(@(x) x.TBDIS,obj.TwinObj,'UniformOutput',0)');

                EffCorr = ones(obj.TwinObj{i}.nqU,148);
                WGTS_CD_MolPerCm2 = obj.WGTS_CD_MolPerCm2(i);
                %ISXsection        = obj.TwinObj{i}.ISXsection;
                WGTS_CD_MolPerCm2_SubRun= obj.WGTS_CD_MolPerCm2_SubRun(:,i);
                WGTS_MolFrac_TT=obj.WGTS_MolFrac_TT(i);
                WGTS_MolFrac_DT=obj.WGTS_MolFrac_DT(i);
                WGTS_MolFrac_HT=obj.WGTS_MolFrac_HT(i);
                WGTS_MolFrac_TT_SubRun=obj.WGTS_MolFrac_TT_SubRun(:,i);
                WGTS_MolFrac_DT_SubRun=obj.WGTS_MolFrac_DT_SubRun(:,i);
                WGTS_MolFrac_HT_SubRun=obj.WGTS_MolFrac_HT_SubRun(:,i);
                
                if contains(obj.RunObj.DataSet,'FirstTritium') % e.g. First Tritium -> no rings
                    BkgPix = (repmat(FitBkg(i),[1,148])./numel(obj.TwinObj{i}.FPD_PixList))';
                else
                    PixLogic     = ismember(1:148,obj.TwinObj{i}.FPD_PixList);
                    BkgRelPix    = obj.GetRingWiseBkg;
                    BkgPix       = BkgRelPix.*FitBkg(i); % Background per pixel. Non-active pixels have BkgPix = 0
                end
                BkgPixSubRun = obj.TwinObj{i}.qUfrac.*obj.TwinObj{i}.TimeSec.*BkgPix';
                
               
                switch obj.RunObj.AnaFlag
                    case 'StackPixel'
                        nqU =  numel(WGTS_MolFrac_HT_SubRun);
                        PixLogic              = repmat(ismember(1:148,obj.TwinObj{i}.FPD_PixList),[nqU,1]); % logical array, says whether pixel is active or non-active
                        TBDIS_NoBkg           = PixLogic.*repmat(TBDIS_NoBkg,[1,148])./numel(obj.TwinObj{i}.FPD_PixList); % calculate spectrum per pixel
                        TBDIS                 = TBDIS_NoBkg+BkgPixSubRun;
                        TBDISE                = sqrt(TBDIS_NoBkg);
                        qU                    = repmat(obj.qU(:,i),[1,148]);
                        qUfrac                = repmat(obj.qUfrac(:,i),[1,148]);
                        TimeSec               = repmat(obj.TimeSec(i),[148,1]);%obj.RunData.TimeSec.*(obj.TwinBias_Time/obj.RunData.TimeSec);
                        TimeperSubRunperPixel = qUfrac.*obj.TimeSec(i); 
                end
                
                % copy real data file to twin
                copyfile(savename_real,savename_twin);
                
                save(savename_twin,...% overwrite copied file with new parameters
                    'TBDarg',... % keep track what's inside twins
                    'TBDIS','TBDISE','EffCorr','TimeSec',...
                    'TimeperSubRunperPixel',...%'TimeperPixel',...
                    'WGTS_CD_MolPerCm2','WGTS_CD_MolPerCm2_SubRun',...
                    'WGTS_MolFrac_TT','WGTS_MolFrac_DT','WGTS_MolFrac_HT',...
                    'WGTS_MolFrac_DT_SubRun','WGTS_MolFrac_HT_SubRun','WGTS_MolFrac_TT_SubRun',...
                    'qU','qUfrac',...
                    '-append'); % do not overwrite entire file, only these variables 
            end

            obj.RunObj.chi2 = chi2_prev;  
        end
    end
    methods
        function ReadRealData(obj)
            % read real data
            tmpDataType = obj.RunObj.DataType; % tmp save here, set back at end of function
            obj.RunObj.DataType = 'Real';
            if isa(obj.RunObj,'MultiRunAnalysis')
                obj.RunObj.StackRuns;
                obj.RealData = obj.RunObj.SingleRunData;
                obj.RealData.matFilePath = obj.RunObj.RunData.matFilePath;
            else
                obj.RunObj.ReadData;
                obj.RealData = obj.RunObj.RunData;
            end
            
            obj.RunObj.DataType = tmpDataType;
        end
        function [FitBkg, FitN, FitE0] = FitRealData(obj)
            % simulate and fit real data
            tmpDataType      = obj.RunObj.DataType;
            tmpexclDataStart = obj.RunObj.exclDataStart;
            tmpfixPar        = obj.RunObj.fixPar;
            tmpKTF = obj.RunObj.KTFFlag;
            obj.RunObj.KTFFlag = 'WGTSMACE';
            
            %set range and fixed parameter to something reliable
            if strcmp(obj.RunObj.DataSet,'FirstTritium.katrin')
                obj.RunObj.exclDataStart = 12; % 125eV range
            elseif strcmp(obj.RunObj.DataSet,'Knm1')
                obj.RunObj.exclDataStart = 13; %27 subruns -> 40 eV range
            else% 40eV range
                [~, obj.RunObj.exclDataStart] = min(abs(obj.RunObj.RunData.qU-18574+40));
            end
            
            nPixels =numel(obj.RunObj.RunData.MACE_Ba_T);
            nPar = 4*nPixels+12; %avaibale fit parameter
            obj.RunObj.DataType = 'Real';
            obj.RunObj.fixPar   = ConvertFixPar('freePar','E0 Bkg Norm',...
                'nPar',nPar,'nPixels',nPixels); %fix everything, except endpoint normalization background 
            
            if isa(obj.RunObj,'MultiRunAnalysis')
                obj.RunObj.SimulateStackRuns;
                obj.RunObj.FitRunList('Recompute','ON');
                FitBkg = obj.RunObj.SingleRun_FitResults.chi2Stat.B;
                FitN   = obj.RunObj.SingleRun_FitResults.chi2Stat.N;
                FitE0  = obj.RunObj.SingleRun_FitResults.chi2Stat.E0;
            else
                obj.RunObj.SimulateRun;
                obj.RunObj.Fit;
                FitBkg = obj.RunObj.FitResult.par(3)+obj.RunObj.ModelObj.BKG_RateSec_i;
                FitN   = obj.RunObj.FitResult.par(4);
                FitE0  = obj.RunObj.FitResult.par(2)+obj.RunObj.ModelObj.Q_i;
            end
            
            %reset to initial values
            obj.RunObj.DataType      = tmpDataType;
            obj.RunObj.exclDataStart = tmpexclDataStart; 
            obj.RunObj.fixPar        = tmpfixPar;
            obj.RunObj.KTFFlag       = tmpKTF;
        end
        function BkgRelPix = GetRingWiseBkg(obj)
            % Get (relative) Background for 12 rings
            RingList = 1:12;
            
            % only 1 per dat set
            switch obj.RunObj.DataSet
                case 'Knm1'
                    RunList = 'KNM1';
                case 'Knm2'
                    RunList = 'KNM2_Prompt';
            end
            % label
            savepath = [getenv('SamakPath'),sprintf('tritium-data/fit/%s/',obj.RunObj.DataSet)];
            MakeDir(savepath);
            savename = sprintf('RelBkgPerPix_%s_%.0fnPix.mat',upper(obj.RunObj.DataSet),numel(obj.RunObj.PixList));
            savefile = [savepath,savename];
            
            if exist(savefile,'file') %first load to check if pixels are identical to current pixel selection
                d = importdata(savefile);
            end
            
            if exist(savefile,'file') && all(d.PixList==obj.RunObj.PixList) 
                load(savefile,'BkgRelPix');
            else
                % free parameter = E0 Bkg Norm
                freePar_local = ConvertFixPar('freePar','E0 Bkg Norm',...
                    'nPar',obj.RunObj.nPar,'nPixels',numel(obj.RunObj.RunData.MACE_Ba_T));
                
                % fit only region close to endpoint (40eV)
                [~, exclDataStart_local] = min(abs(obj.RunObj.RunData.qU-obj.RunObj.ModelObj.Q_i+40));
                
                RunArg = {'RunList',RunList,'RingList',RingList,...
                    'exclDataStart',exclDataStart_local,'fixPar',freePar_local,...
                    'RingMerge','None',...
                    'DataType','Real'};
                
                RunAnaObj = MultiRunAnalysis(RunArg{:});
                R = RingAnalysis('RunAnaObj',RunAnaObj,'RingList',RingList);
                R.FitRings('SaveResult','ON','RecomputeFlag','OFF','AsymErr','OFF');
                
                % background from fit per ring
                par = R.FitResult.par;
                BkgRing = (par(:,3)+cell2mat(arrayfun(@(x) x.ModelObj.BKG_RateSec_i,R.MultiObj,'UniformOutput',0))); 
                
                %number of pixels per ring in fit
                PixList = RunAnaObj.PixList;
                nPixFit = cell2mat(arrayfun(@(x) numel(x.PixList),R.MultiObj,'UniformOutput',0))';
                
                % background from fit: convert to pixelwise bkg
                % asumption: every pixel within a ring sees the same background
                BkgPix = zeros(148,1);
                for i=1:numel(RingList)
                    BkgPix(RunAnaObj.RingPixList{i}) = repmat(BkgRing(i),numel(RunAnaObj.RingPixList{i}),1)./nPixFit(i); % background per pixel
                end
                
                % relative background from fit per pixel
                BkgRelPix = BkgPix./sum(BkgPix);
                
                % save
                MakeDir(savepath);  
                save(savefile,'BkgRelPix','BkgPix','BkgRing','nPixFit','PixList','RunArg');
            end
        end
        function SetSameTwinBias(obj)
            % Set Twin Bias as such that all single runs end up with same qU/ rho d/,... values as the average (stacked) run
            % Only applied for MultiRunAnalysis Object!
           
            if ~isa(obj.RunObj,'MultiRunAnalysis')
                return
                
            elseif isa(obj.RunObj,'MultiRunAnalysis') && (...
                    strcmp(obj.RunObj.Twin_SameqUFlag,'ON') ||strcmp(obj.RunObj.Twin_SameqUfracFlag,'ON') ...
                    || strcmp(obj.RunObj.Twin_SameCDFlag,'ON') || strcmp(obj.RunObj.Twin_SameIsotopFlag,'ON')) 
                
                % set all qU values to average
                if strcmp(obj.RunObj.Twin_SameqUFlag,'ON')
                    obj.RunObj.TwinBias_qU     = obj.RunObj.RunData.qU-obj.RealData.qU;
                end
                
                % set qUfrac values to average
                if strcmp(obj.RunObj.Twin_SameqUfracFlag,'ON')
                    obj.RunObj.TwinBias_qUfrac = obj.RunObj.RunData.qUfrac-obj.RealData.qUfrac;
                     
                end
                
                % set column density to average value
                if strcmp(obj.RunObj.Twin_SameCDFlag,'ON')
                    obj.RunObj.TwinBias_WGTS_CD_MolPerCm2 = obj.RunObj.RunData.WGTS_CD_MolPerCm2./obj.RealData.WGTS_CD_MolPerCm2;
                end
                
                % set isotopologue concentrations to average value
                if strcmp(obj.RunObj.Twin_SameIsotopFlag,'ON')
                    obj.RunObj.TwinBias_WGTS_MolFrac_TT = ...
                        obj.RunObj.RunData.WGTS_MolFrac_TT./obj.RealData.WGTS_MolFrac_TT;
                    obj.RunObj.TwinBias_WGTS_MolFrac_DT = ...
                        obj.RunObj.RunData.WGTS_MolFrac_DT./obj.RealData.WGTS_MolFrac_DT;
                    obj.RunObj.TwinBias_WGTS_MolFrac_HT = ...
                        obj.RunObj.RunData.WGTS_MolFrac_HT./obj.RealData.WGTS_MolFrac_HT;
                end
            end
        end
        function SetTwinBias(obj)
             %% bias SC
            obj.WGTS_CD_MolPerCm2        = obj.RealData.WGTS_CD_MolPerCm2       .*obj.RunObj.TwinBias_WGTS_CD_MolPerCm2;
            obj.WGTS_CD_MolPerCm2_SubRun = obj.RealData.WGTS_CD_MolPerCm2_SubRun.*obj.RunObj.TwinBias_WGTS_CD_MolPerCm2;
            obj.WGTS_MolFrac_TT          = obj.RealData.WGTS_MolFrac_TT         .*obj.RunObj.TwinBias_WGTS_MolFrac_TT;
            obj.WGTS_MolFrac_TT_SubRun   = obj.RealData.WGTS_MolFrac_TT_SubRun  .*obj.RunObj.TwinBias_WGTS_MolFrac_TT;
            obj.WGTS_MolFrac_DT          = obj.RealData.WGTS_MolFrac_DT         .*obj.RunObj.TwinBias_WGTS_MolFrac_DT;
            obj.WGTS_MolFrac_DT_SubRun   = obj.RealData.WGTS_MolFrac_DT_SubRun  .*obj.RunObj.TwinBias_WGTS_MolFrac_DT;
            obj.WGTS_MolFrac_HT          = obj.RealData.WGTS_MolFrac_HT         .*obj.RunObj.TwinBias_WGTS_MolFrac_HT;
            obj.WGTS_MolFrac_HT_SubRun   = obj.RealData.WGTS_MolFrac_HT_SubRun  .*obj.RunObj.TwinBias_WGTS_MolFrac_HT;
            obj.qU                       = obj.RealData.qU+obj.RunObj.TwinBias_qU;
            obj.qUfrac                   = obj.RealData.qUfrac+obj.RunObj.TwinBias_qUfrac;
                                                        
            % bias time
            if ~isa(obj.RunObj,'MultiRunAnalysis') || strcmp(obj.RunObj.Twin_SameqUfracFlag,'OFF')
                
                if isempty(obj.RunObj.TwinBias_Time)
                    obj.TimeSec = obj.RealData.TimeSec;
                elseif numel(obj.RunObj.TwinBias_Time)==1
                    obj.TimeSec =  repmat(obj.RunObj.TwinBias_Time,size(obj.RealData.TimeSec));
                else
                    obj.TimeSec = obj.RunObj.TwinBias_Time;
                end
                
            elseif strcmp(obj.RunObj.Twin_SameqUfracFlag,'ON')
                if isempty(obj.RunObj.TwinBias_Time)
                    obj.TimeSec = repmat(mean(obj.RealData.TimeSec),size(obj.RealData.TimeSec)); % set to average time
                elseif numel(obj.RunObj.TwinBias_Time)==1
                    obj.TimeSec =  repmat(obj.RunObj.TwinBias_Time,size(obj.RealData.TimeSec));
                else
                    fprintf(2,'multiple time given, but SameqUfracFlag ON! give only 1 time for all twins or switch off flag \n')
                    return
                end
            end
        end    
    end

    methods % fake methods
        function SetFakeBias(obj,varargin)
            p=inputParser;
            p.addParameter('qU_Err',obj.RunObj.TwinBias_qU,@(x)isfloat(x));
            p.addParameter('qUfrac_RelErr',obj.RunObj.TwinBias_qUfrac,@(x)isfloat(x));
            p.addParameter('CD_RelErr',obj.RunObj.TwinBias_WGTS_CD_MolPerCm2,@(x)isfloat(x));
            p.parse(varargin{:});
            qU_Err         = p.Results.qU_Err;
            qUfrac_RelErr = p.Results.qUfrac_RelErr;
            CD_RelErr     = p.Results.CD_RelErr;
            % create first TwinRunObj to access the default parameters
            TmpObj = obj.InitFile();
            
            if ~isa(obj.RunObj,'MultiRunAnalysis')
                nRuns = 1;
            else
                nRuns = obj.RunObj.nRuns;
            end
            
            if isempty(obj.RunObj.TwinBias_qUfrac)
               qUfrac_RelErr = 0;
            end
            % vary input parameters according to relative error (RelErr) and distribution (VarDist)
            obj.qU                = TmpObj.qU+qU_Err.*randn(TmpObj.nqU,nRuns);
            obj.qUfrac            = TmpObj.qUfrac.*(1+qUfrac_RelErr.*randn(TmpObj.nqU,nRuns));
            
            obj.qUfrac            = obj.VaryParameter(TmpObj.qUfrac,qUfrac_RelErr,'Rel');
            NormFactor            =  sum(TmpObj.qUfrac)./sum(obj.qUfrac);
            obj.qUfrac            = obj.qUfrac.*NormFactor; %renormalize
            
            if CD_RelErr~=1
                obj.WGTS_CD_MolPerCm2 = obj.VaryParameter(TmpObj.WGTS_CD_MolPerCm2,CD_RelErr,'Rel');
            else
                obj.WGTS_CD_MolPerCm2 =  TmpObj.WGTS_CD_MolPerCm2;
            end
%             obj.TimeSec           = obj.VaryParameter(obj.TwinObj{1}.TimeSec,obj.TimeSec_RelErr);
%             obj.WGTS_MolFrac_TT   = obj.VaryParameter(obj.TwinObj{1}.WGTS_MolFrac_TT,obj.WGTS_MolFrac_TT_RelErr);
%             obj.WGTS_MolFrac_DT   = obj.VaryParameter(obj.TwinObj{1}.WGTS_MolFrac_DT,obj.WGTS_MolFrac_DT_RelErr);
%             obj.WGTS_MolFrac_HT   = obj.VaryParameter(obj.TwinObj{1}.WGTS_MolFrac_HT,obj.WGTS_MolFrac_HT_RelErr);
        end
        function ComputeFakeRun(obj,varargin)
             obj.SetFakeBias; %draw random parameters according to their TwinBias (for qU + qUfrac)
            
            % calculate fake runs
            fprintf('calculate fake run %s \n',func2str(obj.InitFile))
            
            if ~isa(obj.RunObj,'MultiRunAnalysis')
                RunListAll = obj.RunObj.RunNr;
            else
                RunListAll = obj.RunObj.RunList;
            end
            
            if numel(obj.RunObj.TwinBias_Q)==1
                TwinQ = repmat(obj.RunObj.TwinBias_Q,[numel(RunListAll),1]);
            elseif numel(obj.RunObj.TwinBias_Q)==numel(RunListAll)
                TwinQ = obj.RunObj.TwinBias_Q;
            else
                fprintf(2,'Twin bias Q invalied \n');
                return
            end
            
            progressbar('Compute Fake MC runs')
            for i=1:numel(RunListAll)
                progressbar(i/numel(RunListAll));
                TBDArg = {'PixList',obj.RunObj.PixList,...
                        'Q_i',TwinQ(i),...
                        'qU',obj.qU(:,i),...
                        'qUfrac',obj.qUfrac(:,i),...
                        'WGTS_CD_MolPerCm2',obj.WGTS_CD_MolPerCm2(i)};
                    
                if ~isempty(obj.RunObj.TwinBias_Time)
                    obj.TwinObj = obj.InitFile(TBDArg{:},'TimeSec',obj.RunObj.TwinBias_Time);
                else
                    obj.TwinObj = obj.InitFile(TBDArg{:});
                end
                %% Calculate Asimov pixel-wise spectra
                obj.TwinObj.ComputeTBDDS;
                obj.TwinObj.ComputeTBDIS;
                
                % save into file
                TBDIS                    = obj.TwinObj.TBDIS;
                EffCorr                  = ones(obj.TwinObj.nqU,148);
                WGTS_CD_MolPerCm2        = obj.TwinObj.WGTS_CD_MolPerCm2;
                WGTS_CD_MolPerCm2_SubRun = obj.TwinObj.WGTS_CD_MolPerCm2_SubRun;
                WGTS_MolFrac_TT          = obj.TwinObj.WGTS_MolFrac_TT;
                WGTS_MolFrac_DT          = obj.TwinObj.WGTS_MolFrac_DT;
                WGTS_MolFrac_HT          = obj.TwinObj.WGTS_MolFrac_HT;
                WGTS_MolFrac_TT_SubRun   = obj.TwinObj.WGTS_MolFrac_TT_SubRun;
                WGTS_MolFrac_DT_SubRun   = obj.TwinObj.WGTS_MolFrac_DT_SubRun;
                WGTS_MolFrac_HT_SubRun   = obj.TwinObj.WGTS_MolFrac_HT_SubRun;
                ISXsection               = obj.TwinObj.ISXsection;
                
                if strcmp(obj.TwinObj.FPD_Segmentation,'OFF')
                    TBDIS                 = repmat(TBDIS,[1,148])./numel(obj.TwinObj.FPD_PixList);
                    qU                    = repmat(obj.TwinObj.qU,[1,148]);
                    qUfrac                = repmat(obj.TwinObj.qUfrac,[1,148]);
                    TimeSec               = repmat(obj.TwinObj.TimeSec,[148,1]);%obj.RunData.TimeSec.*(obj.TwinBias_Time/obj.RunData.TimeSec);
                    MACE_Ba_T             = repmat(obj.TwinObj.MACE_Ba_T,[148,1]);
                    MACE_Bmax_T           = repmat(obj.TwinObj.MACE_Bmax_T,[148,1]);
                else
                    TBDIS    = zeros(obj.TwinObj.nqU,148);
                    qU        = zeros(obj.TwinObj.nqU,148);
                    qUfrac    = zeros(obj.TwinObj.nqU,148);
                    TimeSec   = zeros(148,1);
                    MACE_Ba_T = zeros(148,1);
                    MACE_Bmax_T = zeros(148,1);
                    nPixRing = cell2mat(cellfun(@(x) numel(x),obj.RunObj.RingPixList,'UniformOutput',0));
                    switch obj.RunObj.RingMerge
                        case 'Full'
                            [~,AllRingPixList] = Ring2PixelCombi(1:4,1:148);
                        case 'Half'
                            [~,AllRingPixList] = Ring2PixelHalfCombi(1:2,1:148);
                    end
                    AllnPixRing = cell2mat(cellfun(@(x) numel(x),AllRingPixList,'UniformOutput',0));
                    for i=1:obj.RunObj.nRings
                        TBDIS(:,AllRingPixList{i})   = repmat(obj.TwinObj.TBDIS(:,i)./nPixRing(i),[1,AllnPixRing(i)]);
                        qU(:,AllRingPixList{i})      = repmat(obj.TwinObj.qU(:,i),[1,AllnPixRing(i)]);
                        qUfrac(:,AllRingPixList{i})  = repmat(obj.TwinObj.qUfrac(:,i),[1,AllnPixRing(i)]);
                        TimeSec(AllRingPixList{i})   = repmat(obj.TwinObj.TimeSec(i),[AllnPixRing(i),1]);
                        MACE_Ba_T(AllRingPixList{i}) = repmat(obj.TwinObj.MACE_Ba_T(i),[AllnPixRing(i),1]);
                        MACE_Bmax_T(AllRingPixList{i}) = repmat(obj.TwinObj.MACE_Bmax_T(i),[AllnPixRing(i),1]);
                    end
                end
                TBDISE                = sqrt(TBDIS);
                TimeperSubRunperPixel = qUfrac.*TimeSec(1);
                StartTimeStamp = datetime('today');
                
                if contains(func2str(obj.InitFile),'KNM1')
                    savepath = [getenv('SamakPath'),'tritium-data/mat/FakeKnm1/'];
                elseif contains(func2str(obj.InitFile),'KNM2')
                    savepath = [getenv('SamakPath'),'tritium-data/mat/FakeKnm2/'];
                else
                    savepath = [getenv('SamakPath'),'tritium-data/mat/Fake/'];
                end
                MakeDir(savepath);
                
                savename = [extractAfter(func2str(obj.InitFile),'ref_'),obj.RunObj.SetTwinOrFakeFileName,...
                    sprintf('_Run%.0f',RunListAll(i)),'.mat'];
                fprintf('saving fake run %s \n',[savepath,savename]);
                
                if isa(ISXsection,'function_handle')
                    ISXsection = ISXsection(18575);
                end
                
                RW_BiasVoltage = 0;

                save([savepath,savename],...
                    'TBDIS','TBDISE','EffCorr','TimeSec',...
                    'TimeperSubRunperPixel',...%'TimeperPixel',...
                    'WGTS_CD_MolPerCm2','WGTS_CD_MolPerCm2_SubRun',...
                    'WGTS_MolFrac_TT','WGTS_MolFrac_DT','WGTS_MolFrac_HT',...
                    'WGTS_MolFrac_DT_SubRun','WGTS_MolFrac_HT_SubRun','WGTS_MolFrac_TT_SubRun',...
                    'ISXsection','MACE_Bmax_T',...
                    'MACE_Ba_T','qU','qUfrac',...
                    'StartTimeStamp',...
                    'RW_BiasVoltage',...
                    '-v7.3','-nocompression');
            end
        end
        function CleanUpFakeRuns(obj)
            % deletes all fake MC runs, created with this InitFile
            savedir = [getenv('SamakPath'),'tritium-data/mat/',extractAfter(func2str(obj.InitFile),'ref_')];
            system(['rm ',savedir,'/*']);
            fprintf(2,'Deleting all fake MC runs in directory: tritium-data/mat/%s ! \n',extractAfter(func2str(obj.InitFile),'ref_'))
        end
        function out = VaryParameter(obj,Parameter,Err,Mode)
            if ~isa(obj.RunObj,'MultiRunAnalysis')
                nRuns = 1;
            else
                nRuns = obj.RunObj.nRuns;
            end
            
            if strcmp(Mode,'Rel')
                out = Parameter.*(1+randn(numel(Parameter),nRuns).*Err);
            elseif strcmp(Mode,'Abs')
                out = Parameter+randn(numel(Parameter),nRuns).*Err;
            end
            %                 case 'Uniform'
            %                     if strcmp(Mode,'Rel')
            %                         out = Parameter.*(1+(2.*rand(numel(Parameter),obj.nRuns)-1).*RelErr);
            %                     elseif strcmp(Mode,'Abs')
%                         out = Parameter+(2.*rand(numel(Parameter),obj.nRuns)-1).*RelErr;
%                     end
%             end         
        end
    end

end