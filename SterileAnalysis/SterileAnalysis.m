classdef SterileAnalysis < handle
    % class to collect and organise sterile analysis routines
    %%
    properties (Access=public)
        RunAnaObj; % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        nGridSteps;   
        SmartGrid;
        RecomputeFlag;
        SysEffect;
        RandMC; 
        range; % in eV
        
        % grid parameters
        mNu4Sq;
        sin2T4;   
        chi2;
        
        % contour
        mNu4Sq_contour;
        sin2T4_contour;
        
        % best fit
        chi2_ref;
        mNu4Sq_bf;
        sin2T4_bf; 
        chi2_bf;
        
        % null hypothesis
        chi2_Null; 
        
        DeltaChi2;
        ConfLevel;
        dof;
        
        
        % plot
        PlotColors;
        PlotLines;
        InterpMode; %lin or spline
    end
    
    methods % constructor
        function obj = SterileAnalysis(varargin)
            p = inputParser;
            p.addParameter('RunAnaObj','', @(x)  isa(x,'RunAnalysis') || isa(x,'MultiRunAnalysis'));
            p.addParameter('nGridSteps',50,@(x)isfloat(x)); % time on server: 25 points ~ 7-10 minutes, 50 points ~ 40 minutes on csltr server
            p.addParameter('SmartGrid','OFF',@(x)ismember(x,{'ON','OFF'})); % work in progress
            p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SysEffect','all',@(x)ischar(x)); % if chi2CMShape: all or only 1
            p.addParameter('RandMC','OFF',@(x)ischar(x) || isfloat(x)); % randomize twins if RandMC is float
            p.addParameter('range',65,@(x)isfloat(x));
            p.addParameter('ConfLevel',95,@(x)isfloat(x));
            p.addParameter('InterpMode','spline',@(x)ismember(x,{'lin','spline'}));
          
            p.parse(varargin{:});
            
            obj.RunAnaObj = p.Results.RunAnaObj;
            obj.nGridSteps    = p.Results.nGridSteps;
            obj.SmartGrid     = p.Results.SmartGrid;
            obj.RecomputeFlag = p.Results.RecomputeFlag;
            obj.SysEffect     = p.Results.SysEffect;
            obj.RandMC        = p.Results.RandMC;
            obj.range         = p.Results.range;
            obj.ConfLevel     = p.Results.ConfLevel;
            obj.InterpMode    = p.Results.InterpMode; 
            
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
        function GridSearch(varargin) %

        end        
    end
    
    methods
        function Interp1Grid(obj,varargin)
            % mesh grid interpolation
            % Finer binning, nicer appearance!
            p = inputParser;
            p.addParameter('nInter',1e3,@(x)isfloat(x));
            p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'})); % force recompute
            p.parse(varargin{:});
            nInter        = p.Results.nInter;
           
            RecomputeFlag = p.Results.RecomputeFlag;
            
            if size(obj.mNu4Sq,1)>=nInter && strcmp(RecomputeFlag,'OFF')
                fprintf('Interp1 stoped - mNuSq size already large than interpolation \n')
                return
            end
            
            %% define maximum m4:
            if obj.range==65
                Maxm4Sq = 59^2;
            elseif obj.range==95 && strcmp(obj.RunAnaObj.DataType,'Real')
                freePar = ConvertFixPar('freePar',obj.RunAnaObj.fixPar,'Mode','Reverse');
                if contains(freePar,'mNu')
                    Maxm4Sq = 83^2;
                else
                    Maxm4Sq =  84.5^2;
                end
            elseif obj.range==40 && strcmp(obj.RunAnaObj.DataType,'Real')
                Maxm4Sq =  36^2;%(obj.range-3)^2;
            else
                Maxm4Sq =  (obj.range-5)^2;
            end
            
            [X,Y] = meshgrid(obj.mNu4Sq(:,1),obj.sin2T4(1,:));
            mNu4tmp = logspace(log10(min(min(obj.mNu4Sq))),log10(Maxm4Sq),nInter);
            obj.mNu4Sq = repmat(mNu4tmp,nInter,1);
            obj.sin2T4 = repmat(logspace(log10(min(min(obj.sin2T4))),log10(max(max(obj.sin2T4))),nInter),nInter,1)';
            obj.chi2   = reshape(interp2(X,Y,obj.chi2,obj.mNu4Sq,obj.sin2T4,obj.InterpMode),nInter,nInter);
            
            obj.chi2(obj.chi2<0) = NaN;
            %             if strcmp(InterpMode,'spline')
            %                 [row, col]    = find(obj.chi2 < 0);
            %                 obj.chi2(row,:) = NaN;
            %                 obj.chi2(:,col) = NaN;
%                 obj.sin2T4(row,col) = NaN;
%                 obj.mNu4Sq(row,col) = NaN;
%             end
            
             obj.chi2_ref = min(min(obj.chi2));
            fprintf('Grid interpolation sucessfull \n');
        end
        function FindBestFit(obj,varargin)
            % best fit parameters of m4 and sin2t4
            % based on interpolated grid
            if size(obj.mNu4Sq,1)<1e3
                obj.Interp1Grid;
            end
            [row, col]    = find(obj.chi2 == min(obj.chi2(:)));
            
            obj.mNu4Sq_bf =   obj.mNu4Sq(row,col);%obj.mNu4Sq(col,row);
            obj.sin2T4_bf = obj.sin2T4(row,col);    
            obj.chi2_bf   = obj.chi2_ref;    
        end
        function CompareBestFitNull(obj,varargin)
            if isempty(obj.chi2_bf)
                obj.FindBestFit;
            end
            
            %  best fit
            fprintf('Best fit sinTsq = %.3f and m4sq = %.1f eV^2 \n',obj.sin2T4_bf,obj.mNu4Sq_bf);
            fprintf('Best fit:        chi2 = %.3f (%.0f dof) -> p-value = %.2f\n',obj.chi2_bf,obj.dof,1-chi2cdf(obj.chi2_bf,obj.dof));
            
            % null 
            fprintf('Null hypothesis: chi2 = %.3f (%.0f dof) -> p-value = %.2f\n',obj.chi2_Null,obj.dof+2,1-chi2cdf(obj.chi2_Null,obj.dof));
            x = linspace(10,99,1e2);
            y =GetDeltaChi2(x,2);
            SignificanceBF = interp1(y,x,obj.chi2_Null-obj.chi2_bf,'spline');
            fprintf('Delta chi2 = %.2f -> %.1f%% C.L. significance \n',obj.chi2_Null-obj.chi2_bf,SignificanceBF);
          
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
            
            obj.nGridSteps = 25;
            obj.InterpMode  = 'lin';
            
            for i=1:numel(Ranges)
                %  progressbar(i/numel(Ranges));
                obj.range = Ranges(i);
                
                obj.RunAnaObj.chi2 = 'chi2Stat';
                obj.LoadGridFile('CheckSmallerN','OFF','CheckLargerN','OFF');
                obj.Interp1Grid('RecomputeFlag','ON','nInter',1e3);
                [M,c]= contour(obj.sin2T4,obj.mNu4Sq,obj.chi2-obj.chi2_ref,...
                    [DeltaChi2_1Par DeltaChi2_1Par]);
                sinStattmp = M(1,2:end);
                mNu4Stattmp = M(2,2:end);
                c.delete;
                
                obj.RunAnaObj.chi2 = 'chi2CMShape';
                obj.LoadGridFile('CheckSmallerN','OFF','CheckLargerN','OFF');
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
    end
    
    %% Plotting
    methods
        function pHandle = ContourPlot(obj,varargin)
            p = inputParser;
            p.addParameter('BestFit','OFF',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('CL',obj.ConfLevel,@(x)isfloat(x));   
            p.addParameter('HoldOn','OFF',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('Color',obj.PlotColors{1},@(x) isfloat(x));
            p.addParameter('LineStyle',obj.PlotLines{1},@(x) ischar(x));
            p.addParameter('PlotSplines','OFF',@(x) ismember(x,{'ON','OFF'}));
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','png'}));
            
            p.parse(varargin{:});  
            CL          = p.Results.CL;      % also works with vector
            HoldOn      = p.Results.HoldOn;
            myColor     = p.Results.Color;
            myLineStyle = p.Results.LineStyle;
            BestFit     = p.Results.BestFit;
            PlotSplines = p.Results.PlotSplines;
            SavePlot    = p.Results.SavePlot;
            
            if strcmp(HoldOn,'ON')
                hold on;
            elseif strcmp(HoldOn,'OFF')
                GetFigure;
            end
            
            obj.DeltaChi2 = GetDeltaChi2(CL,2);
            
             % contour plot
            PlotArg = {'LineWidth',2.5,'LineStyle',myLineStyle};
            if numel(CL)==1
                PlotArg = [PlotArg,{'LineColor',myColor}];
            end
            
            [M,pHandle]= contour(obj.sin2T4,obj.mNu4Sq,obj.chi2-obj.chi2_ref,...
                [obj.DeltaChi2 obj.DeltaChi2],...
                PlotArg{:});
            
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
                colormap('cool');
            end
            
            % best fit
            if strcmp(BestFit,'ON')
                obj.FindBestFit;
                hold on;
                plot(obj.sin2T4_bf,obj.mNu4Sq_bf,'x','MarkerSize',9,'Color',pHandle.LineColor,'LineWidth',pHandle.LineWidth);
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
            legStr = sprintf(' %.0f%% C.L. -',CL);  
            legend(legStr(1:end-1),'EdgeColor',rgb('Silver'),'Location','southwest');
            
              if ~strcmp(SavePlot,'OFF')
                 if strcmp(SavePlot,'ON')
                     plotname = sprintf('%s_Contour_%.2gCL.pdf',obj.DefPlotName,obj.ConfLevel);
                     export_fig(gcf,plotname);
                 elseif strcmp(SavePlot,'png')
                     plotname = sprintf('%s_Contour_%.2gCL.png',obj.DefPlotName,obj.ConfLevel);
                     print(gcf,plotname,'-dpng','-r450');
                 end
                 fprintf('save plot to %s \n',plotname);
             end
             
        end
        function pHandle = ContourPlotOsci(obj,varargin)
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
            p.addParameter('Stereo','ON',@(x) ismember(x,{'ON','OFF'}));

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
            Stereo      = p.Results.Stereo;
            
            [DeltamNu41Sq,sin2T4Sq] = obj.Convert2Osci;
            
            if strcmp(HoldOn,'ON')
                hold on;
            elseif strcmp(HoldOn,'OFF')
                figure('Units','normalized','Position',[0.1,0.1,0.382,0.618]);
            end
            
            obj.DeltaChi2 = GetDeltaChi2(CL,2);
            legHandle = cell(0,0);
            legStr = '';
            savedirOther = [getenv('SamakPath'),'SterileAnalysis/GridSearchFiles/Knm1/Others/'];
            %% Mainz
            if strcmp(Mainz,'ON')
                filenameMainz = sprintf('%scoord_Mainz_95CL.mat',savedirOther);
                dMainz = importdata(filenameMainz);
                pMainz = plot(dMainz.SinSquare2Theta_X,dMainz.DmSquare41_Y,'-.','LineWidth',1.5,'Color',rgb('SeaGreen'));
                legHandle{numel(legHandle)+1} = pMainz;
                legStr = [legStr,{sprintf('Mainz 95%% C.L.')}];
                hold on;
            end
            %% Troitsk
            if strcmp(Troitsk,'ON')
                filenameTroitsk = sprintf('%scoord_Troitsk_95CL.mat',savedirOther);
                dTroitsk = importdata(filenameTroitsk);
                pTroitsk = plot(dTroitsk.SinSquare2Theta_X,dTroitsk.DmSquare41_Y,'-.','LineWidth',1.5,'Color',rgb('Gold'));
                legHandle{numel(legHandle)+1} = pTroitsk;
                legStr = [legStr,{sprintf('Troitsk 95%% C.L.')}];
                hold on;
            end
            %% Prospect
            if strcmp(Prospect,'ON')
                filenameProspect = sprintf('%scoord_Prospect_95CL.mat',savedirOther);
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
                pDANSS = plot(dDANSS.SinSquare2Theta_X,dDANSS.DmSquare41_Y,'-','LineWidth',1,'Color',rgb('LightGreen'));
                legHandle{numel(legHandle)+1} = pDANSS;
                legStr = [legStr,{sprintf('DANSS 95%% C.L.')}];
                hold on;
            end
            %% Stereo
            if strcmp(Stereo,'ON')
                filenameStereo = sprintf('%scoord_Stereo_95CL.mat',savedirOther);
                dStereo = importdata(filenameStereo);
                pStereo = plot(dStereo.SinSquare2Theta_X,dStereo.DmSquare41_Y,'-','LineWidth',1,'Color',rgb('LightGray'));
                legHandle{numel(legHandle)+1} = pStereo;
                legStr = [legStr,{sprintf('Stéréo 95%% C.L.')}];
                hold on;
            end
            %% RAA
            if strcmp(RAA,'ON')
                filenameRAA1 = sprintf('%scoord_RAA_95_A.mat',savedirOther);
                filenameRAA2 = sprintf('%scoord_RAA_95_B.mat',savedirOther);
                dRAA1 = importdata(filenameRAA1,'file');
                hold on;
                dRAA2 = importdata(filenameRAA2,'file');
                pRAA = plot(dRAA1.sith4_X,dRAA1.m4_Y,'-','LineWidth',2,'Color',rgb('Orange'));
                plot(dRAA2.sith4_X,dRAA2.m4_Y,'-','LineWidth',pRAA.LineWidth,'Color',pRAA.Color);
                legHandle{numel(legHandle)+1} = pRAA;
                legStr = [legStr,{sprintf('RAA+GA 95%% CL')}];%-PRD 83, 073006 (2011) -
            end
            %% Neutrino 4
            if strcmp(Neutrino4,'ON')
                filenameN4 = sprintf('%scoord_Neutrino4_123sigma.mat',savedirOther);
                dN4 = importdata(filenameN4);
                pN4 = plot(dN4.SinSquare2Theta_X_2sigma,dN4.DmSquare41_Y_2sigma,'-','LineWidth',1.5,'Color',rgb('Red'));
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
            set(gca,'YScale','log');
            set(gca,'XScale','log');
            xlabel(sprintf('{\\itsin}^2(2\\theta_{ee})'));
            ylabel(sprintf('\\Delta{\\itm}_{41}^2 (eV^2)'));
            PrettyFigureFormat('FontSize',22);
            legStr = [legStr,{sprintf('KATRIN KSN1 %.0f%% C.L. (%s)',obj.ConfLevel,obj.GetPlotTitle('Mode','chi2'))}];
            
            %%
          leg = legend([legHandle{:}],legStr{:},'EdgeColor','none','Location','northoutside',...
              'FontSize',12);
          xlim([1.2e-02 1]);
          
          if numel(legStr)>4
              leg.NumColumns = 2;
          end
          if ~strcmp(SavePlot,'OFF')
              if strcmp(SavePlot,'ON')
                  plotname = sprintf('%s_OsciContour_%.2gCL.pdf',obj.DefPlotName,obj.ConfLevel);
                  export_fig(gcf,plotname);
              elseif strcmp(SavePlot,'png')
                  plotname = sprintf('%s_OsciContour_%.2gCL.png',obj.DefPlotName,obj.ConfLevel);
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
            p.parse(varargin{:});
            CL       = p.Results.CL;
            HoldOn   = p.Results.HoldOn;
            BestFit  = p.Results.BestFit;
            Contour  = p.Results.Contour;
            SavePlot = p.Results.SavePlot;
                 
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
            [M,pContour] = contour(obj.sin2T4,obj.mNu4Sq,obj.chi2-obj.chi2_ref,...
                [obj.DeltaChi2 obj.DeltaChi2],...
                'LineWidth',3,'LineStyle','-','Color',rgb('Black'));
            end
               
            if strcmp(Contour,'ON')  && strcmp(BestFit,'ON')
                leg = legend([pContour,pbf],sprintf('%.0f%% C.L.',obj.ConfLevel),'Best fit',...
                    'EdgeColor','none','Location','southwest','Color',rgb('White'),...
                    'FontSize',get(gca,'FontSize'));
                set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.3]));
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
         
            
            ylim([1,max(max(obj.mNu4Sq))]);
            xlim([min(min(obj.sin2T4)), max(max(obj.sin2T4))]);
            
            title(sprintf('%s',obj.GetPlotTitle),'FontWeight','normal','FontSize',get(gca,'FontSize'));
            
         if obj.range == 65
             xlim([2e-03,0.5]);
         end
             %% save
           if ~strcmp(SavePlot,'OFF')
               name_i = obj.DefPlotName;
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
             obj.LoadGridFile('CheckSmallerN','ON');
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
        
        function PlotmNuSqOverview(obj,varargin)
            p = inputParser;
            p.addParameter('BestFit','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF','png'}));
            p.parse(varargin{:});
            BestFit  = p.Results.BestFit;
            SavePlot = p.Results.SavePlot;
            fixPar_i = obj.RunAnaObj.fixPar;
            pull_i = obj.RunAnaObj.pullFlag;

            %% 1. nuissance nu-mass without pull
            obj.RunAnaObj.fixPar = 'mNu E0 Norm Bkg'; obj.RunAnaObj.InitFitPar;
            obj.RunAnaObj.pullFlag = 99;
            obj.LoadGridFile('CheckSmallerN','ON');
            obj.Interp1Grid('RecomputeFlag','ON');
            pFree = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','OFF',...
                'Color',rgb('ForestGreen'),'LineStyle','-','BestFit',BestFit);

            %%  nuissance nu-mass + pull
            obj.RunAnaObj.pullFlag = 12;
                        obj.LoadGridFile('CheckSmallerN','ON');
            obj.Interp1Grid('RecomputeFlag','ON');
            pPull = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                'Color',rgb('Orange'),'LineStyle','-.','BestFit',BestFit);
            
            %% fixed nu-mass
            obj.RunAnaObj.fixPar = 'E0 Norm Bkg'; obj.RunAnaObj.InitFitPar;
            obj.RunAnaObj.pullFlag = 99;
            obj.LoadGridFile('CheckSmallerN','ON');
            obj.Interp1Grid('RecomputeFlag','ON');
            pFix = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                'Color',rgb('DodgerBlue'),'LineStyle',':','BestFit',BestFit);
            
            PrettyFigureFormat('FontSize',22);
                legend([pFree,pPull,pFix],...
                    sprintf('Free {\\itm}_\\nu^2 without pull term'),...
                    sprintf('Free {\\itm}_\\nu^2 with pull term \\sigma({\\itm}_\\nu^2) = 1.94 eV^2'),...
                    sprintf('Fixed {\\itm}_\\nu^2 = 0 eV^2'),...
                    'EdgeColor',rgb('Silver'),'Location','southwest');
                obj.RunAnaObj.fixPar = fixPar_i;
                obj.RunAnaObj.pullFlag = pull_i ;
               
                if obj.range==65
                    ylim([1 6e3]);
                    xlim([2e-03 0.5]);
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
            obj.LoadGridFile('CheckSmallerN','ON');
            obj.Interp1Grid('RecomputeFlag','ON');
            pFix1 = obj.ContourPlot('CL',95,'HoldOn','ON',...
                'Color',rgb('Orange'),'LineStyle','-.','BestFit',BestFit);
            
            %% 95CL - Wilks Theorem Corrected wia MC simulation
            obj.RunAnaObj.fixPar = 'E0 Norm Bkg'; obj.RunAnaObj.InitFitPar;
            obj.RunAnaObj.pullFlag = 99;
            obj.LoadGridFile('CheckSmallerN','ON');
            obj.Interp1Grid('RecomputeFlag','ON');
            pFix2 = obj.ContourPlot('CL',96.45,'HoldOn','ON',...
                'Color',rgb('DodgerBlue'),'LineStyle',':','BestFit',BestFit);
            
            PrettyFigureFormat('FontSize',22);
                legend([pFix1,pFix2],...
                    sprintf('Fixed {\\itm}_\\nu^2 = 0 eV^2 - \\Delta\\chi^2 = 5.99'),...
                    sprintf('Fixed {\\itm}_\\nu^2 = 0 eV^2 - \\Delta\\chi^2 = 6.68'),...
                    'EdgeColor',rgb('Silver'),'Location','southwest');
                obj.RunAnaObj.fixPar = fixPar_i;
                obj.RunAnaObj.pullFlag = pull_i ;
               
                if obj.range==65
                    ylim([1 6e3]);
                    xlim([2e-03 0.5]);
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
            obj.LoadGridFile('CheckSmallerN','ON');
            obj.Interp1Grid('RecomputeFlag','ON');
            pTwin = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','OFF',...
                'Color',rgb('Orange'),'LineStyle','-.','BestFit','OFF');
            
            % data
            obj.RunAnaObj.DataType = 'Real';
            obj.LoadGridFile('CheckSmallerN','ON');
            obj.Interp1Grid('RecomputeFlag','ON');
            pData = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                'Color',rgb('DodgerBlue'),'LineStyle','-','BestFit',BestFit);
            
            titleStr = sprintf('%.0f%% C.L. , %.0feV range (%s)',obj.ConfLevel,obj.range,obj.GetPlotTitle('Mode','chi2'));
            title(titleStr,'FontWeight','normal','FontSize',get(gca,'FontSize'));
            legend([pData,pTwin],'Data','Twin','EdgeColor',rgb('Silver'),'Location','southwest');
            
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
            
            p.parse(varargin{:});
            SavePlot = p.Results.SavePlot;
            PlotStat = p.Results.PlotStat;
            PlotTot  = p.Results.PlotTot;
            chi2_i   = obj.RunAnaObj.chi2;
        
            if strcmp(obj.RunAnaObj.DataType,'Real')
                BestFit = 'ON';
            else
                BestFit = 'OFF';
            end
        
            LineWidth = 2.5;
            %% load samak
            if strcmp(PlotStat,'ON')
                obj.RunAnaObj.chi2 = 'chi2Stat';
                obj.LoadGridFile('CheckSmallerN','ON');
                obj.Interp1Grid('RecomputeFlag','ON');
                pStat = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','OFF',...
                    'Color',rgb('DodgerBlue'),'LineStyle','-','BestFit',BestFit,'PlotSplines','OFF');
            end
            
            if strcmp(PlotTot,'ON')
                obj.RunAnaObj.chi2 = 'chi2CMShape';
                obj.LoadGridFile('CheckSmallerN','ON');
                obj.Interp1Grid('RecomputeFlag','ON');
                pSys = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                    'Color',rgb('FireBrick'),'LineStyle','-','BestFit',BestFit,'PlotSplines','OFF');
            end
           %% load fitrium
           savedirF = [getenv('SamakPath'),'SterileAnalysis/GridSearchFiles/Knm1/Others/'];
           fstat = sprintf('%scontour_KSN1_Fitrium_%s_%.0feV_stat_95CL_0.txt',savedirF,obj.RunAnaObj.DataType,obj.range);
           dfStat = importdata(fstat);
           
           fsys = sprintf('%scontour_KSN1_Fitrium_%s_%.0feV_total_95CL_0.txt',savedirF,obj.RunAnaObj.DataType,obj.range);
           dfSys = importdata(fsys);
           
           if strcmp(PlotStat,'ON')
               pFStat = plot(dfStat.data(:,1),dfStat.data(:,2),'LineStyle','-.','Color',rgb('PowderBlue'),'LineWidth',LineWidth);
           end
           
           if strcmp(PlotTot,'ON')
               pFSys  = plot(dfSys.data(:,1),dfSys.data(:,2),'LineStyle','-.','Color',rgb('Orange'),'LineWidth',LineWidth);
               if obj.range==95 &&strcmp(obj.RunAnaObj.DataType,'Real')
                   fsys1 = sprintf('%scontour_KSN1_Fitrium_%s_%.0feV_total_95CL_1.txt',savedirF,obj.RunAnaObj.DataType,obj.range);
                   dfSys1 = importdata(fsys1);
                   fsys2 = sprintf('%scontour_KSN1_Fitrium_%s_%.0feV_total_95CL_2.txt',savedirF,obj.RunAnaObj.DataType,obj.range);
                   dfSys2 = importdata(fsys2);
                   plot(dfSys1.data(:,1),dfSys1.data(:,2),'LineStyle','-.','Color',rgb('Orange'),'LineWidth',LineWidth);
                   plot(dfSys2.data(:,1),dfSys2.data(:,2),'LineStyle','-.','Color',rgb('Orange'),'LineWidth',LineWidth);
               end
           end
           
           if strcmp(obj.RunAnaObj.DataType,'Real') && strcmp(BestFit,'ON')
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
                       pF_bfStat = plot(3.676e-02, 7.218e+01,'o','MarkerSize',8,'Color',pFStat.Color,'LineWidth',pFStat.LineWidth);
                   end
                   if strcmp(PlotTot,'ON')
                       pF_bfSys = plot(3.247e-02, 7.218e+01,'o','MarkerSize',8,'Color',pFSys.Color,'LineWidth',pFSys.LineWidth);
                   end
               end
           end
           
           if strcmp(PlotStat,'ON') && strcmp(PlotTot,'ON')
               legStr = {'Samak (stat. only)','Fitrium (stat. only)','Samak (stat. and syst.)','Fitrium (stat. and syst.)'};
               legend([pStat,pFStat,pSys,pFSys],legStr,'EdgeColor',rgb('Silver'),'Location','southwest');
               extraStr = '';
           elseif strcmp(PlotStat,'ON')
                legStr = {'Samak (stat. only)','Fitrium (stat. only)'};
               legend([pStat,pFStat],legStr,'EdgeColor',rgb('Silver'),'Location','southwest');
               extraStr = '_StatOnly';
           elseif strcmp(PlotTot,'ON')
                legStr = {'Samak (stat. and syst.)','Fitrium (stat. and syst.)'};
               legend([pSys,pFSys],legStr,'EdgeColor',rgb('Silver'),'Location','southwest');
               extraStr = '_Tot';
           end
           
           obj.RunAnaObj.chi2 = chi2_i;
            if obj.range==65
                xlim([4e-03 0.5])
                ylim([1 1e4])
            elseif obj.range==40
                xlim([1e-02 0.5])
                ylim([1 3e3])
            elseif obj.range==95
                xlim([3e-03 0.5])
                ylim([1 2e4]) 
            end
            
            title(sprintf('%s , %.0f eV range , %.0f%% C.L.',obj.GetPlotTitle('Mode','data'),obj.range,obj.ConfLevel),'FontWeight','normal','FontSize',get(gca,'FontSize'));
          
           %% save
           if ~strcmp(SavePlot,'OFF')
               name_i = strrep(obj.DefPlotName,sprintf('_%s',chi2_i),'');
               if strcmp(SavePlot,'ON')
                   plotname = sprintf('%s_Fitrium_%.2gCL%s.pdf',name_i,obj.ConfLevel,extraStr);
                   export_fig(gcf,plotname);
               elseif strcmp(SavePlot,'png')
                   plotname = sprintf('%s_Fitrium_%.2gCL%s.png',name_i,obj.ConfLevel,extraStr);
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
            obj.LoadGridFile('CheckSmallerN','ON');
            obj.Interp1Grid('RecomputeFlag','ON');
            pStat = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','OFF',...
                'Color',rgb('DodgerBlue'),'LineStyle','-','BestFit',BestFit,'PlotSplines','OFF');
            
            obj.RunAnaObj.chi2 = 'chi2CMShape';
            obj.LoadGridFile('CheckSmallerN','ON');
            obj.Interp1Grid('RecomputeFlag','ON');
            pSys = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                'Color',rgb('Orange'),'LineStyle','-','BestFit',BestFit,'PlotSplines','OFF');

        %% legend
        
           
               legStr = {'Stat. only','All syst. combined'};
               legend([pStat,pSys],legStr,'EdgeColor',rgb('Silver'),'Location','southwest');
               extraStr = '';
           
           obj.RunAnaObj.chi2 = chi2_i;
            if obj.range==65
                xlim([4e-03 0.5])
                ylim([1 1e4])
            elseif obj.range==40
                xlim([1e-02 0.5])
                ylim([1 3e3])
            elseif obj.range==95
                xlim([3e-03 0.5])
                ylim([1 2e4]) 
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
            
            p.parse(varargin{:});
            BestFit  = p.Results.BestFit;
            SavePlot = p.Results.SavePlot;
            Troitsk  = p.Results.Troitsk;
            Mainz    = p.Results.Mainz;
            
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
                legStr = [legStr,{sprintf('Mainz 95%% C.L.  - {\\itm}_\\nu^2 = 0 eV^2')}];
                hold on;
            end
            %% Troitsk
            if strcmp(Troitsk,'ON')
                filenameTroitsk = sprintf('%scoord_Troitsk_95CL.mat',savedirOther);
                dTroitsk = importdata(filenameTroitsk);
                pTroitsk = plot(dTroitsk.SinSquareTheta_X,dTroitsk.m4Square_Y,'--','LineWidth',1.5,...
                    'Color',rgb('DarkSlateGrey')); % Orange
                legHandle{numel(legHandle)+1} = pTroitsk;
                legStr = [legStr,{sprintf('Troitsk 95%% C.L. - {\\itm}_\\nu^2 = 0 eV^2')}];
                hold on;
            end
            %             %% 1. nuissance nu-mass without pull
            %             obj.RunAnaObj.fixPar = 'mNu E0 Norm Bkg'; obj.RunAnaObj.InitFitPar;
            %             obj.RunAnaObj.pullFlag = 99;
            %             obj.LoadGridFile('CheckSmallerN','ON');
            %             obj.Interp1Grid('RecomputeFlag','ON');
            %             pFree = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','OFF',...
            %                 'Color',rgb('ForestGreen'),'LineStyle','-','BestFit',BestFit);
            %% fixed nu-mass
            obj.RunAnaObj.fixPar = 'E0 Norm Bkg'; obj.RunAnaObj.InitFitPar;
            obj.RunAnaObj.pullFlag = 99;
            obj.LoadGridFile('CheckSmallerN','ON');
            obj.Interp1Grid('RecomputeFlag','ON');
            pFix = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                'Color',rgb('DodgerBlue'),'LineStyle','-','BestFit',BestFit);
            legHandle{numel(legHandle)+1} = pFix;
            %%  nuissance nu-mass + pull
            obj.RunAnaObj.fixPar = 'mNu E0 Norm Bkg'; obj.RunAnaObj.InitFitPar;
            obj.RunAnaObj.pullFlag = 12;
            obj.LoadGridFile('CheckSmallerN','ON');
            obj.Interp1Grid('RecomputeFlag','ON');
            pPull = obj.ContourPlot('CL',obj.ConfLevel,'HoldOn','ON',...
                'Color',rgb('PowderBlue'),'LineStyle',':','BestFit',BestFit);
            legHandle{numel(legHandle)+1} = pPull;
            
            %% appearance + legend
            legStr = [legStr,{sprintf('KATRIN KSN1 %.0f%% C.L. - {\\itm}_\\nu^2 = 0 eV^2',obj.ConfLevel),...
                sprintf('KATRIN KSN1 %.0f%% C.L. - {\\itm}_\\nu^2 free - \\sigma({\\itm}_\\nu^2) = 1.94 eV^2',obj.ConfLevel)}];
            PrettyFigureFormat('FontSize',20);
            leg =  legend([legHandle{:}],legStr{:},...
                'EdgeColor',rgb('Silver'),'Location','southwest');
            set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
            leg.FontSize = get(gca,'FontSize')-2;
            legend boxoff
            if obj.range==65
                ylim([1 6e3]);
                xlim([2e-03 0.5]);
            elseif obj.range==95
                ylim([1 1e4]);
                xlim([9e-04 0.5]);
            elseif obj.range==40
                ylim([1 1.6e3]);
                xlim([6e-03 0.5]);
            end
            title('');%remove title
            %grid on
            %% save
            if ~strcmp(SavePlot,'OFF')
                name_i = strrep(obj.DefPlotName,'_mNuE0BkgNorm','');
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
            PlotStatDom = 'ON';
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
            ylim([0 6.5])
            PrettyFigureFormat;
            leg = legend([pref,pl{:}],legStr{:},'EdgeColor',rgb('Silver'),'Location','northwest');
            
            
        end
        
        function titleStr = GetPlotTitle(obj,varargin)
            p=inputParser;
            p.addParameter('Mode','all',@(x)ismember(x,{'all','chi2','data'}));
            p.parse(varargin{:});
            Mode = p.Results.Mode;
            
            % some useful default title
            if strcmp(obj.RunAnaObj.DataType,'Real')
                DataStr = 'Data';
            else
                 DataStr = 'Twin';
            end
            
            if strcmp(obj.RunAnaObj.chi2,'chi2Stat')
                chi2Str = 'stat. only';
            else
                chi2Str = 'stat. and syst.';
            end
            
            switch Mode
                case 'all'
            titleStr = sprintf('%s , %.0f eV range (%s)',DataStr,obj.range,chi2Str);
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
    
    % Data stream: labels, loading, saving
    methods
        function f = LoadGridFile(obj,varargin)
            p = inputParser;
            p.addParameter('CheckLargerN','ON',@(x)ismember(x,{'ON','OFF'})); % if specified ngrids do not exist - look also for larger n grids
            p.addParameter('CheckSmallerN','OFF',@(x)ismember(x,{'ON','OFF'})); % if specified ngrids  +larger do not exist - look also for smaller n grids
            p.parse(varargin{:});
            CheckLargerN  = p.Results.CheckLargerN;
            CheckSmallerN = p.Results.CheckSmallerN;
            
            filename = obj.GridFilename;
            
            loadSuccess = 0;
            
            if exist(filename,'file')
                f = importdata(filename);
                fprintf('load grid from file %s \n',filename)
                loadSuccess = 1;
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
                fprintf('Cannot find grid %s \n',filename);
                f = 0;
            else
                obj.mNu4Sq = f.mnu4Sq;
                obj.sin2T4 = f.sin2T4;
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
        function filename = GridFilename(obj)
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
            
                 if isfloat(obj.RandMC) && strcmp(obj.RunAnaObj.DataType,'Twin')
                     extraStr = sprintf('%s_RandMC%.0f',extraStr,obj.RandMC);
                     savedir = strrep(savedir,'Twin/','TwinRandomizedMC/');
                 end
                 
                 if ~strcmp(obj.RunAnaObj.ELossFlag,'KatrinT2')
                     extraStr = [extraStr,sprintf('_%s',obj.RunAnaObj.ELossFlag)];
                 end
                 
                 if ~strcmp(obj.RunAnaObj.AngularTFFlag,'OFF')
                     extraStr = [extraStr,'_AngTF'];
                 end
                 
                 if obj.RunAnaObj.pullFlag<=12
                     extraStr = sprintf('%s_pull%.0f',extraStr,obj.RunAnaObj.pullFlag);
                 end
                 
                 
                 MakeDir(savedir);

            % get runlist-name
            RunList = extractBefore(obj.RunAnaObj.RunData.RunName,'_E0');
            if isempty(RunList)
                RunList = obj.RunAnaObj.RunData.RunName;  
            end
            
            freeParStr =  ConvertFixPar('freePar',obj.RunAnaObj.fixPar,'Mode','Reverse');
            filename = sprintf('%sKSN1_GridSearch_%s_%s_%s_%.0feVrange_%s_%.0fnGrid%s.mat',...
                savedir,RunList,obj.RunAnaObj.DataType,strrep(freeParStr,' ',''),...
                obj.range,obj.RunAnaObj.chi2,obj.nGridSteps,extraStr);
  
        end  
        function plotname = DefPlotName(obj)
            % generic plot name
            filename = obj.GridFilename;
            filename = extractBetween(filename,'KSN1_','.mat');
            plotdir = sprintf('%sSterileAnalysis/plots/%s/%s/',...
                     getenv('SamakPath'),obj.RunAnaObj.DataSet,obj.RunAnaObj.DataType);
            MakeDir(plotdir);
            if isfloat(obj.RandMC) && strcmp(obj.RunAnaObj.DataType,'Twin')
                plotdir = strrep(savedir,'Twin/','TwinRandomizedMC/'); 
            end
            
            
             plotname = sprintf('%s%s',plotdir,filename{:});
        end
    end
    
    % small auxillary methods
    methods
        function SetNPfactor(obj)
            if ~strcmp(obj.RunAnaObj.chi2,'chi2Stat')
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
            end 
        end
    end
    
end

