classdef RingAnalysis2 < handle
    %% 
    % class to facilitate ring fits using the MultiRunAnalysis/RunAnalysis class
    % FPD segmentation is off, ringwise fit achieved with PixList
    %
    % neccessary input:
    % 1. RunAnaObj --> 1 (Multi)RunAnalysis Object
    % 2. RingList  --> Vector e.g. 1:12 (13 doesnt work because of zero entries)
    %
    % L. Schlueter - April 2019
    %%
    properties (Access=public)
        RunAnaObj; % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        MultiObj;  % n*RunAnaObj, 1 for each ring
        nRings;    
        RingList;
        FitResult;
    end
    methods % constructor
        function obj = RingAnalysis2(varargin) % constructor
            p=inputParser;
            p.addParameter('RunAnaObj','', @(x) isa(x,'RunAnalysis') || isa(x,'MultiRunAnalysis'));
            p.addParameter('RingList',1:10,@(x)isfloat(x));
            p.parse(varargin{:});
            obj.RunAnaObj = p.Results.RunAnaObj;
            obj.RingList  = p.Results.RingList;
            
            obj.nRings = numel(obj.RingList);
            
            GetSamakPath; %sets current samak path as enviromental variable
            obj.InitMultiObj;
        end
    end
    methods
        function InitMultiObj(obj)
            % Initialize nRing * MultiRunAnalysis/RunAnalysis Objects
            CommonArg = {'DataType','Real','AnaFlag','StackPixel','exclDataStart',obj.RunAnaObj.exclDataStart,...
                'exclDataStart',obj.RunAnaObj.exclDataStart,'PullFlag',obj.RunAnaObj.pullFlag};
            
            if ~isempty(obj.RunAnaObj.RunNr) %RunAnalysis Object
                obj.MultiObj = arrayfun(@(x) RunAnalysis(CommonArg{:},'RunNr',obj.RunAnaObj.RunNr,...
                    'RingList',x),obj.RingList);
            else %MultiRunAnalysis Object
                obj.MultiObj = arrayfun(@(x) MultiRunAnalysis(CommonArg{:},'RunList',obj.RunAnaObj.StackFileName,...
                    'RingList',x),obj.RingList);
            end
        end
        function FitRings(obj,varargin)
            p=inputParser;
            p.addParameter('SaveResult','OFF',@(x)ismember(x,{'ON','OFF'}))
            p.addParameter('RecomputeFlag','ON',@(x)ismember(x,{'ON','OFF'}))
            p.parse(varargin{:});
            SaveResult = p.Results.SaveResult;
            RecomputeFlag = p.Results.RecomputeFlag;
                        %% label
            range = round(abs(obj.RunAnaObj.ModelObj.qU(obj.RunAnaObj.exclDataStart)-obj.RunAnaObj.ModelObj.Q_i),0);
            savename = sprintf('./results/FitRingResult_%s_%.0fruns_%.0feVrange.mat',obj.RunAnaObj.ModelObj.TD,numel(obj.RunAnaObj.RunList),range);
            if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
                obj.FitResult = importdata(savename);
            else
                par     = zeros(obj.nRings,obj.RunAnaObj.nPar);
                err     =  zeros(obj.nRings,obj.RunAnaObj.nPar);
                chi2min =  zeros(obj.nRings,1);
                dof     = zeros(obj.nRings,1);
                
                parfor i=1:obj.nRings
                    obj.MultiObj(i).Fit;
                    par(i,:) = obj.MultiObj(i).FitResult.par;
                    err(i,:) = obj.MultiObj(i).FitResult.err;
                    chi2min(i) = obj.MultiObj(i).FitResult.chi2min;
                    dof(i) = obj.MultiObj(i).FitResult.dof;
                end
                
                obj.FitResult = struct('par',par,'err',err,'chi2min',chi2min,'dof',dof);
            end
            if strcmp(SaveResult,'ON') 
                FitResult = obj.FitResult;
                save(savename,'FitResult');
                
            end
        end
        function PlotFits(obj,varargin)
            p=inputParser;
            p.addParameter('PlotPar',2,@(x)isfloat(x));
            p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('linFit','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('PlotMode','Rel',@(x)ismember(x,{'Rel','Abs'}));
            p.parse(varargin{:});
            PlotPar  = p.Results.PlotPar;
            SavePlot = p.Results.SavePlot;
            linFitFlag   = p.Results.linFit;
            PlotMode   = p.Results.PlotMode;
            %% linear fit
            if strcmp(linFitFlag,'ON')
                [linFitpar, linFiterr, linFitchi2min,linFitdof] =...
                    linFit(obj.RingList',obj.FitResult.par(:,PlotPar),obj.FitResult.err(:,PlotPar));
            end
            %% label
            range = round(abs(obj.RunAnaObj.ModelObj.qU(obj.RunAnaObj.exclDataStart)-obj.RunAnaObj.ModelObj.Q_i),0);
            %% plot fit result
            fig2 = figure('Renderer','opengl');
            set(fig2,'units','normalized','pos',[0.1, 0.1,0.7,0.7]);  
            if PlotPar==3
                ScaleFactor  = cell2mat(arrayfun(@(x) numel(x.PixList),obj.MultiObj,'UniformOutput',0))'; %number of pixels per ring
            else
                ScaleFactor = ones(numel(obj.RingList),1);
            end
            if strcmp(PlotMode,'Rel')
                 meanPar =wmean(obj.FitResult.par(:,PlotPar),1./obj.FitResult.err(:,PlotPar).^2);
            elseif strcmp(PlotMode,'Abs')
                if PlotPar==3
                    meanPar = -cell2mat(arrayfun(@(x) x.ModelObj.BKG_RateSec_i,obj.MultiObj,'UniformOutput',0))';
                else
                    meanPar = zeros(obj.nRings,1);
                end
                
            end
            
            e = errorbar(obj.RingList,(obj.FitResult.par(:,PlotPar)-meanPar)./ScaleFactor,obj.FitResult.err(:,PlotPar),'o',...
                'LineWidth',2,'LineStyle','none','MarkerSize',8,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('SlateGray'));
            
            hold on;
            plot(linspace(0.5,12.5,obj.nRings),zeros(obj.nRings,1),'--','Color',rgb('SlateGray'),'LineWidth',2);
            if strcmp(linFitFlag,'ON')
                l = plot(obj.RingList,(linFitpar(1).*obj.RingList+linFitpar(2)-meanPar)./ScaleFactor','-','Color',rgb('SkyBlue'),'LineWidth',3);
                linFitleg =  sprintf('linear fit slope: (%.1f \\pm %.1f) meV @ \\chi2 = %.1f / %.0f dof',linFitpar(1)*1e3,linFiterr(1)*1e3,linFitchi2min,linFitdof);
                leg = legend([e,l],sprintf('%.0f eV range',range),linFitleg);
            else
                leg = legend([e],sprintf('%.0f eV range',range));
            end
            hold off;
            PrettyFigureFormat;
            if PlotPar==1
                ylabel('m^2_\nu - <m^2_\nu> (eV^2)');
            elseif PlotPar==2
                ylabel('E_0 - <E_0> (eV)');
            elseif PlotPar==3
                ylabel('(B - <B>) / npixel   (cps)');
            elseif  PlotPar==4
                ylabel('N - <N> ');
            end
            xlabel('ring')
            set(gca,'FontSize',20);
            title(sprintf('%.0f stacked runs (%.0f - %.0f) - %.0f pixels total',numel(obj.RunAnaObj.StackedRuns),obj.RunAnaObj.RunList(1),obj.RunAnaObj.RunList(end),numel(obj.RunAnaObj.PixList)));
            legend boxoff
            xlim([min(obj.RingList)-0.5,max(obj.RingList)+0.5])
            leg.Location = 'northwest';
            %ylim([-0.6,0.6]);
            if strcmp(SavePlot,'ON')
                if ~exist('./plots/','dir')
                    system('mkdir ./plots/');
                end
                
                savename = sprintf('./plots/E0ringwise_%sRuns_%.0frange.png',obj.RunAnaObj.ModelObj.TD,range);
                print(savename,'-dpng','-r500');
            end
        end
    end
end