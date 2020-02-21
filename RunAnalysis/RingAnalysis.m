classdef RingAnalysis < handle
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
        function obj = RingAnalysis(varargin) % constructor
            p=inputParser;
            p.addParameter('RunAnaObj','', @(x) isa(x,'RunAnalysis') || isa(x,'MultiRunAnalysis'));
            p.addParameter('RingList',1:12,@(x)isfloat(x));
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
            CommonArg = {'DataType',obj.RunAnaObj.DataType,'AnaFlag','StackPixel','exclDataStart',obj.RunAnaObj.exclDataStart,...
                'PullFlag',obj.RunAnaObj.pullFlag,'RingMerge',obj.RunAnaObj.RingMerge,...
                'fixPar',obj.RunAnaObj.fixPar,'i_Q',obj.RunAnaObj.i_Q,'fitter',obj.RunAnaObj.fitter,...
                'minuitOpt',obj.RunAnaObj.minuitOpt,'FSDFlag',obj.RunAnaObj.FSDFlag,...
                'NonPoissonScaleFactor',obj.RunAnaObj.NonPoissonScaleFactor,'ELossFlag',obj.RunAnaObj.ELossFlag,...
                'FakeInitFile',obj.RunAnaObj.FakeInitFile,'ROIFlag',obj.RunAnaObj.ROIFlag,...
                'MosCorrFlag',obj.RunAnaObj.MosCorrFlag};
            %             if ~isempty(obj.RunAnaObj.RunNr) %RunAnalysis Object
            %                 obj.MultiObj = arrayfun(@(x) RunAnalysis(CommonArg{:},'RunNr',obj.RunAnaObj.RunNr,...
            %                     'RingList',x),obj.RingList);
            %
            %             else %MultiRunAnalysis Object
            %                 obj.MultiObj = arrayfun(@(x) MultiRunAnalysis(CommonArg{:},'RunList',obj.RunAnaObj.StackFileName,...
            %                     'RingList',x),obj.RingList);
            %
            %             end
            if ~isempty(obj.RunAnaObj.RunNr) %RunAnalysis Object
                obj.MultiObj = arrayfun(@(x) RunAnalysis(CommonArg{:},'RunNr',obj.RunAnaObj.RunNr,...
                    'PixList',cell2mat(x)),obj.RunAnaObj.RingPixList);
            else %MultiRunAnalysis Object
                obj.MultiObj = arrayfun(@(x) MultiRunAnalysis(CommonArg{:},'RunList',obj.RunAnaObj.StackFileName,...
                    'PixList',cell2mat(x)),(obj.RunAnaObj.RingPixList));
            end
        end
        
        
        function FitRings(obj,varargin)
            p=inputParser;
            p.addParameter('SaveResult','OFF',@(x)ismember(x,{'ON','OFF'}))
            p.addParameter('RecomputeFlag','ON',@(x)ismember(x,{'ON','OFF'}))
            p.addParameter('AsymErr','OFF',@(x)ismember(x,{'ON','OFF'}))
            p.parse(varargin{:});
            SaveResult    = p.Results.SaveResult;
            RecomputeFlag = p.Results.RecomputeFlag;
            AsymErr       = p.Results.AsymErr;
            %% label
            range = round(abs(obj.RunAnaObj.ModelObj.qU(obj.RunAnaObj.exclDataStart)-obj.RunAnaObj.ModelObj.Q_i),0);
            savedir = [getenv('SamakPath'),sprintf('tritium-data/fit/%s/SingleRingFit/',obj.RunAnaObj.DataSet)];
            MakeDir(savedir);
            freePar = ConvertFixPar('freePar',obj.RunAnaObj.fixPar,'Mode','Reverse'); %convert fixed parameters (numbers) back to readable string with free parameters
            switch obj.RunAnaObj.ROIFlag
                case 'Default'
                    RoiStr = '';
                case '14keV'
                    RoiStr = '_14keVROI';
            end
            if strcmp(obj.RunAnaObj.MosCorrFlag,'ON')
                MosStr = '_MosCorr';
            else
                MosStr = '';
            end
            savename = [savedir,sprintf('SingleRingFitResult_%s_%s_%.0fruns_fix%s_%s_%.0feVrange%s.mat',...
                obj.RunAnaObj.RingMerge,obj.RunAnaObj.RunData.RunName,...
                numel(obj.RunAnaObj.RunList),freePar,obj.RunAnaObj.chi2,range,RoiStr,MosStr)];
            if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
                obj.FitResult = importdata(savename);
            else
                par         = zeros(obj.nRings,obj.RunAnaObj.nPar);
                err         =  zeros(obj.nRings,obj.RunAnaObj.nPar);
                errNeg      =  zeros(obj.nRings,obj.RunAnaObj.nPar);
                errPos      =  zeros(obj.nRings,obj.RunAnaObj.nPar);
                chi2min     =  zeros(obj.nRings,1);
                dof         = zeros(obj.nRings,1);
                ScanResults = cell(obj.nRings,1); %only filled when AsymErr is ON
                
                progressbar('Fitting single rings');
                for i=1:obj.nRings
                    progressbar(i/obj.nRings);
                    obj.MultiObj(i).Fit;
                    %obj.MultiObj(i).Fit;
                    par(i,:) = obj.MultiObj(i).FitResult.par;
                    err(i,:) = obj.MultiObj(i).FitResult.err;
                    chi2min(i) = obj.MultiObj(i).FitResult.chi2min;                 
                    dof(i) = obj.MultiObj(i).FitResult.dof;
                    
                    if isfield(obj.MultiObj(i).FitResult,'errNeg')
                       errNeg(i,1:4) = obj.MultiObj(i).FitResult.errNeg;
                       errPos(i,1:4) = obj.MultiObj(i).FitResult.errPos;
                    end
                    
                    if strcmp(AsymErr,'ON') % asymmetric uncertainties on neutrino mass
                        if ~contains(obj.RunAnaObj.fixPar,'fix 1 ')% if error calculation failed err(i,1)<0.1 &&
                            ScanResults{i} =  obj.MultiObj(i).GetAsymFitError('Mode','Smart','Parameter','mNu','ParScanMax',5);
                            err(i,1) = (abs(ScanResults{i}.AsymErr(1))+abs(ScanResults{i}.AsymErr(2)))/2;
                        end
                    else
                        ScanResults{i} = 0;
                    end
                end
                
                obj.FitResult = struct('par',par,'err',err,'chi2min',chi2min,'dof',dof,...
                    'ScanResults',ScanResults,'errNeg',errNeg,'errPos',errPos);
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
            p.addParameter('Blind','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('YLim','',@(x)isfloat(x));
            p.parse(varargin{:});
            PlotPar      = p.Results.PlotPar;
            SavePlot     = p.Results.SavePlot;
            linFitFlag   = p.Results.linFit;
            PlotMode     = p.Results.PlotMode;
            Blind        = p.Results.Blind;
            YLim         = p.Results.YLim;
          
            y      = obj.FitResult.par;%(:,PlotPar);
            y      = y(:,PlotPar);
            yErr   = obj.FitResult.err;%(:,PlotPar);
            yErr   = yErr(:,PlotPar);

            %% linear fit
            if strcmp(linFitFlag,'ON')
                [linFitpar, linFiterr, linFitchi2min,linFitdof] =...
                    linFit(obj.RingList',y,yErr);
            end
            %% label
            range = round(abs(obj.RunAnaObj.ModelObj.qU(obj.RunAnaObj.exclDataStart)-obj.RunAnaObj.ModelObj.Q_i),0);
            %% plot fit result
            fig2 = figure('Renderer','painters');
            set(fig2,'units','normalized','pos',[0.1, 0.1,0.6,0.5]);  
            plot(linspace(0.5,12.5,obj.nRings),zeros(obj.nRings,1),'-','Color',rgb('SlateGrey'),'LineWidth',2);
            hold on;
            %hold on;
            if PlotPar==3
                ScaleFactor  = cell2mat(arrayfun(@(x) numel(x.PixList),obj.MultiObj,'UniformOutput',0))'; %number of pixels per ring
            else
                ScaleFactor = ones(numel(obj.RingList),1);
            end
            if strcmp(PlotMode,'Rel')
                if PlotPar == 3
                    BKG_i = cell2mat(arrayfun(@(x) x.ModelObj.BKG_RateSec_i,obj.MultiObj,'UniformOutput',0))';
                    %BKG_Tot = sum(BKG_i);
                else
                    meanPar =wmean(y,1./yErr.^2);
                end
                
            elseif strcmp(PlotMode,'Abs')
                if PlotPar==3
                    BKG_i = cell2mat(arrayfun(@(x) x.ModelObj.BKG_RateSec_i,obj.MultiObj,'UniformOutput',0))';
                    meanPar = -cell2mat(arrayfun(@(x) x.ModelObj.BKG_RateSec_i,obj.MultiObj,'UniformOutput',0))';
                else
                    % meanPar = zeros(obj.nRings,1); 
                    meanPar = 0;
                end
                
            end
            if PlotPar == 3            
                e = errorbar(obj.RingList,1e3.*y+BKG_i./ScaleFactor,1e3.*yErr,'o',...
                'LineWidth',2,'LineStyle','none','MarkerSize',8,'Color',obj.RunAnaObj.PlotColor,'MarkerFaceColor',obj.RunAnaObj.PlotColor);
            BKG_Tot = sum(y+BKG_i);
            else
                e = errorbar(obj.RingList,(y-meanPar)./ScaleFactor,yErr,'o',...
                'LineWidth',2,'LineStyle','none','MarkerSize',8,'Color',obj.RunAnaObj.PlotColor,'MarkerFaceColor',obj.RunAnaObj.PlotColor);
            end
            e.CapSize = 0;
             PrettyFigureFormat('FontSize',24);
           
            if strcmp(linFitFlag,'ON') && PlotPar ~= 3 

                runname = strrep(obj.RunAnaObj.RunData.RunName,'_',' ');
                if contains(runname,'afterFix')
                    runname = 'KNM2 after RW fix I';
                elseif contains(runname,'RW3')
                    runname = 'KNM2 after RW fix II';
                end
                l = plot(obj.RingList,(linFitpar(1).*obj.RingList+linFitpar(2)-meanPar)./ScaleFactor','-','Color',obj.RunAnaObj.PlotColor,'LineWidth',2);
                if linFitpar(1)<0.03
                    if PlotPar==1
                        linFitleg =  sprintf('Linear fit %.1f \\pm %.1f meV^2/ring , p-value = %.2f',linFitpar(1)*1e3,linFiterr(1)*1e3,1-chi2cdf(linFitchi2min,linFitdof));
                    elseif PlotPar==2
                        linFitleg =  sprintf('Linear fit %.1f \\pm %.1f meV/ring , p-value = %.2f ',linFitpar(1)*1e3,linFiterr(1)*1e3,1-chi2cdf(linFitchi2min,linFitdof));
                    end
                else
                    if PlotPar==1
                        linFitleg =  sprintf('Linear fit %.2f \\pm %.2f eV^2/ring , p-value = %.2f ',linFitpar(1),linFiterr(1),1-chi2cdf(linFitchi2min,linFitdof));
                  
                    elseif PlotPar==2
                        linFitleg =  sprintf('Linear fit %.2f \\pm %.2f eV/ring , p-value = %.2f',linFitpar(1),linFiterr(1),1-chi2cdf(linFitchi2min,linFitdof));
                    end
                end
                leg = legend([e,l(1)],sprintf('%s , %.0f eV range',runname,range),linFitleg);
            else
                leg = legend(e,sprintf('%s: %.0f eV range',strrep(obj.RunAnaObj.RunData.RunName,'_',' '),range));
            end
            hold off;
            leg.FontSize= get(gca,'FontSize');
            leg.EdgeColor = rgb('Silver');
            if PlotPar==1
                ylabel(sprintf('{\\itm}_\\nu^2 - \\langle{\\itm}_\\nu^2\\rangle (eV^2)'));
            elseif PlotPar==2
                ylabel(sprintf('{\\itE}_0^{fit} - \\langle{\\itE}_0^{fit}\\rangle (eV)'));
            elseif PlotPar==3
                ylabel(sprintf('({\\itB} - \\langle{\\itB}\\rangle) / npixel   (mcps)'));
            elseif  PlotPar==11
                ylabel(sprintf('qU_{offset} - \\langleqU_{offset}\\rangle (V)'));
            elseif  PlotPar==4
                ylabel(sprintf('{\\itN} - \\langle{\\itN}\\rangle '));
            end
            xlabel('Rings')
            ax = gca;
            ax.YAxis.Exponent = 0;
            meanPar =wmean(y,1./yErr.^2);
            if PlotPar==1
                AnaType = 'm2';
                t =title(sprintf('%.0f stacked runs (%.0f - %.0f) - %.0f pixels - \\langle{\\itm}^2\\rangle = %.3f eV^2',numel(obj.RunAnaObj.StackedRuns),obj.RunAnaObj.RunList(1),obj.RunAnaObj.RunList(end),numel(obj.RunAnaObj.PixList), meanPar));
            elseif PlotPar==2
               t = title(sprintf('%.0f stacked runs (%.0f - %.0f) - %.0f pixels - \\langle{\\itE}_0^{fit}\\rangle = %.2f eV',numel(obj.RunAnaObj.StackedRuns),obj.RunAnaObj.RunList(1),obj.RunAnaObj.RunList(end),numel(obj.RunAnaObj.PixList), meanPar + obj.RunAnaObj.ModelObj.Q_i));
                AnaType = 'E0';
            elseif PlotPar==3
               t= title(sprintf('%.0f stacked runs (%.0f - %.0f) - %.0f pixels - B_{Tot} = %.1f mcps',numel(obj.RunAnaObj.StackedRuns),obj.RunAnaObj.RunList(1),obj.RunAnaObj.RunList(end),numel(obj.RunAnaObj.PixList), 1e3*BKG_Tot));
                AnaType = 'BKG';
            elseif PlotPar==4
                AnaType = 'Normalization';
                t =title(sprintf('%.0f stacked runs (%.0f - %.0f) - %.0f pixels - <N> = %.3f',numel(obj.RunAnaObj.StackedRuns),obj.RunAnaObj.RunList(1),obj.RunAnaObj.RunList(end),numel(obj.RunAnaObj.PixList), meanPar));
            elseif PlotPar==11
                AnaType = 'qUoffset';
                t =title(sprintf('%.0f stacked runs (%.0f - %.0f) - %.0f pixels - <qU_{offset}> = %.3f',numel(obj.RunAnaObj.StackedRuns),obj.RunAnaObj.RunList(1),obj.RunAnaObj.RunList(end),numel(obj.RunAnaObj.PixList), meanPar));
                
            end
            
            if strcmp(Blind,'ON')
                t = title(sprintf('%.0f stacked runs (%.0f - %.0f) - %.0f pixels ',...
                    numel(obj.RunAnaObj.StackedRuns),obj.RunAnaObj.RunList(1),obj.RunAnaObj.RunList(end),numel(obj.RunAnaObj.PixList)));
            end
             t.FontWeight = 'normal';
             t.FontSize = get(gca,'FontSize');
            switch obj.RunAnaObj.RingMerge
                case 'InnerPseudo3'
                    xticks([1]);
                    xticklabels({'1,2,3,4,5,6,7,8,9'})
                case 'Full'
                    xticks([1 2 3 4]);
                    xticklabels({'1-3','4-6','7-9','10-12'})
                case 'Default'
                    xticks([1:1:10]);
                    xticklabels({'1-2',3,4,5,6,7,8,9,10,'11-12'})
                case 'None'
                    xticks([1:12])
            end
            set(gca,'XMinorTick','off');
            ymin = min((y-meanPar)./ScaleFactor-yErr);
            ymax = max((y-meanPar)./ScaleFactor+yErr);
            if ~isempty(YLim)
                ylim([min(YLim),max(YLim)]);
            elseif strcmp(PlotMode,'Rel')
                ylim([1.1*ymin,1.7*ymax]);
            else
                ylim([0.5*ymin,1.5*ymax]);
            end
            
            
            xlim([min(obj.RingList)-0.2,max(obj.RingList)+0.2])
            leg.Location = 'northwest';
            
            
            if strcmp(SavePlot,'ON')
                savedir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/SingleRingFit/',obj.RunAnaObj.DataSet)];
                MakeDir(savedir);
                switch obj.RunAnaObj.ROIFlag
                    case 'Default'
                        RoiStr = '';
                    case '14keV'
                        RoiStr = '_14keVROI';
                end
                if strcmp(obj.RunAnaObj.MosCorrFlag,'ON')
                    MosStr = '_MosCorr';
                else
                    MosStr = '';
                end
                freePar = ConvertFixPar('freePar',obj.RunAnaObj.fixPar,'Mode','Reverse');
                savename = sprintf('%s_ringwise_%sRuns_%.0frange_%s_Merge_%s_%s_Blind%s%s%s.pdf',...
                    AnaType,obj.RunAnaObj.ModelObj.TD,range,obj.RunAnaObj.RingMerge,freePar,obj.RunAnaObj.DataType,Blind,RoiStr,MosStr);
                export_fig(gcf,[savedir,savename],'-painters');
                fprintf('Save plot to %s \n',[savedir,savename])
            end
%             N = 10000;
%             e = repmat([obj.FitResult.err(:,2)]',N,1);
%             a = randn(N,numel(obj.FitResult.err(:,2)));
%             d = e.*a;
%             stdsim = std(d,0,2);
%             
%             fighist = figure;
%             h = histogram(stdsim,'Normalization','probability');
%             
%             dataknm1 = [(obj.FitResult.par(:,PlotPar)-meanPar)./ScaleFactor] ;
%             stddata = std(dataknm1);
%             hold on
%             line([stddata stddata], [0 max(h.Values)],'Color','Red');
%             hold off
%             PrettyFigureFormat
         
        end
        function AverageFitResults(obj,varargin)
            % Constant Fit of Fit ParameterAvPar
            
            p      = inputParser;
            p.addParameter('AvPar',1,@(x)isfloat(x));
            p.parse(varargin{:});
            AvPar  = p.Results.AvPar;
            
            disp([ obj.FitResult.par(:,1) obj.FitResult.err(:,1)]);
            %[ round(18573.7+obj.FitResult.par(:,2),3) round(obj.FitResult.err(:,2),3)]
            
             %% Cte fit
            tmparg = sprintf(['set pri -10;  fix 1 ; min ; minos ;'],'');
            Data = [obj.RingList',obj.FitResult.par(:,AvPar),obj.FitResult.err(:,AvPar)];
            parInit = [0 0];
            Args   = {parInit, Data, '-c', tmparg};
            [par, err, chi2min, ~] = fminuit('chi2lin',Args{:});
            dof = numel(obj.FitResult.par(:,AvPar))-numel(parInit)+1;     
            
                        
            fig2 = figure('Renderer','opengl');
            set(fig2,'units','normalized','pos',[0.1, 0.1,0.9,0.7]);  
            
            e = errorbar(obj.RingList,(obj.FitResult.par(:,AvPar)),obj.FitResult.err(:,AvPar),'o',...
                'LineWidth',2,'LineStyle','none','MarkerSize',8,'Color',obj.RunAnaObj.PlotColor,'MarkerFaceColor',obj.RunAnaObj.PlotColor);
            
            hold on
            boundedline(obj.RingList,obj.RingList./obj.RingList.*par(2),obj.RingList./obj.RingList.*err(2),'alpha','cmap','transparency', 0.2,rgb('CadetBlue'));
            hold off
            
            switch AvPar
                case 1
                    mytitle=sprintf('m^2 = %.3f +- %.3f  eV^2 (%.2f/%.0f dof) \n ',par(2),err(2),chi2min,dof);
                    ylabel('m^2 (eV^2)');
                case 2
                    mytitle=sprintf('E_0 = %.3f +- %.3f eV  (%.2f/%.0f dof) \n ',obj.MultiObj(1).ModelObj.Q_i+par(2),err(2),chi2min,dof);
                    ylabel('E_0 (eV)');
            end
            
            switch obj.RunAnaObj.RingMerge
                case 'InnerPseudo3'
                    xticks([1]);
                    xticklabels({'Rings 1,2,3,4,5,6,7,8,9'})
                case 'Full'
                    xticks([1 2 3 4]);
                    xticklabels({'Rings 1,2,3','Rings 4,5,6','Rings 7,8,9','Rings 10,11,12'})
                case 'Default'
                    xticks([1:1:10]);
                    xticklabels({'Rings 1,2',3,4,5,6,7,8,9,10,'11,12'})
                case 'None'
                    xticks([1:12])
            end
            xlim([min(obj.RingList)-0.5,max(obj.RingList)+0.5]);
            title(mytitle);

            PrettyFigureFormat;
            
        end
    end
end