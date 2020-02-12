% ----------------------------------------------------------------------- %
% Plot class for the KATRIN experiment in SAMAK
% ----------------------------------------------------------------------- %
% This class contains options to plot the data produced
% by the KATRIN experiment, with its respective analysis.
%
% It includes:
% - Function for plot of the data-spectra comparison, and residuals
% - Function to automatically create ring fit summary
% - Function to automatically create ring rho d summary
%
% Last update by P. Morales: 06/2018
%
% T. Lasserre
% CEA - Saclay
%
% L. Schlueter
% TUM / MPP
%
% P. I. Morales Guzman
% TUM / MPP
%
% ----------------------------------------------------------------------- %

classdef PLOTC < handle
    
    properties (Constant = true, Hidden = true)
        
    end
    
    properties (Access=protected)
        
    end
    
    properties (Dependent=true,Access=public)
        
    end
    
    properties (Access=public)
        
        
        Xdata; % Data to plot in the Xaxis
        Ydata; % Data to plot in Y axis
        CovMat;  % Covariance Matrix
        
        ModelObj; % Model Object
        FitResult; % Results of fit performed
        startqU; % qU value in which the plotting on the data should start
        RunList; % List of Runs in Analysis
        
        figN; % number of the figure plotted
        plotTitle; % Title of the plot, given directly or selected from a list
        titleFlag; % Flag for predefined selected titles
        
        RingList;
        PixelList;
        StackFileName;
        
        saveplot; % save plot in PDF file
        savename; % name of the file where the plot is stored
        
    end
    
    methods % constructor
        function obj = PLOTC(varargin)
            % BEGIN Parser
            p = inputParser;
            
            p.addParameter('Xdata',[],@(x) ismatrix(x));
            p.addParameter('Ydata',[],@(x) ismatrix(x));
            p.addParameter('CovMat',[],@(x) ismatrix(x));
            
            p.addParameter('ModelObj',[],@(x) isa((x),'TBD'));
            p.addParameter('FitResult',[],@(x) isstruct(x));
            p.addParameter('startqU',1,@(x) isfloat(x) && x>0);
            p.addParameter('RunList',[],@(x) isvector(x) && all(x>0));
            
            p.addParameter('figN',1,@(x) isfloat(x) && x>0);
            p.addParameter('plotTitle',[],@(x) ischar(x));
            p.addParameter('titleFlag','Standard',@(x) ismember(x,{'Standard','Stack','RhoDScanRing','RhoDScan'}));
            
            % multipixel/multiring plots
            p.addParameter('RingList',1:13,@(x) isfloat(x) && all(x>0))
            p.addParameter('PixelList',1:148,@(x) isfloat(x) && LL(x>0))
            
            
            % save options
            p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF','export'}));
            p.addParameter('savename','SamakPlot',@(x) ischar(x));
            
            p.parse(varargin{:});
            obj.Xdata            = p.Results.Xdata;
            obj.Ydata            = p.Results.Ydata;
            obj.CovMat           = p.Results.CovMat;
            
            obj.ModelObj         = p.Results.ModelObj;
            obj.FitResult        = p.Results.FitResult;
            obj.startqU          = p.Results.startqU;
            obj.RunList          = p.Results.RunList;
            
            obj.figN             = p.Results.figN;
            obj.plotTitle        = p.Results.plotTitle;
            obj.titleFlag        = p.Results.titleFlag;
            
            obj.RingList         = p.Results.RingList;
            obj.PixelList        = p.Results.PixelList;
            
            obj.saveplot         = p.Results.saveplot;
            obj.savename         = p.Results.savename;
            
            % END Parser
            
            
            SelectTitle(obj);
        end % constructor
        
        
        
    end % methods
    
    methods % Initialize
        function SelectTitle(obj)
            if isempty(obj.plotTitle)
                switch obj.titleFlag
                    case 'Standard'
                        obj.plotTitle = ...
                            'Cool Samak Fit';
                    case 'Stack'
                        obj.plotTitle = ...
                            sprintf('Stacked Data Run %u-%u - Samak Fit \n(18575 eV - qUmin) = %.0f eV',...
                            obj.RunList(1),obj.RunList(end),obj.ModelObj.Q_i-obj.Xdata(obj.startqU));
                    case 'RhoDScanRing'
                        obj.plotTitle = ...
                            sprintf('Stacked Data Run %u-%u ring-wise \\rho d scan',...
                            obj.RunList(1),obj.RunList(end));
                    case 'RhoDScan'
                        obj.plotTitle = ...
                            sprintf('Stacked Data Run %u-%u \\rho d scan \n(18575 eV - qUmin) = %.0f eV',...
                            obj.RunList(1),obj.RunList(end),18575-obj.Xdata(obj.startqU));
                end
            end
        end
        
    end
    
    methods % Single/Stacked-pixel Plots
        
        %         function SimplePlotSpectrumAndResiduals(obj)
        %
        %         end
        
        function PlotSpectrumAndResiduals(obj)
            
            % Prepare Variables
            CMStat = diag(obj.Ydata);
            if isempty(obj.CovMat); obj.CovMat = CMStat; end
            Xmodel = obj.ModelObj.qU-obj.ModelObj.Q_i;
            Ymodel = obj.ModelObj.TBDIS./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec;
            
            XdataPlot = obj.Xdata - obj.ModelObj.Q_i;
            YdataPlot = obj.Ydata./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec;
            
            xlimmin = XdataPlot(obj.startqU);
            xlimmax = XdataPlot(end);
            ylimmin = 0.75*min(YdataPlot);
            ylimmax = 1.25*max(YdataPlot(obj.startqU:end));
            
            % Spectrum + Fit with Residuals
            figX = figure(obj.figN);
            set(figX, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.8, 1.2]);
            subplot(2,1,1);
            [lfit , ~] = boundedline(Xmodel,Ymodel,diag(sqrt(obj.CovMat)),...
                'alpha','cmap',rgb('CadetBlue'));
            hold on;
            hdata = errorbar(XdataPlot,YdataPlot,sqrt(obj.Ydata)./obj.ModelObj.qUfrac./obj.ModelObj.TimeSec,'ko',...
                'MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
            hdata.CapSize = 0;
            hold off;
            set(gca,'yscale','log');
            % legends
            myfitleg = sprintf(' Fit: \\chi2 / dof=%.2f/%.0f \n m^2=%.2f \t\\pm %.2f eV \n E0=%.2f \t\\pm %.2f eV \n B=%.2f \t\\pm %.2f mcps\n N=%.2f \t\\pm %.2f',...
                obj.FitResult.chi2min,obj.FitResult.dof,((obj.ModelObj.mnuSq_i+obj.FitResult.par(1))),obj.FitResult.err(1),(obj.ModelObj.Q_i+obj.FitResult.par(2)),obj.FitResult.err(2),...
                (obj.ModelObj.BKG_RateSec_i+obj.FitResult.par(3))*1e3,obj.FitResult.err(3)*1e3,obj.FitResult.par(4)+1,obj.FitResult.err(4));
            myleg = legend([hdata lfit],'Data',myfitleg,'Location','SouthWest');
            legend('boxoff')
            
            % labels
            title(obj.plotTitle);
            xlabel('retarding potential - 18575 (V)');
            ylabel('counts per second');
            
            PrettyFigureFormat;
            myleg.FontSize = 16;
            set(gca,'FontSize',16);
            set(gca,'YMinorTick','off');
            set(gca,'TickLength',[0.01 0.01]);
            xlim([xlimmin xlimmax]);
            ylim([ylimmin ylimmax]);
            
            subplot(2,1,2);
            [lstat,pstat]  = boundedline(XdataPlot,zeros(length(XdataPlot),1),sqrt(diag(CMStat))./sqrt(diag(obj.CovMat)),...
                'alpha','transparency', 0.4,'cmap',rgb('SteelBlue'));%[0 0 0]);%0.8*[0 0 0]);
            hold on;
            [lsys,psys] = boundedline(XdataPlot,zeros(length(XdataPlot),1),ones(length(XdataPlot),1),...
                'alpha','transparency',0.4,'cmap',rgb('DarkCyan'));
            lstat.LineStyle= '--'; lsys.LineStyle= '--';
            plot(XdataPlot,(obj.Ydata-obj.ModelObj.TBDIS)./sqrt(diag(obj.CovMat)),'o','MarkerEdgeColor','k');
            legend([pstat psys],'stat','stat+sys','Location','NorthWest'); %hsyst
            PrettyFigureFormat;
            legend('boxoff')
            hold off;
            xlabel('retarding potential - 18575 (V)');
            ylabel('norm. residuals');
            xlim([xlimmin xlimmax]);
            set(gca,'YMinorTick','off');
            set(gca,'XMinorTick','off');
            set(gca,'TickLength',[0.01 0.01]);
            set(gca,'FontSize',16);
            
            
            obj.savePlotFun();
            
        end
    end % Single/Stacked-pixel Plots
    
    
    methods % Ring Plots
        
        function plotRingIterative(obj)
            
            if size(obj.FitResult.par,2) > 4
                obj.FitResult.par(:,[5 6]) = [];
            end
            
            if size(obj.FitResult.err,2) > 4
                obj.FitResult.err(:,[5 6]) = [];
            end
            
            %obj.saveplot = 'export';
            allqU = obj.Xdata;
            allTBDIS = obj.Ydata;
            PlotName = cell(length(obj.RingList),1);
            RingResultsTable = nan(13,10);
            
            RingResultsTable(obj.RingList,:) = [obj.FitResult.par, obj.FitResult.chi2min, ...
                obj.FitResult.dof, obj.FitResult.err];
            
            
            for ri = 1:length(obj.RingList)
                ring = obj.RingList(ri);
                pixelsInRing = obj.ModelObj.ring{ring};
                
                obj.Xdata = mean(allqU(:,pixelsInRing),2);
                obj.Ydata = sum(allTBDIS(:,pixelsInRing),2);
                obj.FitResult.par = RingResultsTable(ring,1:4);
                obj.FitResult.err = RingResultsTable(ring,7:10);
                obj.FitResult.chi2min = RingResultsTable(ring,5);
                obj.FitResult.dof = RingResultsTable(ring,6);
                
                obj.savename = sprintf('DataSpectrumResiduals_Ring%.0d.pdf',ring);
                PlotName{ri} = ['plots/',obj.savename];
                
                obj.ModelObj.qU = obj.Xdata;
                obj.ModelObj = ref_RunAnalysis(obj.StackFileName,'','mpix',...
                    'FPD_Segmentation','RING',...
                    'FPD_Ring',ring,'nTeBinningFactor',20);
                obj.ModelObj.ComputeTBDDS('mSq_bias',obj.FitResult.par(1),...
                    'E0_bias',obj.FitResult.par(2),...
                    'N_bias',obj.FitResult.par(4),...
                    'B_bias',obj.FitResult.par(3));
                obj.ModelObj.ComputeTBDIS();
                obj.CovMat = diag(obj.Ydata);
                fprintf('-----PLOTTING SPECTRUM AND RESIDUALS RING %0d (please be patient) ------\n',ring)
                obj.PlotSpectrumAndResiduals();
                
                close all;
            end
            
            titlelist = {['Neutrino mass squared [eV^2]'],['Endpoint bias [eV]'],...
                ['Background [mcps]'],['Normalization'],...
                ['\chi^2 (DoF = ',num2str(obj.FitResult.dof),')']};
            
            
            RingResultsTable(obj.RingList,2) =  RingResultsTable(obj.RingList,2);% + obj.ModelObj.Q_i;
            RingResultsTable(obj.RingList,3) =  RingResultsTable(obj.RingList,3)*1e3;
            RingResultsTable(obj.RingList,9) =  RingResultsTable(obj.RingList,9)*1e3;
            RingResultsTable(obj.RingList,4) =  RingResultsTable(obj.RingList,4) + 1;
            
            FPDName = cell(5,1);
            for fpd = 1:5
                figfpd = figure(fpd);
                set(figfpd,'color','white');
                FPDViewer(RingResultsTable(:,fpd),'ReDrawSkeleton','OFF');
                titlehandle = title(titlelist{fpd});
                pos = get(figfpd,'position');
                set(figfpd,'position',[pos(1:2)/3.5 pos(3:4)*1.5])
                obj.savename = sprintf('FPD%u.pdf',fpd);
                FPDName{fpd} = ['plots/',obj.savename];
                obj.savePlotFun();
            end
            close all;
            
            BarName = cell(5,1);
            for barc = 1:5
                figbar = figure(barc);
                pos = get(figbar,'position');
                set(figbar,'position',[pos(1:2)/3.5 pos(3:4)*1.5])
                %screensize = (get(0, 'Screensize'));
                %set(figbar, 'Position', [screensize(1:2)*0.95,screensize(4),screensize(4)*0.95]);
                hold on
               bar(obj.RingList,RingResultsTable(obj.RingList,barc),'facecolor',rgb('CadetBlue'));
                if barc < 5
                    errorb(obj.RingList,RingResultsTable(obj.RingList,barc),RingResultsTable(obj.RingList,barc+6),'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1)
                    weightedMean = wmean(RingResultsTable(obj.RingList,barc),(1./(RingResultsTable(obj.RingList,barc+6))).^2);
                    weightedUncertainty = 1/sqrt(sum((1./(RingResultsTable(obj.RingList,barc+6))).^2));
                else
                    weightedMean = mean(RingResultsTable(obj.RingList,barc));
                    weightedUncertainty = 0;
                end
                set(gca,'XTick',obj.RingList);
                hold off
                titlehandle = title(titlelist{barc});
                titlehandle.FontSize = 7;
                PrettyFigureFormat;
                legend(sprintf('Mean = %0.3g +/- %0.3g',weightedMean,weightedUncertainty));
                obj.savename = sprintf('Bar%u.pdf',barc);
                BarName{barc} = ['plots/',obj.savename];
                obj.savePlotFun();
            end
            close all;
            
            if ~strcmp(obj.saveplot,'OFF')
                allNames = [FPDName;BarName;PlotName];
                obj.savename = 'Ring.pdf';
                obj.appendPDFex(allNames);
            end
            
        end % plotRingIterative
        
        function plotRingRhoDScan(obj,saverhoDmin,saverhoDUnc_low,saverhoDUnc_up,plotname)
            standardAveCD = 4.45e17;
            figRhoD = figure;
            screensize = (get(0, 'Screensize'));
            set(figRhoD, 'Position', [screensize(1:2)*0.95,screensize(3)*0.95,screensize(4)*0.95]);
            %set(figRhoD, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.7]);
            left_color = [0 0 0];
            right_color = [0 0 0];
            set(figRhoD,'defaultAxesColorOrder',[left_color; right_color]);
            
            yyaxis left
            ebhandle = errorbar(obj.RingList,saverhoDmin,saverhoDmin-saverhoDUnc_low,saverhoDUnc_up-saverhoDmin,...
                'ks','MarkerSize',12,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
            xlim([0.8 13.2])
            ylim([0 2.25*standardAveCD])
            xlabel('Ring')
            ylabel('\rho d_{min} (mol/cm^2)')
            
            weightedMeanRhoD = wmean(saverhoDmin,(1./(saverhoDUnc_up-saverhoDUnc_low)).^2);
            weightedUncertaintyRhoD = 1/sqrt(sum((1./(saverhoDUnc_up-saverhoDUnc_low)).^2));
            hold on
            plot(obj.RingList,weightedMeanRhoD*ones(1,length(obj.RingList)),'Color','red','LineStyle','--','LineWidth',3);
            
            yyaxis right
            errorbar(obj.RingList,saverhoDmin./standardAveCD*100,...
                (saverhoDmin-saverhoDUnc_low)./standardAveCD*100,...
                (saverhoDUnc_up-saverhoDmin)./standardAveCD*100,...
                'ks','MarkerSize',12,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
            set(gca,'XTick',obj.RingList);
            xlim([0.8 13.2])
            ylim([0 225])
            ylabel('100 % Col. Density = 4.45e17 mol/cm^2')
            title(obj.plotTitle,'interpreter','tex')
            PrettyFigureFormat;
            grid on;
            
            hold off
            legend(ebhandle,sprintf('Mean = %0.3g +/- %0.3g (error on mean)',weightedMeanRhoD,weightedUncertaintyRhoD));
            
            plotname = [plotname;sprintf('plots/RhoDScan_ring%.0d.pdf',obj.RingList(length(plotname))+1)];
            export_fig(plotname{length(plotname)},'-pdf');
            close all
            
            figbar = figure;
            set(figbar,'color','white');
            screensize = (get(0, 'Screensize'));
            set(figbar, 'Position', [screensize(1:2)*0.95,screensize(4)*0.95,screensize(4)*0.95]);
            [mu,sigma] = normfit(saverhoDmin);
            hold on
            [~,Nbins,Wbins] = nhist(saverhoDmin);
            
            binWidth = Wbins(2) - Wbins(1);
            x_N = linspace(min(saverhoDmin),max(saverhoDmin),length(Nbins)*10);
            plot(x_N,length(obj.RingList)*binWidth*normpdf(x_N,mu,sigma),'linewidth',2);
            hold off
            
            
            textString = {['\sigma = ',num2str(round(sigma,3,'significant')),' mol/cm^2'],['mean = ',num2str(round(mu,3,'significant')),' mol/cm^2']};
            l = legend(textString);
            l.FontSize = 16;
            t = annotation('textbox');
            t.String = textString;
            t.FontSize = 16;
            t.Position = l.Position;
            t.FitBoxToText = 'on';
            delete(l);
            
            plotname = [plotname;sprintf('plots/RhoDScan_ring%.0d.pdf',obj.RingList(length(plotname)-1)+2)];
            export_fig(plotname{length(plotname)},'-pdf');
            close all
            
            obj.savename = 'RhoDScanRing.pdf';
            obj.appendPDFex(plotname)
            
        end
        
    end % Ring Plots
    
    
    
    methods % Multipixel Plots
        
    end % Multipixel Plots
    
    methods % Save plots
        function savePlotFun(obj)
            switch obj.saveplot
                case 'ON'
                    try
                        publish_figurePDF(obj.figN,['plots/',obj.savename]);
                    catch
                        mkdir plots;
                        publish_figurePDF(obj.figN,['plots/',obj.savename]);
                    end
                case 'export'
                    try
                        export_fig(['plots/',obj.savename],'-pdf');
                    catch
                        mkdir plots;
                        export_fig(['plots/',obj.savename],'-pdf');
                    end
            end
            
        end %savePlotFun
        
        function appendPDFex(obj,FileNames)
            
            
            if exist(['plots/Summary',obj.savename],'file') == 2
                summaryNr = 1;
                savename_tmp = [num2str(summaryNr),obj.savename];
                while exist(['plots/Summary',savename_tmp],'file') == 2
                    summaryNr = summaryNr + 1;
                    savename_tmp = [num2str(summaryNr),obj.savename];
                end
                obj.savename = savename_tmp;
            end
            
            
            if isunix || ismac
                unix(['touch plots/Summary',obj.savename])
            end
            
            append_pdfs(['plots/Summary',obj.savename],FileNames{:});
            delete(FileNames{:});
            
        end % appendPDFex
        
    end % Save plots
    
    methods %Secret Methods (don't use)
        
        function savepdfiterative(obj,pixlist)
            par = obj.RESULTS{1}; err = obj.RESULTS{2};
            chi2min = obj.RESULTS{3}; DoF = obj.RESULTS{5};
            
            mnu_fit = par(pixlist,1);
            e0_fit = par(pixlist,2);
            norms_fit = par(pixlist,4);
            bcks_fit = par(pixlist,3);
            chi2min = chi2min(pixlist);
            %             mnu_fit = par(:,1);
            %             e0_fit = par(:,2);
            %             norms_fit = par(:,4);
            %             bcks_fit = par(:,3);
            %             chi2min = chi2min(:);
            
            fprintf('Setting figure visibility to OFF. \n')
            % Plots from the focal plane detector
            
            % neutrino mass
            fprintf('Plot mnu per pixel in FPD \n')
            pixel_mnu = NaN(148,1);
            pixel_mnu(pixlist) = mnu_fit;
            
            figure(99);
            set(gcf,'Visible','off');
            set(gcf,'color','white');
            FPDViewer(pixel_mnu,'ReDrawSkeleton','ON');
            titlehandle = title(['Neutrino mass squared [eV^2] ',obj.SO.TD],'interpreter','none');
            titlehandle.FontSize = 7;
            if ispc; export_fig('plots/FPDmnu.pdf','-pdf');
            elseif isunix || ismac; publish_figure(99,'plots/ftmultipix00001.eps'); end
            fprintf('Plot mnu done \n')
            
            % endpoint
            fprintf('Plot E0 per pixel in FPD \n')
            pixel_e0 = NaN(148,1);
            pixel_e0(pixlist) = e0_fit;
            
            figure(100);
            set(gcf,'Visible','off');
            set(gcf,'color','white');
            FPDViewer(pixel_e0,'ReDrawSkeleton','ON');
            titlehandle = title(['Endpoint bias [eV] ',obj.SO.TD],'interpreter','none');
            titlehandle.FontSize = 7;
            if ispc; export_fig('plots/FPDE0.pdf','-pdf');
            elseif isunix || ismac; publish_figure(100,'plots/ftmultipix0001.eps'); end
            fprintf('Plot E0 done \n')
            
            % normalization
            fprintf('Plot NORM per pixel in FPD \n')
            pixel_norms = NaN(148,1);
            pixel_norms(pixlist) = norms_fit;
            
            figure(1);
            set(gcf,'Visible','off');
            set(gcf,'color','white');
            FPDViewer(pixel_norms+1,'ReDrawSkeleton','ON');
            titlehandle = title(['Normalization ',obj.SO.TD],'interpreter','none');
            titlehandle.FontSize = 7;
            if ispc; export_fig('plots/FPDNorm.pdf','-pdf');
            elseif isunix || ismac; publish_figure(1,'plots/ftmultipix001.eps'); end
            fprintf('Plot NORM done \n')
            
            % background
            fprintf('Plot BCK per pixel in FPD \n')
            pixel_bcks = NaN(148,1);
            pixel_bcks(pixlist) = bcks_fit*1e3;
            
            figure(2);
            set(gcf,'Visible','off');
            set(gcf,'color','white');
            FPDViewer(pixel_bcks,'ReDrawSkeleton','ON');
            titlehandle = title(['Background [mcps] ',obj.SO.TD],'interpreter','none');
            titlehandle.FontSize = 7;
            if ispc; export_fig('plots/FPDBCK.pdf','-pdf');
            elseif isunix || ismac; publish_figure(2,'plots/ftmultipix002.eps'); end
            fprintf('Plot BCK done \n')
            
            % chi2
            fprintf('Plot chi2 per pixel in FPD \n')
            pixel_chi2 = NaN(148,1);
            pixel_chi2(pixlist) = chi2min;
            
            figure(98);
            set(gcf,'Visible','off');
            set(gcf,'color','white');
            FPDViewer(pixel_chi2,'ReDrawSkeleton','ON');
            titlehandle = title(['Chi2 ',obj.SO.TD, 'DoF = ',num2str(DoF)],'interpreter','none');
            titlehandle.FontSize = 7;
            if ispc; export_fig('plots/FPDchi2.pdf','-pdf');
            elseif isunix || ismac; publish_figure(98,'plots/ftmultipix0025.eps'); end
            fprintf('Plot chi2 done \n')
            
            % Distributions of the parameters per pixel
            nlength = length(pixlist);
            
            % neutrino mass distribution
            figure(101);
            set(gcf,'Visible','off');
            set(gcf,'color','white');
            [mu_N,sigma_N] = normfit(mnu_fit);
            hold on
            [~,Nbins,Wbins] = nhist(mnu_fit);
            x_N = linspace(min(mnu_fit),max(mnu_fit),length(Nbins)*obj.SO.nTeBinningFactor);
            binWidth = Wbins(2) - Wbins(1);
            plot(x_N,nlength*binWidth*normpdf(x_N,mu_N,sigma_N));
            hold off
            annotation('textbox',[0.6 0.5 0.2 0.2],'String',...
                {['\sigma = ',num2str(sigma_N),' eV^2'],['mean = ',num2str(mu_N),' eV^2']},'FitBoxToText','on');
            title('neutrino mass squared');
            xlabel('neutrino mass squared');
            if ispc; export_fig('plots/MNUdist.pdf','-pdf');
            elseif isunix || ismac; publish_figure(101,'plots/ftmultipix003.eps'); end
            fprintf('Plot mnu dist done \n')
            
            % endpoint distribution
            figure(200);
            set(gcf,'Visible','off');
            set(gcf,'color','white');
            [mu_N,sigma_N] = normfit(e0_fit);
            hold on
            [~,Nbins,Wbins] = nhist(e0_fit);
            x_N = linspace(min(e0_fit),max(e0_fit),length(Nbins)*obj.SO.nTeBinningFactor);
            binWidth = Wbins(2) - Wbins(1);
            plot(x_N,nlength*binWidth*normpdf(x_N,mu_N,sigma_N));
            hold off
            annotation('textbox',[0.6 0.5 0.2 0.2],'String',...
                {['\sigma = ',num2str(sigma_N),' eV'],['mean = ',num2str(mu_N),' eV']},'FitBoxToText','on');
            title('endpoint bias from 18575 [eV]');
            xlabel('endpoint bias [eV]');
            if ispc; export_fig('plots/E0dist.pdf','-pdf');
            elseif isunix || ismac; publish_figure(200,'plots/ftmultipix004.eps'); end
            fprintf('Plot e0 dist done \n')
            
            % background
            figure(4);
            set(gcf,'Visible','off');
            set(gcf,'color','white');
            bcks_residuals = bcks_fit*1e3;
            [mu_b,sigma_b] = normfit(bcks_residuals);
            hold on
            [~,Nbins,Wbins] = nhist(bcks_residuals);
            x_b = linspace(min(bcks_residuals),max(bcks_residuals),length(Nbins)*obj.SO.nTeBinningFactor);
            binWidth = Wbins(2) - Wbins(1);
            plot(x_b,nlength*binWidth*normpdf(x_b,mu_b,sigma_b));
            hold off
            annotation('textbox',[0.6 0.5 0.2 0.2],'String',...
                {['\sigma = ',num2str(sigma_b),' mcps'],['mean = ',num2str(mu_b),' mcps']},'FitBoxToText','on');
            title('background [mcps]');
            xlabel('background [mcps]');
            if ispc; export_fig('plots/BCKdist.pdf','-pdf');
            elseif isunix || ismac; publish_figure(4,'plots/ftmultipix005.eps'); end
            fprintf('Plot BCK dist done \n')
            
            % normalization
            figure(5);
            set(gcf,'Visible','off');
            set(gcf,'color','white');
            [mu_N,sigma_N] = normfit(norms_fit);
            hold on
            [~,Nbins,Wbins] = nhist(norms_fit);
            x_N = linspace(min(norms_fit),max(norms_fit),length(Nbins)*obj.SO.nTeBinningFactor);
            binWidth = Wbins(2) - Wbins(1);
            plot(x_N,nlength*binWidth*normpdf(x_N,mu_N,sigma_N));
            hold off
            annotation('textbox',[0.6 0.5 0.2 0.2],'String',...
                {['\sigma = ',num2str(sigma_N),''],['mean = ',num2str(mu_N),'']},'FitBoxToText','on');
            title('normalization');
            xlabel('normalization');
            if ispc; export_fig('plots/NORMdist.pdf','-pdf');
            elseif isunix || ismac; publish_figure(5,'plots/ftmultipix0055.eps'); end
            fprintf('Plot norm dist done \n')
            
            % chi2min
            figure(309);
            set(gcf,'Visible','off');
            set(gcf,'color','white');
            [mu_N,sigma_N] = normfit(chi2min);
            hold on
            [~,Nbins,Wbins] = nhist(chi2min);
            x_N = linspace(min(chi2min),max(chi2min),length(Nbins)*obj.SO.nTeBinningFactor);
            binWidth = Wbins(2) - Wbins(1);
            plot(x_N,nlength*binWidth*normpdf(x_N,mu_N,sigma_N));
            hold off
            annotation('textbox',[0.6 0.5 0.2 0.2],'String',...
                {['DoF = ',num2str(DoF),'']},'FitBoxToText','on');
            title('\chi^2');
            xlabel('\chi^2');
            if ispc; export_fig('plots/CHI2dist.pdf','-pdf');
            elseif isunix || ismac; publish_figure(309,'plots/ftmultipix0056.eps'); end
            fprintf('Plot chi2 dist done \n')
            
            
            
            % Plots of correlation from background and normalization
            figure(6);
            set(gcf,'Visible','off');
            set(gcf,'color','white');
            corrplotm([mnu_fit,e0_fit,bcks_fit*1e3,norms_fit],'varNames',{'MNU','E0','BCK','NORM'});
            set(gcf,'Visible','off');
            ax = axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
            set(get(ax,'Title'),'Visible','on');
            ax.Title.String = {' ',['Correlation Multipixel Fit']};
            if ispc; export_fig('plots/CORRALL.pdf','-pdf');
            elseif isunix || ismac; publish_figure(6,'plots/ftmultipix006.eps'); end
            fprintf('Plot correlation done \n')
            
            
            % Plots from each of the spectra of each pixel
            
            plotNameSpectraPixel = cell(1,nlength);
            D_IS = obj.counts; D_ISE = obj.c_error;
            rangefit = obj.exclDataStart:obj.exclDataStop;
            M_IS = zeros(obj.SO.nqU,nlength);
            for p = 1:nlength
                obj.SO.ComputeTBDDS(...
                    'mSq_bias',mnu_fit(p),...
                    'E0_bias',e0_fit(p),...
                    'B_bias',bcks_fit(p),...
                    'N_bias',norms_fit(p));
                
                obj.SO.ComputeTBDIS();
                M_IS(:,p) = obj.SO.TBDIS;
                
                
                fprintf('Plotting data, spectrum and residuals pixel %0.d \n',pixlist(p))
                figure(p+2000);
                set(gcf,'Visible','off');
                set(gcf,'color','white');
                % Plot of the data points, spectra and error bars
                subplot(2,1,1)
                hold on
                hline = plot(obj.qUdata(rangefit,p) - obj.SO.Q, M_IS(rangefit,p),...
                    'LineWidth',1,'LineStyle','-','Color','Black');
                switch obj.chi2name
                    case {'chi2Stat','chi2P'}
                        hdata = errorbar(obj.qUdata(rangefit,p) - obj.SO.Q, D_IS(rangefit,p), D_ISE(rangefit,p),...
                            'o','Color','black','MarkerEdgeColor','black','MarkerFaceColor','black','MarkerSize',5,'LineWidth',1);
                        
                    case {'chi2CM','chi2CMFrac'}
                        hdata = errorbar(obj.qUdata(rangefit,p) - obj.SO.Q, D_IS(rangefit,p), sqrt(diag(covmat{pixlist(p)})),...
                            'o','Color','black','MarkerEdgeColor','black','MarkerFaceColor','black','MarkerSize',5,'LineWidth',1);
                end
                mydata = ['Data: E_0eff = ??? eV'];
                myline = sprintf('Fit: \\delta E_0eff = %.3f \\pm %.3f eV',e0_fit(p),err(p,2));
                
                if strcmp(obj.SO.FPD_Segmentation,'OFF') || strcmp(obj.SO.FPD_Segmentation,'SINGLEPIXEL')
                    hchi2 = plot(obj.qUdata(rangefit,p) - obj.SO.Q, M_IS(rangefit,p),...
                        'LineWidth',1,'LineStyle','-','Color','Black');
                    mychi2 = sprintf('\\chi^2/DoF = %0.3f / %0.d',chi2min(p),DoF);
                    lgd_fit = legend([hline hchi2 hdata],myline,mychi2,mydata,'Location','Best');
                    set(gca, 'yscale', 'log')
                    
                else
                    lgd_fit = legend([hdata hline],'Data','Model',...
                        'Location','Best');
                end
                
                hold off
                lgd_fit.FontSize = 8;
                
                xlabel('qU-E_0 [eV]','FontSize',14);
                ylabel('Counts','FontSize',14);
                title(sprintf('KATRIN - %s - pix %g - %0.0f s ',obj.SO.TD,pixlist(p),obj.SO.TimeSec),'interpreter','none');
                set(gca,'FontSize',12);
                PrettyFigureFormat;
                
                % Plot of the residuals according to uncertainty of each data point
                subplot(2,1,2)
                hold on
                switch obj.chi2name
                    case {'chi2Stat','chi2P'}
                        hresiduals = scatter(obj.qUdata(rangefit,p)-obj.SO.Q,...
                            ((D_IS(rangefit,p)-M_IS(rangefit,p))./D_ISE(rangefit,p)),...
                            'o','filled','LineWidth',1);
                        hresidualserr = errorbar(obj.qUdata(rangefit,p)-obj.SO.Q,...
                            (D_IS(rangefit,p)-M_IS(rangefit,p))./D_ISE(rangefit,p),...
                            D_ISE(rangefit,p)./D_ISE(rangefit,p),'k','LineWidth',1);
                        plot(linspace(min(obj.qUdata(rangefit,p))-obj.SO.Q,...
                            max(obj.qUdata(rangefit,p))-obj.SO.Q,length(obj.qUdata(rangefit,p))),...
                            zeros(length(obj.qUdata(rangefit,p)),1),'k');
                        lgd_res = legend([hresiduals],'(data-fit)/uncertainty',...
                            'Location','Best');
                        
                    case {'chi2CM','chi2CMFrac'}
                        hresiduals = scatter(obj.qUdata(rangefit,p)-obj.SO.Q,...
                            ((D_IS(rangefit,p)-M_IS(rangefit,p))./sqrt(diag(covmat{pixlist(p)}))),...
                            'o','filled','k','LineWidth',1);
                        hresidualsline = plot(obj.qUdata(rangefit,p)-obj.SO.Q,...
                            ((D_IS(rangefit,p)-M_IS(rangefit,p))./sqrt(diag(covmat{pixlist(p)}))),...
                            'k','LineWidth',1);
                        plot(linspace(min(obj.qUdata(rangefit,p))-obj.SO.Q,...
                            max(obj.qUdata(rangefit,p))-obj.SO.Q,length(obj.qUdata(rangefit,p))),...
                            zeros(length(obj.qUdata(rangefit,p)),1),'k');
                        lgd_res = legend([hresiduals],'(data-fit)/(total uncertainty)',...
                            'Location','Best');
                end
                hold off
                grid on
                xlabel('qU-E_0 (eV)','FontSize',14);
                ylabel('Residuals','FontSize',14);
                set(gca,'FontSize',12);
                
                lgd_res.FontSize = 8;
                PrettyFigureFormat;
                
                
                % Saving the generated plot
                if ispc
                    plotNameSpectraPixel(1,p) = {sprintf('plots/spectrumpix%g.pdf',p)};
                    export_fig(plotNameSpectraPixel{1,p},'-pdf');
                elseif isunix || ismac
                    spectra_figname = sprintf('plots/ftmultipix%04d.eps',p+2000);
                    publish_figure(p+2000,spectra_figname);
                end
            end
            
            % Fit results per pixel
            if false
                latex_preamble = {'\documentclass[10pt,a4paper]{article}';...
                    '\usepackage[utf8]{inputenc}';...
                    '\usepackage{amsmath}';...
                    '\usepackage{amsfonts}';...
                    '\usepackage{amssymb}';...
                    '\begin{document}';...
                    '\begin{small}'};...
                    
                fid = fopen('latex_preamble.tex','w');
                fprintf(fid,'%s\n', latex_preamble{:});
                fclose(fid);
                if strcmp(obj.SO.FPD_Segmentation,'MULTIPIXEL')
                    results_size = length(bcks_fit);
                    index1 = ceil(results_size/3);
                    index2 = 2*index1;
                    index3 = results_size;
                    space = cell(index1,1);
                    space(:) = {' '};
                    extraentries = cell(3*index1-index3,3);
                    horizontallabels = {'Pix.','Bck.','Norm.',' ','Pix.','Bck.','Norm.',' ','Pix.','Bck.','Norm.';...
                        ' ','[mcps]',' ',' ',' ','[mcps]',' ',' ',' ','[mcps]',' '};
                    numformatlatex = {'%0.d','%0.4g','%0.3g','%0.d','%0.d','%0.4g','%0.3g','%0.d','%0.d','%0.4g','%0.3g','%0.d'};
                    bcks_fit = bcks_fit*1e3;
                    celltableforlatex = [num2cell([(1:index1)',bcks_fit(1:index1)',norms_fit(1:index1)']),space,...
                        num2cell([(index1+1:index2)',bcks_fit(index1+1:index2)',norms_fit(index1+1:index2)']),space,...
                        [num2cell([(index2+1:index3)',bcks_fit(index2+1:index3)',norms_fit(index2+1:index3)']);extraentries]];
                    latextable(celltableforlatex,'name','fitresultstable.tex','Horiz',horizontallabels,'Hline',[0,2,NaN],...
                        'Vline',[1,2,3,5,6,7,9,10],'format',numformatlatex);
                    bcks_fit = bcks_fit/1e3;
                else
                    horizontallabels = {'$\nu^2$','unc.','$E_0$','unc.','Bck.','unc.','Norm.','unc.';...
                        '(ev$^2$)','(ev$^2$)','(eV)','(eV)','(mcps)','(mcps)',' ',' '};
                    numformatlatex = {'%0.3g','%0.3g','%0.3g','%0.3g','%0.d','%0.3g','%0.4g','%0.3g'};
                    bcks_fit = bcks_fit*1e3;
                    celltableforlatex = num2cell([par(1),err(1),par(2),err(2),bcks_fit,err(3),norms_fit,err(4)]);
                    latextable(celltableforlatex,'name','fitresultstable.tex','Horiz',horizontallabels,'Hline',[0,NaN],...
                        'Vline',[0,1,2,3,5,4,6,7,8],'format',numformatlatex);
                    bcks_fit = bcks_fit/1e3;
                end
                
                latex_end = {'\end{small}';...
                    '\end{document}'};
                
                fid = fopen('latex_end.tex','w');
                fprintf(fid,'%s\n', latex_end{:});
                fclose(fid);
                
                system('copy latex_preamble.tex+fitresultstable.tex+latex_end.tex fitresultslatex.tex')
                
                system('pdflatex fitresultslatex.tex')
                
                movefile fitresultslatex.pdf plots
            end
            % Save all plots in one PDF and erase the individual PDFs
            
            if ispc
                delete('plots/ResultsSummary.pdf');
                append_pdfs('plots/ResultsSummary.pdf',...
                    'plots/FPDmnu.pdf','plots/FPDe0.pdf','plots/FPDBCK.pdf','plots/FPDNORM.pdf','plots/FPDchi2.pdf',...
                    'plots/MNUdist.pdf','plots/E0dist.pdf','plots/BCKdist.pdf','plots/NORMdist.pdf','plots/CHI2dist.pdf','plots/CORRALL.pdf',...
                    plotNameSpectraPixel{:})
                delete('plots/spectrumpix*.pdf','plots/FPDNorm.pdf','plots/FPDBCK.pdf','plots/FitResults.pdf',...
                    'plots/BCKdist.pdf','plots/NORMdist.pdf','plots/CORRBCKNorm.pdf','plots/FitResultsPix.pdf',...
                    'plots/FPDmnu.pdf','plots/FPDe0.pdf','plots/FPDchi2.pdf','plots/MNUdist.pdf','plots/E0dist.pdf',...
                    'plots/CHI2dist.pdf','plots/CORRALL.pdf');
            elseif isunix || ismac
                PATH = getenv('PATH');
                setenv('PATH', [PATH,':/usr/local/bin/']);
                cd plots
                command = 'gs -sDEVICE=pdfwrite -sOutputFile="ResultsSummaryTMP.pdf" -dNOPAUSE -dEPSCrop -c "<</Orientation 0>> setpagedevice" -f ftmultipix*.eps -c quit';
                unix(command);
                unix('rm ftmultipix*.eps');
                
                unix('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="ResultsSummary.pdf" ResultsSummaryTMP.pdf fitresultslatex.pdf')
                cd ..
            end
            
            fprintf('Saving PDFs done. \n')
            close all
        end
        
        function savepdffunction(obj)
            par = obj.RESULTS{1}; err = obj.RESULTS{2};
            chi2min = obj.RESULTS{3};
            DoF = obj.SO.nPixels*(obj.exclDataStop-obj.exclDataStart) - ...
                length(obj.parinit);
            norms_fit = (3+obj.SO.nPixels:3+2*obj.SO.nPixels-1);
            bcks_fit = par(3:2+obj.SO.nPixels);
            
            if strcmp(obj.SO.FPD_Segmentation,'MULTIPIXEL')
                fprintf('Setting figure visibility to OFF. \n')
                % Plots from the focal plane detector
                fprintf('Plot NORM per pixel in FPD \n')
                pixel_norms = NaN(148,1);
                pixel_norms(obj.SO.FPD_Pixel) = norms_fit;
                
                figure(1);
                set(gcf,'Visible','off');
                set(gcf,'color','white');
                FPDViewer(pixel_norms,'ReDrawSkeleton','ON');
                title(['Normalization bias ',obj.SO.TD,' ',num2str(obj.SO.TimeSec), 's']);
                if ispc; export_fig('plots/FPDNorm.pdf','-pdf');
                elseif isunix || ismac; publish_figure(1,'plots/ftmultipix001.eps'); end
                fprintf('Plot NORM done \n')
                
                
                fprintf('Plot BCK per pixel in FPD \n')
                pixel_bcks = NaN(148,1);
                pixel_bcks(obj.SO.FPD_Pixel) = bcks_fit;
                figure(2);
                set(gcf,'Visible','off');
                set(gcf,'color','white');
                FPDViewer(pixel_bcks,'ReDrawSkeleton','ON');
                title(['Background ',obj.SO.TD,' ',num2str(obj.SO.TimeSec), 's']);
                
                if ispc; export_fig('plots/FPDBCK.pdf','-pdf');
                elseif isunix || ismac; publish_figure(2,'plots/ftmultipix002.eps'); end
                fprintf('Plot BCK done \n')
                
                
                % Distributions of the parameters per pixel
                figure(4);
                set(gcf,'Visible','off');
                set(gcf,'color','white');
                bcks_residuals = bcks_fit - obj.counts(end,:);
                [mu_b,sigma_b] = normfit(bcks_residuals);
                hold on
                [~,Nbins,Wbins] = nhist(bcks_residuals);
                x_b = linspace(min(bcks_residuals),max(bcks_residuals),length(Nbins)*obj.SO.nTeBinningFactor);
                binWidth = Wbins(2) - Wbins(1);
                plot(x_b,obj.SO.nPixels*binWidth*normpdf(x_b,mu_b,sigma_b));
                hold off
                annotation('textbox',[0.6 0.5 0.2 0.2],'String',...
                    {['\sigma = ',num2str(sigma_b),' cps'],['mean = ',num2str(mu_b),' cps']},'FitBoxToText','on');
                title('\delta background [cps]');
                xlabel('\delta background [cps]');
                if ispc; export_fig('plots/BCKdist.pdf','-pdf');
                elseif isunix || ismac; publish_figure(4,'plots/ftmultipix004.eps'); end
                fprintf('Plot BCK dist done \n')
                
                figure(5);
                set(gcf,'Visible','off');
                set(gcf,'color','white');
                norms_fit = norms_fit + 1;
                [mu_N,sigma_N] = normfit(norms_fit);
                hold on
                [~,Nbins,Wbins] = nhist(norms_fit);
                x_N = linspace(min(norms_fit),max(norms_fit),length(Nbins)*obj.SO.nTeBinningFactor);
                binWidth = Wbins(2) - Wbins(1);
                plot(x_N,obj.SO.nPixels*binWidth*normpdf(x_N,mu_N,sigma_N));
                hold off
                annotation('textbox',[0.6 0.5 0.2 0.2],'String',...
                    {['\sigma = ',num2str(sigma_N),''],['mean = ',num2str(mu_N),'']},'FitBoxToText','on');
                title('normalization');
                xlabel('normalization');
                if ispc; export_fig('plots/NORMdist.pdf','-pdf');
                elseif isunix || ismac; publish_figure(5,'plots/ftmultipix005.eps'); end
                fprintf('Plot norm dist done \n')
                
                % Plots of correlation from background and normalization
                figure(6);
                set(gcf,'Visible','off');
                set(gcf,'color','white');
                corrplotm([bcks_fit',norms_fit'],'varNames',{'BCK','NORM'});
                set(gcf,'Visible','off');
                ax = axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
                set(get(ax,'Title'),'Visible','on');
                ax.Title.String = {' ',['Correlation Multipixel Fit']};
                if ispc; export_fig('plots/CORRBCKNorm.pdf','-pdf');
                elseif isunix || ismac; publish_figure(6,'plots/ftmultipix006.eps'); end
                fprintf('Plot correlation done \n')
            end
            
            % Plots from each of the spectra of each pixel
            
            M_IS = obj.SO.TBDIS;
            D_IS = obj.counts;
            D_ISE = obj.c_error;
            
            plotNameSpectraPixel = cell(1,obj.SO.nPixels);
            for p = 1:obj.SO.nPixels
                fprintf('Plotting data, spectrum and residuals pixel %0.d \n',p)
                figure(p+6);
                set(gcf,'Visible','off');
                set(gcf,'color','white');
                % Plot of the data points, spectra and error bars
                subplot(2,1,1)
                hold on
                hline = plot(obj.qUdata(:,p) - obj.SO.Q, M_IS(:,p),...
                    'LineWidth',1,'LineStyle','-','Color','Black');
                hdata = errorbar(obj.qUdata(:,p) - obj.SO.Q, D_IS(:,p), D_ISE(:,p),...
                    'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
                mydata = ['Data: E_0eff = ??? eV'];
                myline = sprintf('Fit: \\delta E_0eff = %.3f \\pm %.3f eV',par(2),err(2));
                
                if strcmp(obj.SO.FPD_Segmentation,'OFF') || strcmp(obj.SO.FPD_Segmentation,'SINGLEPIXEL')
                    hchi2 = plot(obj.qUdata(:,p) - obj.SO.Q, M_IS(:,p),...
                        'LineWidth',1,'LineStyle','-','Color','Black');
                    mychi2 = sprintf('\\chi^2/DoF = %0.3f / %0.d',chi2min,DoF);
                    lgd_fit = legend([hline hchi2 hdata],myline,mychi2,mydata,'Location','Best') ;
                else
                    lgd_fit = legend([hdata hline],'Data','Model',...
                        'Location','Best');
                end
                hold off
                lgd_fit.FontSize = 8;
                
                xlabel('qU-E_0 [eV]','FontSize',14);
                ylabel('Counts','FontSize',14);
                title(sprintf('KATRIN - %s - pix %g - %0.0f s ',obj.SO.TD,p,obj.SO.TimeSec));
                set(gca,'FontSize',12);
                PrettyFigureFormat;
                
                % Plot of the residuals according to uncertainty of each data point
                subplot(2,1,2)
                hold on
                hresiduals = scatter(obj.qUdata(:,p)-obj.SO.Q,...
                    ((D_IS(:,p)-M_IS(:,p))./D_ISE(:,p)),...
                    'ks','LineWidth',1);
                hresidualserr = errorbar(obj.qUdata(:,p)-obj.SO.Q,...
                    (D_IS(:,p)-M_IS(:,p))./D_ISE(:,p),...
                    D_ISE(:,p)./D_ISE(:,p),'ks','LineWidth',1);
                plot(linspace(min(obj.qUdata(:,p))-obj.SO.Q,...
                    max(obj.qUdata(:,p))-obj.SO.Q,length(obj.qUdata(:,p))),...
                    zeros(length(obj.qUdata(:,p)),1),'k')
                hold off
                grid on
                xlabel('qU-E_0 (eV)','FontSize',14);
                ylabel('Residuals','FontSize',14);
                set(gca,'FontSize',12);
                lgd_res = legend([hresidualserr],'(data-fit)/uncertainty',...
                    'Location','Best');
                lgd_res.FontSize = 8;
                PrettyFigureFormat;
                
                % Saving the generated plot
                if ispc
                    plotNameSpectraPixel(1,p) = {sprintf('plots/spectrumpix%g.pdf',p)};
                    export_fig(plotNameSpectraPixel{1,p},'-pdf');
                elseif isunix || ismac
                    spectra_figname = sprintf('plots/ftmultipix%03d.eps',p+6);
                    publish_figure(p+6,spectra_figname);
                end
            end
            
            % Fit results per pixel
            
            latex_preamble = {'\documentclass[10pt,a4paper]{article}';...
                '\usepackage[utf8]{inputenc}';...
                '\usepackage{amsmath}';...
                '\usepackage{amsfonts}';...
                '\usepackage{amssymb}';...
                '\begin{document}';...
                '\begin{small}'};...
                
            fid = fopen('latex_preamble.tex','w');
            fprintf(fid,'%s\n', latex_preamble{:});
            fclose(fid);
            if strcmp(obj.SO.FPD_Segmentation,'MULTIPIXEL')
                results_size = length(bcks_fit);
                index1 = ceil(results_size/3);
                index2 = 2*index1;
                index3 = results_size;
                space = cell(index1,1);
                space(:) = {' '};
                extraentries = cell(3*index1-index3,3);
                horizontallabels = {'Pix.','Bck.','Norm.',' ','Pix.','Bck.','Norm.',' ','Pix.','Bck.','Norm.';...
                    ' ','[mcps]',' ',' ',' ','[mcps]',' ',' ',' ','[mcps]',' '};
                numformatlatex = {'%0.d','%0.4g','%0.3g','%0.d','%0.d','%0.4g','%0.3g','%0.d','%0.d','%0.4g','%0.3g','%0.d'};
                bcks_fit = bcks_fit*1e3;
                celltableforlatex = [num2cell([(1:index1)',bcks_fit(1:index1)',norms_fit(1:index1)']),space,...
                    num2cell([(index1+1:index2)',bcks_fit(index1+1:index2)',norms_fit(index1+1:index2)']),space,...
                    [num2cell([(index2+1:index3)',bcks_fit(index2+1:index3)',norms_fit(index2+1:index3)']);extraentries]];
                latextable(celltableforlatex,'name','fitresultstable.tex','Horiz',horizontallabels,'Hline',[0,2,NaN],...
                    'Vline',[1,2,3,5,6,7,9,10],'format',numformatlatex);
                bcks_fit = bcks_fit/1e3;
            else
                horizontallabels = {'$\nu^2$','unc.','$E_0$','unc.','Bck.','unc.','Norm.','unc.';...
                    '(ev$^2$)','(ev$^2$)','(eV)','(eV)','(mcps)','(mcps)',' ',' '};
                numformatlatex = {'%0.3g','%0.3g','%0.3g','%0.3g','%0.d','%0.3g','%0.4g','%0.3g'};
                bcks_fit = bcks_fit*1e3;
                celltableforlatex = num2cell([par(1),err(1),par(2),err(2),bcks_fit,err(3),norms_fit,err(4)]);
                latextable(celltableforlatex,'name','fitresultstable.tex','Horiz',horizontallabels,'Hline',[0,NaN],...
                    'Vline',[0,1,2,3,5,4,6,7,8],'format',numformatlatex);
                bcks_fit = bcks_fit/1e3;
            end
            
            latex_end = {'\end{small}';...
                '\end{document}'};
            
            fid = fopen('latex_end.tex','w');
            fprintf(fid,'%s\n', latex_end{:});
            fclose(fid);
            
            system('copy latex_preamble.tex+fitresultstable.tex+latex_end.tex fitresultslatex.tex')
            
            system('pdflatex fitresultslatex.tex')
            
            %movefile fitresultslatex.pdf plots
            
            % Save all plots in one PDF and erase the individual PDFs
            
            if ispc
                delete('plots/ResultsSummary.pdf');
                if strcmp(obj.SO.FPD_Segmentation,'MULTIPIXEL')
                    append_pdfs('plots/ResultsSummary.pdf',...
                        'plots/FitResults.pdf','plots/fitresultslatex.pdf','plots/FPDNorm.pdf','plots/FPDBCK.pdf',...
                        'plots/BCKdist.pdf','plots/NORMdist.pdf','plots/CORRBCKNorm.pdf',...
                        plotNameSpectraPixel{:})
                else
                    append_pdfs('plots/ResultsSummary.pdf',...
                        'plots/FitResults.pdf','fitresultslatex.pdf',...
                        plotNameSpectraPixel{:})
                end
                delete('plots/spectrumpix*.pdf','plots/FPDNorm.pdf','plots/FPDBCK.pdf','plots/FitResults.pdf',...
                    'plots/BCKdist.pdf','plots/NORMdist.pdf','plots/CORRBCKNorm.pdf','plots/FitResultsPix.pdf');
            elseif isunix || ismac
                PATH = getenv('PATH');
                setenv('PATH', [PATH,':/usr/local/bin/']);
                cd plots
                command = 'gs -sDEVICE=pdfwrite -sOutputFile="ResultsSummaryTMP.pdf" -dNOPAUSE -dEPSCrop -c "<</Orientation 0>> setpagedevice" -f ftmultipix*.eps -c quit';
                unix(command);
                unix('rm ftmultipix*.eps');
                
                unix('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="ResultsSummary.pdf" ResultsSummaryTMP.pdf fitresultslatex.pdf')
                cd ..
            end
            
            fprintf('Saving PDFs done. \n')
            close all
        end
           
    end %Secret Methods (don't use)
    
    
    
end % class
