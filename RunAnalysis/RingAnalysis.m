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
                'MosCorrFlag',obj.RunAnaObj.MosCorrFlag,...
                'TwinBias_Q',obj.RunAnaObj.TwinBias_Q,...
                'SynchrotronFlag',obj.RunAnaObj.SynchrotronFlag,'AngularTFFlag',obj.RunAnaObj.AngularTFFlag,...
                'FSD_Sigma',obj.RunAnaObj.FSD_Sigma};
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
                RunList = strrep(strrep(obj.RunAnaObj.StackFileName,obj.RunAnaObj.SetTwinOrFakeFileName,''),'Twin','');
                obj.MultiObj = arrayfun(@(x) MultiRunAnalysis(CommonArg{:},'RunList',RunList,...
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
            savename = [savedir,sprintf('SingleRingFitResult_%s_%s_%.0fruns_fix%s_%s_%.0feVrange%s%s.mat',...
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
                
                progressbar('Fit single rings');
                for i=1:obj.nRings
                    progressbar(i/obj.nRings);
                    obj.MultiObj(i).Fit;
                    obj.MultiObj(i).Fit; % irregularities observed... repeat fit to double check
                    %                     if obj.MultiObj(i).FitResult.err(1)<0.05
                    %                         % sometime fit doesn work properly -> nu-mass error tiny, reason unknown
                    %                         % repeat fit often works
                    %                         obj.MultiObj(i).Fit;
                    %                     end
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
                
                FitResult.par = par;
                FitResult.err = err;
                FitResult.chi2min = chi2min;
                FitResult.dof = dof;
                FitResult.ScanResults = ScanResults;
                FitResult.errNeg = errNeg;
                FitResult.errPos = errPos;
                
                obj.FitResult = FitResult;
                
                if strcmp(SaveResult,'ON')
                    save(savename,'FitResult');
                end
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
            
            if isfield(obj.FitResult,'errNeg') && PlotPar==1
                fprintf('switch to asymmetric errors for PlotPar %.0f \n',PlotPar)
                yErr = 0.5*(abs(obj.FitResult.errNeg(:,PlotPar))+obj.FitResult.errPos(:,PlotPar));
                %                 if any(yErr<1e-2)
                %                     Index = yErr<1e-2;
                %                     meanYerr = mean(yErr(~Index));
                %                      yErr(Index) = meanYerr;
                %                 end
            elseif PlotPar==2
                y = y+obj.RunAnaObj.ModelObj.Q_i;
            elseif PlotPar==3
                BKG_i = cell2mat(arrayfun(@(x) x.ModelObj.BKG_RateSec_i,obj.MultiObj,'UniformOutput',0));
                nPixperPseudoRing = cell2mat(arrayfun(@(x) numel(x.PixList),obj.MultiObj,'UniformOutput',0)); %number of pixels per ring
                y = 1e3.*(y+BKG_i)./nPixperPseudoRing; % background per pixel within each pseudo-ring (mcps)
                yErr = yErr./nPixperPseudoRing.*1e3; % in mcps
                BKG_Tot = sum(y.*nPixperPseudoRing);
            elseif PlotPar==4
                y = 1+y;
            end
            
            %             if any(yErr<1e-2) && PlotPar==2
            %                 % temporary solution for plot:
            %                 % scale weird error up to average error (for lin fit)
            %                 Index = yErr<1e-2;
            %                 meanYerr = mean(yErr(~Index));
            %                 yErr(Index) = meanYerr;
            %             end
            
            if strcmp(PlotMode,'Rel')
                meanPar = wmean(y,1./yErr.^2);
            elseif strcmp(PlotMode,'Abs')
                meanPar = 0;
            end
            
            %% linear fit
            if strcmp(linFitFlag,'ON')
                [linFitpar, linFiterr, linFitchi2min,linFitdof] =...
                    linFit(obj.RingList',y-meanPar,yErr);
                %                 % tmp
                %                 corrcoeff = -0.94872; % from minuit terminal output...
                %                 linFitCovMat = [1,corrcoeff;corrcoeff,1].*linFiterr.*linFiterr';
                %                 linFitParSamples = mvnrnd(linFitpar,linFitCovMat,1e3)';
                %                 linFitSamples = linFitParSamples(1,:).*obj.RingList'-linFitParSamples(2,:);
                %                 linFitStd = std(linFitSamples,0,2);
            end
            %% label
            range = round(abs(obj.RunAnaObj.ModelObj.qU(obj.RunAnaObj.exclDataStart)-obj.RunAnaObj.ModelObj.Q_i),0);
            %% plot fit result
            fig2 = figure('Renderer','painters');
            set(fig2,'units','normalized','pos',[0.1, 0.1,0.5,0.5]);
            if strcmp(PlotMode,'Rel')
                pm = plot(linspace(0.5,12.5,obj.nRings),zeros(obj.nRings,1),':','Color',rgb('SlateGrey'),'LineWidth',2);
            else
                pm = plot(linspace(0.5,12.5,obj.nRings),wmean(y,1./yErr.^2).*ones(obj.nRings,1),':','Color',rgb('Silver'),'LineWidth',2);
            end
            hold on;
            
            
            if strcmp(linFitFlag,'ON')
                l = plot(obj.RingList,(linFitpar(1).*obj.RingList+linFitpar(2)),'-','Color',obj.RunAnaObj.PlotColor,'LineWidth',2);
                %                 ylinFit = (linFitpar(1).*obj.RingList+linFitpar(2))./ScaleFactor';
                %                 ylinErrFit =((linFitpar(1)+linFiterr(1)).*obj.RingList+linFitpar(2))./ScaleFactor';
                %                 [l,a] =  boundedline(obj.RingList,ylinFit,ylinErrFit);
                %                 l.Color = obj.RunAnaObj.PlotColor; l.LineWidth=2;
                %                 a.FaceColor = rgb('PowderBlue'); a.FaceAlpha = 0.5;
            end
            
            %  if PlotPar == 3
            %                 e = errorbar(obj.RingList,1e3.*((y+BKG_i')./ScaleFactor'),1e3.*yErr,'.',...
            %                 'LineWidth',2,'LineStyle','none','MarkerSize',22,...
            %                 'Color',rgb('Orange'),'MarkerFaceColor',rgb('Orange'));
            %'Color',obj.RunAnaObj.PlotColor,'MarkerFaceColor',obj.RunAnaObj.PlotColor);
            
            %             else
            e = errorbar(obj.RingList,y-meanPar,yErr,'.',...
                'LineWidth',2,'LineStyle','none','MarkerSize',22,...
                'Color',obj.RunAnaObj.PlotColor,'MarkerFaceColor',obj.RunAnaObj.PlotColor);
            %  'Color',rgb('Orange'),'MarkerFaceColor',rgb('Orange'));
            %           end
            e.CapSize = 0;
            PrettyFigureFormat('FontSize',25);
            
            if strcmp(obj.RunAnaObj.RunData.RunName,'KNM2_Prompt')
                obj.RunAnaObj.RunData.RunName = 'KNM2 data';
            end
            if strcmp(linFitFlag,'ON') %&& PlotPar ~= 3
                runname = strrep(obj.RunAnaObj.RunData.RunName,'_',' ');
                if contains(runname,'afterFix')
                    runname = 'KNM2 after RW fix I';
                elseif contains(runname,'RW3')
                    runname = 'KNM2 after RW fix II';
                end
                
                %                 if linFitpar(1)<0.03
                %                     if PlotPar==1
                %                         linFitleg =  sprintf('Linear fit %.1f \\pm %.1f meV^2/ring , p-value = %.2f',linFitpar(1)*1e3,linFiterr(1)*1e3,1-chi2cdf(linFitchi2min,linFitdof));
                %                     elseif PlotPar==2
                %                         linFitleg =  sprintf('Linear fit %.1f \\pm %.1f meV/ring , p-value = %.2f ',linFitpar(1)*1e3,linFiterr(1)*1e3,1-chi2cdf(linFitchi2min,linFitdof));
                %                     end
                %                 else
                if PlotPar==1
                    linFitleg =  sprintf('Linear slope: (%.2f \\pm %.2f) eV^{ 2} per pseudo-ring ',linFitpar(1),linFiterr(1));%,1-chi2cdf(linFitchi2min,linFitdof));%, {\\it p} = %.2f
                    
                elseif PlotPar==2
                    if round(linFitpar(1),2) <0.1 || round(linFiterr(1),2)<0.1
                        linFitleg =  sprintf('Linear slope: (%.0f \\pm %.0f) meV per pseudo-ring',1e3.*linFitpar(1),1e3.*linFiterr(1));%,1-chi2cdf(linFitchi2min,linFitdof));
                    else
                        linFitleg =  sprintf('Linear slope: (%.2f \\pm %.2f) eV per pseudo-ring',linFitpar(1),linFiterr(1));%,1-chi2cdf(linFitchi2min,linFitdof));
                    end
                elseif PlotPar==3
                    if round(linFiterr(1),2)>=0.01
                        linFitleg =  sprintf('Linear slope: (%.2f \\pm %.2f) mcps per pseudo-ring',linFitpar(1),linFiterr(1));%,1-chi2cdf(linFitchi2min,linFitdof));
                    else
                        linFitleg =  sprintf('Linear slope: (%.3f \\pm %.3f) mcps per pseudo-ring',linFitpar(1),linFiterr(1));%,1-chi2cdf(linFitchi2min,linFitdof));
                        
                    end
                    
                elseif PlotPar==4
                    linFitleg =  sprintf('Linear slope: (%.1g \\pm %.1g) per pseudo-ring',linFitpar(1),linFiterr(1));%,1-chi2cdf(linFitchi2min,linFitdof));
                end
                %                 end
                leg = legend([e,l],sprintf('Pseudo-ring-wise fits'),linFitleg);% (%s),runname
            else
                leg = legend(e,sprintf('%s',strrep(obj.RunAnaObj.RunData.RunName,'_',' ')));
            end
            hold off;
            legend box off;
            leg.FontSize= get(gca,'FontSize');
            PrettyLegendFormat(leg);
            % leg.EdgeColor = rgb('Silver');
            if PlotPar==1
                ystr = sprintf('{\\itm}_\\nu^2');
                yUnit = sprintf('(eV^{ 2})');
            elseif PlotPar==2
                ystr = sprintf('{\\itE}_0^{fit}');
                yUnit = sprintf('(eV)');
            elseif PlotPar==3
                ystr = sprintf('{\\itB}_{base}');
                yUnit = sprintf(' per pixel (mcps)');
            elseif  PlotPar==11
                ystr = sprintf('{\\itqU}_{offset}');
                yUnit = sprintf('(eV)');
            elseif  PlotPar==4
                ystr = sprintf('{\\itN}');
                yUnit = '';
            end
            if strcmp(PlotMode,'Rel')
                ylabel(sprintf('%s - \\langle%s\\rangle %s',ystr,ystr,yUnit));
            else
                ylabel(sprintf('%s %s',ystr,yUnit));
            end
            xlabel('Detector rings')
            ax = gca;
            ax.YAxis.Exponent = 0;
            meanPar =wmean(y,1./yErr.^2);
            if PlotPar==1
                AnaType = 'm2';
                t =title(sprintf('%.0f stacked runs (%.0f - %.0f) - %.0f pixels - %.0f eV range - \\langle{\\itm}_\\nu^2\\rangle = %.2f eV^2 ',...
                    numel(obj.RunAnaObj.StackedRuns),obj.RunAnaObj.RunList(1),obj.RunAnaObj.RunList(end),numel(obj.RunAnaObj.PixList),range, mean(y)));
            elseif PlotPar==2
                t = title(sprintf('%.0f stacked runs (%.0f - %.0f) - %.0f pixels - %.0f eV range - \\langle{\\itE}_0^{fit}\\rangle = %.2f eV',...
                    numel(obj.RunAnaObj.StackedRuns),obj.RunAnaObj.RunList(1),obj.RunAnaObj.RunList(end),numel(obj.RunAnaObj.PixList),range, meanPar + obj.RunAnaObj.ModelObj.Q_i));
                AnaType = 'E0';
            elseif PlotPar==3
                t= title(sprintf('%.0f stacked runs (%.0f - %.0f) - %.0f pixels - %.0f eV range - B_{Tot} = %.1f mcps',...
                    numel(obj.RunAnaObj.StackedRuns),obj.RunAnaObj.RunList(1),obj.RunAnaObj.RunList(end),numel(obj.RunAnaObj.PixList),range, 1e3*BKG_Tot));
                AnaType = 'BKG';
            elseif PlotPar==4
                AnaType = 'Normalization';
                t =title(sprintf('%.0f stacked runs (%.0f - %.0f) - %.0f pixels - %.0f eV range - <N> = %.3f',...
                    numel(obj.RunAnaObj.StackedRuns),obj.RunAnaObj.RunList(1),obj.RunAnaObj.RunList(end),numel(obj.RunAnaObj.PixList), range, meanPar));
            elseif PlotPar==11
                AnaType = 'qUoffset';
                t =title(sprintf('%.0f stacked runs (%.0f - %.0f) - %.0f pixels - %.0f eV range - <qU_{offset}> = %.3f',...
                    numel(obj.RunAnaObj.StackedRuns),obj.RunAnaObj.RunList(1),obj.RunAnaObj.RunList(end),numel(obj.RunAnaObj.PixList),range, meanPar));
                
            end
            
            if strcmp(Blind,'ON')
                t = title(sprintf('%.0f stacked runs (%.0f - %.0f) - %.0f pixels - %.0f eV range',...
                    numel(obj.RunAnaObj.StackedRuns),obj.RunAnaObj.RunList(1),obj.RunAnaObj.RunList(end),numel(obj.RunAnaObj.PixList),range));
            end
            t.FontWeight = 'normal';
            t.FontSize = get(gca,'FontSize');
            t.delete;
            
            switch obj.RunAnaObj.RingMerge
                case 'Half'
                    xticks([1 2])
                    %  xticklabels({'1-6','7-12'})
                    xticklabels({'Bullseye & 1-5','6-11'});
                case 'InnerPseudo3'
                    xticks([1]);
                    %xticklabels({'1,2,3,4,5,6,7,8,9'})
                    xticklabels({'Bullseye',1,2,3,4,5,6,7,8});
                case 'Full'
                    xticks([1 2 3 4]);
                    %xticklabels({'1-3','4-6','7-9','10-12'})
                    xticklabels({'Bullseye & 1-2','3-5','6-8','9-11'})
                case 'Default'
                    xticks([1:1:10]);
                    %xticklabels({'1-2',3,4,5,6,7,8,9,10,'11-12'})
                    xticklabels({'Bullseye & 1',2,3,4,5,6,7,8,9,'10-11'});
                case 'None'
                    xticks([1:12])
                    xticklabels({'Bullseye',1,2,3,4,5,6,7,8,9,10,11});
                case 'Slice'
                    xticks([1:28])
                    %xticklabels({'Bullseye',1,2,3,4,5,6,7,8,9,10,11});
            end
            set(gca,'XMinorTick','off');
            
            if PlotPar==2 && strcmp(PlotMode,'Abs')
                %yticks([18572.1:0.2:18575])
            elseif PlotPar==4 && strcmp(PlotMode,'Abs')
                
            end
            
            if ~isempty(YLim)
                ylim([min(YLim),max(YLim)]);
            elseif strcmp(PlotMode,'Rel')
                ymin = min((y-meanPar)-yErr);
                ymax = max((y-meanPar)+yErr);
                ylim([1.1*ymin,1.7*ymax]);
            else
                ymin = min((y)-yErr);
                ymax = max((y)+yErr);
                if ymin<0
                    ylim([1.1*ymin,1.5*ymax]);
                else
                    ylim([0.5*ymin,1.5*ymax]);
                end
            end
            
            if obj.nRings<10
                xlim([min(obj.RingList)-0.2,max(obj.RingList)+0.2])
            else
                xlim([min(obj.RingList)-0.5,max(obj.RingList)+0.5])
            end
            leg.Location = 'northwest';
            % set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.5]));
            
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
                export_fig(gcf,[savedir,savename]);
                fprintf('Save plot to %s \n',[savedir,savename])
                
                if strcmp(linFitFlag,'ON')
                    linFitp = 1-chi2cdf(linFitchi2min,linFitdof);
                    linFitData  =  [linFitpar,linFitchi2min,linFitp;linFiterr,linFitdof,0];
                    linFitfile = strrep([savedir,savename],'.pdf','_LinFitResult');
                    Write2Txt('filename',linFitfile,'Format','txt','variable',linFitData',...
                        'variableName','slope   offset   chi2min/dof  pvalue','nCol',size(linFitData,2));
                end
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