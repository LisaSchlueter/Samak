function RingwiseBestFit(varargin)
    p=inputParser;
    p.addParameter('Nrings',4,@(x)isfloat(x));
    p.addParameter('mnuSq',0,@(x)isfloat(x));
    p.addParameter('pullFlag',3);
    p.addParameter('DataType','Real',@(x)ismember(x,{'Twin','Real'}));
    p.addParameter('Syst','ON',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('RunList','KNM1',@(x)ischar(x));
    p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.parse(varargin{:});
    Nrings   = p.Results.Nrings;
    mnuSq    = p.Results.mnuSq;
    pullFlag = p.Results.pullFlag;
    DataType = p.Results.DataType;
    Syst     = p.Results.Syst;
    RunList  = p.Results.RunList;
    saveplot = p.Results.saveplot;

    
    if exist(sprintf('./EtaFitResult_PSR%.1f_real.mat',Nrings),'file')
        load(sprintf('./EtaFitResult_PSR%.1f_real.mat',Nrings));
    else
        if strcmp(Syst,'OFF')
            Chi2opt='chi2Stat';
            NP=1;
        else
            Chi2opt='chi2CMShape';
            NP=1.064;
        end
        D = MultiRunAnalysis('RunList',RunList,...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                    'chi2',Chi2opt,...              % uncertainties: statistical or stat + systematic uncertainties
                    'DataType',DataType,...                 % can be 'Real' or 'Twin' -> Monte Carlo
                    'fixPar','mNu E0 Norm Bkg',...        % free Parameter!!
                    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                    'NonPoissonScaleFactor',NP,...     % background uncertainty are enhanced
                    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                    'pullFlag',pullFlag,...
                    'FSDFlag','SibilleFull',...           % final state distribution
                    'ELossFlag','KatrinT2',...            % energy loss function
                    'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
                    'DopplerEffectFlag','FSD',...
                    'Twin_SameCDFlag','OFF',...
                    'Twin_SameIsotopFlag','OFF',...
                    'SynchrotronFlag','ON',...
                    'AngularTFFlag','OFF',...
                    'TwinBias_Q',18573.73,...
                    'TwinBias_mnuSq',mnuSq,...
                    'RingMerge','Full');
        A=RingAnalysis('RunAnaObj',D,'RingList',1:Nrings);
        
        M = MultiRunAnalysis('RunList',RunList,...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                    'chi2',Chi2opt,...              % uncertainties: statistical or stat + systematic uncertainties
                    'DataType',DataType,...                 % can be 'Real' or 'Twin' -> Monte Carlo
                    'fixPar','mNu E0 Norm Bkg eta',...    % free Parameter!!
                    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                    'NonPoissonScaleFactor',NP,...     % background uncertainty are enhanced
                    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                    'pullFlag',pullFlag,...
                    'FSDFlag','SibilleFull',...           % final state distribution
                    'ELossFlag','KatrinT2',...            % energy loss function
                    'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
                    'DopplerEffectFlag','FSD',...
                    'Twin_SameCDFlag','OFF',...
                    'Twin_SameIsotopFlag','OFF',...
                    'SynchrotronFlag','ON',...
                    'AngularTFFlag','OFF',...
                    'TwinBias_Q',18573.73,...
                    'TwinBias_mnuSq',mnuSq,...
                    'RingMerge','Full');
        B=RingAnalysis('RunAnaObj',M,'RingList',1:Nrings);

        fitresults = zeros(22,Nrings);
        for j=1:Nrings
            D=A.MultiObj(j);
            M=B.MultiObj(j);
            M.exclDataStart = M.GetexclDataStart(40);
            M.Fit;
            D.exclDataStart = D.GetexclDataStart(40);
            D.Fit;
            fitresults(1,j)= M.FitResult.par(1);
            fitresults(2,j)= M.FitResult.err(1);
            fitresults(3,j)= M.ModelObj.Q_i+M.FitResult.par(2);
            fitresults(4,j)= M.FitResult.err(2);
            fitresults(5,j)= M.ModelObj.BKG_RateSec_i+M.FitResult.par(3);
            fitresults(6,j)= M.FitResult.err(3);
            fitresults(7,j)= M.FitResult.par(3+M.ModelObj.nPixels:3+2*M.ModelObj.nPixels-1) + 1;
            fitresults(8,j)= M.FitResult.err(3+M.ModelObj.nPixels:3+2*M.ModelObj.nPixels-1);
            fitresults(9,j)= M.FitResult.par(17).*1e10;
            fitresults(10,j)=M.FitResult.err(17).*1e10;
            fitresults(11,j)=M.FitResult.chi2min;
            fitresults(12,j)=D.FitResult.par(1);
            fitresults(13,j)=D.FitResult.err(1);
            fitresults(14,j)=D.ModelObj.Q_i+M.FitResult.par(2);
            fitresults(15,j)=D.FitResult.err(2);
            fitresults(16,j)=D.ModelObj.BKG_RateSec_i+D.FitResult.par(3);
            fitresults(17,j)=D.FitResult.err(3);
            fitresults(18,j)=D.FitResult.par(3+D.ModelObj.nPixels:3+2*D.ModelObj.nPixels-1) + 1;
            fitresults(19,j)=D.FitResult.err(3+D.ModelObj.nPixels:3+2*D.ModelObj.nPixels-1);
            fitresults(20,j)=D.FitResult.par(17).*1e10;
            fitresults(21,j)=D.FitResult.err(17).*1e10;
            fitresults(22,j)=D.FitResult.chi2min;
            save(sprintf('./EtaFitResult_PSR%.1f_real.mat',Nrings),'fitresults');

            %% Global variables
            times = M.ModelObj.qUfrac*M.ModelObj.TimeSec;
            qU    = M.ModelObj.qU; qU    = qU-M.TwinBias_Q; % Energy axis

            % Spectrum no relics 
            IS  = D.ModelObj.TBDIS; 
            YIs = IS./times;
            DIS = D.RunData.TBDIS;
            DIS = DIS./times;

            % Spectrum relics
            YI = M.ModelObj.TBDIS; 
            YI = YI./times;

            % Error bar
            err  = (diag(sqrt(D.FitCMShape)));
            err  = err./times;
            err  = err./YI;

            %% Constraining everything to qULimiteV
            qULimit = -40;
            YIs=YIs(qU>qULimit);
            DIS=DIS(qU>qULimit);
            YI=YI(qU>qULimit);
            err=err(qU>qULimit);
            times=times(qU>qULimit);
            qU=qU(qU>qULimit);

            LocalFontSize = 20;

            fig = figure('Renderer','painters');
            set(fig, 'Units', 'normalized', 'Position', [0.001, 0.001,0.45, 0.8]);

            plot_title = sprintf('Relic neutrino model: \\eta = %.2g',D.ModelObj.eta);
            prlG = [81 126 102]/255;
            prlB = [50 148 216]/255;
            FitStyleArg = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',4,'MarkerEdgeColor',rgb('Black')};

            % Spectra  - First Subplot

            s1=subplot(4,1,[1 2]);

            % Plot

            pfit = plot(M.RunData.qU(M.exclDataStart:end)-M.TwinBias_Q,...
                M.RunData.TBDIS(M.exclDataStart:end)./...
                (M.ModelObj.qUfrac(M.exclDataStart:end)*M.ModelObj.TimeSec),...
                'DisplayName','No Relics','color',prlB,'LineWidth',3,'LineStyle','-')
            hold on

            pdata = errorbar(D.RunData.qU(D.exclDataStart:end)-M.TwinBias_Q,...
                D.RunData.TBDIS(D.exclDataStart:end)./...
                (D.ModelObj.qUfrac(D.exclDataStart:end)*D.ModelObj.TimeSec),...
                sqrt(D.RunData.TBDIS(D.exclDataStart:end))./(D.ModelObj.qUfrac(D.exclDataStart:end)*D.ModelObj.TimeSec)*20,FitStyleArg{:},'CapSize',0)

            yl1 = ylabel('Count rate (cps)');
            legend([pdata,pfit],{'Twin data with 1\sigma error bars \times 20','\beta-decay model'},'Location','northeast','box','off');
            lgd=legend;
            lgd.FontSize = LocalFontSize-2;

            xlim([min(qU-5) 10]);
            %ylim([0.007 2*max(YI_N)]);
            PRLFormat;
            set(gca,'FontSize',LocalFontSize);
            set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
            set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);
            set(gca, 'YScale', 'log');
            ylim([0.05 10])

            % Ratio    - Second Subplot
            s2=subplot(4,1,3);

            % Ratio
            RSP  = (YI./YIs);
            RSPd = RSP;

            % Plot
            hr1 = plot(qU,RSP,'color',rgb('Salmon'),'LineWidth',3,'LineStyle','-');
            hold on;
            hr2 = plot(qU,ones(1,numel(qU)),'color',prlB,'LineWidth',3,'LineStyle',':');
            hr3 = errorbar(qU,(DIS./YIs),err,FitStyleArg{:},'CapSize',0);
            yl2=ylabel('Ratio');
            katrinsim   = sprintf('\\eta=0');
            sterilemod  = sprintf('Best fit: \\eta=%.2g',M.ModelObj.eta);
            hl=legend([hr2 hr1],{katrinsim,sterilemod},'Location','southeast','box','off');
            hl.NumColumns=2;
            hl.FontSize = LocalFontSize-2;

            xlim([min(qU-5) 10]);
            %ylim([min((DIS./YI-1)./err) max((DIS./YI-1)./err)]);

            PRLFormat;
            set(gca,'FontSize',LocalFontSize);
            set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
            set(get(gca,'YLabel'),'FontSize',LocalFontSize+4); 
            hl.Position(2) = 0.333;
            ax2 = gca;

            % MTD      -  Third Subplot
            s3=subplot(4,1,4);

            bar(qU,times./(60*60*24),0.5,...
            'FaceColor',prlB,'EdgeColor',prlB);

            xlabel('Retarding energy - 18575 (eV)');
            ylh = ylabel('Time (d)');
            ylh.Position(1) = ax2.YLabel.Position(1)-4;%-range-6.8;
            yl1.Position(1) = ax2.YLabel.Position(1)-4;
            yl2.Position(1) = ax2.YLabel.Position(1)-4;
            PRLFormat;
            set(gca,'FontSize',LocalFontSize);
            set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
            set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);     

            linkaxes([s1,s2,s3],'x');
            hold off;
        end
    end
end