function PlotBestFit(varargin)
    p=inputParser;
    p.addParameter('mnuSq',0,@(x)isfloat(x));
    p.addParameter('Nfit',500,@(x)sfloat(x));
    p.parse(varargin{:});
    mnuSq = p.Results.mnuSq;
    Nfit  = p.Results.Nfit;

    
    if exist(sprintf('./EtaFitResult_AllParams_mnuSq%g_Nfit%g.mat',mnuSq,Nfit),'file')
        load(sprintf('./EtaFitResult_AllParams_mnuSq%g_Nfit%g.mat',mnuSq,Nfit));
    else
        D = MultiRunAnalysis('RunList','KNM1',...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                    'chi2','chi2Stat',...              % uncertainties: statistical or stat + systematic uncertainties
                    'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                    'fixPar','mNu E0 Norm Bkg',...        % free Parameter!!
                    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                    'NonPoissonScaleFactor',1,...     % background uncertainty are enhanced
                    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                    'FSDFlag','SibilleFull',...           % final state distribution
                    'ELossFlag','KatrinT2',...            % energy loss function
                    'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
                    'DopplerEffectFlag','FSD',...
                    'Twin_SameCDFlag','OFF',...
                    'Twin_SameIsotopFlag','OFF',...
                    'SynchrotronFlag','ON',...
                    'AngularTFFlag','OFF',...
                    'TwinBias_Q',18573.73,...
                    'TwinBias_mnuSq',mnuSq);

        fitresults = zeros(11,Nfit);
        for j=1:Nfit
            M = MultiRunAnalysis('RunList','KNM1',...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                        'chi2','chi2Stat',...              % uncertainties: statistical or stat + systematic uncertainties
                        'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                        'fixPar','mNu E0 Norm Bkg eta',...    % free Parameter!!
                        'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                        'NonPoissonScaleFactor',1,...     % background uncertainty are enhanced
                        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
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
                        'RelicPeakPosition','');


            M.RunData.TBDIS(end-5:end)=mean(M.RunData.TBDIS(end-5:end)./(M.RunData.qUfrac(end-5:end).*M.RunData.TimeSec)).*M.RunData.qUfrac(end-5:end).*M.RunData.TimeSec;
            statfluct = zeros(numel(D.RunData.qU),1);
            for i=1:numel(D.RunData.qU)
                gm=gmdistribution(D.RunData.TBDIS(i),D.RunData.TBDIS(i));
                statfluct(i) = random(gm)-D.RunData.TBDIS(i);
            end
            M.RunData.TBDIS = M.RunData.TBDIS+statfluct;
            M.exclDataStart = M.GetexclDataStart(40);
            M.Fit;
            fitresults(1,j)=M.FitResult.par(1);
            fitresults(2,j)=M.FitResult.err(1);
            fitresults(3,j)= M.ModelObj.Q_i+M.FitResult.par(2);
            fitresults(4,j)= M.FitResult.err(2);
            fitresults(5,j)= M.ModelObj.BKG_RateSec_i+M.FitResult.par(3);
            fitresults(6,j)= M.FitResult.err(3);
            fitresults(7,j)= M.FitResult.par(3+M.ModelObj.nPixels:3+2*M.ModelObj.nPixels-1) + 1;
            fitresults(8,j)= M.FitResult.err(3+M.ModelObj.nPixels:3+2*M.ModelObj.nPixels-1);
            fitresults(9,j)=M.FitResult.par(17).*1e10;
            fitresults(10,j)=M.FitResult.err(17).*1e10;
            fitresults(11,j)=M.FitResult.chi2min;
        end
        save(sprintf('./EtaFitResult_AllParams_mnuSq%g_Nfit%g.mat',mnuSq,Nfit),'fitresults');
        D.exclDataStart = D.GetexclDataStart(40);
        D.RunData.TBDIS = D.RunData.TBDIS+statfluct;
        D.Fit;

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
            sqrt(D.RunData.TBDIS(D.exclDataStart:end))./(D.ModelObj.qUfrac(D.exclDataStart:end)*D.ModelObj.TimeSec)*50,FitStyleArg{:},'CapSize',0)

        yl1 = ylabel('Count rate (cps)');
        legend([pdata,pfit],{'Twin data with 1\sigma error bars \times 50','\beta-decay model'},'Location','northeast','box','off');
        lgd=legend;
        lgd.FontSize = LocalFontSize-2;

        xlim([min(qU-5) 10]);
        %ylim([0.007 2*max(YI_N)]);
        PRLFormat;
        set(gca,'FontSize',LocalFontSize);
        set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
        set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);
        set(gca, 'YScale', 'log');
        ylim([0.2 25])

        % Ratio    - Second Subplot
        s2=subplot(4,1,3);

        % Ratio
        RSP  = (YI./YIs);
        RSPd = RSP;

        % Plot
        hr1 = plot(qU,RSP,'color',rgb('Salmon'),'LineWidth',3,'LineStyle','-');
        hold on;
        hr2 = plot(qU,ones(1,numel(qU)),'color',prlB,'LineWidth',3,'LineStyle',':');
        hr3 = errorbar(qU,(YI./DIS),err,FitStyleArg{:},'CapSize',0);
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

%     for i=1:983
%         if fitresults(1,i)<0.6
%             fitresults(:,i)=[];
%         end
%     end
    fig2 = figure('Renderer','painters');
    set(fig2, 'Units', 'normalized', 'Position', [0.001, 0.001,0.45, 0.8]);
    PlotStyle = { 'o','MarkerSize',2,'MarkerFaceColor',rgb('IndianRed'),...%rgb('IndianRed'),...%
                'LineWidth',2,'Color',rgb('IndianRed')};
    sub1=subplot(4,4,1);
    plot(fitresults(1,:),fitresults(3,:)-18573.73,PlotStyle{:});
    ylabel('\Delta E_{0} (eV)');
    set(gca,'xtick',[]);
    PrettyFigureFormat;
    sub2=subplot(4,4,5);
    plot(fitresults(1,:),fitresults(7,:)-1,PlotStyle{:});
    ylabel('\Delta N');
    set(gca,'xtick',[]);
    PrettyFigureFormat;
    sub2=subplot(4,4,6);
    plot(fitresults(3,:)-18573.73,fitresults(7,:)-1,PlotStyle{:});
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    PrettyFigureFormat;
    sub2=subplot(4,4,9);
    plot(fitresults(1,:),fitresults(5,:),PlotStyle{:});
    ylabel('Bkg (cps)');
    set(gca,'xtick',[]);
    PrettyFigureFormat;
    sub2=subplot(4,4,10);
    plot(fitresults(3,:)-18573.73,fitresults(5,:),PlotStyle{:});
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    PrettyFigureFormat;
    sub2=subplot(4,4,11);
    plot(fitresults(7,:)-1,fitresults(5,:),PlotStyle{:});
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    PrettyFigureFormat;
    sub2=subplot(4,4,13);
    plot(fitresults(1,:),fitresults(9,:),PlotStyle{:});
    ylabel('\eta');
    xlabel('m_{\nu}^2 (eV^2)');
    PrettyFigureFormat;
    sub2=subplot(4,4,14);
    plot(fitresults(3,:)-18573.73,fitresults(9,:),PlotStyle{:});
    xlabel('\Delta E_0 (eV)');
    set(gca,'ytick',[]);
    PrettyFigureFormat;
    sub2=subplot(4,4,15);
    plot(fitresults(7,:)-1,fitresults(9,:),PlotStyle{:});
    xlabel('\Delta N');
    set(gca,'ytick',[]);
    PrettyFigureFormat;
    sub2=subplot(4,4,16);
    plot(fitresults(5,:),fitresults(9,:),PlotStyle{:});
    xlabel('Bkg (cps)');
    set(gca,'ytick',[]);
    PrettyFigureFormat;
    sub2=subplot(4,4,[3 4 7 8]);
    histfit(fitresults(11,:));
    xlabel('\chi^2');
    PrettyFigureFormat;
end


% mnuSqVec=linspace(0,2,6);
% etaBestFitVec=mnuSqVec;
% mnufitvec=mnuSqVec;
% mnufitveceta=mnuSqVec;
% e0fitvec=mnuSqVec;
% e0fitveceta=mnuSqVec;
% NormFitVec=mnuSqVec;
% normfitveceta=mnuSqVec;
% bkgfitvec=mnuSqVec;
% bkgfitveceta=mnuSqVec;
% for i=1:numel(mnuSqVec)
%     R = MultiRunAnalysis('RunList','KNM1',...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
%             'chi2','chi2CMShape',...              % uncertainties: statistical or stat + systematic uncertainties
%             'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
%             'fixPar','mNu E0 Norm Bkg eta',...    % free Parameter!!
%             'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
%             'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
%             'minuitOpt','min ; minos',...         % technical fitting options (minuit)
%             'FSDFlag','SibilleFull',...           % final state distribution
%             'ELossFlag','KatrinT2',...            % energy loss function
%             'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
%             'DopplerEffectFlag','FSD',...
%             'Twin_SameCDFlag','OFF',...
%             'Twin_SameIsotopFlag','OFF',...
%             'SynchrotronFlag','ON',...
%             'AngularTFFlag','OFF',...
%             'TwinBias_Q',18573.73,...
%             'TwinBias_mnuSq',mnuSqVec(i));
%         D = MultiRunAnalysis('RunList','KNM1',...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
%             'chi2','chi2CMShape',...              % uncertainties: statistical or stat + systematic uncertainties
%             'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
%             'fixPar','mNu E0 Norm Bkg',...        % free Parameter!!
%             'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
%             'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
%             'minuitOpt','min ; minos',...         % technical fitting options (minuit)
%             'FSDFlag','SibilleFull',...           % final state distribution
%             'ELossFlag','KatrinT2',...            % energy loss function
%             'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
%             'DopplerEffectFlag','FSD',...
%             'Twin_SameCDFlag','OFF',...
%             'Twin_SameIsotopFlag','OFF',...
%             'SynchrotronFlag','ON',...
%             'AngularTFFlag','OFF',...
%             'TwinBias_Q',18573.73,...
%             'TwinBias_mnuSq',mnuSqVec(i));
%         R.exclDataStart = R.GetexclDataStart(40);
%         D.exclDataStart = D.GetexclDataStart(40);
%         R.Fit;
%         D.Fit;
%         etaBestFitVec(i)=R.FitResult.par(17).*1e10;
%         mnufitvec(i)=D.FitResult.par(1);
%         mnufitveceta(i)=R.FitResult.par(1);
%         e0fitvec(i)=D.ModelObj.Q_i+D.FitResult.par(2);
%         e0fitveceta(i)=R.ModelObj.Q_i+R.FitResult.par(2);
%         NormFitVec(i)=D.FitResult.par(3+D.ModelObj.nPixels:3+2*D.ModelObj.nPixels-1) + 1;
%         normfitveceta(i)=R.FitResult.par(3+R.ModelObj.nPixels:3+2*R.ModelObj.nPixels-1) + 1;
%         bkgfitvec(i)=D.ModelObj.BKG_RateSec_i+D.FitResult.par(3);
%         bkgfitveceta(i)=R.ModelObj.BKG_RateSec_i+R.FitResult.par(3);
% end
% save('./RelicNuBkg/Misc/BestFitEta_KNM1.mat','etaBestFitVec','mnufitvec','mnufitveceta','e0fitvec','e0fitveceta','NormFitVec','normfitveceta','bkgfitvec','bkgfitveceta');

% fig2=figure(2);
% load('./EtaFitResultsFluctuated.mat');
% fitresults=fitresults(1:312);
% pdG = fitdist(fitresults','Normal');
% pdfG = @(b) 1/sqrt(2*pi*pdG.sigma.^2) .* exp(-(b-pdG.mu).^2/(2*pdG.sigma.^2));
% b=linspace(min(fitresults),max(fitresults),100);
% histogram(fitresults,20,'Normalization','pdf');
% hold on;
% plot(b,pdfG(b),'Color',rgb('IndianRed'),'LineWidth',4);