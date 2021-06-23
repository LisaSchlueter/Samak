 %% understand stat over syst plot
 % Non poisson 
%% configure RunAnalysis object
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar','E0 Norm Bkg',...%free par
        'SysBudget',40,...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'RadiativeFlag','ON',...
        'FSDFlag','KNM2_0p5eV',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'NonPoissonScaleFactor',1,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'TwinBias_Q',18573.7,...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',3*1e-06,...
        'TwinBias_BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
    A = MultiRunAnalysis(RunAnaArg{:});
    %%
    savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
    savename = sprintf('%sksn2_SystBreakdown_StatOverSyst_%s_%.0feV_RasterScan%s.mat',savedir,'Twin',40,'ON');
    if exist(savename,'file')
        d = importdata(savename);
        fprintf('load file %s \n',savename);
    else
        return
        
    end
    
    %% relatitive statistical uncertainty
   
    Time = A.ModelObj.qUfrac(11:end-5).*A.ModelObj.TimeSec;
    qU = abs(A.ModelObj.qU(11:end-5)-18574);
    TBDIS = A.ModelObj.TBDIS(11:end-5);
    TBDISE = sqrt(A.ModelObj.TBDIS(11:end-5));
    
    SigmaSqRel_Stat = 1./TBDISE;
    SigmaSqAbs_Stat = SigmaSqRel_Stat.*TBDIS;
    
    %% siganl to background ratio
    A.InitModelObj_Norm_BKG;
    TimeSubRun   = A.ModelObj.qUfrac.*A.ModelObj.TimeSec;
    TimeAvSubrun = TimeSubRun./A.nRuns;
    BkgRate_PngSlope = 0.5.*A.ModelObj.BKG_PtSlope.*TimeAvSubrun;
    Bkg_PtSlope      = BkgRate_PngSlope.*TimeSubRun;
    BkgCounts = (BkgRate_PngSlope+A.ModelObj.BKG_RateSec).*TimeSubRun;
    TBDIS_Signal = A.ModelObj.TBDIS-BkgCounts;
    TBDIS_Signal(TBDIS_Signal<0) = 0;
    %%
    
    f1 = figure('Units','normalized','Position',[-0.1,0.1,0.5,0.7]);  
    s1 = subplot(3,1,1);
    
    % plot 1
    Ratio = d.sin2t4_Sys(4,:).^2./d.sin2t4_Tot(1,:).^2;
    x = -sqrt(d.mNu4Sq(1,:));
    plot(x,Ratio,'-k','LineWidth',2)
    xlabel(sprintf('- {\\itm}_4 (eV)'));
    ylabel(sprintf('\\sigma_{syst.}^2 / \\sigma_{total}^2'));
    PrettyFigureFormat
    leg = legend('Syst. effect: Non-Poisson background','Location','northwest');
    PrettyLegendFormat(leg);
    ylim([0 0.15])
    ax1 = gca;
   
    s2 = subplot(3,1,2);
   e1 = errorbar(A.ModelObj.qU-18574,BkgCounts,sqrt(BkgCounts),...
        '.-','Color',rgb('DodgerBlue'),'LineWidth',2,'MarkerSize',15,'CapSize',1);
    PrettyFigureFormat
    xlabel(sprintf('{\\itqU} - 18574 (eV)'))
    ylabel(sprintf('{\\itB} counts')); 
    ax2 = gca;
     
    s3 = subplot(3,1,3);
     bar(A.ModelObj.qU-18574,(A.ModelObj.qUfrac.*A.ModelObj.TimeSec)./(60*60));
    PrettyFigureFormat
    xlabel(sprintf('{\\itqU} - 18574 (eV^2)'))
    ylabel(sprintf('Time (hours)'));
     ylim([0 70]);
    ax3 = gca;

    linkaxes([s1,s2,s3],'x');
    xlim([-37.5,0]);
    
    ax1.Position(4) = 0.2;
    ax2.Position(4) = 0.2;
    ax3.Position(4) = 0.2;
    ax1.YLabel.Position(1) = -40.5;
    ax2.YLabel.Position(1) =   ax1.YLabel.Position(1);
    ax3.YLabel.Position(1) =   ax1.YLabel.Position(1);
 
    pltname = strrep(strrep(savename,'results','plots'),'.mat','_NonPois.png');
    print(gcf,pltname,'-dpng','-r350')