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
    
    f1 = figure('Units','normalized','Position',[-0.1,0.1,0.5,0.7]);
    s1 = subplot(3,1,1);
    % plot(qU.^2,SigmaSqRel_Stat.*1e2,'.:','MarkerSize',15);
    errorbar(qU.^2,TBDIS,TBDISE,'.','MarkerSize',15,'CapSize',0);
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    PrettyFigureFormat
    xlabel(sprintf('({\\itqU} - 18574)^2'))
    % ylabel('Rel. statistic uncertainty (%)')
    %  ylabel(sprintf('Rel. \\sigma_{stat.}^2 (%%)'));
    ylabel(sprintf('Counts'));
    %  ylim([0 0.6]);
    ylim([3*10^4 10^6]);
    leg = legend('KNM-2 integral tritium spectrum','Location','northwest');
    PrettyLegendFormat(leg);
    
    s2 = subplot(3,1,2);
%    SigmaSqStat = d.sin2t4_Stat(1,:).^2;   
%     % syst: constant with m4^2
%     SigmaSqSys_const = median(d.sin2t4_Sys(1,:).^2);
%     SigmaSqTot_const = SigmaSqStat+SigmaSqSys_const;
%     Ratio_const = SigmaSqSys_const./SigmaSqTot_const;
%     
%     %syst: linear increase with m4^2
%      SigmaSqSys_lin = d.mNu4Sq(1,:).*mean(d.sin2t4_Sys(1,:).^2)./max(d.mNu4Sq(1,:));
%      SigmaSqTot_lin = SigmaSqStat+SigmaSqSys_lin;
%      Ratio_lin = SigmaSqSys_lin./SigmaSqTot_lin;
   
    % plot 2
    x = d.mNu4Sq(1,:);%sqrt(d.mNu4Sq(1,:));
    plot(x,d.Ratio(1,:),'-k','LineWidth',2)
    set(gca,'XScale','log');
    xlabel(sprintf('{\\itm}_4^2 (eV^2)'));
    ylabel(sprintf('\\sigma_{syst.}^2 / \\sigma_{total}^2'));
    PrettyFigureFormat
     leg = legend('All syst. effects','Location','northwest');
     PrettyLegendFormat(leg);
     ylim([0 0.5])
     
    s3 = subplot(3,1,3);
    psys =  plot(x,d.sin2t4_Sys(1,:).^2,'-k','LineWidth',2);
    hold on;
    pstat=  plot(x,d.sin2t4_Stat(1,:).^2,'Color',rgb('DodgerBlue'),'LineWidth',2);
    ptot = plot(x,d.sin2t4_Tot(1,:).^2,'Color',rgb('Orange'),'LineWidth',2);
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    xlabel(sprintf('{\\itm}_4^2 (eV^2)'));
    ylabel(sprintf('\\sigma^2'));
     PrettyFigureFormat
     leg = legend([ptot,pstat,psys],'Total','Stat. only','Syst. only','Location','southwest');
     PrettyLegendFormat(leg);
     
     linkaxes([s1,s2,s3],'x');
%    %%
%    SigmaSqRel_Syst_const = 0.0028;
%    SigmaSqRel_Syst_lin = qU*0.0005;%.*1e-06;
%    
%    SigmaSqAbs_Syst_const = SigmaSqRel_Syst_const.*TBDIS;
%    SigmaSqAbs_Syst_lin = SigmaSqRel_Syst_lin.*TBDIS;
%    
%    SigmaSqAbs_Tot_const = SigmaSqAbs_Stat+SigmaSqAbs_Syst_const;
%     SigmaSqAbs_Tot_lin = SigmaSqAbs_Stat+SigmaSqAbs_Syst_lin;
%      
%      plot(d.mNu4Sq(1,:),d.Ratio(1,:))
%      hold on;
%      p1 =  plot(qU.^2,SigmaSqAbs_Syst_const./SigmaSqAbs_Tot_const);
%      hold on;
%      p2 =  plot(qU.^2,SigmaSqAbs_Syst_lin./SigmaSqAbs_Tot_lin);
% 
%     leg = legend(sprintf('\\sigma_{syst.}^2: constant'),sprintf('\\sigma_{syst.}^2: linear increase'));
%     PrettyLegendFormat(leg);
%     set(gca,'XScale','log');
%     xlabel(sprintf('({\\itqU} - 18574)^2'));
%    ylabel(sprintf('\\sigma_{syst.}^2 / \\sigma_{total}^2'));
%    
%    
%      