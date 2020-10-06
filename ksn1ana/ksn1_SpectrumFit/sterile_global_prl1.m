% spectrum and fit for sterile neutrino

% set path
savedirFile = [getenv('SamakPath'),'ksn1ana/ksn1_SpectrumFit/results/'];
savedirPlot = [getenv('SamakPath'),'ksn1ana/ksn1_SpectrumFit/plots/'];
MakeDir(savedirFile);
MakeDir(savedirPlot);

%% Settings
E0 = 18573.7;                                         % Endpoint in eV for display
sterile_mass  = 10;                                   % Sterile neutrino mass in eV
mixing_angle = 0.01;                                  % sin2(th4)
range = 40;
ErrorBarScaling = 50;
SavePlot = 'ON';
Subplot2 = 'Ratio'; % Ratio or Residuals
%% Data
datafile = sprintf('%sksn1_sterile_global_datafit_%.0feV.mat',savedirFile,range);
if exist(datafile,'file') 
    d = importdata(datafile);
    qU                 = d.qU;
    RateData           = d.RateData;
    RateErrData        = d.RateErrData;
    RateFitData3mNu    = d.RateFitData3mNu;
    RateFitErrData3mNu = d.RateFitErrData3mNu;
    TimeSubrun         = d.TimeSubrun;
    FitCMShape         = d.FitCMShape;
else
    RunAnaArg = {'RunList','KNM1',...
        'chi2','chi2CMShape',...
        'DataType','Real',...
        'fixPar','mNu E0 Bkg Norm',...
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',1.064,...
        'minuitOpt','min ; migrad',...
        'FSDFlag','SibilleFull',...
        'ELossFlag','KatrinT2',...
        'SysBudget',22,...
        'SynchrotronFlag','OFF',...
        'AngularTFFlag','OFF'};
    
    D = MultiRunAnalysis(RunAnaArg{:});
    
    D.exclDataStart = D.GetexclDataStart(range);
    D.ComputeCM('BkgMode','SlopeFit','nTrials',1000);%('SysEffect',struct('FSD','ON','TASR','ON','Stacking','ON'),'BkgCM','ON');
    D.Fit; 
    qU = D.RunData.qU(D.exclDataStart:end)-E0;
    TimeSubrun = (D.ModelObj.qUfrac(D.exclDataStart:end)*D.ModelObj.TimeSec);
    RateData = D.RunData.TBDIS(D.exclDataStart:end)./TimeSubrun;
    RateErrData = sqrt(D.RunData.TBDIS(D.exclDataStart:end))./TimeSubrun;
    
    RateFitData3mNu    = D.ModelObj.TBDIS(D.exclDataStart:end)./TimeSubrun;
    RateFitErrData3mNu = diag(sqrt(FitCMShape(D.exclDataStart:end,D.exclDataStart:end)))./TimeSubrun;
    %sqrt(D.ModelObj.TBDIS(D.exclDataStart:end))./TimeSubrun; 
    FitCMShape = D.FitCMShape;
    save(datafile,'qU','RateData','RateErrData','RateFitData3mNu','RateFitErrData3mNu','TimeSubrun','RunAnaArg','FitCMShape');
end
%% Simulation
simfile = sprintf('%sksn1_sterile_global_sim_m4%.1feV_sint4Sq%.2f_%.0feV.mat',savedirFile,sterile_mass,mixing_angle,range);
if exist(simfile,'file')
    d = importdata(simfile);
    RateSim3mNu    = d.RateSim3mNu;
    RateErrSim3mNu = d.RateErrSim3mNu;
    RateSim4mNu    = d.RateSim4mNu;
    RateErrSim4mNu = d.RateErrSim4mNu;
else
    RunAnaArg = {'RunList','KNM1',...
        'chi2','chi2Stat',...
        'DataType','Twin',...
        'fixPar','Norm',...
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',1.064,...
        'minuitOpt','min ; migrad',...
        'FSDFlag','SibilleFull',...
        'ELossFlag','KatrinT2',...
        'SysBudget',22,...
        'exclDataStart',13,...
        'SynchrotronFlag','OFF',...
        'AngularTFFlag','OFF'};
    R = MultiRunAnalysis(RunAnaArg{:});
    
    R.exclDataStart = R.GetexclDataStart(range);
    
    % 3 neutrino model
    RateSim3mNu     = R.ModelObj.TBDIS(D.exclDataStart:end)./TimeSubrun;
    RateErrSim3mNu  = sqrt(R.ModelObj.TBDIS(D.exclDataStart:end))./TimeSubrun;
    
    % 3+1 neutrino model
    R.ModelObj.SetFitBiasSterile(sterile_mass^2,mixing_angle);
    R.ModelObj.ComputeTBDDS;
    R.ModelObj.ComputeTBDIS;
    
    RateSim4mNu     = R.ModelObj.TBDIS(D.exclDataStart:end)./TimeSubrun;
    RateErrSim4mNu  = sqrt(R.ModelObj.TBDIS(D.exclDataStart:end))./TimeSubrun;
    save(simfile,'RateSim3mNu','RateErrSim3mNu','RateSim4mNu','RateErrSim4mNu','RunAnaArg');
end

%% ===== PLOTTING =====
LocalFontSize = 20;

fig = figure('Renderer','painters');
set(fig, 'Units', 'normalized', 'Position', [0.001, 0.001,0.45, 0.8]);

plot_title = sprintf('3+1 model: m_{4} = %.1f eV   |U_{e4}|^2 = %.2f',sterile_mass,mixing_angle);
prlG = [81 126 102]/255;
prlB = [50 148 216]/255;
FitStyleArg = {'Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerEdgeColor',rgb('Black')};

% Spectra  - First Subplot
s1=subplot(4,1,[1 2]);
pfit = plot(qU,RateFitData3mNu,'color',prlB,'LineWidth',3,'LineStyle','-');
hold on
pdata = errorbar(qU,RateData,RateErrData*ErrorBarScaling,'.','MarkerSize',15,FitStyleArg{:},'CapSize',0);

% apperance: legend, labels etc.
yl1 = ylabel('Count rate (cps)');
legend([pdata,pfit],{sprintf('KATRIN data with 1 \\sigma error bars \\times 50'),'3-\nu model'},'Location','northeast','box','off');
lgd=legend;
lgd.FontSize = LocalFontSize-2;
xlim([min(qU-5) max(qU+5)]);
ylim([0.18 2*max(RateData+RateErrData)]);
PRLFormat;
set(gca,'FontSize',LocalFontSize);
set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);
set(gca, 'YScale', 'log');
yticks([1 10 100])
ylim([0.15 45])
t1 = text(-41.243,31,'a)','FontSize',LocalFontSize,'FontName',get(gca,'FontName'));
ax1 = gca;
%% Ratio    - Second Subplot
s2=subplot(4,1,3);
if strcmp(Subplot2,'Ratio')
    Ratio    = RateSim4mNu./RateSim3mNu;
    RatioErr = RateErrSim4mNu./RateSim3mNu ;
    
    hr1 = plot(qU,RateSim4mNu./RateSim3mNu,'color',rgb('Salmon'),'LineWidth',3,'LineStyle','-');
    hold on;
    hr2 = plot(qU,ones(1,numel(qU)),'color',prlB,'LineWidth',3,'LineStyle',':');
    hr3 = errorbar(qU,Ratio,RatioErr,'d','MarkerSize',4,FitStyleArg{:},'CapSize',0);
    pnone = plot(qU,zeros(numel(qU),1),'LineStyle','none');
    hold off;
    ylabel('Ratio');
    katrinsim   = sprintf('3+1 simulation {\\itm}_{4} = %.1f eV   |{\\itU}_{e4}|^2 = %.2f',sterile_mass,mixing_angle);
    sterilemod  = sprintf('3+1 model');
    hl=legend([hr2 hr1 pnone hr3],{'3-\nu model',sterilemod,'',katrinsim},'Location','southwest','box','off');
    ylim([0.975 1.012]);
    hl.FontSize = LocalFontSize-2;
    
    saveStr = 'Ratio';
else
    Residuals = (RateData-RateFitData3mNu)./RateFitErrData3mNu;
    hr2 = plot(qU,zeros(1,numel(qU)),'color',rgb('DodgerBlue'),'LineWidth',3,'LineStyle',':');
    hold on;
    hr1 = plot(qU,(RateSim4mNu-RateSim3mNu)./RateErrSim4mNu,'color',rgb('Salmon'),'LineWidth',3,'LineStyle','-.');
    hr3 = plot(qU,Residuals,'.','MarkerSize',15,FitStyleArg{:});
    ylabel(sprintf('Residuals (\\sigma)'));
    PRLFormat
    katrinsim1   = sprintf('3-\\nu simulation');
    katrinsim2   = sprintf('3+1 simulation {\\itm}_{4} = %.1f eV   |{\\itU}_{e4}|^2 = %.2f',sterile_mass,mixing_angle);
    hl=legend([hr2 hr1],{katrinsim1,katrinsim2},'Location','southwest','box','off');
    hl.FontSize = LocalFontSize-2.5;
   
    ylim([-4 3])
    hold off;
    
     saveStr = 'Residuals';
end


%% apperance: legend, labels etc.
hl.NumColumns=2;
xlim([min(qU-5) max(qU+5)]);
ax2 = gca;
PRLFormat;
set(gca,'FontSize',LocalFontSize);
set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);

hl.Position(2) = 0.333;
if strcmp(Subplot2,'Residuals')
    hl.Position(1)= 0.135;
    text(t1.Position(1),2,'b)','FontSize',LocalFontSize,'FontName',get(gca,'FontName'));
    
            ax1.YLabel.Position(1) = -51;
else
    text(ax2.YLabel.Position(1)+10,1.006,'b)','FontSize',LocalFontSize,'FontName',get(gca,'FontName'));
end
%% MTD   -  Third Subplot
s3=subplot(4,1,4);
bar(qU,TimeSubrun./(60*60),0.5,'FaceColor',prlB,'EdgeColor',prlB);
% apperance: legend, labels etc.
xlabel('Retarding energy - 18574 (eV)');
ylh = ylabel('Time (h)');
ylh.Position(1) = ax2.YLabel.Position(1);%
yl1.Position(1) = ax2.YLabel.Position(1);
ylim([0 50])
yticks([0 25 50])
xlim([-range-2 50])
PRLFormat;
set(gca,'FontSize',LocalFontSize);
set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);
text(t1.Position(1),42,'c)','FontSize',get(gca,'FontSize'),'FontName',get(gca,'FontName'));
linkaxes([s1,s2,s3],'x');

ax3 = gca;
if strcmp(Subplot2,'Residuals')  
    ax1.YLabel.Position(1) = -51;
    ax2.YLabel.Position(1) = -51;
    ax3.YLabel.Position(1) = -51;
end
%% save plot
if strcmp(SavePlot,'ON')
    plotname = sprintf('%sksn1_spectrum_prl1_%s.pdf',savedirPlot,saveStr);
    export_fig(fig,plotname);
    fprintf('save plot to %s \n',plotname);
    
    print(fig,strrep(plotname,'.pdf','.png'),'-dpng','-r250');
end