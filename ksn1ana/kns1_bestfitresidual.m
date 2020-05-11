%
% Build KSN1 Plot of:
% - Spectrum (Data)
% - Residuals Best Fit
% - MTD
%
% T. Lasserre
% May 2020
%

rangeeV    = '95eV';
Type       = 'Residual'; % 'Residual' or 'Ratio'

% Inputs
E0         = 18573.7;       % Endpoint in eV

switch rangeeV
    case '65eV'
        m4Sq       = 78.1921;       % Sterile mass in eV^2
        Ue4Sq      = 0.0251;        % sin2(th4)
        rangeIndex = 7;             % 65eV;
    case '95eV'
        m4Sq       = 3.913419e+03;        % Sterile mass in eV^2
        Ue4Sq      = 0.01433;         % sin2(th4)
        rangeIndex = 1;             % 95eV;
end

%% Data
kns1Data   = MultiRunAnalysis('RunList','KNM1',...
    'chi2','chi2CMShape',...
    'DataType','Real',...
    'fixPar','E0 Norm Bkg',...
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1.064,...
    'minuitOpt','min ; minos',...
    'FSDFlag','SibilleFull',...
    'ELossFlag','KatrinT2',...
    'SysBudget',22,...
    'exclDataStart',rangeIndex,...
    'SynchrotronFlag','OFF',...
    'AngularTFFlag','OFF');
% Data
KNS1hv        = kns1Data.RunData.qU-E0;
KNS1dataRate  = kns1Data.RunData.TBDIS./(kns1Data.ModelObj.qUfrac*kns1Data.ModelObj.TimeSec);
KNS1dataRateE = sqrt(kns1Data.RunData.TBDIS)./(kns1Data.ModelObj.qUfrac*kns1Data.ModelObj.TimeSec)*50;
KNS1hv        = KNS1hv(rangeIndex:end);
KNS1dataRate  = KNS1dataRate(rangeIndex:end);
KNS1dataRateE = KNS1dataRateE(rangeIndex:end);

% NullHypothesis Fit
kns1Data.Fit;
kns1FitNHhv        = kns1Data.RunData.qU-E0;
kns1FitNHRate      = kns1Data.ModelObj.TBDIS./(kns1Data.RunData.qUfrac*kns1Data.RunData.TimeSec);
kns1FitNHRateE     = diag(sqrt(kns1Data.FitCMShape))./(kns1Data.RunData.qUfrac*kns1Data.RunData.TimeSec);
kns1FitNHhv        = kns1FitNHhv(rangeIndex:end);
kns1FitNHRate      = kns1FitNHRate(rangeIndex:end);
kns1FitNHRateE     = kns1FitNHRateE(rangeIndex:end);

% Sterile Fit
kns1Fit   = MultiRunAnalysis('RunList','KNM1',...
    'chi2','chi2CMShape',...
    'DataType','Real',...
    'fixPar','E0 Norm Bkg',...
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1.064,...
    'minuitOpt','min ; minos',...
    'FSDFlag','SibilleFull',...
    'ELossFlag','KatrinT2',...
    'SysBudget',22,...
    'exclDataStart',rangeIndex,...
    'SynchrotronFlag','OFF',...
    'AngularTFFlag','OFF');

kns1Fit.ModelObj.mnu4Sq_i   = m4Sq;
kns1Fit.ModelObj.sin2T4_i   = Ue4Sq;
kns1Fit.Fit;

%% Fit Model
kns1Fithv        = kns1Fit.RunData.qU-E0;
kns1FitRate      = kns1Fit.ModelObj.TBDIS./(kns1Fit.RunData.qUfrac*kns1Fit.RunData.TimeSec);
kns1FitRateE     = diag(sqrt(kns1Fit.FitCMShape))./(kns1Fit.RunData.qUfrac*kns1Fit.RunData.TimeSec);
kns1Fithv        = kns1Fithv(rangeIndex:end);
kns1FitRate      = kns1FitRate(rangeIndex:end);
kns1FitRateE     = kns1FitRateE(rangeIndex:end);

%% Plot
fig = figure('Renderer','painters');
set(fig,'Units','normalized','Position',[0.001, 0.001,0.45, 0.8]);
FitStyleArg = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',4,'MarkerEdgeColor',rgb('Black')};
prlG        = [81 126 102]/255;
prlB        = [50 148 216]/255;


% First Subplot - Spectral Data
s1=subplot(4,1,[1 2]);
plot_title  = sprintf('KATRIN 3+1 model: m_{4} = %.1f eV   |U_{e4}|^2 = %.3f',sqrt(m4Sq),Ue4Sq);
hm = plot(KNS1hv,KNS1dataRate,...
    'color',prlB,'LineWidth',3,'LineStyle','-');
hold on
hd = errorbar(KNS1hv,KNS1dataRate,KNS1dataRateE,...
    FitStyleArg{:},'CapSize',0);

title(plot_title);
PRLFormat;

ylabel('Rate (cps)'); set(gca, 'YScale', 'log');
lgd=legend([hm hd],...
    {'NH model','KATRIN data with errobars \times50'},...
    'Location','northeast','box','off');

yticks([1 10 100])
ylim([0.18 2*max(KNS1dataRate)]);


% - Second Subplot - Residuals Best Fit Sterile
s2=subplot(4,1,3);
switch Type
    case 'Residual'
        hbfr = errorbar(KNS1hv,...
            (KNS1dataRate-kns1FitRate)./kns1FitRateE,...
            (KNS1dataRate-kns1FitRate)./kns1FitRateE.*0,...
            FitStyleArg{:},'CapSize',0);
        ylim([-2 2]);ylabel('Residuals');
        hold on
        plot(KNS1hv,KNS1hv.*0,'--','Color',prlG);
        hold off
    case 'Ratio'
        hbfr = errorbar(KNS1hv,...
            (KNS1dataRate./kns1FitRate),...
            (KNS1dataRate./kns1FitRate).*0,...
            FitStyleArg{:},'CapSize',0);
        ylabel('Ratio');
        hold on
        plot(KNS1hv,KNS1hv./KNS1hv,'--','Color',prlG);
        hold off
end


lgd=legend([hbfr],...
    {sprintf('3+1 \\chi^2=%.1f/%.0f',kns1Fit.FitResult.chi2min,kns1Fit.FitResult.dof)},...
    'Location','northeast','box','off');
PRLFormat;


% - Third Subplot - Residuals Best Fit NH
s3=subplot(4,1,4);
switch Type
    case 'Residual'
        hbfr = errorbar(kns1FitNHhv,...
            (KNS1dataRate-kns1FitNHRate)./kns1FitNHRateE,...
            (KNS1dataRate-kns1FitNHRate)./kns1FitNHRateE.*0,...
            FitStyleArg{:},'CapSize',0);
        ylim([-2 2]);ylabel('Residuals');
        hold on
        plot(KNS1hv,KNS1hv.*0,'--','Color',prlG);
        hold off
    case 'Ratio'
        hbfr = errorbar(KNS1hv,...
            (KNS1dataRate./kns1FitNHRate),...
            (KNS1dataRate./kns1FitNHRate).*0,...
            FitStyleArg{:},'CapSize',0);
        ylabel('Ratio');
        hold on
        plot(KNS1hv,KNS1hv./KNS1hv,'--','Color',prlB);
        hold off
end


lgd=legend([hbfr],...
    {sprintf('NH \\chi^2=%.1f/%.0f',kns1Data.FitResult.chi2min,kns1Data.FitResult.dof)},...
    'Location','northeast','box','off');
PRLFormat;

linkaxes([s1,s2,s3],'x');
xlim([min(KNS1hv-5) max(KNS1hv+5)]);

% Save File
figtitle = sprintf('./plots/ksn1_spectrumresiduals_%s_m%0.feV2_%s.pdf',rangeeV,floor(m4Sq),Type);
export_fig(fig,figtitle);
