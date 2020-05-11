%% ===== DATA =====

%% Settings
E0 = 18573.7;                                         % Endpoint in eV
%sterile_mass = sqrt(70.3);                                    % Sterile neutrino mass in eV
sterile_mass  = 10;                              % Sterile neutrino mass in eV

mixing_angle_1 = 0.01;                                % sin2(th4)
%mixing_angle_1 = 0.0224;                                % sin2(th4)
% mixing_angle_2 = 1-(1-2*mixing_angle_1)^2;          % sin2(2th4)

% mixing_angle_2 = 0.1;                               % sin2(2th4)
%mixing_angle_1 = (1-sqrt(1-mixing_angle_2))/2;       % sin2(th4)
%mixing_angle_1  = 0.5;

%% Data
D = MultiRunAnalysis('RunList','KNM1',...
            'chi2','chi2CMShape',...
            'DataType','Real',...
            'fixPar','',...
            'RadiativeFlag','ON',...
            'NonPoissonScaleFactor',1.064,...
            'minuitOpt','min ; migrad',...
            'FSDFlag','SibilleFull',...
            'ELossFlag','KatrinT2',...
            'SysBudget',22,...
            'exclDataStart',7,...
            'SynchrotronFlag','OFF',...
            'AngularTFFlag','OFF');

%% No sterile

R = MultiRunAnalysis('RunList','KNM1',...
            'chi2','chi2Stat',...
            'DataType','Twin',...
            'fixPar','Norm',...
            'RadiativeFlag','ON',...
            'NonPoissonScaleFactor',1.064,...
            'minuitOpt','min ; migrad',...
            'FSDFlag','SibilleFull',...
            'ELossFlag','KatrinT2',...
            'SysBudget',22,...
            'exclDataStart',7,...
            'SynchrotronFlag','OFF',...
            'AngularTFFlag','OFF');

% Global variables
times = R.ModelObj.qUfrac*R.ModelObj.TimeSec;
qU    = R.ModelObj.qU; qU    = qU-E0; % Energy axis

% Spectrum
R.ModelObj.ComputeTBDDS(); YD = R.ModelObj.TBDDS;
R.ModelObj.ComputeTBDIS(); YI = R.ModelObj.TBDIS; YI = YI./times;

%% Sterile theoretical

Rs = MultiRunAnalysis('RunList','KNM1',...
            'chi2','chi2Stat',...
            'DataType','Twin',...
            'fixPar','Norm',...
            'RadiativeFlag','ON',...
            'NonPoissonScaleFactor',1.064,...
            'minuitOpt','min ; migrad',...
            'FSDFlag','SibilleFull',...
            'ELossFlag','KatrinT2',...
            'SysBudget',22,...
            'exclDataStart',7,...
            'SynchrotronFlag','OFF',...
            'AngularTFFlag','OFF');
        
Rs.ModelObj.mnu4Sq_i = sterile_mass^2;
Rs.ModelObj.sin2T4_i = mixing_angle_1;

% Spectrum sterile
Rs.Fit;
Rs.ModelObj.ComputeTBDDS(); YDs = Rs.ModelObj.TBDDS;
Rs.ModelObj.ComputeTBDIS(); IS  = Rs.ModelObj.TBDIS; YIs = IS./times;

% Spectrum no-sterile
R.Fit;
R.ModelObj.ComputeTBDDS(); YD = R.ModelObj.TBDDS;
R.ModelObj.ComputeTBDIS(); YI = R.ModelObj.TBDIS; YI = YI./times;

%% Sterile "data"

YIsd = IS;
% Error - stat - 
% err  = sqrt(YIsd) ;
% Error - stat +syst
err  = (diag(sqrt(D.FitCMShape))) ;

% Fluctuations (data sim)
%YIsd = YIsd + err.*randn(length(YIsd),1);
YIsd = YIsd./times;

% Error bar
err  = err./times;
err  = err./YI;

% Constraining everything to qULimiteV
qULimit = -65;
YIsd=YIsd(qU>qULimit);
YIs=YIs(qU>qULimit);
YI=YI(qU>qULimit);
sum(YI);
err=err(qU>qULimit);
qUc=qU(qU>qULimit);

% ===== PLOTTING =====

fig = figure('Renderer','painters');
set(fig, 'Units', 'normalized', 'Position', [0.001, 0.001,0.45, 0.8]);

plot_title = sprintf('3+1 model: m_{4} = %.1f eV   |U_{e4}|^2 = %.2f',sterile_mass,mixing_angle_1);
%plot_title = sprintf('\\Deltam^{2}_{14} = %.1f eV^2 sin^2(2\\theta_{ee}) = %.2f',sterile_mass^2,mixing_angle_2);
% plot_title = sprintf('\\Deltam^{2}_{14} = %.1f eV^2',sterile_mass^2);
% Colors
prlG = [81 126 102]/255;
prlB = [50 148 216]/255;
FitStyleArg = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',4,'MarkerEdgeColor',rgb('Black')};

% Spectra  - First Subplot

s1=subplot(4,1,[1 2]);

% Normalisation
% YI_N = YIs; bkg=YI_N(length(YI_N));
% YI_N = (YI_N-bkg).*1.3 + bkg;
% 
% % Sterile Branch somehow
% YS = YI_N-YI+bkg;
YI_N = YIs; bkg = YI_N(length(YI_N));
YI_N = (YI_N-bkg).*(YI(1)/YI_N(1)) + bkg;

% Plot
%plot(qUc,YI,'DisplayName','No Sterile','color',prlB,'LineWidth',3,'LineStyle','-')
plot(D.RunData.qU-E0,D.RunData.TBDIS./(D.ModelObj.qUfrac*D.ModelObj.TimeSec),'DisplayName','No Sterile','color',prlB,'LineWidth',3,'LineStyle','-')
hold on
%plot(qUc,YI_N,'--','DisplayName','With Sterile','color',prlB,'LineWidth',3,'LineStyle','-')
%hold on
%errorbar(qUc,YIsd,err.*50,FitStyleArg{:},'CapSize',0)
errorbar(D.RunData.qU-E0,D.RunData.TBDIS./(D.ModelObj.qUfrac*D.ModelObj.TimeSec),sqrt(D.RunData.TBDIS)./(D.ModelObj.qUfrac*D.ModelObj.TimeSec)*50,FitStyleArg{:},'CapSize',0)

% hold on;
% plot(qUc,YS,'DisplayName','Sterile branch')

% Plot style
% xlabel('Retarding energy - 18574 (eV)');
ylabel('Rate (cps)');
%legend({'3-\nu model',plot_title,'KATRIN data with errobars \times50'},'Location','northeast','box','off');
legend({'3-\nu model','KATRIN data with errobars \times50'},'Location','northeast','box','off');
lgd=legend;

xlim([min(qUc-5) max(qUc+5)]);
ylim([0.18 2*max(YI_N)]);
PRLFormat;
set(gca, 'YScale', 'log');
yticks([1 10 100])

% Title
% plot_title = sprintf('\\Deltam^{2}_{14} = %.1f eV^2 and sin^2(2\\theta_{ee}) = %.2f',sterile_mass^2,mixing_angle_2);
% % ('Trtitum beta decay spectra comparison\nwith and without sterile neutrino\n\
% title(plot_title)


% Ratio    - Second Subplot

s2=subplot(4,1,3);

% Ratio
RSP  = YIs./YI;
RSPd = YIsd./YI;

% Plot
hr1 = plot(qUc,RSP,'color',rgb('Salmon'),'LineWidth',3,'LineStyle','-');
hold on;
hr2 = plot(qUc,ones(1,numel(qUc)),'color',prlB,'LineWidth',3,'LineStyle',':');
hr3 = errorbar(qUc,RSPd,err,FitStyleArg{:},'CapSize',0);
hold on;
hr4 = plot(qUc,RSPd,'o','MarkerSize',2,'MarkerFaceColor',...
                rgb('Black'),'MarkerEdgeColor',rgb('Black'),'LineWidth',3);

% Plot style
%ylabel('Ratio \nu_4/\nu_{\beta}');
ylabel('Ratio');
katrinsim   = sprintf('3+1 simulation m_{4}=%.1f eV   |U_{e4}|^2=%.2f',sterile_mass,mixing_angle_1);
sterilemod  = sprintf('3+1 model');
hl=legend([hr2 hr1 hr3],{'3-\nu model',sterilemod,katrinsim},'Location','southwest','box','off');
hl.NumColumns=1;

xlim([min(qUc-5) max(qUc+5)]);
% ylim([min([min(RSP) min(RSPd)])*0.99 max([max(RSP) max(RSPd)])*1.01])
ylim([min(RSP)*0.985 1.01])
PRLFormat;
% yticks([0.9 1])

% MTD      -  Third Subplot

s3=subplot(4,1,4)

bar(qUc,times(qU>qULimit)./(60*60),0.5,...
    'FaceColor',prlB,'EdgeColor',prlB);

xlabel('Retarding energy - 18574 (eV)');
ylh = ylabel('Time (h)');
ylh.Position(1) = -105;
ylim([0 50])
yticks([0 25 50])
PRLFormat;
linkaxes([s1,s2,s3],'x');

export_fig(fig,'./plots/ksn1_spectrum.pdf');