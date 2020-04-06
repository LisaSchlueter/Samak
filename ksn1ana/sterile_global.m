%% ===== DATA =====

%% Settings
E0 = 18570;                                         % Endpoint in eV
sterile_mass = sqrt(2.3);                                  % Sterile neutrino mass in eV

% mixing_angle_1 = 0.1;                               % sin2(th4)
% mixing_angle_2 = 1-(1-2*mixing_angle_1)^2;          % sin2(2th4)

mixing_angle_2 = 0.15;                              % sin2(2th4)
mixing_angle_1 = (1-sqrt(1-mixing_angle_2))/2;        % sin2(th4)

%% No sterile

R = MultiRunAnalysis('RunList','KNM1',...
            'chi2','chi2Stat',...
            'DataType','Twin',...
            'fixPar','',...
            'RadiativeFlag','ON',...
            'NonPoissonScaleFactor',1.064,...
            'minuitOpt','min ; migrad',...
            'FSDFlag','Sibille0p5eV',...
            'ELossFlag','KatrinT2',...
            'SysBudget',22);

% Global variables
times = (R.ModelObj.qUfrac*R.ModelObj.TimeSec);

qU = R.ModelObj.qU;     % Energy axis
qU = qU-E0;

% Spectrum
R.ModelObj.ComputeTBDDS();
YD=R.ModelObj.TBDDS;
R.ModelObj.ComputeTBDIS();

YI = R.ModelObj.TBDIS;
sum(YI)
YI = YI./times;

%% Sterile theoretical

Rs = MultiRunAnalysis('RunList','KNM1',...
            'chi2','chi2Stat',...
            'DataType','Twin',...
            'fixPar','',...
            'RadiativeFlag','ON',...
            'NonPoissonScaleFactor',1.064,...
            'minuitOpt','min ; migrad',...
            'FSDFlag','Sibille0p5eV',...
            'ELossFlag','KatrinT2',...
            'SysBudget',22);

Rs.ModelObj.mnu4Sq_i = sterile_mass^2;
Rs.ModelObj.sin2T4_i = mixing_angle_1;

% Spectrum
Rs.ModelObj.ComputeTBDDS();
YDs=Rs.ModelObj.TBDDS;
Rs.ModelObj.ComputeTBDIS();

IS = Rs.ModelObj.TBDIS;     % Counts
YIs = IS./times;

%% Sterile "data"

YIsd = IS;

% Error
err  = sqrt(YIsd);

% Fluctuations (data sim)
YIsd = YIsd + err.*randn(length(YIsd),1);

YIsd = YIsd./times;

% Error bar
err  = err./times;
err  = err./YI;


%% Constraining everything to -90eV
YIsd=YIsd(qU>-90);
YIs=YIs(qU>-90);
YI=YI(qU>-90);
sum(YI)
err=err(qU>-90);
qUc=qU(qU>-90);



%% ===== PLOTTING =====

fig = figure('Renderer','painters');
set(fig, 'Units', 'normalized', 'Position', [0.001, 0.001,0.45, 0.8]);

plot_title = sprintf('\\Deltam^{2}_{14} = %.1f eV^2 sin^2(2\\theta_{ee}) = %.2f',sterile_mass^2,mixing_angle_2);
% plot_title = sprintf('\\Deltam^{2}_{14} = %.1f eV^2',sterile_mass^2);
%% Colors
prlG = [81 126 102]/255;
prlB = [50 148 216]/255;
FitStyleArg = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',4,'MarkerEdgeColor',rgb('Black')};

%% Spectra  - First Subplot

subplot(4,1,[1 2]);

% Normalisation
% YI_N = YIs; bkg=YI_N(length(YI_N));
% YI_N = (YI_N-bkg).*1.3 + bkg;
% 
% % Sterile Branch somehow
% YS = YI_N-YI+bkg;
YI_N = YIs; bkg=YI_N(length(YI_N));
YI_N = (YI_N-bkg).*(YI(1)/YI_N(1)) + bkg;

% Plot
plot(qUc,YI,'DisplayName','No Sterile','color',prlB,'LineWidth',3)
hold on
plot(qUc,YI_N,'--','DisplayName','With Sterile','color',prlG,'LineWidth',3)
hold on
errorbar(qUc,YIsd,err.*50,FitStyleArg{:},'CapSize',0)

% hold on;
% plot(qUc,YS,'DisplayName','Sterile branch')

% Plot style
% xlabel('Retarding energy - 18574 (eV)');
ylabel('Rate (cps)');
legend({'KATRIN model with \Deltam^{2}_{14} = 0 eV^2',plot_title,'KATRIN MC data with errobars \times50'},'Location','northeast','box','off');
lgd=legend;

xlim([min(qUc-5) max(qUc+5)]);
ylim([0.18 2*max(YI_N)]);
PRLFormat;
set(gca, 'YScale', 'log');
yticks([1 10 100])

%% Title
% plot_title = sprintf('\\Deltam^{2}_{14} = %.1f eV^2 and sin^2(2\\theta_{ee}) = %.2f',sterile_mass^2,mixing_angle_2);
% % ('Trtitum beta decay spectra comparison\nwith and without sterile neutrino\n\
% title(plot_title)


%% Ratio    - Second Subplot

subplot(4,1,3);

% Ratio
RSP  = YIs./YI;
RSPd = YIsd./YI;

% Plot
plot(qUc,RSP,'color',prlB,'LineWidth',3)
hold on;
errorbar(qUc,RSPd,err,FitStyleArg{:},'CapSize',0)
hold on;
plot(qUc,RSPd,'o','MarkerSize',2,'MarkerFaceColor',...
                rgb('Black'),'MarkerEdgeColor',rgb('Black'),'LineWidth',3)

% Plot style
ylabel('Ratio \nu_4/\nu_{\beta}');
legend({'Model','KATRIN MC Data'},'Location','southeast','box','off');

xlim([min(qUc-5) max(qUc+5)]);
% ylim([min([min(RSP) min(RSPd)])*0.99 max([max(RSP) max(RSPd)])*1.01])
ylim([min(RSP)*0.99 1.01])
PRLFormat;
% yticks([0.9 1])

%% MTD      -  Third Subplot

subplot(4,1,4)

bar(qUc,times(qU>-90)./(60*60),0.5,...
    'FaceColor',prlB,'EdgeColor',prlB);

xlabel('Retarding energy - 18574 (eV)');
ylh = ylabel('Time (h)');
ylh.Position(1) = -105;
ylim([0 50])
yticks([0 25 50])
PRLFormat;