%% ===== DATA =====

%% Settings
E0 = 18570;                                         % Endpoint in eV
sterile_mass = 20;                                  % Sterile neutrino mass in eV

mixing_angle_1 = 0.1;                               % sin2(th4)
mixing_angle_2 = 1-(1-2*mixing_angle_1)^2;          % sin2(2th4)

% mixing_angle_2 = 0.15;                              % sin2(2th4)
% mixing_angle_1 = (1-sqrt(1-mixing_angle))/2;        % sin2(th4)

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
err=err(qU>-90);
qUc=qU(qU>-90);



%% ===== PLOTTING =====

figure;

%% Colors
prlG = [81 126 102]/255;
prlB = [50 148 216]/255;
FitStyleArg = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',4,'MarkerEdgeColor',rgb('Black')};

%% Spectra  - First Subplot

subplot(2,1,1);

% Normalisation
YI_N = YIs; bkg=YI_N(length(YI_N));
YI_N = (YI_N-bkg).*1.3 + bkg;
% 
% % Sterile Branch somehow
% YS = YI_N-YI+bkg;

% Plot
plot(qUc,YI,'DisplayName','No Sterile','color',prlB,'LineWidth',3)
hold on;
plot(qUc,YI_N,'--','DisplayName','With Sterile','color',prlG,'LineWidth',3)
% hold on;
% plot(qUc,YS,'DisplayName','Sterile branch')

% Plot style
xlabel('Retarding energy - 18574 (eV)');
ylabel('Rate (cps)');
legend({'Without sterile','With sterile'},'Location','southwest','box','off');
lgd=legend;

ylim([0.18 2*max(YI_N)]);
PRLFormat;
set(gca, 'YScale', 'log');

%% Title
plot_title = sprintf('Trtitum beta decay spectra comparison with and without sterile neutrino\nwith \\Deltam_{14} = %.1f eV, sin^2(\\theta_{ee}) = %.2f and sin^2(2\\theta_{ee}) = %.3f',sterile_mass,mixing_angle_1,mixing_angle_2);
title(plot_title)


%% Ratio    - Second Subplot

subplot(2,1,2);

% Ratio
RSP  = YIs./YI;
RSPd = YIsd./YI;

% Plot
plot(qUc,RSP,'color',prlB,'LineWidth',3)
hold on;
errorbar(qUc,RSPd,err,FitStyleArg{:},'CapSize',0)
hold on;
plot(qUc,RSPd,'o','MarkerSize',4,'MarkerFaceColor',...
                rgb('Black'),'MarkerEdgeColor',rgb('Black'),'LineWidth',3)

% Plot style
xlabel('Retarding energy - 18574 (eV)');
ylabel('Spectra ratio sterile/no sterile');
legend({'Model','Fake Data'},'Location','southwest','box','off');

PRLFormat;
set(gca, 'YScale', 'log');
