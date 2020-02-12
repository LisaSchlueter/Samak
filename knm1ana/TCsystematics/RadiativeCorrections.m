%% settings
RunList               = 'KNM1';
exclDataStart         = 14;
NonPoissonScaleFactor = 1.1;
nTrials               = 1000;
RecomputeFlag         = 'OFF';
SysEffects            = struct('RF_RX','OFF','RF_EL','OFF','RF_BF','OFF','BkgShape','OFF');

%% Init Model Object and covariance matrix object
TwinRadON_T3 = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Twin',...
    'exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11',...
    'RadiativeFlag','ON');
TwinRadOFF = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Twin',...
    'exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11',...
    'RadiativeFlag','OFF');
%%
qUmin = round(TwinRadON.ModelObj.qU(TwinRadON.exclDataStart));
qUmax = round(TwinRadON.ModelObj.qU(end));
range = round(TwinRadON.ModelObj.qU(TwinRadON.exclDataStart)-TwinRadON.ModelObj.Q_i);

%% Plot The Differential Spectrum - 1
myMainTitle = sprintf('Tritium Beta Decay Differential Spectrum: Radiative Corrections');
maintitle   = myMainTitle;
savefile    = sprintf('plots/KNM1_TBDDS_RadiativeCorrections_1.png');
fig1      = figure('Name','Tritium Beta Decay Differential Spectrum: Radiative Corrections','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
a=annotation('textbox', [0 0.91 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=22;a.FontWeight='bold';
Norm1=((TwinRadOFF.ModelObj.TBDDS+TwinRadON_T1.ModelObj.TBDDS)/2);
h1 = stairs(TwinRadOFF.ModelObj.Te-TwinRadOFF.ModelObj.Q,...
    (TwinRadOFF.ModelObj.TBDDS-TwinRadON_T1.ModelObj.TBDDS)./Norm1,'LineWidth',5);
hold on
Norm2=((TwinRadOFF.ModelObj.TBDDS+TwinRadON_T2.ModelObj.TBDDS)/2);
h2 = stairs(TwinRadOFF.ModelObj.Te-TwinRadOFF.ModelObj.Q,...
    (TwinRadOFF.ModelObj.TBDDS-TwinRadON_T2.ModelObj.TBDDS)./Norm2,'LineWidth',5);
Norm3=((TwinRadOFF.ModelObj.TBDDS+TwinRadON_T3.ModelObj.TBDDS)/2);
h3 = stairs(TwinRadOFF.ModelObj.Te-TwinRadOFF.ModelObj.Q,...
    (TwinRadOFF.ModelObj.TBDDS-TwinRadON_T3.ModelObj.TBDDS)./Norm3,'LineWidth',5);
hold off
set(gca, 'YScale', 'lin');
set(gca, 'XScale', 'log');
grid on
xlabel(sprintf('E-%.1f (eV)',TwinRadOFF.ModelObj.Q_i),'FontSize',22);
ylabel('ON-OFF / ON+OFF','FontSize',22);
a = legend([h1 h2 h3],'Repko Approximation, Phys. Rev. C28 (1983) 2433','Repko Generic, Phys. Rev. C28 (1983) 2433','Sirlin, Phys. Rev. 164, 1767 (1967)'); 
PrettyFigureFormat;
xlim([-205 1e-1]);
set(gca,'FontSize',24);
print(gcf,savefile,'-dpng','-r400');

 %% Plot The Differential Spectrum - 2
myMainTitle = sprintf('Tritium Beta Decay Differential Spectrum: Radiative Corrections');
maintitle   = myMainTitle;
savefile    = sprintf('plots/KNM1_TBDDS_RadiativeCorrections_2.png');
fig1      = figure('Name','Tritium Beta Decay Differential Spectrum: Radiative Corrections','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
a=annotation('textbox', [0 0.91 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=22;a.FontWeight='bold';
Norm1=((TwinRadOFF.ModelObj.TBDDS+TwinRadON_T1.ModelObj.TBDDS)/2);
h1 = stairs(TwinRadOFF.ModelObj.Te-TwinRadOFF.ModelObj.Q,...
    (TwinRadOFF.ModelObj.TBDDS-TwinRadON_T1.ModelObj.TBDDS)./Norm1,'LineWidth',5);
hold on
Norm2=((TwinRadOFF.ModelObj.TBDDS+TwinRadON_T2.ModelObj.TBDDS)/2);
h2 = stairs(TwinRadOFF.ModelObj.Te-TwinRadOFF.ModelObj.Q,...
    (TwinRadOFF.ModelObj.TBDDS-TwinRadON_T2.ModelObj.TBDDS)./Norm2.*min(h1.YData)./min(h2.YData),'LineWidth',5);
Norm3=((TwinRadOFF.ModelObj.TBDDS+TwinRadON_T3.ModelObj.TBDDS)/2);
h3 = stairs(TwinRadOFF.ModelObj.Te-TwinRadOFF.ModelObj.Q,...
    (TwinRadOFF.ModelObj.TBDDS-TwinRadON_T3.ModelObj.TBDDS)./Norm3.*min(h1.YData)./min(h3.YData),'LineWidth',5);
hold off
set(gca, 'YScale', 'lin');
set(gca, 'XScale', 'log');
grid on
xlabel(sprintf('E-%.1f (eV)',TwinRadOFF.ModelObj.Q_i),'FontSize',22);
ylabel('ON-OFF / ON+OFF','FontSize',22);
a = legend([h1 h2 h3],'Repko Approximation, Phys. Rev. C28 (1983) 2433','Repko Generic, Phys. Rev. C28 (1983) 2433','Sirlin, Phys. Rev. 164, 1767 (1967)'); 
PrettyFigureFormat;
xlim([-205 1e-1]);
 set(gca,'FontSize',24);
print(gcf,savefile,'-dpng','-r400');
