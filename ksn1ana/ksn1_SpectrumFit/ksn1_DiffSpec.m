%% plot diff spec w and w/o sterile neutrino
mNu4Sq = 20^2;
sin2T4  = 0.3;
%% settings for runanalysis
DataType = 'Real';
RunAnaArg = {'RunList','KNM1',...
    'fixPar','E0 Norm Bkg',...
    'DataType',DataType,...
    'FSDFlag','SibilleFull',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'ROIFlag','Default',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','OFF',...
    'ISCSFlag','Edep',...
    'TwinBias_Q',18573.73,...
    'SysBudget',24,...
    'pullFlag',99,...
    'NonPoissonScaleFactor',1};
T = MultiRunAnalysis(RunAnaArg{:});
T.InitModelObj_Norm_BKG;
%% diff spec w/o sterile neutrino
T.ModelObj.SetFitBiasSterile(0,0);
T.ModelObj.ComputeTBDDS;

Esur = T.ModelObj.Te-T.ModelObj.Q;  % surplus Energie
TBDDS  = T.ModelObj.TBDDS;

%% %% diff spec w sterile neutrino
T.ModelObj.SetFitBiasSterile(mNu4Sq,sin2T4);
T.ModelObj.ComputeTBDDS;
TBDDS_s  = T.ModelObj.TBDDS;

% renormalize
NormI = find(Esur>=-1,1);
TBDDS_s = TBDDS_s.*TBDDS(NormI)./TBDDS_s(NormI);
%% plot

PlotSterile = 'OFF';

LineWidth = 3;
GetFigure;
if strcmp(PlotSterile,'ON')
    ptot = plot(Esur,TBDDS_s,'-','LineWidth',LineWidth,'Color',rgb('Black'));
    hold on;
end
p = plot(Esur,TBDDS,'-.','LineWidth',LineWidth,'Color',rgb('Orange'));


if strcmp(PlotSterile,'ON')
   % hold on;
   % ptot = plot(Esur,TBDDS_s,'-.','LineWidth',LineWidth,'Color',rgb('Black'));
    ps = plot(Esur,TBDDS_s-TBDDS,':','LineWidth',LineWidth,'Color',rgb('DodgerBlue'));
    
    leg = legend([p,ps,ptot],sprintf('Active neutrino branch'),...
        sprintf('Sterile neutrino branch'),...% ({\\itm}_4 = %.0f eV , |{\\itU}_{e4}|^2 = %.1f)',sqrt(mNu4Sq),sin2T4),...
        'Total');
    SaveStr = '_sterile';
else
    leg = legend(p,'active neutrino branch');
    SaveStr = '';
end
leg.EdgeColor = rgb('Silver');
PrettyFigureFormat('FontSize',24);
xlabel(sprintf('{\\itE} - {\\itE}_0 (eV)'));
ylabel('Rate (arb. units)');
xlim([-40 1]);
ylim([-0.25 3.5])
savepath = [getenv('SamakPath'),'ksn1ana/ksn1_SpectrumFit/plots/'];
savename = sprintf('%sksn1_DiffSpec%s.png',savepath,SaveStr);
print(savename,'-dpng','-r350');
