% KSNX KATRIN final sensitivity
% Lisa, May 2020
% 750 days
% 130 mcps background
% TDR like systematics
% simple plot
%% settings for RunAnalysis
nGridSteps = 50;
DataType   = 'Fake';
chi2       = 'chi2Stat';


savedir = [getenv('SamakPath'),'ksnxana/ksnx_sensitivity/results/'];
MakeDir(savedir);
savename1 = sprintf('%sksnx_finalsensitivity_contour_%s_mNuE0NormBkg_nGridSteps%.0f.mat',...
    savedir,chi2,nGridSteps);
d1 = importdata(savename1);
savename2 = sprintf('%sksnx_finalsensitivity_contour_%s_E0NormBkg_nGridSteps%.0f.mat',...
    savedir,chi2,nGridSteps);
d2 = importdata(savename2);

%% 
GetFigure;
pNu = plot(d1.sin2T4,d1.mNu4Sq,'-.','LineWidth',2);
hold on;
pFix = plot(d2.sin2T4,d2.mNu4Sq,'LineWidth',2);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
PrettyFigureFormat('FontSize',22);
leg = legend([pFix,pNu],sprintf('{\\itm}_\\nu^2 = 0 eV^2'),sprintf('{\\itm}_\\nu^2 free'),...
    'Location','southwest');
PrettyLegendFormat(leg);
leg.Title.String = sprintf('KATRIN Final 95%% C.L.');
leg.Title.FontWeight = 'normal'; leg.Title.FontSize= 16;
xlim([1e-03 0.5]);
ylim([0.1 1600]);


