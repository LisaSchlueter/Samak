RunList = 'KNM1';
[Twin, ~, ~, ~, ~, ~,~] = ComputeLoadTwinObjects('RunList',RunList);
savedir = [getenv('SamakPath'),'knm1ana/knm1Twins/results/'];

FitRange40 = round(-Twin.ModelObj.qU(14)+18573.7); %for labeling
save_file = [savedir,'FitResult_Twins_all',RunList,sprintf('%.0feVrange',FitRange40),'.mat'];
d40 = importdata(save_file);
mNuSq40     = d40.mNuSq;
mNuSqErr40  = d40.mNuSqErr;
E040        = d40.E0;
E0Err40     = d40.E0Err;
chi2min40   = d40.chi2min;
dof40 = d40.dof;

FitRange30 = round(-Twin.ModelObj.qU(17)+18573.7); %for labeling
save_file = [savedir,'FitResult_Twins_all',RunList,sprintf('%.0feVrange',FitRange30),'.mat'];
d30 = importdata(save_file);
mNuSq30     = d30.mNuSq;
mNuSqErr30  = d30.mNuSqErr;
E030        = d30.E0;
E0Err30     = d30.E0Err;
chi2min30   = d30.chi2min;
dof30 = d30.dof;
%% plot
Twin.GetPlotColor;
twinLabel = {'Twin',sprintf('\\langle\\rhod\\rangle'),...
    sprintf('\\langlec_{T_2}\\rangle'),...
    sprintf('\\langlet_{{\\itqU}}\\rangle'),...
    sprintf('\\langle{\\itqU}\\rangle'),...
    sprintf('\\langlet_{{\\itqU}}, {\\itqU}\\rangle'),...
    sprintf('\\langleTwin\\rangle')};%,'all same -qUf'
lineArg = {'.','LineWidth',3,'MarkerSize',27,'Color',rgb('Red')};
x = 1:numel(twinLabel);

% get mNuSq40 from presentation in analysis workshop
mNuSq40 = [-4.5,-3.52,-3.56,-4.05,-1.02,-0.35,0];
% neutrino mass shift
%figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
GetFigure;
%  plot(x,mNuSq30,lineArg{:});
%  hold on;
plot(x,mNuSq40,lineArg{:});
hold off
xticks(1:numel(x));
xticklabels(twinLabel);
xlim([0.8 numel(x)+0.2])
%xlabel('MC Twin slow-control configuration');
PrettyFigureFormat('FontSize',24);
ylabel(sprintf('\\Delta{\\itm}^2_\\nu \\times 10^{-3} (eV^{ 2})'));
xlabel('Simulated spectrum')
set(gca,'XMinorTick','off');

 xticks(1:numel(x));
 xticklabels(twinLabel)
 %xtickangle(35)
%leg = legend('30 eV','40 eV'); legend boxoff
%leg.Title.String = 'range'; leg.Location = 'best';

plotdir =[getenv('SamakPath'),'knm1ana/knm1_Stacking/plots/'];
plotname = [plotdir,sprintf('knm1_Stacking_FitTwins.pdf')];
export_fig(plotname);


% % % endpoint
% fig2 = figure(2);
% set(fig2, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
% plot(x,E030,lineArg{:});
% hold on;
% plot(x,E040,lineArg{:},'LineStyle','--','Color',rgb('DarkGoldenRod'),'MarkerFaceColor',rgb('GoldenRod'));
% hold on;
% xticks(1:numel(x));
% xticklabels(twinLabel)
% xtickangle(35)
% ylabel(sprintf('E0 shift (eV)'));
% leg = legend('30 eV','40 eV'); legend boxoff
% leg.Title.String = 'range'; leg.Location = 'best';
% PrettyFigureFormat;
% set(gca,'FontSize',22);
% grid on;

% % chi2 /dof
% fig3 = figure(3);
% set(fig3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
% plot(x,chi2min30./dof30,lineArg{:});
% hold on;
% plot(x,chi2min40./dof40,lineArg{:},'LineStyle','--','Color',rgb('DarkGoldenRod'),'MarkerFaceColor',rgb('GoldenRod'));
% hold off;
% xticks(1:numel(x));
% xticklabels(twinLabel)
% xtickangle(35)
% ylabel(sprintf('\\chi2 /  dof'));
% leg = legend('30 eV','40 eV'); legend boxoff
% leg.Title.String = 'range'; leg.Location = 'best';
% PrettyFigureFormat;
% xlim([0.5,numel(x)+0.5])
% grid on;
% set(gca,'FontSize',22);


%%
% save plot

savedirplots = strrep(savedir,'results','plots');
savedirplots = strrep(savedirplots,'knm1Twins','knm1_Stacking');
if ~exist(savedirplots,'dir')
    system(['mkdir ',savedirplots]);
end
figname = [savedirplots,sprintf('TwinStackingEffect_%s',RunList)];

print(fig1,[figname,'mNu3040.png'],'-dpng','-r450');
% print(fig2,[figname,'E03040.png'],'-dpng','-r450');
% print(fig3,[figname,'chi23040.png'],'-dpng','-r450');
