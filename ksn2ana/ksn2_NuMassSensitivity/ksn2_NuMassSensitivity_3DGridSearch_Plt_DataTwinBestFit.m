% plot best fits
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];

DataType = 'Real';

%% data
savefileCombi = sprintf('%sksn2_%s_NuMassSensitivityGridSearch_CombiSTEREO-%s.mat',...
    savedir,DataType,'OFF');

if exist(savefileCombi,'file') 
    d = importdata(savefileCombi);
else
    fprintf('file not available \n')
end

% load contour with free mNu
ContourFile = sprintf('%sSterileAnalysis/MiniFiles/KSN2/KSN2contour_%s_mNuE0NormBkg_chi2CMShape_40eV.mat',getenv('SamakPath'),DataType);
dA = importdata(ContourFile);
%%
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.45]);
sp = scatter3(d.sin2T4_bf,d.mNu4Sq_bf,d.FixmNuSq_all,50,d.FixmNuSq_all,'filled');%,30,cmp,'filled');%,'MarkerSize',10)
view(2)
grid off;
c = colorbar; 
c.Limits =([min(d.FixmNuSq_all) max(d.FixmNuSq_all)]);
hold on;
pD = plot(dA.sin2T4_contour,dA.mNu4Sq_contour,'-k','LineWidth',2);

set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_{4}^2 (eV^2)'));
c.Label.String = sprintf('{\\itm}_\\nu^2 (eV^2)');
PrettyFigureFormat('FontSize',19);
c.Label.FontSize = get(gca,'FontSize')+4;

%%
if strcmp(DataType,'Real')
    pDbf = plot3(dA.sin2T4_bf,dA.mNu4Sq_bf,99,'pk','LineWidth',2,'MarkerSize',14,'MarkerIndices',[1,1],'MarkerFaceColor','k');
    leg = legend([pD,pDbf,sp],...
    sprintf('KNM2 exclusion contour at 95%% C.L. with free {\\itm}_\\nu^2'),...
    sprintf('Best fit for {\\itm}_\\nu^2 free'),...
    sprintf('Best fits for different fixed {\\itm}_\\nu^2'),...
     'Location','southwest');
 %leg.Title.String = sprintf('Data');
else
    leg = legend([sp,pD,pDbf],...
    sprintf('Best fits for different fixed {\\itm}_\\nu^2'),...
    sprintf('Sensitivity contour at 95%% C.L. with free {\\itm}_\\nu^2'),...
     sprintf('Best fit for {\\itm}_\\nu^2 free'),...
     'Location','southwest');
 leg.Title.String = sprintf('MC twins');
end


leg.Title.FontSize = get(gca,'FontSize'); leg.Title.FontWeight = 'normal';
leg.FontSize = get(gca,'FontSize')+2;
PrettyLegendFormat(leg);
xlim([4e-03 0.5]);
ylim([0.3 2e3]);
%%
pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);
pltname = strrep(strrep(savefileCombi,'results','plots'),'_CombiSTEREO-OFF.mat','_BestFitPos.pdf');
%print(gcf,pltname,'-dpng','-r350');
export_fig(pltname);
fprintf('save plot to %s \n',pltname);