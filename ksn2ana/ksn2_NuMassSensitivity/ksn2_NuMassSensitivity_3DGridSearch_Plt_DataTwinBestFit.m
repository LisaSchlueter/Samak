% plot best fits
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];

DataType = 'Real';

%% data
savefileCombi = sprintf('%sksn2_%s_NuMassSensitivityGridSearch_CombiSTEREO-%s.mat',...
    savedir,DataType,'ON');

if exist(savefileCombiR,'file') 
    d = importdata(savefileCombi);
else
    fprintf('file not available \n')
end

% load contour with free mNu
ContourFile = sprintf('%sSterileAnalysis/MiniFiles/KSN2/KSN2contour_%s_mNuE0NormBkg_chi2CMShape_40eV.mat',getenv('SamakPath'),DataType);
dA = importdata(ContourFile);
%%
GetFigure;
sp = scatter3(d.sin2T4_bf,d.mNu4Sq_bf,d.FixmNuSq_all,40,d.FixmNuSq_all,'filled');%,30,cmp,'filled');%,'MarkerSize',10)
view(2)
grid off;
c = colorbar; 
c.Limits =([min(d.FixmNuSq_all) max(d.FixmNuSq_all)]);
hold on;
pD = plot(dA.sin2T4_contour,dA.mNu4Sq_contour,'-k','LineWidth',2);

PrettyFigureFormat('FontSize',22);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_{4}^2 (eV^2)'));
c.Label.String = sprintf('{\\itm}_\\nu^2 (eV^2)');
c.Label.FontSize = get(gca,'FontSize');


if strcmp(DataType,'Real')
    pDbf = plot3(dA.sin2T4_bf,dA.mNu4Sq_bf,99,'xk','LineWidth',2,'MarkerSize',8);
    leg = legend([sp,pD,pDbf],...
    sprintf('Best fits for different fixed {\\itm}_\\nu^2'),...
    sprintf('Exclusion contour at 95%% C.L. with free {\\itm}_\\nu^2'),...
     sprintf('Best fit for {\\itm}_\\nu^2 free'),...
     'Location','southwest');
 leg.Title.String = sprintf('Data');
else
    leg = legend([sp,pD,pDbf],...
    sprintf('Best fits for different fixed {\\itm}_\\nu^2'),...
    sprintf('Sensitivity contour at 95%% C.L. with free {\\itm}_\\nu^2'),...
     sprintf('Best fit for {\\itm}_\\nu^2 free'),...
     'Location','southwest');
 leg.Title.String = sprintf('MC twins');
end


leg.Title.FontSize = get(gca,'FontSize'); leg.Title.FontWeight = 'normal';
PrettyLegendFormat(leg);
xlim([4e-03 0.5]);
ylim([0.1 40^2]);

pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);
pltname = strrep(strrep(savefileCombi,'results','plots'),'_CombiSTEREO-OFF.mat','_BestFitPos.png');
print(gcf,pltname,'-dpng','-r350');
fprintf('save plot to %s \n',pltname);