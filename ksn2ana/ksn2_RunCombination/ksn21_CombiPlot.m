% ksn2+1 compare results ksn1 and ksn2 stand-alone

% load results fro MiniFiles
chi2 = 'chi2CMShape';
DataType = 'Real';
range = 40;
BF = 'OFF';
%% load ksn-2 only results
savedir = [getenv('SamakPath'),'SterileAnalysis/MiniFiles/KSN2/'];
savefilemNuFix = sprintf('%sKSN2contour_%s_%s_%s_%.0feV.mat',savedir,DataType,'E0NormBkg',chi2,range);
d = importdata(savefilemNuFix);
fprintf('load %s \n',savefilemNuFix);

savefilemNuFree = sprintf('%sKSN2contour_%s_%s_%s_%.0feV.mat',savedir,DataType,'mNuE0NormBkg',chi2,range);
dmNu = importdata(savefilemNuFree);
fprintf('load %s \n',savefilemNuFree);

%% load combined results
savedirCombi = sprintf('%sksn2ana/ksn2_RunCombination/results/',getenv('SamakPath'));

savefileCombimNuFix = sprintf('%sksn21_Combination_ReAna_%s.mat',savedirCombi,DataType);                            
dCombi = importdata(savefileCombimNuFix);
fprintf('load %s \n',savefileCombimNuFix);

savefileCombimNuFree = sprintf('%sksn21_Combi_freemNuSq_ReAna_%s.mat',savedirCombi,DataType);
dCombimNu = importdata(savefileCombimNuFree);
fprintf('load %s \n',savefileCombimNuFree);

legHandle = cell(0,0);
legStr = '';
GetFigure;

% ksn1
p1mNuFix  = plot(dCombi.sin2T4_contour_1,dCombi.mNu4Sq_contour_1,'-.','LineWidth',3,'Color',rgb('SkyBlue'));
hold on;
p1mNuFree = plot(dCombimNu.sin2T4_k1_contour,dCombimNu.mNu4Sq_1_contour,'-.','LineWidth',3,'Color',rgb('Salmon'));

% ksn2
p2mNuFix = plot(d.sin2T4_contour,d.mNu4Sq_contour,'-','LineWidth',3,'Color',rgb('DodgerBlue')); hold on;
hold on;
p2mNuFree = plot(dmNu.sin2T4_contour,dmNu.mNu4Sq_contour,'-','LineWidth',3,'Color',rgb('Crimson'));

% combi
pCombimNuFix = plot(dCombi.sin2T4_contour_12,dCombi.mNu4Sq_contour_12,':','LineWidth',3,'Color',rgb('Navy'));
pCombimNuFree = plot(dCombimNu.sin2T4_k12_contour,dCombimNu.mNu4Sq_k12_contour,':','LineWidth',3,'Color',rgb('DarkRed'));


xlim([3e-03 0.5]);
ylim([0.3 2200]);
yticks([0.1,1,10,100,1e3]);

if strcmp(BF,'ON')
   pmNuFix_bf2 = plot(d.sin2T4_bf,d.mNu4Sq_bf,'x','LineWidth',3,'Color',p2mNuFix.Color,'MarkerSize',11);
   pmNuFree_bf2 = plot(dmNu.sin2T4_bf,dmNu.mNu4Sq_bf,'x','LineWidth',3,'Color',p2mNuFree.Color,'MarkerSize',11);
   pCombimNuFix_bf = plot(dCombi.sin2T4_bf_12,dCombi.mNu4Sq_bf_12,'x','LineWidth',2.5,'Color',pCombimNuFix.Color,'MarkerSize',11);
   pCombimNuFree_bf = plot(dCombimNu.sin2T4bf_k12,dCombimNu.mNu4Sqbf_k12,'x','LineWidth',2.5,'Color',pCombimNuFree.Color,'MarkerSize',11);
   pmNuFix_bf1 = plot(dCombi.sin2T4_bf_1,dCombi.mNu4Sq_bf_1,'x','LineWidth',2.5,'Color',p1mNuFix.Color,'MarkerSize',11);
  pmNuFree_bf1 = plot(dCombimNu.sin2T4bf_k1,dCombimNu.mNu4Sqbf_k1,'x','LineWidth',2.5,'Color',p1mNuFree.Color,'MarkerSize',11);

end

pNone = plot(NaN,NaN,'Color','none');
pNone2 = plot(NaN,NaN,'Color','none');
set(gca,'YScale','log');
set(gca,'XScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^{2})'));

% legend
leg = legend([p1mNuFix,p2mNuFix,pCombimNuFix,p1mNuFree,p2mNuFree,pCombimNuFree],...
  sprintf('KSN1'),...%:     {\\itm}_\\nu^2 = 0 eV^2'),...
   sprintf('KSN2'),...%:     {\\itm}_\\nu^2 = 0 eV^2'),...
      sprintf('KSN1+2'),...%: {\\itm}_\\nu^2 = 0 eV^2'),...
    sprintf('KSN1'),...%: {\\itm}_\\nu^2 free'),...
    sprintf('KSN2'),...%: {\\itm}_\\nu^2 free'),...
    sprintf('KSN1+2'),...%: {\\itm}_\\nu^2 free'),...
    'Location','southwest','box','off');

PrettyLegendFormat(leg);
leg.NumColumns = 2;
%PrettyFigureFormat('FontSize',22);
PRLFormat;
set(gca,'FontSize',30);

%ax = gca;
%ax.FontSize =  26;
leg.FontSize = 26;
% ax.XLabel.FontSize = 30;
% ax.YLabel.FontSize = 30;

axleg=axes('Position',get(gca,'Position'),'Visible','Off');
hl2 = legend(axleg,[pNone,pNone2],...
    sprintf('{\\itm}_\\nu^2 = 0 eV^2'),...
    sprintf('{\\itm}_\\nu^2 free'),...
    'Location','southwest','box','off');
hl2.FontSize = 26;
hl2.NumColumns = 2;
PRLFormat;
hl2.Position(2) = 0.46;
hl2.Position(1)  = 0.125;
%%
pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_RunCombination/plots/'];
MakeDir(pltdir);
pltname = sprintf('%sksn21_CombiPlot_%s.png',pltdir,DataType);
print(gcf,pltname,'-dpng','-r350');
pltname = sprintf('%sksn21_CombiPlot_%s.pdf',pltdir,DataType);
export_fig(gcf,pltname);


