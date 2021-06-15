% ksn2 compare results in Mainz/Troitsk
% also with ksn1+2 combination
% load results fro MiniFiles
chi2 = 'chi2CMShape';
DataType = 'Real';
range = 40;
Mainz = 'ON';
Troitsk = 'ON';
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
%% load others
  savedirOther = [getenv('SamakPath'),'SterileAnalysis/GridSearchFiles/Knm1/Others/'];   
if strcmp(Mainz,'ON')
    filenameMainz = sprintf('%scoord_Mainz_95CL.mat',savedirOther);
    dMainz = importdata(filenameMainz);
    sinTsq = 0.5*(1-sqrt(1-dMainz.SinSquare2Theta_X));
    pMainz = plot(sinTsq,dMainz.DmSquare41_Y,'-.','LineWidth',1.5,'Color',rgb('LimeGreen'));
    hold on;
end
%% Troitsk
if strcmp(Troitsk,'ON')
    filenameTroitsk = sprintf('%scoord_Troitsk_95CL.mat',savedirOther);
    dTroitsk = importdata(filenameTroitsk);
    pTroitsk = plot(dTroitsk.SinSquareTheta_X,dTroitsk.m4Square_Y,'--','LineWidth',1.5,...
        'Color',rgb('DarkGreen')); % Orange
    hold on;
end

set(gca,'YScale','log');
set(gca,'XScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^{ 2})'));
PrettyFigureFormat('FontSize',22);

%% plot KATRIN
pmNuFix = plot(d.sin2T4_contour,d.mNu4Sq_contour,'-','LineWidth',3,'Color',rgb('DodgerBlue')); hold on;
pmNuFree = plot(dmNu.sin2T4_contour,dmNu.mNu4Sq_contour,'-','LineWidth',3,'Color',rgb('Salmon'));
pCombimNuFix = plot(dCombi.sin2T4_contour_12,dCombi.mNu4Sq_contour_12,':','LineWidth',2.5,'Color',rgb('Navy'));
pCombimNuFree = plot(dCombimNu.sin2T4_k12_contour,dCombimNu.mNu4Sq_k12_contour,'-.','LineWidth',2.5,'Color',rgb('FireBrick'));
xlim([3e-03 0.5]);
ylim([0.2 2200]);

%% legend
leg = legend([pMainz,pTroitsk,pmNuFix,pmNuFree,pCombimNuFix,pCombimNuFree],...
    sprintf('Mainz: {\\itm}_\\nu^2 = 0 eV^2'),...
    sprintf('Troitsk: {\\itm}_\\nu^2 = 0 eV^2'),...
    sprintf('KSN2: {\\itm}_\\nu^2 = 0 eV^2'),...
    sprintf('KSN2: {\\itm}_\\nu^2 free'),...
     sprintf('KSN1+2: {\\itm}_\\nu^2 = 0 eV^2'),...
    sprintf('KSN1+2: {\\itm}_\\nu^2 free'),...
    'Location','south');
PrettyLegendFormat(leg);
leg.NumColumns = 3;
leg.FontSize = 15;

pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_OtherExp/plots/'];
MakeDir(pltdir);
% pltname = sprintf('%sksn2_CompareMainzTroitskCombi_%s.png',pltdir,DataType);
% print(gcf,pltname,'-dpng','-r350');
pltname = sprintf('%sksn2_CompareMainzTroitskCombi_%s.pdf',pltdir,DataType);
%export_fig(gcf,pltname);


%% find differences in sin^2 between ksn-2 and combination
Case = 'FreemNu';
if strcmp(Case,'FixmNu')
    xmin = max([min(dCombi.mNu4Sq_contour_12),min(d.mNu4Sq_contour)]);
    xmax = min([max(dCombi.mNu4Sq_contour_12),max(d.mNu4Sq_contour)]);
elseif strcmp(Case,'FreemNu')
    xmin = max([min(dCombimNu.mNu4Sq_k12_contour),min(dmNu.mNu4Sq_contour)]);
    xmax = min([max(dCombimNu.mNu4Sq_k12_contour),max(dmNu.mNu4Sq_contour)]);
end
mNu4Sq = logspace(log10(xmin),log10(xmax),1e3);
if strcmp(Case,'FixmNu')
    sin2T4  = interp1(d.mNu4Sq_contour,d.sin2T4_contour,mNu4Sq,'spline');
    sin2T4C = interp1(dCombi.mNu4Sq_contour_12,dCombi.sin2T4_contour_12,mNu4Sq,'spline');
elseif strcmp(Case,'FreemNu')
    sin2T4  = interp1(dmNu.mNu4Sq_contour,dmNu.sin2T4_contour,mNu4Sq,'spline');
    sin2T4C = interp1(dCombimNu.mNu4Sq_k12_contour,dCombimNu.sin2T4_k12_contour,mNu4Sq,'spline');
end
GetFigure;
plot(sin2T4-sin2T4C,mNu4Sq,'-','LineWidth',2);
set(gca,'YScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2(KSN2) - |{\\itU}_{e4}|^2(Combi)'));
ylabel(sprintf('{\\itm}_4^2 (eV^{ 2})'));
PrettyFigureFormat('FontSize',22);
