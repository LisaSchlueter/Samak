% plot asimov sensitivity vs exclusion

savedirCombi = sprintf('%sksn2ana/ksn2_RunCombination/results/',getenv('SamakPath'));

filenameD= sprintf('%sksn21_Combi_freemNuSq_ReAna_%s.mat',savedirCombi,'Real');
dD = importdata(filenameD);
fprintf('load %s \n',filenameD);

filenameT= sprintf('%sksn21_Combi_freemNuSq_ReAna_%s.mat',savedirCombi,'Twin');
dT = importdata(filenameT);
fprintf('load %s \n',filenameT);


%% 
DataSet = 'KNM1+2';
LocalFontSize = 19;
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.45]);

pNone = plot(NaN,NaN,'.w');
hold on;
switch DataSet
    case 'KNM1'
        pD = plot(dD.sin2T4_k1_contour,dD.mNu4Sq_1_contour,'LineWidth',3,'Color',rgb('Salmon'));
        pT = plot(dT.sin2T4_k1_contour,dT.mNu4Sq_1_contour,'-.','LineWidth',2.5,'Color',rgb('DimGray'));
        plot(dD.sin2T4bf_k1,dD.mNu4Sqbf_k1,'p','MarkerIndices',[1,1],'MarkerFaceColor',pD.Color,...
       'LineWidth',1,'Color',pD.Color,'MarkerSize',17);
    case 'KNM2'
        pD = plot(dD.sin2T4_k2_contour,dD.mNu4Sq_k2_contour,'LineWidth',3,'Color',rgb('Crimson'));
        pT = plot(dT.sin2T4_k2_contour,dT.mNu4Sq_k2_contour,'-.','LineWidth',2.5,'Color',rgb('DimGray'));
         plot(dD.sin2T4bf_k2,dD.mNu4Sqbf_k2,'p','MarkerIndices',[1,1],'MarkerFaceColor',pD.Color,...
       'LineWidth',1,'Color',pD.Color,'MarkerSize',17);
    case 'KNM1+2'
        pD = plot(dD.sin2T4_k12_contour,dD.mNu4Sq_k12_contour,'LineWidth',3,'Color',rgb('DarkRed'));
        pT = plot(dT.sin2T4_k12_contour,dT.mNu4Sq_k12_contour,'-.','LineWidth',2.5,'Color',rgb('DimGray'));
         plot(dD.sin2T4bf_k12,dD.mNu4Sqbf_k12,'p','MarkerIndices',[1,1],'MarkerFaceColor',pD.Color,...
       'LineWidth',1,'Color',pD.Color,'MarkerSize',17);
end
leg = legend([pNone,pT,pD],...
    sprintf('Case II) {\\itm}_\\nu^2 free'),...
    'Sensitivity (Asimov)',...
    sprintf('%s data exclusion',DataSet),...
    'Location','southwest');

set(gca,'YScale','log');
set(gca,'XScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^{2})'));
ylim([5 2000]);
xlim([4e-03 0.5]);

PrettyFigureFormat('FontSize',LocalFontSize)
PrettyLegendFormat(leg);
leg.FontSize = LocalFontSize-2;
%
CombiPltDir = strrep(savedirCombi,'results','plots');
MakeDir(CombiPltDir);
pltname = sprintf('%s%s_DataTwin_mNuSq.pdf',CombiPltDir,DataSet);
export_fig(pltname);
fprintf('save plot to %s \n',pltname);

