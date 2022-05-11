range = 40;%
SavePlt = 'ON';
chi2Str = 'chi2CMShape';
InterpMode = 'lin';
savedir = [getenv('SamakPath'),'ksn1ana/ksn1_WilksTheorem/results/'];
MakeDir(savedir);
savename = sprintf('%sksn1_WilksTheorem_%.0frange_%s_%s.mat',savedir,range,chi2Str,InterpMode);

if exist(savename,'file')
    d = importdata(savename);
    fprintf('load %s\n',savename);
else
    fprintf(2,'file not found, run ksn1_WT_MergeFiles.m \n');
    return
end
mNu4Sq_bf = d.mnu4Sq_bf;
sin2T4_bf = d.sin2T4_bf;

%% load asimov contour
savedirCombi = sprintf('%sksn2ana/ksn2_RunCombination/results/',getenv('SamakPath'));
savefileCombimNuFix = sprintf('%sksn21_Combination_ReAna_%s.mat',savedirCombi,'Twin');                            
dCombi = importdata(savefileCombimNuFix);
fprintf('load %s \n',savefileCombimNuFix);
mNu4Sq_contour_Asimov = dCombi.mNu4Sq_contour_1;
sin2T4_contour_Asimov = dCombi.sin2T4_contour_1;


yedge = sort(mNu4Sq_bf);
xedge = sort(sin2T4_bf);
GetFigure;
h = histogram2(sin2T4_bf,mNu4Sq_bf,xedge,yedge,'FaceColor','flat','Normalization','probability');
hold on
pAsimov = plot3(sin2T4_contour_Asimov',mNu4Sq_contour_Asimov',ones(numel(mNu4Sq_contour_Asimov),1),'k-','LineWidth',2);
view([0 0 1])
grid off
c = colorbar;
colormap('cool')
PrettyFigureFormat('FontSize',22);
c.Label.String = 'Best fit probability';
c.Label.FontSize = get(gca,'FontSize');
set(gca,'YScale','log');
set(gca,'XScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));

xlim([1e-03,0.5]);
ylim([0.1,2000]);
leg = legend([h,pAsimov],sprintf('Best fits %.0f pseudo-experiments',numel(sin2T4_bf)),sprintf('KNM1 sensitivity at %.0f%% C.L.',95),'EdgeColor',rgb('Silver'),'Location','southwest');
PrettyLegendFormat(leg);
%%
if strcmp(SavePlt,'ON')
    pltdir = strrep(savedir,'results','plots');
    MakeDir(pltdir);
    plotnameContourBf = [pltdir,'ksn1_WT_BestFitsH0'];%strrep(,'.mat','_BestFits.png');
    print(plotnameContourBf,'-dpng','-r450');
    fprintf('save plot to %s \n',plotnameContourBf);

end