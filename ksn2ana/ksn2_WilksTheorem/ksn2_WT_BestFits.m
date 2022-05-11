% Test of Wilk's theorem (coverage)
% plot best fits on contour plot
Hypothesis = 'H0';
InterpMode = 'lin';
SavePlt = 'ON';
MergeNew = 'ON';
RmDuplicates = 'ON';

switch Hypothesis
    case 'H0'
        randMC =[1001:1260,1294:1300,1349:1500];%11:1e3;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        chi2 = 'chi2CMShape';
        randMC_new  = [1:1250];
    case 'H1'
        randMC = [1:1500];
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
        chi2 = 'chi2CMShape';
        MergeNew = 'OFF'; % nothing new
        NrandMC = numel(randMC);
    case 'H2'
        randMC      = 1:1500;
        Twin_sin2T4 = 0.07;
        Twin_mNu4Sq = 20;
        chi2        = 'chi2CMShape';
        MergeNew    = 'OFF'; % nothing new
          NrandMC = numel(randMC);
end

if strcmp(MergeNew,'ON')
    MergeStr = sprintf('_MergeNew%.0f',numel(randMC_new));
else
    MergeStr = '';
end

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];

if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_Interp%s_%.0fsamples%s_RmDouble%s.mat',...
        savedir,InterpMode,numel(randMC),MergeStr,RmDuplicates);
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_Interp%s_%.0fsamples.mat',...
        savedir,Twin_mNu4Sq,Twin_sin2T4,InterpMode,NrandMC);
end

if exist(savefile,'file')
    load(savefile);
    fprintf('load file from %s \n',savefile);
else
     fprintf('file does not exist: %s \n',savefile);
     return
end

   NrandMC = numel(chi2_bf);
%% plot with best fits
yedge = sort(mNu4Sq_bf);
xedge = sort(sin2T4_bf);
GetFigure;
h = histogram2(sin2T4_bf,mNu4Sq_bf,xedge,yedge,'FaceColor','flat','Normalization','probability');
hold on
%h = scatter(sin2T4_bf(ClosedLog95),mNu4Sq_bf(ClosedLog95),'MarkerFaceColor','none');%rgb('Orange'));
%hold on;
%h = scatter(sin2T4_bf(chi2_bf<50),mNu4Sq_bf(chi2_bf<50),'MarkerFaceColor',rgb('DodgerBlue'));
%hold on;
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
leg = legend([h,pAsimov],sprintf('Best fits %.0f pseudo-experiments',NrandMC),sprintf('KNM2 sensitivity at %.0f%% C.L.',95),'EdgeColor',rgb('Silver'),'Location','southwest');
PrettyLegendFormat(leg);

if strcmp(SavePlt,'ON')
    %plotnameContourBf = strrep(strrep(savefile,'results','plots'),'.mat','_BestFits.png');
    pltdir= strrep(savedir,'results','plots');
    plotnameContourBf = sprintf('%sksn2_WT_BestFits%s',pltdir,Hypothesis);
    print(gcf,plotnameContourBf,'-dpng','-r450');
    fprintf('save plot to %s \n',plotnameContourBf);
end
