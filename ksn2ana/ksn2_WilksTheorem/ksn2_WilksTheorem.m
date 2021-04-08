% Test of Wilk's theorem (coverage)
% plot best fits on contour plot
NrandMC = 50;%752;
Twin_sin2T4 = 0;
Twin_mNu4Sq = 0;
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_%.0fsamples.mat',savedir,NrandMC);
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_%.0fsamples.mat',savedir,Twin_mNu4Sq,Twin_sin2T4,NrandMC);
end

if exist(savefile,'file')
    load(savefile);
    fprintf('load file from %s \n',savefile);
else
     fprintf('file does not exist: %s \n',savefile);
     return
end


%% plot with best fits
yedge = sort(mNu4Sq_bf);
xedge = sort(sin2T4_bf);
GetFigure;
%h = histogram2(sin2T4_bf,mNu4Sq_bf,xedge,yedge,'FaceColor','flat','Normalization','probability');
hold on
%h = scatter(sin2T4_bf(chi2_bf>50),mNu4Sq_bf(chi2_bf>50),'MarkerFaceColor',rgb('Orange'));
hold on;
h = scatter(sin2T4_bf(chi2_bf<50),mNu4Sq_bf(chi2_bf<50),'MarkerFaceColor',rgb('DodgerBlue'));
hold on;
pAsimov = plot3(sin2T4_contour_Asimov',mNu4Sq_contour_Asimov',ones(numel(mNu4Sq_contour_Asimov),1),'k-','LineWidth',2);
view([0 0 1])
grid off
c = colorbar;
colormap('cool')
c.Label.String = 'Best fit probability';
c.Label.FontSize = get(gca,'FontSize');
set(gca,'YScale','log');
set(gca,'XScale','log');
xlabel('|U_{e4}|^2');
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));

xlim([5e-03,0.5]);
PrettyFigureFormat;
leg = legend([h,pAsimov],sprintf('Best fits %.0f pseudo-experiments',NrandMC),sprintf('Sensitivity (%.0f%% C.L.)',95),'EdgeColor',rgb('Silver'),'Location','southwest');
PrettyLegendFormat(leg);
