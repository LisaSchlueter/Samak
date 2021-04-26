% Test of Wilk's theorem (coverage)
% chi2 distribution of best fits
Hypothesis = 'H1';
switch Hypothesis
    case 'H0'
        NrandMC = 1e3;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
    case 'H1'
        NrandMC = 751;
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
end
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

%% plot
% tmp
dof = 2;

plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);

GetFigure;
deltachi2 = chi2_null-chi2_bf;
% deltachi2(deltachi2<0) = 0;
hchi2 = histogram(deltachi2,'BinWidth',0.3,...
    'FaceAlpha',1,'FaceColor',rgb('DeepSkyBlue'),'EdgeColor',rgb('SteelBlue'),'Normalization','probability');
hold on;
x = linspace(0,dof*10,1e3);
y = chi2pdf(x,dof);
pchi2 = plot(x,y*hchi2.BinWidth,'Color',rgb('Black'),'LineWidth',2);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\chi^2_{null} - \\chi^2_{min} (%.0f dof)',dof));
ylabel('Frequency');
%title(resultsStr,'FontWeight','normal','FontSize',get(gca,'FontSize'))
leg = legend([hchi2,pchi2],sprintf('%.0f pseudo-experiments',numel(chi2_bf)),...
                 sprintf('\\chi^2 distribution for %.0f dof',dof),...
                 'EdgeColor',rgb('Silver'));
%xlim([0 15]);

plotnameChi2 = strrep(strrep(savefile,'results','plots'),'.mat','_DeltaChi2Dist.png');
%print(gcf,plotnameChi2,'-dpng','-r450');
fprintf('save plot to %s \n',plotnameChi2);
%export_fig(gcf,strrep(plotnameChi2,'.png','.pdf'));
