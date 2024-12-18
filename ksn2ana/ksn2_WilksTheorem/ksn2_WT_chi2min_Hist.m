% Test of Wilk's theorem (coverage)
% chi2 distribution of best fits
Hypothesis = 'H2';
MergeNew = 'ON';
RmDuplicates = 'ON';
InterpMode = 'lin';
switch Hypothesis
    case 'H0'
        NrandMC = 419;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        randMC_new  = 1:1250;
        MergeStr = sprintf('_MergeNew%.0f',numel(randMC_new));
    case 'H1'
         randMC = [1:1500];
       % excl = [1:139,577:757];
       % randMC = randMC(~ismember(randMC,excl));
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
end

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_Interp%s_%.0fsamples%s_RmDouble%s.mat',...
        savedir,InterpMode,NrandMC,MergeStr,RmDuplicates);
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

%% plot
% tmp
dof = 23;
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);

GetFigure;
hchi2 = histogram(chi2_bf,'BinWidth',3,...
    'FaceAlpha',1,'FaceColor',rgb('DeepSkyBlue'),'EdgeColor',rgb('SteelBlue'),'Normalization','probability');
hold on;
x = linspace(0,dof*3,1e3);
y = chi2pdf(x,dof);
pchi2 = plot(x,y*hchi2.BinWidth,'Color',rgb('Black'),'LineWidth',2);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\chi^2_{min} (%.0f dof)',dof));
ylabel('Frequency');
%title(resultsStr,'FontWeight','normal','FontSize',get(gca,'FontSize'))
leg = legend([hchi2,pchi2],sprintf('%.0f pseudo-experiments',numel(chi2_bf)),...
                 sprintf('\\chi^2 distribution for %.0f dof',dof),...
                 'EdgeColor',rgb('Silver'));
%xlim([0 70]);

plotnameChi2 = strrep(strrep(savefile,'results','plots'),'.mat','_Chi2minDist.png');
print(gcf,plotnameChi2,'-dpng','-r450');
fprintf('save plot to %s \n',plotnameChi2);
%export_fig(gcf,strrep(plotnameChi2,'.png','.pdf'));
