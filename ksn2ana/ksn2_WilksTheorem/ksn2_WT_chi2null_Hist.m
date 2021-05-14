% Test of Wilk's theorem (coverage)
% chi2 distribution of null hypothesis
Hypothesis = 'H1';
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
         randMC_new  = 1:1251;
    case 'H1' 
        randMC = [1:129,578:748];
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
        chi2 = 'chi2CMShape';
        MergeNew = 'OFF'; % nothing new
end

if strcmp(MergeNew,'ON')
    MergeStr = sprintf('_MergeNew%.0f',numel(randMC_new));
      NrandMC = numel(randMC)+numel(randMC_new);
else
    NrandMC = numel(randMC);
end

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];

if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_Interp%s_%.0fsamples%s_RmDouble%s.mat',...
        savedir,InterpMode,numel(randMC),MergeStr,RmDuplicates);
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_Interp%s_%.0fsamples.mat',...
        savedir,Twin_mNu4Sq,Twin_sin2T4,InterpMode,numel(randMC));
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
dof = 25;
% if strcmp(Hypothesis,'H1')
%     chi2_null =  ReCalc_chi2Null_i;
% end
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);

GetFigure;
hchi2 = histogram(chi2_null,'BinWidth',3,...
    'FaceAlpha',1,'FaceColor',rgb('DeepSkyBlue'),...
    'EdgeColor',rgb('SteelBlue'),'Normalization','probability');
hold on;
x = linspace(0,dof*3,1e3);
y = chi2pdf(x,dof);
pchi2 = plot(x,y*hchi2.BinWidth,'Color',rgb('Black'),'LineWidth',2);

PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\chi^2_{NH} (%.0f dof)',dof));
ylabel('Frequency');
%title(resultsStr,'FontWeight','normal','FontSize',get(gca,'FontSize'))
leg = legend([hchi2,pchi2],sprintf('%.0f pseudo-experiments',numel(chi2_bf)-1),...
                 sprintf('\\chi^2 distribution for %.0f dof',dof),...
                 'EdgeColor',rgb('Silver'));
%xlim([0 70]);

if strcmp(SavePlt,'ON')
plotnameChi2 = strrep(strrep(savefile,'results','plots'),'.mat','_Chi2NullDist.png');
print(gcf,plotnameChi2,'-dpng','-r450');
fprintf('save plot to %s \n',plotnameChi2);
%export_fig(gcf,strrep(plotnameChi2,'.png','.pdf'));
end
