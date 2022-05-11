% Test of Wilk's theorem (coverage)
% chi2 distribution of best fits
Hypothesis = 'H2';
InterpMode = 'Mix';
SavePlt = 'ON';
MergeNew = 'ON';
RmDuplicates = 'ON';

switch Hypothesis
    case 'H0'
        randMC =[1001:1260,1294:1300,1349:1500];%11:1e3;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        chi2 = 'chi2CMShape';
        randMC_new  = 1:1250;
    case 'H1' 
       randMC = 1:1500;
       Twin_sin2T4 = 0.0240;
       Twin_mNu4Sq = 92.7;
       chi2 = 'chi2CMShape';
       MergeNew = 'OFF'; % nothing new
    case 'H2'
        randMC      = 1:1500;
        Twin_sin2T4 = 0.07;
        Twin_mNu4Sq = 20;
        chi2        = 'chi2CMShape';
        MergeNew    = 'OFF'; % nothing new
end

if strcmp(MergeNew,'ON')
    MergeStr = sprintf('_MergeNew%.0f',numel(randMC_new));
      NrandMC = numel(randMC)+numel(randMC_new);
else
    MergeStr = '';
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
dof = 2;
chi2_delta(chi2_delta<0) = 0;
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
%chi2_delta = ReCalc_chi2Null_i- chi2_bf;
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.45]);
% chi2_delta(chi2_delta<0) = 0;
hchi2 = histogram(chi2_delta,'BinWidth',0.4,...
    'FaceAlpha',1,'FaceColor',rgb('DeepSkyBlue'),'EdgeColor',rgb('SteelBlue'),'Normalization','probability');
hold on;

x = linspace(0,dof*10,1e3);
y = chi2pdf(x,dof);
pchi2 = plot(x,y*hchi2.BinWidth,'Color',rgb('Black'),'LineWidth',2);

xlabel(sprintf('\\Delta\\chi^2'));
ylabel('Frequency');
PrettyFigureFormat('FontSize',24);
%title(resultsStr,'FontWeight','normal','FontSize',get(gca,'FontSize'))
pNone = plot(NaN,NaN,'w','LineStyle','none');
leg = legend([hchi2,pNone,pchi2],sprintf('%.0f KNM2 pseudo-experiments',numel(chi2_bf)),...
    sprintf('with MC truth: {\\itm}_4^2 = %.3g eV^2, |{\\itU}_{e4}|^2 = %.2g',Twin_mNu4Sq,Twin_sin2T4),...
                 sprintf('Chi-squared distribution for %.0f dof',dof),...
                 'EdgeColor',rgb('Silver'));
PrettyLegendFormat(leg);
leg.FontSize = get(gca,'FontSize')+2;
xlim([0 15]);
%%
%pltname = strrep(strrep(savefile,'results','plots'),'.mat','_DeltaChi2Dist.png');
%print(gcf,pltname,'-dpng','-r450');
 pltdir= strrep(savedir,'results','plots');
 pltname = sprintf('%sksn2_WT_DeltaChi2Hist%s.pdf',pltdir,Hypothesis);
   export_fig(pltname);
    

fprintf('save plot to %s \n',pltname);
%export_fig(gcf,strrep(plotnameChi2,'.png','.pdf'));
%%

 