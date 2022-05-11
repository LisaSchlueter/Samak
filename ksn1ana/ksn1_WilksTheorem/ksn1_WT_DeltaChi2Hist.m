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

chi2_delta = d.chi2min_null-d.chi2min_bf;
%%
dof = 2;
chi2_delta(chi2_delta<0) = 0;
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
%chi2_delta = ReCalc_chi2Null_i- chi2_bf;
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.45]);
% chi2_delta(chi2_delta<0) = 0;
hchi2 = histogram(chi2_delta,'BinWidth',0.3,...
    'FaceAlpha',1,'FaceColor',rgb('DeepSkyBlue'),'EdgeColor',rgb('SteelBlue'),'Normalization','probability');
hold on;
x = linspace(0,dof*10,1e3);
y = chi2pdf(x,dof);
pchi2 = plot(x,y*hchi2.BinWidth,'Color',rgb('Black'),'LineWidth',2);

xlabel(sprintf('\\Delta\\chi^2'));
ylabel('Frequency');
PrettyFigureFormat('FontSize',22);
%title(resultsStr,'FontWeight','normal','FontSize',get(gca,'FontSize'))

             pNone = plot(NaN,NaN,'w','LineStyle','none');
 leg = legend([hchi2,pNone,pchi2],sprintf('%.0f KNM1 pseudo-experiments',numel(chi2_delta)),...
    sprintf('with MC truth: {\\itm}_4^2 = %.3g eV^2, |{\\itU}_{e4}|^2 = %.2g',0,0),...
                 sprintf('Chi-squared distribution for %.0f dof',dof),...
                 'EdgeColor',rgb('Silver'));
             
leg.FontSize = get(gca,'FontSize')+2;
PrettyLegendFormat(leg);
xlim([0 15]);


pltdir= strrep(savedir,'results','plots');
 pltname = sprintf('%sksn1_WT_DeltaChi2HistH0.pdf',pltdir);
   export_fig(pltname);
fprintf('save plot to %s \n',pltname);