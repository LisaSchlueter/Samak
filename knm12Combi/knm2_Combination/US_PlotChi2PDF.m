% plot chi2 distribution for illustration

dof = 1;
chi2min = 1.92;% 21.7;

x = linspace(0,5,1e3);
y = chi2pdf(x,dof);
%%
close all
GetFigure
a = area(x(x>chi2min),y(x>chi2min),'FaceColor',rgb('LightGray'),'EdgeColor','none');
hold on;
plot(x,y,'LineWidth',2,'Color',rgb('Black'));
pchi2min = plot([chi2min,chi2min],[0,chi2pdf(chi2min,dof)],':','LineWidth',2,'Color',rgb('Gray'));
xlabel(sprintf('\\chi^2 (%.0f dof)',dof));
ylabel('Probability density');
leg = legend([pchi2min,a],sprintf('\\chi^2_{min} = %.1f (%.0f dof) ',chi2min,dof),...
    sprintf('p = %.2f',1-chi2cdf(chi2min,dof)));
PrettyFigureFormat('FontSize',22);
PrettyLegendFormat(leg);
savedir = [getenv('SamakPath'),'knm2ana/knm2_Combination/plots/'];
MakeDir(savedir);
savename = sprintf('%sPlotChi2PDF_%.0fdof_chi2min%.1f.png',savedir,dof,chi2min);
print(gcf,savename,'-dpng','-r300');

