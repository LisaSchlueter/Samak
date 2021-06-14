% compute confidence level over DeltaChi^2
nPar = 2;
x = linspace(0,8,1e3);
y = chi2pdf(x,nPar);
xT = 2;

GetFigure
a = area(x(x<=xT),y(x<=xT),'FaceColor',rgb('Silver'),'EdgeColor','none');
hold on;
plot(x,y,'-k','LineWidth',2);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\chi^2 (%.0f dof)',nPar));
ylabel('Probability density');
leg = legend(a,sprintf('Confidence level = %.2f',chi2cdf(xT,nPar)),'Location','northeast');
PrettyLegendFormat(leg);

pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_Significance/plots/'];
MakeDir(pltdir);
pltname = [pltdir,'ksn2_ConfLevelChi2.png'];
print(pltname,'-dpng','-r350');