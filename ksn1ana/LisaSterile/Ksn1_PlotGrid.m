SavePlot = 'ON';
PlotContour = 'OFF';
range =95;%
nGridSteps =50;%25;%50;
chi2Str = 'chi2CMShape';
DataType = 'Real';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'OFF';
CL = 0.82;
[mnu4Sq,sin2T4,chi2,chi2_ref,savefile] = KSN1GridSearch('range',range,...
    'nGridSteps',nGridSteps,...
    'chi2',chi2Str,...
    'DataType',DataType,...
    'freePar',freePar,...
    'RunList',RunList,...
    'SmartGrid',SmartGrid,...
    'RecomputeFlag','OFF');

%

GetFigure;
DeltaChi2 = GetDeltaChi2(CL,2);
chi2_ref = min(min(chi2));
%chi2_ref = chi2(1,1); 
chi2((chi2-chi2_ref)>DeltaChi2) =  DeltaChi2+chi2_ref;% NaN;
%chi2((chi2-chi2_ref)<0) =  NaN;% -1+chi2_ref;% NaN;
surf(sin2T4,mnu4Sq,chi2'-chi2_ref,'EdgeColor','interp','FaceColor','interp');
zlim([0 DeltaChi2])
set(gca,'XScale','log')
set(gca,'YScale','log')
c =colorbar;
c.Label.String = sprintf('\\Delta\\chi^2');
PrettyFigureFormat('FontSize',24)
c.Label.FontSize = get(gca,'FontSize')+2;
c.Limits=[0 DeltaChi2];
xlabel('|U_{e4}|^2');
ylabel(sprintf('{\\itm}_4 (eV^2)'));
zlabel(sprintf('\\Delta\\chi^2'))
view([0 0 1])
 grid off
if range==40
    xlim([1e-02 0.5])
elseif range>=90
    xlim([1e-03 0.5])
end
ylim([1 (range+5)^2]);
if strcmp(chi2Str,'chi2Stat')
    chi2Label = 'stat. only';
else
    chi2Label = 'stat. and syst.';
end
if strcmp(DataType,'Real')
    DataLabel = 'Data';
else
    DataLabel = 'Twin';
end
 title(sprintf('%s %.0f eV range (%s)',DataLabel,range,chi2Label),...
     'FontWeight','normal','FontSize',get(gca,'FontSize'));
 
 if strcmp(PlotContour,'ON')
     hold on;
      [mnu4Sq_contour, sin2T4_contour] = ...
        KSN1Grid2Contour(mnu4Sq,sin2T4,chi2,chi2_ref,CL);
    pHandle  = plot(sin2T4_contour,mnu4Sq_contour,'k-','LineWidth',2,'MarkerSize',20);
 end
 
 if strcmp(SavePlot,'ON')
     plotname = strrep(strrep(savefile,'results','plots'),'.mat','_Grid.png');
     print(gcf,plotname,'-dpng','-r450');
 end

 fprintf('save plot to %s \n',plotname)