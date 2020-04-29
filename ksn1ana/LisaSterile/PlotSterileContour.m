% plot script for KSN1 sterile contour plot
% Lisa, April 2020
Nathan = 'ON'; % compare with nathan
CL = 90;
% load contour
savedir = [getenv('SamakPath'),'ksn1ana/LisaSterile/results/'];
MakeDir(savedir);
range = 95;%
nGridSteps = 50;
NonPoissonScaleFactor=1;
chi2 = 'chi2CMShape';
DataType = 'Twin';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
savefile = sprintf('%sSterileTestPar_%s_%s_%s_%.0feVrange_%s_%.0fGridSteps.mat',...
    savedir,RunList,DataType,strrep(freePar,' ',''),range,chi2,nGridSteps);
if exist(savefile,'file')
    load(savefile)
else
    fprintf('file doesnt exist \n');
    return
end

[mnu4Sq_contour, sin2T4_contour] = GetSterileContour(mnu4Sq,sin2T4,chi2Grid,chi2_ref,CL);

% plot contour
GetFigure;
p1 = plot(sin2T4_contour,mnu4Sq_contour,'-','LineWidth',2,'MarkerSize',20);
set(gca,'YScale','log');
set(gca,'XScale','log');
if range==40
   xlim([1e-02 0.5])
elseif range==90
    xlim([5e-03 0.5])
end
PrettyFigureFormat;
xlabel('|U_{e4}|^2');
ylabel(sprintf('{\\itm}_4 (eV^2)'));
if strcmp(chi2,'chi2Stat')
    chi2Str = 'stat. only';
else
     chi2Str = 'stat. and syst.';
end

plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = strrep(strrep(savefile,'results','plots'),'.mat','.png');
if strcmp(Nathan,'ON')
    d = importdata(sprintf('NathanContour_%.0feV_%s_%s.mat',range,DataType,chi2));
    hold on;
    pN = plot(d.sith4_X,d.m4_Y,'-.','LineWidth',2);
    leg = legend('Lisa','Nathan','EdgeColor',rgb('Silver'),'Location','southwest');
    plotname = strrep(plotname,'.png','_CompareN.png');
end
title(sprintf('%s , %s , %.0f eV at %.0f%% C.L.',DataType,chi2Str,range,CL),'FontWeight','normal');

print(plotname,'-dpng','-r450');
fprintf('save plot to %s \n',plotname);