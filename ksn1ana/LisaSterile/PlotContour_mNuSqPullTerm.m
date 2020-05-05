% plot script for KSN1 sterile contour plot
% free neutrino mass: with and without pull term
% Lisa, April 2020
%% settings
CL = 95;
SavePlot = 'ON';
PlotSplines = 'ON'; % plot smooth spline instead of points
range = 95;%
nGridSteps = 50;
chi2Str = 'chi2CMShape';
DataType = 'Real';
freePar = 'mNu E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'OFF';
pullFlag = {99, 12}; %%everything above 12 means no pull
BestFit = 'ON';
%% load grid (or calculate if doesn't exist)
GridArg = {'range',range,...
    'nGridSteps',nGridSteps,...
    'chi2',chi2Str,...
    'RunList',RunList,...
    'SmartGrid',SmartGrid,...
    'RecomputeFlag','OFF',...
    'DataType',DataType};

%% chi2 grids
[mnu4Sq_p1,sin2T4_p1,chi2_p1,chi2_ref_p1,savefile_p1] = KSN1GridSearch(...
    'pullFlag',pullFlag{1},GridArg{:},'freePar',freePar);
[mnu4Sq_p2,sin2T4_p2,chi2_p2,chi2_ref_p2,savefile_p2] = KSN1GridSearch(...
    'pullFlag',pullFlag{2},GridArg{:},'freePar',freePar);
    
[mnu4Sq_free,sin2T4_free,chi2_free,chi2_ref_free,savefile_free] = KSN1GridSearch(...
    'pullFlag',99,GridArg{:},'freePar','E0 Bkg Norm');
%% plot
[p1,~] = KSN1ContourPlot('mnu4Sq',mnu4Sq_p1,'sin2T4',sin2T4_p1,...
    'chi2',chi2_p1,'chi2_ref',chi2_ref_p1,'CL',CL,...
    'HoldOn','OFF','PlotSplines','ON','Color','DodgerBlue','LineStyle','-','BestFit',BestFit);
hold on;
[p2,~] = KSN1ContourPlot('mnu4Sq',mnu4Sq_p2,'sin2T4',sin2T4_p2,...
    'chi2',chi2_p2,'chi2_ref',chi2_ref_p2,'CL',CL,...
    'HoldOn','ON','PlotSplines','ON','Color','Orange','LineStyle','-.','BestFit',BestFit);
[pfree,~] = KSN1ContourPlot('mnu4Sq',mnu4Sq_free,'sin2T4',sin2T4_free,...
    'chi2',chi2_free,'chi2_ref',chi2_ref_free,'CL',CL,...
    'HoldOn','ON','PlotSplines','ON','Color','ForestGreen','LineStyle',':','BestFit',BestFit);

%% find best fits
[row, col] = find(chi2_p1 == min(chi2_p1(:)));
mnu4Sq_p1BF =  mnu4Sq_p1(col,row);
sin2T4_p1BF = sin2T4_p1(col,row);
[row, col] = find(chi2_p2 == min(chi2_p2(:)));
mnu4Sq_p2BF =  mnu4Sq_p2(col,row);
sin2T4_p2BF = sin2T4_p2(col,row);

[row, col] = find(chi2_free == min(chi2_free(:)));
mnu4Sq_freeBF =  mnu4Sq_free(col,row);
sin2T4_freeBF = sin2T4_free(col,row);
%% plot contour
if range==40
    xlim([1e-02 0.4])
elseif range>=90
    xlim([4e-03 0.4])
end
ylim([1 (range+5)^2])

if strcmp(chi2Str,'chi2Stat')
    chi2Label = 'stat. only';
else
    chi2Label = 'stat. and syst.';
end

%% plot label
savedir = [getenv('SamakPath'),'ksn1ana/LisaSterile/results/'];
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = strrep(strrep(savefile_p1,'results','plots'),'.mat','_mNuSqPullTerm.png');

%% legend and title
if strcmp(DataType,'Real')
    DataStr = 'Data';
else
    DataStr = 'Twin';
end
if strcmp(BestFit,'ON')
    legend([p1,p2,pfree],...
        sprintf('Free {\\itm}_\\nu^2 without pull term                       - best fit {\\itm}_4 = %.0f eV, |U_{e4}|^2 = %.3f',sqrt(mnu4Sq_p1BF),sin2T4_p1BF),...
        sprintf('Free {\\itm}_\\nu^2 with pull term \\sigma({\\itm}_\\nu^2) = 1.94 eV^2 - best fit {\\itm}_4 = %.0f eV, |U_{e4}|^2 = %.3f',sqrt(mnu4Sq_p2BF),sin2T4_p2BF),...
        sprintf('Fixed {\\itm}_\\nu^2 = 0 eV^2                                     - best fit {\\itm}_4 = %.0f eV, |U_{e4}|^2 = %.3f',sqrt(mnu4Sq_freeBF),sin2T4_freeBF),...
       'EdgeColor',rgb('Silver'),'Location','southwest');
else
     legend([p1,p2,pfree],...
        sprintf('Free {\\itm}_\\nu^2 without pull term'),...
        sprintf('Free {\\itm}_\\nu^2 with pull term \\sigma({\\itm}_\\nu^2) = 1.94 eV^2'),...
        sprintf('Fixed {\\itm}_\\nu^2 = 0 eV^2'),...
       'EdgeColor',rgb('Silver'),'Location','southwest');
end
title(sprintf('%s %s , %.0f eV at %.0f%% C.L.',DataStr,chi2Label,range,CL),'FontWeight','normal');

%% save plot
if strcmp(SavePlot,'ON')
    print(plotname,'-dpng','-r450');
    fprintf('save plot to %s \n',plotname);
end
