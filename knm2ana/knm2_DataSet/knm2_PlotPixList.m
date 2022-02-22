% FPD viewer of golden pixels

% get Data
savedir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savename = sprintf('%sknm2ubfinal_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s_SysBudget40.mat',...
    savedir,3,'Real',40,'mNuE0BkgNorm','chi2CMShape','StackPixel','KNM2_0p1eV');
d = importdata(savename);


Pixels = NaN.*zeros(148,1);
Pixels(d.A.PixList) = 1;

[plotHandle, cbHandle] = FPDViewer(Pixels,'Label','ON','ReDrawSkeleton','ON');
cbHandle.delete;
%colormap(summer);
%colormap(spring);
ax = gca;

xlabel('Golden pixels');

ax.XLabel.FontSize = 24;
pltdir = [getenv('SamakPath'),'knm2ana/knm2_DataSet/plots/'];
export_fig(gcf,[pltdir,'FPDView_GoldenPixels.pdf']);
