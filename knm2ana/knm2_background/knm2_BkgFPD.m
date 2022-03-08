
savedir = [getenv('SamakPath'),'knm2ana/knm2_background/results/'];
savename = sprintf('%sknm2_BkgFPD.mat',savedir);

if exist(savename,'file')
    load(savename)
else
    savedirf = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
    savenamef = sprintf('%sknm2ubfinal_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s_SysBudget40.mat',...
        savedirf,3,'Real',40,'mNuE0BkgNorm','chi2CMShape','StackPixel','KNM2_0p1eV');
    d = importdata(savenamef);
    
    d.A.ReadSingleRunData;
    BkgRate_subrun_run_pixel = d.A.SingleRunData.TBDIS(d.A.RunData.qU>18574,:,:)./(d.A.SingleRunData.qUfrac(d.A.RunData.qU>18574,:,:).*d.A.SingleRunData.TimeSec);
    BkgRate_subrun_pixel = squeeze(mean(BkgRate_subrun_run_pixel,2));
    BkgRate_pixel = squeeze(mean(BkgRate_subrun_pixel));
    
    PixList = d.A.PixList;
    save(savename,'PixList','BkgRate_pixel','BkgRate_subrun_pixel','BkgRate_subrun_run_pixel');
end

BkgRate_mean = sum(BkgRate_pixel(PixList));

% display inactive pixel in grey (-> NaN)
BkgRate_plt = NaN.*zeros(148,1);
BkgRate_plt(PixList) = BkgRate_pixel(PixList);

%% plot
close all
[plotHandle, cbHandle] = FPDViewer(1e3.*BkgRate_plt,'Label','OFF','ReDrawSkeleton','ON');
cbHandle.Label.String = 'Background rate (mcps)';
cbHandle.Label.FontSize = get(gca,'FontSize')+4;

pltdir= [getenv('SamakPath'),'knm2ana/knm2_background/plots/'];
MakeDir(pltdir);
pltname = sprintf('%sknm2_FPD_background.pdf',pltdir);
export_fig(pltname);
fprintf('save plot to %s \n',pltname);
