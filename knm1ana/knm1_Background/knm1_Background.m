savedir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/results/'];
savefile = sprintf('%sknm1_RunwiseFits.mat',savedir);

if exist(savefile,'file')
    load(savefile);
else
    return
end

BkgIdx = find(R.RunData.qU>18573.7,1);


BkgRate = R.RunData.TBDIS(BkgIdx:end)./(R.RunData.qUfrac(BkgIdx:end).*R.RunData.TimeSec);
%BkgRate       = R.SingleRunData.TBDIS(BkgIdx:end,:)./(R.SingleRunData.qUfrac(BkgIdx:end,:).*R.SingleRunData.TimeSec);
%GetNPfactor(R,'SlopeCorr','OFF','RecomputeFlag','OFF','SanityPlot','ON')

fprintf('Mean background %.1f mcps \n',mean(BkgRate)*1e3)

%% get single pixel
   R.ReadSingleRunData; 
   %%
   BkgRate_subrun_pixel =  R.SingleRunData.TBDIS(BkgIdx:end,:,:)./(R.SingleRunData.qUfrac(BkgIdx:end,:,:).*R.SingleRunData.TimeSec);
   BkgRate_pixel_all =   squeeze(mean(squeeze(mean(BkgRate_subrun_pixel))));
   
   BkgRate_pixel = NaN.*zeros(148,1);
   BkgRate_pixel(R.PixList) = BkgRate_pixel_all(R.PixList);
   
   TimeSec_Bkg = mean(squeeze(sum(sum((R.SingleRunData.qUfrac(BkgIdx:end,:,:).*R.SingleRunData.TimeSec)))));

%% fpd
close all

[p1, cb] = FPDViewer(1e3.*BkgRate_pixel,'ReDrawSkeleton','ON','Label','OFF');
ax = gca;
ax.Position(1) = 0.05;
ax.FontSize = 20;
cb.Label.String = sprintf('Background rate (mcps)');
cb.Label.FontSize = ax.FontSize;
cb.Position(3) = 0.025;
cb.Position(1) = 0.76;
cb.Position(1)

%% group into rings and fit
rStart = [0 0.7398 1.4796 1.9573 2.3394 2.6674 2.9592 3.2247 ...
    3.4699 3.699 3.9146 4.119 4.3137]';
rEnd = [0.7398 1.4796 1.9573 2.3394 2.6674 2.9592 3.2247 ...
    3.4699 3.699 3.9146 4.119 4.3137 4.5]';
Radius_Ring_cm = mean([rStart';rEnd']);

[PixList,RingPixList] = Ring2Pixel(1:12,R.PixList);

BkgRate_Ring    = cellfun(@(x) sum(BkgRate_pixel_all(x))./numel(x),RingPixList);
BkgRateErr_Ring = sqrt(BkgRate_Ring.*TimeSec_Bkg)./TimeSec_Bkg;

[linFitpar, linFiterr, linFitchi2min,linFitdof] = linFit(Radius_Ring_cm(1:12)',BkgRate_Ring,BkgRateErr_Ring);
fprintf('Linear fit Background rate(radius): slope = (%.1f +- %.1f) mcps/cm -> %.1f sigma \n',1e3.*linFitpar(1),1e3.*linFiterr(1),abs(linFitpar(1)./linFiterr(1)));

GetFigure;
e1 = errorbar(Radius_Ring_cm(1:12),1e3.*BkgRate_Ring,1e3.*BkgRateErr_Ring,'.','CapSize',0,'LineWidth',2,'MarkerSize',20);
hold on;
x = linspace(min(Radius_Ring_cm(1:12)),max(Radius_Ring_cm(1:12)),10);
p1=plot(x,1e3*(x.*linFitpar(1)+linFitpar(2)),'LineWidth',2);
xlabel('Radial position (cm)');
ylabel('Background rate per pixel (mcps)');
PrettyFigureFormat;
leg = legend([e1,p1],'Average background rate',sprintf('Linear fit slope = (%.2g \\pm %.2g) mcps per pixel / cm , %.1f sigma \n',1e3.*linFitpar(1),1e3.*linFiterr(1),abs(linFitpar(1)./linFiterr(1))));
PrettyLegendFormat(leg);


%%
pltdir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/plots/'];
pltname = sprintf('%sknm1_FPD_Background.pdf',pltdir);
export_fig(pltname);
