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
%%
pltdir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/plots/'];
pltname = sprintf('%sknm1_FPD_Background.pdf',pltdir);
export_fig(pltname);
