% plot analyzing plance mangetic field strength
savedir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/results/'];
savefile = sprintf('%sknm1_AnaPotentials.mat',savedir);
pltdir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/plots/'];
if exist(savefile,'file')
    load(savefile);
else
    savefileModel = sprintf('%sknm1_RunwiseFits.mat',savedir);
    load(savefileModel);
    R.ReadSingleRunData;
    
    PixList   = R.PixList; % active pixel
    Bana      = mean(R.SingleRunData.MACE_Ba_T,2); % pixelwise field
    
    qU        = squeeze(R.SingleRunData.qU(25,1,:)); % pixelwise qU-offset
    qU         = qU-mean(qU);
    save(savefile,'qU','Bana','PixList');
    return
end

Bana_plt = NaN.*ones(148,1);
Bana_plt(PixList) = Bana(PixList);

qU_plt = NaN.*ones(148,1);
qU_plt(PixList) = qU(PixList);


fprintf('Bana std = %.2g T , max. variation = %.2g T\n',nanstd(Bana_plt),max(Bana_plt)-min(Bana_plt));
fprintf('qU   std = %.2g mV, max. variation = %.1f mV\n \n',1e3*nanstd(qU_plt),1e3*(max(qU_plt)-min(qU_plt)));
%%
close all
[p1, cb] = FPDViewer(Bana_plt,'ReDrawSkeleton','ON','Label','OFF');
ax = gca;
ax.Position(1) = 0.05;
ax.FontSize = 20;
cb.Label.String = sprintf('{\\itB}_{ana.} (T)');
cb.Label.FontSize = ax.FontSize;
cb.Position(3) = 0.025;
cb.Position(1) = 0.76;
cb.Position(1)
pltname = sprintf('%sknm1_FPD_Bana.pdf',pltdir);
%export_fig(pltname);
%%
close all

[p1, cb] = FPDViewer(1e3*qU_plt,'ReDrawSkeleton','ON','Label','OFF');
ax = gca;
ax.Position(1) = 0.05;
ax.FontSize = 20;
cb.Label.String = sprintf('\\Delta{\\itqU} (mV)');
cb.Label.FontSize = ax.FontSize;
cb.Position(3) = 0.025;
cb.Position(1) = 0.76;
cb.Position(1)
pltname = sprintf('%sknm1_FPD_qU.pdf',pltdir);
export_fig(pltname);

