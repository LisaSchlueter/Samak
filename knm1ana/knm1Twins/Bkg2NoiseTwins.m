M = MultiRunAnalysis('RunList','KNM1','exclDataStart',17,'RingList',1:12);
R = RingAnalysis('RunAnaObj',M);
M.Fit;
BkgMeanSubRun = (M.FitResult.par(3)+M.ModelObj.BKG_RateSec_i).*M.ModelObj.qUfrac.*M.ModelObj.TimeSec;
Bkg_i = M.ModelObj.BKG_RateSec_i;
%% Get deviation of background and normalization
R.FitRings;
%%
Nmax = max(R.FitResult.par(:,4));
Nmin = min(R.FitResult.par(:,4));

nPixFit = cell2mat(arrayfun(@(x) numel(x.PixList),R.MultiObj,'UniformOutput',0))'; %number of pixels per ring in fit
Bkg_i_Ring = cell2mat(arrayfun(@(x) x.ModelObj.BKG_RateSec_i,R.MultiObj,'UniformOutput',0))'; %number of pixels per ring in fit
BmaxPix = max((R.FitResult.par(:,3)+Bkg_i_Ring)./nPixFit);
BminPix = min((R.FitResult.par(:,3)+Bkg_i_Ring)./nPixFit);

Bmax = BmaxPix*numel(M.PixList);
Bmin = BminPix*numel(M.PixList);
%% signal to background ratio if normalization changes
M.SimulateStackRuns('BKG_RateAllFPDSec',0); %no background + reset to init value
M.ModelObj.ComputeTBDDS('N_bias',Nmax); M.ModelObj.ComputeTBDIS;
TBDIS_Nmax = M.ModelObj.TBDIS;
PlotArg = {'LineWidth',2};
M.SimulateStackRuns('BKG_RateAllFPDSec',0); %reset to init value
M.ModelObj.ComputeTBDDS('N_bias',Nmin); M.ModelObj.ComputeTBDIS;
TBDIS_Nmin = M.ModelObj.TBDIS;

f32 = figure('Renderer','opengl');
set(f32, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.5, 0.5]);

plot(M.ModelObj.qU(14:end)-18575,TBDIS_Nmax(14:end)./BkgMeanSubRun(14:end),PlotArg{:});
hold on;
plot(M.ModelObj.qU(14:end)-18575,TBDIS_Nmin(14:end)./BkgMeanSubRun(14:end),PlotArg{:});
leg = legend(sprintf('N = %.4f',Nmax+1),sprintf('N = %.4f',Nmin+1)); legend boxoff
PrettyFigureFormat;
xlabel('retarding energy -18575 (eV)');
ylabel('signal / background');
MaxDiff = max(abs(TBDIS_Nmax(14:end)./BkgMeanSubRun(14:end)-TBDIS_Nmin(14:end)./BkgMeanSubRun(14:end)));
title(['background fixed - ',sprintf('max \\Delta(S/B) =%.3g',MaxDiff)]);
hold off;


% save plot
savepath = [getenv('SamakPath'),'knm1ana/knm1Twins/plots/'];
savename = [savepath,'Signal2BkgRatio_normalization'];

if ~exist(savepath,'dir')
    system(['mkdir ',savepath]);
end
print(gcf,savename,'-dpng','-r400');


%% signal to background ratio if background 
M.SimulateStackRuns('BKG_RateAllFPDSec',0); %no background + reset to init value
TBDIS = M.ModelObj.TBDIS;
BkgMaxSubRun = Bmax.*M.ModelObj.qUfrac.*M.ModelObj.TimeSec;
BkgMinSubRun = Bmin.*M.ModelObj.qUfrac.*M.ModelObj.TimeSec;

f33 = figure('Renderer','opengl');
set(f33, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.5, 0.5]);

PlotArg = {'LineWidth',2};
plot(M.ModelObj.qU(14:end)-18575,TBDIS(14:end)./BkgMaxSubRun(14:end),PlotArg{:});
hold on;
plot(M.ModelObj.qU(14:end)-18575,TBDIS(14:end)./BkgMinSubRun(14:end),PlotArg{:});
leg = legend(sprintf('B_{max} = %.4f cps',Bmax),sprintf('B_{min} = %.4f cps',Bmin)); legend boxoff
PrettyFigureFormat;
xlabel('retarding energy -18575 (eV)');
ylabel('signal/background');
MaxDiff2 = max(abs(TBDIS(14:end)./BkgMaxSubRun(14:end)-TBDIS(14:end)./BkgMinSubRun(14:end)));
title(['normalization fixed - ',sprintf('max \\Delta(S/B) =%.3g',MaxDiff2)]);
hold off;

% save plot
savepath = [getenv('SamakPath'),'knm1ana/knm1Twins/plots/'];
savename = [savepath,'Signal2BkgRatio_background'];

if ~exist(savepath,'dir')
    system(['mkdir ',savepath]);
end
print(gcf,savename,'-dpng','-r400');

