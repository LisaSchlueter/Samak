% nu-mass bias as a functio of BKG_PtSlope


% unblinded fit with penning track background slope
range     = 40;
freePar   = 'mNu E0 Bkg Norm';
chi2      = 'chi2Stat';
DataType  = 'Twin';
AnaFlag   = 'StackPixel';
RingMerge = 'Full';%'None';
BKG_PtSlope = (0:6).*1e-06;
TwinBias_BKG_PtSlope = 3*1e-06;
FSDFlag   = 'KNM2';
SysBudget = 40;

% init
mNuSq_Twin    = zeros(numel(BKG_PtSlope),1);
mNuSqErr_Twin = zeros(numel(BKG_PtSlope),1);
chi2min_Twin   = zeros(numel(BKG_PtSlope),1);
mNuSq_Real    = zeros(numel(BKG_PtSlope),1);
mNuSqErr_Real    = zeros(numel(BKG_PtSlope),1);
chi2min_Real   = zeros(numel(BKG_PtSlope),1);

for i=1:numel(BKG_PtSlope)
savedir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];

savenameT = sprintf('%sknm2ubfinal_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s.mat',...
    savedir,BKG_PtSlope(i)*1e6,'Twin',range,strrep(freePar,' ',''),chi2,AnaFlag,FSDFlag);
savenameR = sprintf('%sknm2ubfinal_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s.mat',...
    savedir,BKG_PtSlope(i)*1e6,'Real',range,strrep(freePar,' ',''),chi2,AnaFlag,FSDFlag);

if ~strcmp(chi2,'chi2Stat')
    savenameT = strrep(savenameT,'.mat',sprintf('_SysBudget%.0f.mat',SysBudget));
    savenameR = strrep(savenameR,'.mat',sprintf('_SysBudget%.0f.mat',SysBudget));
end
savenameT = strrep(savenameT,'.mat',sprintf('_TwinBpng-%.1fmucpsPers.mat',1e6*TwinBias_BKG_PtSlope));


dT = importdata(savenameT);
mNuSq_Twin(i) = dT.FitResult.par(1);
mNuSqErr_Twin(i) = 0.5.*(dT.FitResult.errPos(1)-dT.FitResult.errNeg(1));
chi2min_Twin(i)   = dT.FitResult.chi2min;

dR = importdata(savenameR);
mNuSq_Real(i) = dR.FitResult.par(1);
mNuSqErr_Real(i) = 0.5.*(dR.FitResult.errPos(1)-dR.FitResult.errNeg(1));
chi2min_Real(i)   = dR.FitResult.chi2min;

end

%%
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.8]);

subplot(3,1,1:2);
pT = plot(BKG_PtSlope*1e6,mNuSq_Twin,'-.','LineWidth',2,'Color',rgb('DodgerBlue'));
hold on;
pR = plot(BKG_PtSlope*1e6,mNuSq_Real,'LineWidth',2,'Color',rgb('Orange'));
PrettyFigureFormat;
grid on
xlabel(sprintf('Model background time slope (\\mucps/s)'))
ylabel(sprintf('{\\itm}_\\nu^2 (eV^2)'))
leg = legend([pT,pR],sprintf('KNM-2 Twins (\\alpha = %.0f \\mucps/s)',TwinBias_BKG_PtSlope*1e06),'KNM-2 data');
leg.EdgeColor = rgb('Silver');
leg.Location = 'northwest';

subplot(3,1,3);
pT = plot(BKG_PtSlope*1e6,chi2min_Twin+mean(chi2min_Real),'-.','LineWidth',2,'Color',rgb('DodgerBlue'));
hold on;
pR = plot(BKG_PtSlope*1e6,chi2min_Real,'LineWidth',2,'Color',rgb('Orange'));
PrettyFigureFormat;
xlabel(sprintf('Model background time slope (\\mucps/s)'))
ylabel(sprintf('\\chi2_{min}'))
leg = legend(pT,sprintf('\\chi2_{min} + %.1f (KNM-2 Twins)',mean(chi2min_Real)));
legend boxoff
leg.Location = 'southwest';

pltname = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/plots/PngBkg_mnuPlot.png'];
print(pltname,'-dpng','-r300');