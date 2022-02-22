% load new results and plot
NonPoissonScaleFactor = 1;
range = 40;
DataType = 'Real';
savedir = [getenv('SamakPath'),'knm12Combi/knm1ReAnalysis/results/'];
savefile_Common = sprintf('%sknm1_UniformFit_ReAna_%s_chi2Stat_NP%.4g_mNuE0NormBkg_%.0feV',savedir,DataType,NonPoissonScaleFactor,range);

nFiles = 6;
mNuSq       = zeros(nFiles,1);
mNuSqErr    = zeros(nFiles,1);
mNuSqTot    = zeros(2,1);
mNuSqErrToT = zeros(2,1);

% knm1: prl config: stat
f_prl = [savefile_Common,'_SibilleFull_KatrinT2_AngTFOFF.mat'];
d_prl = importdata(f_prl);
mNuSq(1) = d_prl.FitResult.par(1);
mNuSqErr(1) = 0.5*(d_prl.FitResult.errPos(1)-d_prl.FitResult.errNeg(1));

% knm1: non-isotropic angular transmission
f_AngTf = [savefile_Common,'_Sibille0p5eV_KatrinT2_AngTFON.mat'];
d_AngTf = importdata(f_AngTf);
mNuSq(2) = d_AngTf.FitResult.par(1);
mNuSqErr(2) = 0.5*(d_AngTf.FitResult.errPos(1)-d_AngTf.FitResult.errNeg(1));

% knm1: KNM-2 FSD
f_fsd = [savefile_Common,'_KNM2_0p1eV_KatrinT2_AngTFOFF.mat'];
d_fsd = importdata(f_fsd);
mNuSq(3) = d_fsd.FitResult.par(1);
mNuSqErr(3) = 0.5*(d_fsd.FitResult.errPos(1)-d_fsd.FitResult.errNeg(1));

% knm1: new energy-loss
f_el = [savefile_Common,'_Sibille0p5eV_KatrinT2A20_AngTFOFF.mat'];
d_el = importdata(f_el);
mNuSq(4) = d_el.FitResult.par(1);
mNuSqErr(4) = 0.5*(d_el.FitResult.errPos(1)-d_el.FitResult.errNeg(1));

% knm1: penning background
f_pt = [savefile_Common,'_Sibille0p5eV_KatrinT2_AngTFOFF_BkgPtSlope-2.2muCpsS.mat'];
d_pt = importdata(f_pt);
mNuSq(5) = d_pt.FitResult.par(1);
mNuSqErr(5) = 0.5*(d_pt.FitResult.errPos(1)-d_pt.FitResult.errNeg(1));

% knm1: all new in 
f_new = [savefile_Common,'_KNM2_0p1eV_KatrinT2A20_AngTFON_BkgPtSlope-2.2muCpsS.mat'];
d_new = importdata(f_new);
mNuSq(6) = d_new.FitResult.par(1);
mNuSqErr(6) = 0.5*(d_new.FitResult.errPos(1)-d_new.FitResult.errNeg(1));

%% stat & syst
% knm1: prl config stat & syst.
f_prl_Tot = [savefile_Common,'_Sibille0p5eV_KatrinT2_AngTFOFF.mat'];
f_prl_Tot = strrep(strrep(f_prl_Tot,'chi2Stat','chi2CMShape_SysBudget24'),'NP1','NP1.064');
d_prl_tot = importdata(f_prl_Tot);
mNuSqTot(1) = d_prl_tot.FitResult.par(1);
mNuSqErrTot(1) = 0.5*(d_prl_tot.FitResult.errPos(1)-d_prl_tot.FitResult.errNeg(1));

% knm1: all new in  stat & syst.
f_new_Tot = strrep(strrep([savefile_Common,'_KNM2_0p1eV_KatrinT2A20_AngTFON_BkgPtSlope-2.2muCpsS.mat'],'chi2Stat','chi2CMShape_SysBudget200'),'NP1','NP1.064');
d_new_Tot = importdata(f_new_Tot);
mNuSqTot(2) = d_new_Tot.FitResult.par(1);
mNuSqErrTot(2) = 0.5*(d_new_Tot.FitResult.errPos(1)-d_new_Tot.FitResult.errNeg(1));


%% plot
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.6]);
plot(mNuSqTot(1).*[1,1],[0 nFiles+1],':','Color',rgb('Silver'),'LineWidth',2);
hold on;
eStat = errorbar(mNuSq,1:nFiles,zeros(nFiles,1),zeros(nFiles,1),mNuSqErr,mNuSqErr,'.','MarkerSize',17,'CapSize',0,'LineWidth',1.5,'Color',rgb('Orange'));
eTot = errorbar(mNuSqTot,[1,nFiles]+0.1,zeros(2,1),zeros(2,1),mNuSqErrTot,mNuSqErrTot,'.','MarkerSize',17,'CapSize',0,'LineWidth',1.5,'Color',rgb('DodgerBlue'));
PrettyFigureFormat;
ylim([0.5 nFiles+0.9])
yticks(1:nFiles);
xlim([-2.2 0.1])
yticklabels({'Original',sprintf('Non-isotropic transmission'),'Final-state distribution','Energy-loss function',sprintf('{\\itB}_{penning}'),'Reanalysis 2021'});
set(gca,'YMinorTick','off');
xlabel(sprintf('{\\itm}_\\nu^2 (eV^{ 2})'));
leg = legend([eStat,eTot],'Stat. only','Stat. and syst');
PrettyLegendFormat(leg);
leg.NumColumns = 2;
leg.Location = 'north';

ax = gca;
ax.Position(3) = 0.5;
% difference with respect to PRL
DeltamNuSqAll  = mNuSqTot(2)-mNuSqTot(1);
DeltamNuSq = mNuSq(2:end)-mNuSq(1);
DeltamNuSq(end) = DeltamNuSqAll;
a.delete;
for i=1:6
    if i==6
        a = text(0.15,i+0.62,sprintf('\\Delta{\\itm}_\\nu^2 (eV^2)'),'FontSize',get(gca,'FontSize'),'Color',rgb('DimGray'));
    else
        if abs(DeltamNuSq(i))<1e-02
            a = text(0.2,i+1.05,sprintf('%+.0f\\times10^{-3}',1e3.*DeltamNuSq(i)),'FontSize',get(gca,'FontSize'),'Color',rgb('DimGray'));
        else
            
            a = text(0.2,i+1.05,sprintf('%+.0f\\times10^{-2}',1e2.*DeltamNuSq(i)),'FontSize',get(gca,'FontSize'),'Color',rgb('DimGray'));
        end
    end
end
%% save
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = sprintf('%sknm1_UniformFit_ReAna_mNuSqOverview.pdf',plotdir);
%print(f1,plotname,'-dpng','-r350');
export_fig(f1,plotname);