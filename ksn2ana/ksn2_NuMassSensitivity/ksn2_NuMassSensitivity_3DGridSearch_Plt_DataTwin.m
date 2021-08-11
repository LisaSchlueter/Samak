% ksn2 calculate chi2-prifle for nu-mass
% perform grid searches for differed nu-masses
% plot data + twin in 1
%% settings that might change

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];
savefileCombiR = sprintf('%sksn2_%s_NuMassSensitivityGridSearch_CombiSTEREO-%s.mat',...
    savedir,'Real','ON');

if exist(savefileCombiR,'file') 
    dR = importdata(savefileCombiR);
else
    fprintf('file not available \n')
end
%% some simple calculations for data
% interpolate
mNuSq_inter = linspace(min(dR.FixmNuSq_all),max(dR.FixmNuSq_all),1e3);
chi2_inter  = interp1(dR.FixmNuSq_all,dR.chi2min,mNuSq_inter,'spline');
% best fit
chi2bf      = min(chi2_inter);
Idx_bf      = find(chi2_inter==chi2bf);
mNuSq_bf    = mNuSq_inter(Idx_bf);
% uncertainties
mNuSqDown = interp1(chi2_inter(mNuSq_inter<mNuSq_bf),mNuSq_inter(mNuSq_inter<mNuSq_bf),chi2bf+1,'spline');
mNuSqUp   = interp1(chi2_inter(mNuSq_inter>mNuSq_bf),mNuSq_inter(mNuSq_inter>mNuSq_bf),chi2bf+1,'spline');
mNuSqErrDown = mNuSq_bf-mNuSqDown;
mNuSqErrUp = mNuSqUp-mNuSq_bf;
mNuSqErr   = 0.5.*(mNuSqErrDown+mNuSqErrUp);
% STEREO combi
chi2_inter_Combi  = interp1(dR.FixmNuSq_all,dR.chi2min_Combi,mNuSq_inter,'spline');
% best fit
chi2bf_Combi      = min(chi2_inter_Combi);
Idx_bf_Combi      = find(chi2_inter_Combi==chi2bf_Combi);
mNuSq_bf_Combi    = mNuSq_inter(Idx_bf_Combi);
% uncertainties
if mNuSq_bf_Combi==min(mNuSq_inter)
    mNuSqDown_Combi = NaN;
else
    mNuSqDown_Combi = interp1(chi2_inter_Combi(mNuSq_inter<mNuSq_bf_Combi),...
        mNuSq_inter(mNuSq_inter<mNuSq_bf_Combi),chi2bf_Combi+1,'spline');
end
mNuSqUp_Combi   = interp1(chi2_inter_Combi(mNuSq_inter>mNuSq_bf_Combi),mNuSq_inter(mNuSq_inter>mNuSq_bf_Combi),chi2bf_Combi+1,'spline');
mNuSqErrDown_Combi = mNuSq_bf_Combi-mNuSqDown_Combi;
mNuSqErrUp_Combi   = mNuSqUp_Combi-mNuSq_bf_Combi;
mNuSqErr_Combi     = 0.5.*(mNuSqErrDown_Combi+mNuSqErrUp_Combi);

%% load twin
savefileCombiT = sprintf('%sksn2_%s_NuMassSensitivityGridSearch_CombiSTEREO-%s.mat',...
    savedir,'Twin','OFF');

if exist(savefileCombiT,'file') 
    dT = importdata(savefileCombiT);
else
    fprintf('file not available \n')
end
%% some simple calculations for twin
% interpolate
mNuSq_interT = linspace(min(dT.FixmNuSq_all),max(dT.FixmNuSq_all),1e3);
chi2_interT  = interp1(dT.FixmNuSq_all,dT.chi2min,mNuSq_interT,'spline');
% best fit
chi2bfT      = min(chi2_interT);
Idx_bf      = find(chi2_interT==chi2bfT);
mNuSq_bfT    = mNuSq_interT(Idx_bf);
% uncertainties
mNuSqDownT = interp1(chi2_interT(mNuSq_interT<mNuSq_bfT),mNuSq_interT(mNuSq_interT<mNuSq_bfT),chi2bfT+1,'spline');
mNuSqUpT   = interp1(chi2_interT(mNuSq_interT>mNuSq_bfT),mNuSq_interT(mNuSq_interT>mNuSq_bfT),chi2bfT+1,'spline');
mNuSqErrDownT = mNuSq_bfT-mNuSqDownT;
mNuSqErrUpT = mNuSqUpT-mNuSq_bfT;
mNuSqErr   = 0.5.*(mNuSqErrDownT+mNuSqErrUpT);
%% plot
CombiSTEREO = 'OFF';
GetFigure;
% load knm2 results
PltKnm2 = 'ON';
if strcmp(PltKnm2,'ON')
    dN = importdata(sprintf('%stritium-data/fit/Knm2/Chi2Profile/Uniform/Chi2Profile_Real_UniformScan_mNu_Knm2_UniformFPD_chi2CMShape_SysBudget40_NP1.112_FitParE0BkgNorm_nFit50_KNM2_0p1eV_min-2.6_max1.mat',getenv('SamakPath')));
    p2 = plot(dN.ScanResults.ParScan(:,1),dN.ScanResults.chi2min(:,1)-min(min(dN.ScanResults.chi2min)),':','Color',rgb('Black'),'LineWidth',3);
   hold on;
   plot(dN.ScanResults.ParScan(2:end,2),dN.ScanResults.chi2min(2:end,2)-min(min(dN.ScanResults.chi2min)),':','Color',p2.Color,'LineWidth',p2.LineWidth);
[l2,a2] = boundedline(dN.ScanResults.BestFit.par.*ones(10,1),linspace(-5,1e2,10),[dN.ScanResults.BestFit.errNeg.*ones(10,1),dN.ScanResults.BestFit.errPos.*ones(10,1)],...
    'orientation','horiz');
l2.delete;
a2.FaceAlpha = 0.3; a2.FaceColor = rgb('DarkGray');
end

[l,a] = boundedline(mNuSq_bf.*ones(10,1),linspace(-5,1e2,10),[mNuSqErrDown.*ones(10,1),mNuSqErrUp.*ones(10,1)],...
    'orientation','horiz');
l.delete; a.FaceAlpha = 0.3; a.FaceColor = rgb('SkyBlue');
hold on;
if strcmp(CombiSTEREO,'ON')
    [lS,aS] = boundedline(mNuSq_bf_Combi.*ones(10,1),linspace(-5,1e2,10),[mNuSqErrUp_Combi.*ones(10,1),mNuSqErrUp_Combi.*ones(10,1)],...
        'orientation','horiz');
    lS.delete; aS.FaceAlpha = 0.3; aS.FaceColor = rgb('BlueViolet');
end
[lT,aT] = boundedline(mNuSq_bfT.*ones(10,1),linspace(-5,1e2,10),[99.*ones(10,1),mNuSqErrUpT.*ones(10,1)],...
    'orientation','horiz');
lT.delete; aT.FaceAlpha = 0.3; aT.FaceColor = rgb('Orange');

plot(mNuSq_inter,zeros(numel(mNuSq_inter),1),'k--','LineWidth',1);
hold on;
plot(mNuSq_inter,ones(numel(mNuSq_inter),1),'k--','LineWidth',1);
yR = smooth(chi2_inter,100);
p1interR = plot(mNuSq_inter,yR-min(yR),'-','LineWidth',3,'Color',rgb('SkyBlue'));

if strcmp(CombiSTEREO,'ON')  
    yS = smooth(chi2_inter_Combi,100);
    p2Stereo = plot(mNuSq_inter,yS-min(yS),':','LineWidth',3,'Color',rgb('BlueViolet'));
end
yT = smooth(chi2_interT,100);
p1interT = plot(mNuSq_interT,yT-min(yT),'-.','LineWidth',3,'Color',rgb('Orange'));


%%
dof = 28-4-2;% nqU-freeFitPar
xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
ylabel(sprintf('\\Delta\\chi^2'));
%PrettyFigureFormat('FontSize',22);
PRLFormat;
set(gca,'FontSize',30);
if strcmp(CombiSTEREO,'OFF') && strcmp(PltKnm2,'OFF')
leg = legend([p1interR,p1interT],...
    sprintf('Data: {\\itm}_\\nu^2 = %.1f^{+%.2f}_{-%.2g} eV^2',mNuSq_bf,mNuSqErrUp,mNuSqErrDown),...
    sprintf('MC  : {\\itm}_\\nu^2 = %.1f^{+%.2f}_{-\\infty} eV^2',0,mNuSqErrUpT),...
    'Location','northwest');
ylim([-0.3 5]);
xlim([-1 2.5]);
elseif strcmp(CombiSTEREO,'OFF') && strcmp(PltKnm2,'ON')
% leg = legend([p1interR,p1interT,p2],...
%     sprintf('3\\nu+1 data: {\\itm}_\\nu^2 = %.1f^{+%.2f}_{-%.2g} eV^2',mNuSq_bf,mNuSqErrUp,mNuSqErrDown),...
%     sprintf('3\\nu+1 twin: {\\itm}_\\nu^2 = %.1f^{+%.2f}_{-\\infty} eV^2',0,mNuSqErrUpT),...
%     sprintf('3\\nu    data: {\\itm}_\\nu^2 = %.1f^{+%.2f}_{-%.2g} eV^2',dN.ScanResults.BestFit.par,dN.ScanResults.BestFit.errPos,dN.ScanResults.BestFit.errNeg),...
%     'Location','northwest');
% mini version
leg = legend([p1interR,p1interT,p2],...
    sprintf('3\\nu+1 data'),...
    sprintf('3\\nu+1 sensitivity'),...
    sprintf('3\\nu     data'),...
    'Location','northeast');
ylim([-0.3 5]);
xlim([-1 2.5]);
    
else
    leg = legend([p1interR,p2Stereo],...
        sprintf('{\\itm}_\\nu^2 = %.1f^{+%.2f}_{-%.2g} eV^2 (KATRIN)',mNuSq_bf,mNuSqErrUp,mNuSqErrDown),...
        sprintf('{\\itm}_\\nu^2 = %.1f^{+%.2f}_{-%.2g} eV^2 (KATRIN + STEREO)',mNuSq_bf_Combi,mNuSqErrUp_Combi,mNuSqErrDown_Combi),...
     'Location','north');
 ylim([-0.5 7]);
 xlim([min(max(dR.FixmNuSq_all)),max(dR.FixmNuSq_all)]);
end
PrettyLegendFormat(leg);

pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);
pltname = sprintf('%sksn2_NuMassSensitivity.png',pltdir);%strrep(strrep(strrep(savefileCombiR,'results','plots'),'.mat','.png'),'Real','');
print(gcf,pltname,'-dpng','-r350');
fprintf('save plot to %s \n',pltname);

export_fig(sprintf('%sksn2_NuMassSensitivity.pdf',pltdir));