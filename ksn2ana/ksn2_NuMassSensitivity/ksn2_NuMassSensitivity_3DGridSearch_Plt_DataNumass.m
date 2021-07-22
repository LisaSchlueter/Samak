% ksn2 calculate chi2-prifle for nu-mass
% perform grid searches for differed nu-masses
% plot data + knm2 result
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

%% plot

GetFigure;

% load knm2 results
% boundedline: 1 sigma band
dN = importdata(sprintf('%stritium-data/fit/Knm2/Chi2Profile/Uniform/Chi2Profile_Real_UniformScan_mNu_Knm2_UniformFPD_chi2CMShape_SysBudget40_NP1.112_FitParE0BkgNorm_nFit50_KNM2_0p1eV_min-2.6_max1.mat',getenv('SamakPath')));
[l2,a2] = boundedline(dN.ScanResults.BestFit.par.*ones(10,1),linspace(-5,1e2,10),[dN.ScanResults.BestFit.errNeg.*ones(10,1),dN.ScanResults.BestFit.errPos.*ones(10,1)],...
    'orientation','horiz');
l2.delete;
a2.FaceAlpha = 0.3; a2.FaceColor = rgb('DarkGray');
hold on;
[l,a] = boundedline(mNuSq_bf.*ones(10,1),linspace(-5,1e2,10),[mNuSqErrDown.*ones(10,1),mNuSqErrUp.*ones(10,1)],...
    'orientation','horiz');
l.delete; a.FaceAlpha = 0.3; a.FaceColor = rgb('SkyBlue');

% ref lines
plot(mNuSq_inter,zeros(numel(mNuSq_inter),1),'k-','LineWidth',1);
plot(mNuSq_inter,ones(numel(mNuSq_inter),1),'k-','LineWidth',1);

% plot smooth chi^2 function
chi2ref2 = 0;%min(min(dN.ScanResults.chi2min); 
pknm2 = plot(dN.ScanResults.ParScan(:,1),dN.ScanResults.chi2min(:,1)-chi2ref2,':','Color',rgb('DimGray'),'LineWidth',3);
plot(dN.ScanResults.ParScan(2:end,2),dN.ScanResults.chi2min(2:end,2)-chi2ref2,':','Color',pknm2.Color,'LineWidth',pknm2.LineWidth);
yR = smooth(chi2_inter,100);
chi2ref1 = 0;%min(yR);
pksn2 = plot(mNuSq_inter,yR-chi2ref1,'-','LineWidth',3,'Color',rgb('SkyBlue'));


%%
xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
ylabel(sprintf('\\chi^2'));
PrettyFigureFormat('FontSize',22);
dofKSN2 = 28-4-2;
dofKNM2 = 28-4;

leg = legend([pksn2,pknm2],...
    sprintf('3\\nu+1: {\\itm}_\\nu^2 = %.2f ^{+ %.2f}_{- %.2g} eV^2 with %.0f dof',...
    mNuSq_bf,mNuSqErrUp,mNuSqErrDown,dofKSN2),...
    sprintf('3\\nu   : {\\itm}_\\nu^2 = %.2f ^{+ %.2f}_{- %.2f} eV^2 with %.0f dof',...
    dN.ScanResults.BestFit.par,...
    dN.ScanResults.BestFit.errPos,dN.ScanResults.BestFit.errNeg,dofKNM2),...
    'Location','northeast');
ylim([25 32]);
xlim([-1 2.5]);
PrettyLegendFormat(leg,'alpha',0.8);
leg.FontSize = get(gca,'FontSize');
t = text(-0.9,1.15,sprintf('68.27%% C.L.'),'FontSize',get(gca,'FontSize'));
% 
return
pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);
pltname = strrep([strrep(extractBefore(savefileCombiR,'_Combi'),'results','plots'),'.png'],'Real','Ksn2Knm2');
%pltname = strrep(strrep(strrep(savefileCombiR,'results','plots'),'.mat','.png'),'Real','Ksn2Knm2');
print(gcf,pltname,'-dpng','-r350');
fprintf('save plot to %s \n',pltname);
