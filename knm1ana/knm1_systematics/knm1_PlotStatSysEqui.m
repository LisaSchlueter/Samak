% plotting script
% fit neutrino mass at different ranges (-> qU Scan)
% once stat. only
% once with systematics
% search for equilibrium of statistical and systematic uncertainties

%% load files
qURange  = [95,20];
savedir = [getenv('SamakPath'),'knm1ana/knm1_qUScan/results/'];
DataType = 'Twin';

% filenames
switch DataType
    case 'Twin'
        savenameSys = sprintf('%sknm1_qUScan_Mini_%s_%s_NP%.4g_%.0feV_to_%.0feV.mat',savedir,'Twin','chi2CMShape',1.064,qURange(1),qURange(2));
        savenameStat = sprintf('%sknm1_qUScan_Mini_%s_%s_NP%.4g_%.0feV_to_%.0feV.mat',savedir,'Twin','chi2Stat',1,qURange(1),qURange(2));
    case 'Real'
        savenameSys = sprintf('%sknm1_qUScan_Mini_%s_%s_NP%.4g_%.0feV_to_%.0feV.mat',savedir,'Real','chi2CMShape',1.064,qURange(1),qURange(2));
        savenameStat = sprintf('%sknm1_qUScan_Mini_%s_%s_NP%.4g_%.0feV_to_%.0feV.mat',savedir,'Real','chi2Stat',1,qURange(1),qURange(2));
end
% import
dSys = importdata(savenameSys);
dStat = importdata(savenameStat);

%%
figureHandle = figure('Units','normalized','Position',[0.1,0.1,0.5,0.4]);

nqU = numel(dSys.chi2qU);
FitRange = dSys.M.RunData.qU(1:nqU)-18574;
Sigma_sys =  sqrt(dSys.mNuSqErr.^2-dStat.mNuSqErr.^2);
hold on;
x = linspace(min(FitRange),max(FitRange),1e3);
Idx_AnaInterval = 13;
% [l,aref]  = boundedline(FitRange(Idx_AnaInterval).*ones(10,1),linspace(-1,2,10),0.8.*ones(10,1),'orientation','horiz');
% l.LineStyle = 'none'; aref.FaceColor = rgb('Silver'); aref.FaceAlpha =0.8;%obj.PlotColorLight;
pref = plot(FitRange(Idx_AnaInterval).*ones(10,1),linspace(-1,2,10),':','LineWidth',2.5,'Color',rgb('DimGray'));
hold on;
pTall = plot(x,smooth(interp1(FitRange,dSys.mNuSqErr,x,'spline'),100),'-','LineWidth',4,'Color',rgb('DodgerBlue'));
pTstat = plot(x,smooth(interp1(FitRange,dStat.mNuSqErr,x,'spline'),100),'-.','LineWidth',3.5,'Color',rgb('Orange'));
pTsys = plot(x,smooth(interp1(FitRange,Sigma_sys,x,'spline'),100),'--','LineWidth',3.5,'Color',rgb('Crimson'));

xlabel('Lower fit boundary below 18574 (eV)');
ylabel(sprintf('1\\sigma sensitivity on {\\itm}_\\nu^2 (eV^2)'));
PrettyFigureFormat('FontSize',18);
xlim([-93,-20]);
ylim([0,1.8]);
leg = legend([pTall,pTstat,pTsys],'Stat. and syst.','Stat. only','Syst. only','Location','northwest');
leg.ItemTokenSize = [40,18];
PrettyLegendFormat(leg);
leg.FontSize = get(gca,'FontSize')+2;
t = text(FitRange(Idx_AnaInterval)-0.2,1,...
    sprintf('   Standard \nanalysis range'),...
    'Rotation',90,'FontSize',get(gca,'FontSize')+2,'Color',rgb('DimGray'));


%%
pltdir = [getenv('SamakPath'),'knm1ana/knm1_qUScan/plots/'];
pltname = sprintf('%sknm1_StatSysEqui_%s.pdf',pltdir,DataType);
export_fig(pltname);
fprintf('save plot to %s \n',pltname)







