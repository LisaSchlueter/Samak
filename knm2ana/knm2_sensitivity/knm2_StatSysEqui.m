% fit neutrino mass at different ranges (-> qU Scan)
% once stat. only
% once with systematics
% search for equilibrium of statistical and systematic uncertainties
% calculate, plot in different script


qURange  = [95,20];
fitter = 'minuit';
Chi2Profile = 'OFF';

% create mini short cut file for plotting
savedir = [getenv('SamakPath'),'knm2ana/knm2_qUScan/results/'];
savenameCM = sprintf('%sknm2_qUScanTwin_Mini_%.0feV_to_%.0feV_%s_NP%.3f_%s_profilechi2%s.mat',...
    savedir,qURange(1),qURange(2),'chi2CMShape',1.112,fitter,Chi2Profile);
savenameStat = sprintf('%sknm2_qUScanTwin_Mini_%.0feV_to_%.0feV_%s_NP%.3f_%s_profilechi2%s.mat',...
    savedir,qURange(1),qURange(2),'chi2Stat',1,fitter,Chi2Profile);

if exist(savenameCM,'file')
    dSys = importdata(savenameCM);
     fprintf('loqe file %s \n',savenameCM)
else
    fprintf(2,'file not found %s \n',savenameCM)
end

if exist(savenameStat,'file')
    dStat = importdata(savenameStat);
     fprintf('loqe file %s \n',savenameStat)
else
    fprintf(2,'file not found %s \n',savenameStat)
end

%% prepare data
FitRange = dSys.D.RunData.qU(1:nqU)-18574;
ErrStat = dStat.errqU(1,:);
ErrSyst = dSys.errqU(1,:);

KeepIdx = dSys.errqU(1,:)>dStat.errqU(1,:);
FitRange = FitRange(KeepIdx);
ErrStat  = ErrStat(KeepIdx);
ErrSyst  = ErrSyst(KeepIdx);
Sigma_sys =  sqrt(ErrSyst.^2-ErrStat.^2);

%
figureHandle = figure('Units','normalized','Position',[0.1,0.1,0.5,0.4]);

nqU = numel(dSys.chi2qU);

hold on;
x = linspace(min(FitRange),max(FitRange),1e3);
Idx_AnaInterval = 11;
% [l,aref]  = boundedline(FitRange(Idx_AnaInterval).*ones(10,1),linspace(-1,2,10),0.8.*ones(10,1),'orientation','horiz');
% l.LineStyle = 'none'; aref.FaceColor = rgb('Silver'); aref.FaceAlpha =0.8;%obj.PlotColorLight;
pref = plot(FitRange(Idx_AnaInterval).*ones(10,1),linspace(-1,2,10),':','LineWidth',2.5,'Color',rgb('DimGray'));
hold on;
pTall = plot(x,smooth(interp1(FitRange,ErrSyst,x,'spline'),100),'-','LineWidth',4,'Color',rgb('DodgerBlue'));
pTstat = plot(x,smooth(interp1(FitRange,ErrStat,x,'spline'),100),'-.','LineWidth',3.5,'Color',rgb('Orange'));
pTsys = plot(x,smooth(interp1(FitRange,Sigma_sys,x,'spline'),100),'--','LineWidth',3.5,'Color',rgb('Crimson'));

%  plot(FitRange,ErrStat,'o','LineWidth',3.5,'Color',rgb('Orange'));
%  plot(FitRange,ErrSyst,'x','LineWidth',3.5,'Color',rgb('Crimson'));



xlabel('Lower fit boundary below 18574 (eV)');
ylabel(sprintf('1\\sigma sensitivity on {\\itm}_\\nu^2 (eV^2)'));
PrettyFigureFormat('FontSize',18);
xlim([-90,-20]);
ylim([0,0.62]);
leg = legend([pTall,pTstat,pTsys],'Stat. and syst.','Stat. only','Syst. only','Location','northwest');
leg.ItemTokenSize = [40,18];
PrettyLegendFormat(leg);
leg.FontSize = get(gca,'FontSize')+2;
t = text(FitRange(Idx_AnaInterval)-0.2,0.375,...
    sprintf('   Standard \nanalysis range'),...
    'Rotation',90,'FontSize',get(gca,'FontSize')+2,'Color',rgb('DimGray'));


%%
return
pltdir = [getenv('SamakPath'),'knm1ana/knm1_qUScan/plots/'];
pltname = sprintf('%sknm1_StatSysEqui_%s.pdf',pltdir,DataType);
export_fig(pltname);
fprintf('save plot to %s \n',pltname)





