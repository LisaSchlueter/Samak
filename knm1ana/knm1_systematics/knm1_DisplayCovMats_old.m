% make nice plots for 1 covariance matrix
SysEffect = 'FSD';
CM = 'CMShape';
RunList = 'KNM1';
if ~exist('M','var')
M = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','DataType','Real');
end
M.ComputeCM;%('SysEffects',struct(SysEffect,'ON'),'BkgCM','OFF')
qUWindowIndexMin = 14;

switch CM
    case 'CM'
        CovMatDisp = M.FitCM_Obj.CovMat;
        str = 'covariance matrix';
    case 'CMFrac'
        CovMatDisp = M.FitCM_Obj.CovMatFrac;
        str = 'fractional covariance matrix';
    case 'CMShape'
        CovMatDisp = M.Convert2ShapeOnlyCM('CMFrac',M.FitCM_Obj.CovMatFrac);
        str = 'fractional shape only covariance matrix';
end
%%
imagesc(CovMatDisp(qUWindowIndexMin:end,qUWindowIndexMin:end));
PrettyFigureFormat('FontSize',24)
c = colorbar;
c.Label.String = 'Covariance';
c.Label.FontSize = get(gca,'FontSize')+4;
c.LineWidth=2;
%title(str);
colormap(parula);
pbaspect([1 1 1])
set(gca,'xtick',[2 40-qUWindowIndexMin]),set(gca,'ytick',[])
qUmin = sprintf('-%.0f eV',abs(M.FitCM_Obj.StudyObject.qU(qUWindowIndexMin)-M.FitCM_Obj.StudyObject.Q_i)+1);
qUmax = sprintf('%.0f eV',abs(M.FitCM_Obj.StudyObject.qU(end)-M.FitCM_Obj.StudyObject.Q_i)+1);
set(gca,'xticklabel',{qUmin,qUmax}); set(gca,'yticklabel',[]);
xlabel(sprintf('Retarding energy below E_{0,eff}'));
set(gca,'XMinorTick','off');
set(gca,'TickLength',[0 0]);
savedir = [getenv('SamakPath'),'knm1ana/knm1_systematics/plots/'];
%print([savedir,CM,SysEffect,'.png'],'-dpng','-r450');
export_fig([savedir,CM,SysEffect,'.pdf']);
%% correlation
corplot(CovMatDisp(qUWindowIndexMin:end,qUWindowIndexMin:end));
c = colorbar;
title('correlation matrix');
PrettyFigureFormat
pbaspect([1 1 1])
set(gca,'xtick',[1 40-qUWindowIndexMin+1]),set(gca,'ytick',[])
qUmin = sprintf('qU_{min} = E_0-%.0fV',abs(M.FitCM_Obj.StudyObject.qU(qUWindowIndexMin)-M.FitCM_Obj.StudyObject.Q_i));
qUmax = sprintf('qU_{max} = E_0+%.0fV',abs(M.FitCM_Obj.StudyObject.qU(end)-M.FitCM_Obj.StudyObject.Q_i));
set(gca,'xticklabel',{qUmin,qUmax}); set(gca,'yticklabel',[]);
set(gca,'FontSize',16)
savedir = [getenv('SamakPath'),'knm1ana/knm1_systematics/plots/'];
print([savedir,CM,SysEffect,'CorrMat.png'],'-dpng','-r450');