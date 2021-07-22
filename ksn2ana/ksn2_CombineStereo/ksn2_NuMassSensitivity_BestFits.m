% position of profile best fits
% based on ksn2_NuMassSensitivity.m
% impact of 3+1 model extionsion on neutrino mass sensitivity

chi2 = 'chi2CMShape';
DataType = 'Real';%
SavePlt = 'ON';

CombiStereo = 'OFF';

if strcmp(DataType,'Real')
mNuSq = -1:0.01:2.5;
else
mNuSq = -1:0.01:2.5;
end

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];
if strcmp(CombiStereo,'OFF')
    CombiStr = '';
else
    CombiStr = 'CombiStereo';
end
savefile = sprintf('%sksn2_NuMassSensitivity%s_%s_%s_mNuSqMin%.2f_mNuSqMax%.2f_mNuSteps%.0f.mat',...
    savedir,CombiStr,DataType,chi2,min(mNuSq),max(mNuSq),numel(mNuSq));%CombiStereo

if exist(savefile,'file') 
    load(savefile);
     fprintf('load file  %s \n',savefile)
else
     fprintf('file does not exist: %s \n',savefile)
    return
end


GetFigure;


mNuSq_leg = [-1:0.25:0,0.3,0.5,1,2];
Mycolors= [autumn(sum(mNuSq_leg<=0));winter(sum(mNuSq_leg>0))];

pHandles = cell(numel(mNuSq_leg),1);
legStr   = cell(numel(mNuSq_leg),1);

for i=1:numel(mNuSq)
    ptmp =  plot(sin2T4_bf(i),mNu4Sq_bf(i),'.','MarkerSize',10,'Color',rgb('Silver'));
    hold on;
    if min(abs(mNuSq(i)-mNuSq_leg))<1e-010%ismember(mNuSq(i),mNuSq_leg)
        Idx = find(abs(mNuSq(i)-mNuSq_leg)==min(abs(mNuSq(i)-mNuSq_leg)));%find(mNuSq(i)==mNuSq_leg);
        
        ptmp.MarkerFaceColor = Mycolors(Idx,:);
        ptmp.Color = Mycolors(Idx,:);
        ptmp.Marker = 'o';
        ptmp.MarkerSize = 8;
        
        pHandles{Idx} = ptmp;
        legStr{Idx} = sprintf('%.1f eV^2',mNuSq(i));
    end
end
%%
plot(sin2T4_contour_bf,mNu4Sq_contour_bf,'kx','LineWidth',2);
plot(sin2T4_contour,mNu4Sq_contour,'k-','LineWidth',2);
PrettyFigureFormat;
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('m_{4}^2 (eV^2)'));
leg = legend([pHandles{:}],legStr,'Location','southwest');
leg.Title.String = sprintf('Local best fit for {\\itm}_\\nu^2 (eV^2)');
leg.Title.FontSize = get(gca,'FontSize'); leg.Title.FontWeight = 'normal';
PrettyLegendFormat(leg);
xlim([8e-04 0.5]);
ylim([0.1 40^2]);
pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/plots/'];
MakeDir(pltdir)
pltfile = sprintf('%sksn2_NuMassSensitivity%s_%s_%s_mNuSqMin%.2f_mNuSqMax%.2f_mNuSteps%.0f.png',...
    pltdir,CombiStr,DataType,chi2,min(mNuSq),max(mNuSq),numel(mNuSq));%CombiStereo
print(pltfile,'-dpng','-r350');
