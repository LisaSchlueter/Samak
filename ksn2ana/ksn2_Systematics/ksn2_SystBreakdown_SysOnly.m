% a la Mainz: (fig. 4) https://arxiv.org/pdf/1210.4194.pdf
% syst. only  with raster scan
% 2 Plot Options: (x,y) = (sin2t4,m4Sq) or vice versa
PltAxism4Sq = 'y';
SavePlt = 'ON';

RasterScan = 'ON';
DataType = 'Real';
range = 40;

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
plotdir = strrep(savedir,'results','plots');
savename = sprintf('%sksn2_SystBreakdown_StatOverSyst_%s_%.0feV_RasterScan%s.mat',savedir,DataType,range,RasterScan);

if exist(savename,'file')
    d = importdata(savename);
    fprintf('load file %s \n',savename);
else
    fprintf('file not found: %s \n',savename);
    return
end

pHandle = cell(numel(d.SysEffectsAll),1);
LineStyle = {'-','-.',':','--','-','-.',':','--','-','-.',':','--','-.','-','--','-.',':','--','-.',':','--'};
Colors = colormap('jet');
fSyst = figure('Units','normalized','Position',[0.1,0.1,0.7,0.5]);

if strcmp(PltAxism4Sq,'y')
    %% plot 2. only systematics contribution in contour plot style
    hold on;
    for i=1:d.nSys
        if strcmp(d.SysEffectsAll{i},'all')
            PltColor = rgb('Black');
            LineWidth = 3;
        else
            PltColor =Colors(floor((i)*256/(d.nSys)),:);
            LineWidth = 2;
        end
        pHandle{i} = plot(d.sin2t4_Sys(i,:),d.mNu4Sq(i,:),'LineWidth',LineWidth,'Color',PltColor,'LineStyle',LineStyle{i});
        
    end
    pStat = plot(d.sin2t4_Stat(1,:),d.mNu4Sq(1,:),'LineWidth',3,'Color',rgb('Silver'));
    xlabel(sprintf('\\sigma_{syst.}'))
    ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
    xlim([2e-04 0.5]);
    ylim([1.3 38^2])
    
    
elseif strcmp(PltAxism4Sq,'x')
    %%
    hold on;
    
    for i=1:d.nSys
        if strcmp(d.SysEffectsAll{i},'all')
            PltColor = rgb('Black');
            LineWidth = 3;
        else
            PltColor =Colors(floor((i)*256/(d.nSys)),:);
            LineWidth = 2;
        end
        pHandle{i} = plot(d.mNu4Sq(i,:),d.sin2t4_Sys(i,:),'LineWidth',LineWidth,'Color',PltColor,'LineStyle',LineStyle{i});
    end
    pStat = plot(d.mNu4Sq(1,:),d.sin2t4_Stat(1,:),'LineWidth',3,'Color',rgb('Silver'));
    
    xlabel(sprintf('{\\itm}_4^2 (eV^2)'));
    ylabel(sprintf('\\sigma_{syst.}'));%\\Delta |{\\itU}_{e4}|^2'));
    ylim([2e-04 0.5]);
    xlim([1.3 38^2]);
end

PrettyFigureFormat('FontSize',22);
set(gca,'XScale','log');
set(gca,'YScale','log');
leg = legend([pHandle{:},pStat],{d.SysEffectLabel{:},'Stat. only'});
PrettyLegendFormat(leg);
leg.Title.String = 'Stystematic effect';
leg.Title.FontWeight ='normal';
leg.Location = 'eastoutside';

if strcmp(SavePlt,'ON')
    plotname = sprintf('%sksn2_SysBreakdown_SystOnly_%s_RasterScan%s_%s.png',plotdir,DataType,RasterScan,PltAxism4Sq);
    print(fSyst,plotname,'-dpng','-r400');
    fprintf('save plot to %s \n',plotname);
    export_fig(fSyst,strrep(plotname,'png','pdf'));
end