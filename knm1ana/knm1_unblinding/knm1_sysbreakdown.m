% plot final knm1 systematics (PRL config) breakdown
% either data or sensitivity
DataType = 'Real';
filedir = [getenv('SamakPath'),'tritium-data/sensitivity/Knm1/'];
twinFile =  sprintf('%sSensitivitySys_Asimov_KNM1_40eV_TCFSDTASRRF_ELRF_BFRF_RXStackFPDeffBkg_chi2CMShape_budget22_MatlabFit_Twin.mat',filedir);
dataFile =  sprintf('%sSensitivitySys_Asimov_KNM1_40eV_TCFSDTASRRF_ELRF_BFRF_RXStackFPDeffBkg_chi2CMShape_budget22_MatlabFit_Real.mat',filedir);

dtwin = importdata(twinFile);
ddata = importdata(dataFile);

%% legend
SysEffectLeg      = {'Theoretical corrections';...
    'Final-state distribution';...
    'Scan fluctuations';...%: Tritium activity';...
    sprintf('Energy-loss function');...
    'Magnetic fields';...
    sprintf('Number of scatterings \\rho{\\itd}\\sigma');...%  
    'Detector efficiency';...
    sprintf('Background {\\itqU} slope')};

PlotColor = {rgb('White'),rgb('DodgerBlue'),rgb('GoldenRod'),rgb('PowderBlue'),...
    rgb('CadetBlue'),rgb('DarkOrange'),rgb('FireBrick'),rgb('DarkSlateGray'),...
    rgb('YellowGreen'),rgb('Navy')};
%% assign variables
nSys = numel(dtwin.SysEffectsAll);
SingleBarY_twin = zeros(1,nSys+1);
SingleBarY_data = zeros(1,nSys+1);
%twin
SingleBarStat_twin = dtwin.MultiLpar.Stat(1);
PlotVarTmp_twin = struct2array(structfun(@(x)x(1),dtwin.MultiLpar,'UniformOutput',0));      
SingleBarY_twin(1)      = sqrt(dtwin.MultiLpar.Stat(1)^2-dtwin.NPcomponent(1)^2);
SingleBarY_twin(2:end)  = sqrt(PlotVarTmp_twin(2:end).^2-SingleBarStat_twin^2);
SingleBarY_twin(4) = sqrt(SingleBarY_twin(4)^2+SingleBarY_twin(8)^2);
SingleBarY_twin(8) = []; % remove HV stacking, because not approved

%data
PlotVarTmp_data = struct2array(structfun(@(x)x(1),ddata.MultiLpar,'UniformOutput',0));      
SingleBarStat_data       = ddata.MultiLpar.Stat(1);
SingleBarY_data(1)      = sqrt(ddata.MultiLpar.Stat(1)^2-ddata.NPcomponent(1)^2);
SingleBarY_data(2:end)  = sqrt(PlotVarTmp_data(2:end).^2-SingleBarStat_data^2);
SingleBarY_data(4) = sqrt(SingleBarY_data(4)^2+SingleBarY_data(8)^2);
SingleBarY_data(8) = []; % remove HV stacking, because not approved

if strcmp(DataType,'Twin')
    y = SingleBarY_twin;
else
    y = SingleBarY_data;
end
%%

f55 = figure('Name','MultiBarPlot','Renderer','painters');
set(f55, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.8]);

leg_str = {'Background rate', SysEffectLeg{:}};
LocalFontSize = 27;
[y,ia] = sort(y,'descend');
leg_str = leg_str(ia);

bsingle = cell(nSys+1,1);
t =cell(nSys+1,1);
SingleBarX = numel(y):-1:1;
for i=1:numel(SingleBarY_twin)
    bsingle{i}  = barh(SingleBarX(i),y(i));
    hold on;
    bsingle{i}.FaceColor = PlotColor{ia(i)}; bsingle{i}.LineStyle ='none';
    bsingle{i}.FaceAlpha=1;
    bsingle{i}.BarWidth = 0.9;
    if i==1
        bsingle{i}.FaceColor = rgb('LightGray');
    end
    
    if y(i)<0.001
        tstr = sprintf('<10^{-3}') ;
    elseif y(i)<0.01
         tstr = sprintf('<10^{-2}') ;
    else  
         tstr = sprintf('%.0f\\cdot10^{-2}',y(i)*1e2) ;
    end
    t{i}= text(0.46,SingleBarX(i),tstr,...
        'HorizontalAlignment','right','FontSize',LocalFontSize-2,...
        'Interpreter','tex',...
        'FontName','Times New Roman');
end

% axis options
xlim([5e-03 0.49])
ylim([min(SingleBarX)-0.9 max(SingleBarX)+0.9])
cl = 0.683;

if strcmp(DataType,'Real')
    xlabel(sprintf('1 \\sigma uncertainty on {\\itm}^2_\\nu at %.0f%% C.L. (eV^{ 2})',cl*100));
else
    xlabel(sprintf('1 \\sigma sensitivity on {\\itm}^2_\\nu at %.0f%% C.L. (eV^{ 2})',cl*100));
end

PRLFormat;
set(gca,'FontSize',LocalFontSize);
% no y-ticks
set(gca,'YMinorTick','off');
set(gca,'TickDir','out');

h = gca;
h.YRuler.TickLength = [0,0];
set(gca,'XScale','log');
xticks([0.0,0.0001,0.001,0.01,0.1,1])
yticks(flip(SingleBarX))
yticklabels(flip(leg_str));


%remove top and right ticks
a = gca;
set(a,'box','off','color','none')% set box property to off and remove background color
b = axes('Position',a.Position,...
    'box','on','xtick',[],'ytick',[],'LineWidth',1.5);% create new, empty axes with box but without ticks
axes(a)% set original axes as active
% linkaxes([a b]) % link axes in case of zooming

savename = sprintf('%s/plots/KNM1_SysBreakdown_%s.pdf',filedir,DataType);
export_fig(savename);