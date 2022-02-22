% plot final knm1 systematics (PRL config) breakdown
% data +  sensitivity as ref
close all
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
    sprintf('Background {\\itqU} slope');...
   };

% PlotColor = {rgb('PowderBlue'),rgb('DodgerBlue'),rgb('GoldenRod'),rgb('ForestGreen'),...
%     rgb('DarkSlateGray'),rgb('DarkOrange'),rgb('FireBrick'),rgb('HotPink'),...
%     rgb('YellowGreen'),rgb('Navy')};

PlotColor = {rgb('DodgerBlue'),rgb('MistyRose'),rgb('YellowGreen'),rgb('PowderBlue'),...
    rgb('HotPink'),rgb('GoldenRod'),rgb('DarkOrange'),rgb('FireBrick'),...
    rgb('ForestGreen'),rgb('MediumBlue'),rgb('DarkSlateGray')};

%% load asymmetric stat. only with MINOS errors
savefileStatR = sprintf('%sknm1finalfitUniform_%s_chi2Stat_NP1_mNuE0NormBkg_40eV_Sibille0p5eV.mat',[getenv('SamakPath'),'knm1ana/knm1_unblinding/results/'],'Real');
dStatR = importdata(savefileStatR);
savefileStatT = sprintf('%sknm1finalfitUniform_%s_chi2Stat_NP1_mNuE0NormBkg_40eV_Sibille0p5eV.mat',[getenv('SamakPath'),'knm1ana/knm1_unblinding/results/'],'Twin');
dStatT = importdata(savefileStatT);

savefileCMR = sprintf('%sknm1finalfitUniform_%s_chi2CMShape_SysBudget22_NP1.064_mNuE0NormBkg_40eV_Sibille0p5eV.mat',[getenv('SamakPath'),'knm1ana/knm1_unblinding/results/'],'Real');
dCmR = importdata(savefileCMR);
savefileCMT = sprintf('%sknm1finalfitUniform_%s_chi2CMShape_SysBudget22_NP1.064_mNuE0NormBkg_40eV_Sibille0p5eV.mat',[getenv('SamakPath'),'knm1ana/knm1_unblinding/results/'],'Twin');
dCmT = importdata(savefileCMT);
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
SingleBarY_twin(end+1) = 0.5*(dStatT.FitResult.errPos(1)-dStatT.FitResult.errNeg(1));
SingleBarY_twin(end+1) = 0.5*(dCmT.FitResult.errPos(1)-dCmT.FitResult.errNeg(1));

%data
PlotVarTmp_data = struct2array(structfun(@(x)x(1),ddata.MultiLpar,'UniformOutput',0));      
SingleBarStat_data       = ddata.MultiLpar.Stat(1);


% add Non-Poisson and Stat only and total
SingleBarY_data(1)      = sqrt(ddata.MultiLpar.Stat(1)^2-ddata.NPcomponent(1)^2);
SingleBarY_data(2:end)  = sqrt(PlotVarTmp_data(2:end).^2-SingleBarStat_data^2);
SingleBarY_data(4) = sqrt(SingleBarY_data(4)^2+SingleBarY_data(8)^2); % merge TASR and HV
SingleBarY_data(8) = []; % remove HV stacking, because not approved on its own

SingleBarY_data(end+1) = 0.5*(dStatR.FitResult.errPos(1)-dStatR.FitResult.errNeg(1));
SingleBarY_data(end+1) = 0.5*(dCmR.FitResult.errPos(1)-dCmR.FitResult.errNeg(1));

y = SingleBarY_data;
%%

f55 = figure('Name','MultiBarPlot','Renderer','painters');
set(f55, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.8]);

leg_str = {'Background rate', SysEffectLeg{:},'Statistics','Total'};
LocalFontSize = 27;
[y,ia] = sort(y,'descend');
leg_str = leg_str(ia);

bsingle = cell(nSys+1,1);
btwin = cell(nSys+1,1);
t =cell(nSys+1,1);
SingleBarX = numel(y):-1:1;
for i=1:numel(y)
    bsingle{i}  = barh(SingleBarX(i),y(i));
    hold on;
    bsingle{i}.FaceColor = PlotColor{ia(i)}; bsingle{i}.LineStyle ='none';
    bsingle{i}.FaceAlpha=1;
    bsingle{i}.BarWidth = 0.9;
%     if i==1
%         bsingle{i}.FaceColor = rgb('LightGray');
%     end
    
    if y(i)<0.001
        tstr = sprintf('<10^{-3}') ;
    elseif y(i)<0.01
         tstr = sprintf('<10^{-2}') ;
    else  
         tstr = sprintf('%.0f\\cdot10^{-2}',y(i)*1e2) ;
    end%0.46
    t{i}= text(2.1,SingleBarX(i),tstr,...
        'HorizontalAlignment','right','FontSize',LocalFontSize-2,...
        'Interpreter','tex',...
        'FontName','Times New Roman');
    
%     % add sensitivity as reference line
%     btwin{i}  = barh(SingleBarX(i),SingleBarY_twin(i));
%     hold on;
%     if SingleBarY_twin(i)>y(i)
%         btwin{i}.EdgeColor = bsingle{i}.FaceColor;
%     else
%         btwin{i}.EdgeColor = 'k';
%     end
%     btwin{i}.FaceColor = 'none';
%     btwin{i}.LineStyle =':';
%     btwin{i}.BarWidth = 0.85;
%     btwin{i}.LineWidth = 2;

      if SingleBarY_twin(ia(i))<0.001
        tstr = sprintf('<10^{-3}') ;
    elseif SingleBarY_twin(ia(i))<0.01
         tstr = sprintf('<10^{-2}') ;
    else  
         tstr = sprintf('%.0f\\cdot10^{-2}',SingleBarY_twin(ia(i))*1e2) ;
      end%0.95
     t{i}= text(5,SingleBarX(i),tstr,...
        'HorizontalAlignment','right','FontSize',LocalFontSize-2,...
        'Interpreter','tex',...
        'FontName','Times New Roman',...
        'Color',rgb('DimGray'));
end

% axis options
%xlim([5e-03 0.49])
xlim([5e-03 2.2]);
ylim([min(SingleBarX)-0.9 max(SingleBarX)+0.9])
cl = 0.683;

    xlabel(sprintf('1\\sigma uncertainty on {\\itm}^2_\\nu at %.1f%% C.L. (eV^{ 2})',cl*100));
 %   xlabel(sprintf('1\\sigma sensitivity on {\\itm}^2_\\nu at %.0f%% C.L. (eV^{ 2})',cl*100));

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

savename = sprintf('%s/plots/KNM1_SysBreakdown_DataTwin.pdf',filedir);
export_fig(savename);