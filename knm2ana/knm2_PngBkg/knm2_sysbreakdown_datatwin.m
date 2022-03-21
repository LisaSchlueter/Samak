% plot final knm1 systematics (PRL config) breakdown
% data +  sensitivity as ref
close all
filedir = [getenv('SamakPath'),'tritium-data/sensitivity/Knm2/'];
twinFile =  sprintf('%sSensitivitySys41_KNM2_Prompt_E018573.70eV_FSDsigma0.122eV_BkgPtSlope3.0muCpsPerS_41eV_TCoff_OTHERFSDTASRRF_ELRF_BFRF_RXStackLongPlasmaFPDeffNPBkgPTBkg_Par1_chi2CMShape_MinuitMinosFit_Twin_ranges11_SysOnly_RingFull.mat',filedir);
dataFile =  sprintf('%sSensitivitySys41_KNM2_Prompt_41eV_TCoff_OTHERFSDTASRRF_ELRF_BFRF_RXStackLongPlasmaFPDeffNPBkgPTBkg_Par1_chi2CMShape_MinuitMinosFit_Real_ranges11_SysOnly_RingFull.mat',filedir);

dtwin = importdata(twinFile);
ddata = importdata(dataFile);

%% legend
SysEffectLeg      = ...
    {'Theoretical corrections';...
    'Final-state distribution';...
    'Scan fluctuations';...%: Tritium activity';...
    sprintf('Energy-loss function');...
    'Magnetic fields';...
    sprintf('Number of scatterings \\rho{\\itd}\\sigma');...%  
    'Source potential';...
    'Detector efficiency';...
    'Non-Poisson background over-dispersion';...
     sprintf('Scan-step-time-dependent background');...
    sprintf('Retarding-potential dependent background');...
   };

% PlotColor = {rgb('PowderBlue'),rgb('DodgerBlue'),rgb('GoldenRod'),rgb('ForestGreen'),...
%     rgb('DarkSlateGray'),rgb('DarkOrange'),rgb('FireBrick'),rgb('HotPink'),...
%     rgb('YellowGreen'),rgb('Navy')};

PlotColor = {rgb('MistyRose'),...%theo
    rgb('YellowGreen'),...% FSD
    rgb('PowderBlue'),...%tasr
    rgb('HotPink'),...%eloss
    rgb('GoldenRod'),...%magnetic fields
    rgb('DarkOrange'),...% num scatterings
    rgb('Salmon'),... % long plasma
    rgb('FireBrick'),... % FPD eff
    rgb('DodgerBlue'),...Non Poiss
    rgb('Indigo'),... % bkg pt slope
    rgb('ForestGreen'),...% qU slope
    rgb('MediumBlue'),...%stat
    rgb('DarkSlateGray')}; % tot

 %% load asymmetric stat. only with MINOS errors
 savedir2 = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
% savefileStatR = sprintf('%sknm2ubfinal_Fit_Bpng-3.0mucpsPers_Real_40eV_mNuE0BkgNormqU_chi2Stat_Ring_KNM2.mat',savedir2);
% dStatR = importdata(savefileStatR);
% savefileStatT = sprintf('%sknm2ubfinal_Fit_Bpng-3.0mucpsPers_Twin_40eV_mNuE0BkgNorm_chi2Stat_Ring_KNM2_TwinBpng-3.0mucpsPers.mat',savedir2);
% dStatT = importdata(savefileStatT);
% 
savefileCMR = sprintf('%sknm2ubfinal_Fit_Bpng-3.0mucpsPers_Real_40eV_mNuE0BkgNormqU_chi2CMShape_Ring_KNM2_SysBudget41.mat',savedir2);
dCmR = importdata(savefileCMR);
% savefileCMT = sprintf('%sknm1finalfitUniform_%s_chi2CMShape_SysBudget22_NP1.064_mNuE0NormBkg_40eV_Sibille0p5eV.mat',savedir2);
% dCmT = importdata(savefileCMT);
%% assign variables
nSys = numel(dtwin.SysEffectsAll);

%twin
PlotVarTmp_twin = struct2array(structfun(@(x)x(1),dtwin.MultiLpar,'UniformOutput',0));  
SingleBarY_twin = sqrt(PlotVarTmp_twin(2:end).^2-PlotVarTmp_twin(1)^2);
SingleBarY_twin(7) = []; % remove HV stacking, because not approved
SingleBarY_twin(end+1) = PlotVarTmp_twin(1); % stat only 
SingleBarY_twin(end+1) = sqrt(sum(SingleBarY_twin.^2));%tmp 0.5*(dCmT.FitResult.errPos(1)-dCmT.FitResult.errNeg(1));

%data
PlotVarTmp_data = struct2array(structfun(@(x)x(1),ddata.MultiLpar,'UniformOutput',0));   
PlotVarTmp_data(PlotVarTmp_data<PlotVarTmp_data(1)) = PlotVarTmp_data(1);
SingleBarY_data = sqrt(PlotVarTmp_data(2:end).^2-PlotVarTmp_data(1)^2);
SingleBarY_data(7) = []; % remove HV stacking, because not approved
SingleBarY_data(end+1) = PlotVarTmp_data(1); % stat only 
SingleBarY_data(end+1) =  0.5*(dCmR.FitResult.errPos(1)-dCmR.FitResult.errNeg(1));
y = SingleBarY_data;
%%

f55 = figure('Name','MultiBarPlot','Renderer','painters');
set(f55, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.9, 0.8]);

leg_str = {SysEffectLeg{:},'Statistics','Total'};
LocalFontSize = 22;
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
        tstr = sprintf('%.0f\\cdot10^{-3}',y(i)*1e3) ;
    else
        tstr = sprintf('%.0f\\cdot10^{-2}',y(i)*1e2) ;
    end%0.46
    t{i}= text(0.94,SingleBarX(i),tstr,...
        'HorizontalAlignment','right','FontSize',LocalFontSize-2,...
        'Interpreter','tex',...
        'FontName','Helvetica');
    
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
          tstr = sprintf('%.0f\\cdot10^{-3}',SingleBarY_twin(ia(i))*1e3) ;
    else  
         tstr = sprintf('%.0f\\cdot10^{-2}',SingleBarY_twin(ia(i))*1e2) ;
      end%0.95
     t{i}= text(2,SingleBarX(i),tstr,...
        'HorizontalAlignment','right','FontSize',LocalFontSize-2,...
        'Interpreter','tex',...
        'FontName','Helvetica',...
        'Color',rgb('DimGray'));
end

% axis options
xlim([1e-03 1])
ylim([min(SingleBarX)-0.9 max(SingleBarX)+0.9])
cl = 0.683;

    xlabel(sprintf('1\\sigma uncertainty on {\\itm}^2_\\nu at %.1f%% C.L. (eV^{ 2})',cl*100));
 %   xlabel(sprintf('1\\sigma sensitivity on {\\itm}^2_\\nu at %.0f%% C.L. (eV^{ 2})',cl*100));

PrettyFigureFormat
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

savename = sprintf('%s/plots/KNM2_SysBreakdown_DataTwin.pdf',filedir);
export_fig(savename);