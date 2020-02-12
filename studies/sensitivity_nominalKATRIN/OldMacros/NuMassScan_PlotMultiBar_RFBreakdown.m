function NuMassScan_PlotMultiBar_RFBreakdown(varargin)
p=inputParser;
p.addParameter('SysBudget','03',@(x)ischar(x));  % Systematics Input Information in Result folder
p.addParameter('range',[30,45,60],@(x)all(isfloat(x)));% MTD energy range: 30,45,60eV below E0 (18575eV)
p.addParameter('MACE_Ba_T',7*1e-04,@(x)isfloat(x));
p.addParameter('WGTS_B_T',0.7*3.6,@(x)isfloat(x));
p.addParameter('BKG_RateSec','',@(x)isfloat(x));%if empty, fill accordning to range and Ba ("optimal")
p.addParameter('TimeSec',3*365*24*60*60,@(x)isfloat(x));
p.addParameter('TD','',@(x)ischar(x)); %if empty, fill accordning to range and Ba ("optimal")
p.addParameter('SavePlot','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('SingleSys','ON',@(x)ismember(x,{'ON','OFF'})); % shows also single systematic contributions
p.addParameter('Anchor6G','OFF',@(x)ismember(x,{'ON','OFF'})); % Background Option
p.addParameter('Anchor6GValue',335e-3,@(x)isfloat(x)); % 335=FT, [26keV, 32 keV]ROI

p.parse(varargin{:});
SysBudget     = p.Results.SysBudget;
MACE_Ba_T     = p.Results.MACE_Ba_T;
WGTS_B_T      = p.Results.WGTS_B_T;
range         = p.Results.range;
TimeSec       = p.Results.TimeSec;
TD            = p.Results.TD;
SavePlot      = p.Results.SavePlot;
SingleSys     = p.Results.SingleSys;
BKG_RateSec   = p.Results.BKG_RateSec;
Anchor6G     = p.Results.Anchor6G;
Anchor6GValue = p.Results.Anchor6GValue;

if isempty(BKG_RateSec)
    BKG_RateSec = GetBackground('MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T);
end

s = GetSysBudget('SysBudget',SysBudget); % struct with information about systematics input
Nranges = numel(range);
nEffects = 4; % 1- stat, 2.rhod, 3.rhod and b-fields, 4.all RF
nSysRF   = 3; % rhod only,  bfields only, eloss only
% Init 
mNu90Nplus1     = zeros(Nranges,nEffects); % 1- stat, 2.rhod, 3.rhod and b-fields, 4.all RF
mNu90Stack      = zeros(Nranges,nEffects); % conversion to stacked format
mNu90Stack_Single     = zeros(Nranges,nSysRF);   % for single bars
rangeLabel      = cell(Nranges,1);
rangePlot       = zeros(Nranges,1);
rangePlotSingle = zeros(Nranges,nSysRF);
%% Import Results

for i=1:Nranges
    if isempty(TD) || contains(TD,'Sensitivity')
        TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range(i),MACE_Ba_T*1e4);
    elseif contains(TD,'Optimize')
        TD  = sprintf('OptimizeMTD%.0f_NuMassFactor_BkgFraction_all_03',range(i));
        if strcmp(SysBudget,'08') || strcmp(SysBudget,'09')
            SysBudget = [SysBudget,'_NoELoss'];
        end
    elseif contains(TD,'MTDcreator')
        TD  = sprintf('MTDcreator_E018575.0_%.0feV_B35_Ba7.0_RedF0.7_NuMF0.59_BkgF0.20',range(i));
        if strcmp(SysBudget,'08') || strcmp(SysBudget,'09')
            SysBudget = [SysBudget,'_NoELoss'];
        end
    end
    TDlabel1 = strrep(TD,sprintf('_Ba%.0fG',MACE_Ba_T*1e4),''); %
    TDlabel = strrep(TDlabel1,'Sensitivity_','');
    results_filename = sprintf('./results/%s_%0.0fd_SensitivityNominal_ResultsNuMassScan_MTD-%s_Ba%.0fG_Bs%.2fT_Systematics.mat',SysBudget,TimeSec/86400,TDlabel,MACE_Ba_T*1e4,WGTS_B_T);
    if exist(results_filename,'file')==2
        d = importdata(results_filename);
    else
            disp(results_filename);
            fprintf('File doesnt exist! \n');
        return
    end
    mNu90Nplus1(i,1) = 1e3*sqrt(d.mNu90(1)); % stat
    mNu90Nplus1(i,2) = 1e3*sqrt(d.mNu90(8)); % RF - rho d
    mNu90Nplus1(i,3) = 1e3*sqrt(d.mNu90(9)); % RF - rho d and b-fields
    mNu90Nplus1(i,4) = 1e3*sqrt(d.mNu90(4)); % RF - all
    mNu90Stack(i,:) = Convert2Stack(mNu90Nplus1(i,:),mNu90Nplus1(i,1));
    
    mNu90Stack_Single(i,1) = 1e3*sqrt(d.mNu90(8))-1e3*sqrt(d.mNu90(1)); % rhod only 
    mNu90Stack_Single(i,2) = 1e3*sqrt(d.mNu90(7))-1e3*sqrt(d.mNu90(1)); % bfields only
    mNu90Stack_Single(i,3) = 1e3*sqrt(d.mNu90(6))-1e3*sqrt(d.mNu90(1)); % eloss only
    
    rangeLabel{i} = sprintf('%.0f eV',range(i));
    rangePlot(i) = i*20;
    if Nranges~=1
    SingleBarWidth = 2.6;
    nOffset = 4.3;
    elseif Nranges==1
      SingleBarWidth=4;
      nOffset = 4.1;
    end
    rangePlotSingle(i,:) = rangePlot(i)-((1:SingleBarWidth:3*SingleBarWidth)+nOffset);
     SysBudget = strrep(SysBudget,'_NoELoss','');
end
%% special case: TDR
if strcmp(TD,'DR30')
    rangePlotSingle(1,1)=rangePlotSingle(1,1)-0.6;
    rangePlot = [rangePlot; 40];
    rangePlotSingle = [rangePlotSingle; 37,38,39];
    mNu90Stack = [mNu90Stack; NaN*zeros(1,4)];
    mNu90Nplus1 = [mNu90Nplus1; NaN*zeros(1,nEffects)];
    mNu90Stack_Single = [mNu90Stack_Single; NaN*zeros(1,nSysRF)];
end
%% MultiBar plot
f5 = figure('Name','MultiBar_RFBreakdown','Renderer','opengl'); 
set(f5, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
if strcmp(SysBudget,'08') || strcmp(SysBudget,'09')
  b = barh(rangePlot,mNu90Stack(:,1:3),'stacked','BaseValue',0);  
else
b = barh(rangePlot,mNu90Stack,'stacked','BaseValue',0);
end
b(1).FaceColor = rgb('White');
b(3).FaceColor = rgb('DarkCyan');
b(2).FaceColor = rgb('PowderBlue');
b(1).LineStyle = '--'; b(2).LineStyle = 'none';b(3).LineStyle = 'none';
% legends
if strcmp(SysBudget,'08') || strcmp(SysBudget,'09')
 legStack = {'Statistical','+ Column Density and Inel. Cross Section','+ Magnetic Fields'};   
else
    b(4).FaceColor = rgb('Orange');b(4).LineStyle = 'none';
legStack = {'Statistical','+ Column Density and Inel. Cross Section','+ Magnetic Fields','+ Energy Loss Functions'};
end
leg = legend(legStack{:},'Location','northwest');

if strcmp(SingleSys,'ON')
    b(1).BarWidth = b(1).BarWidth/2;
    hold on;
    pnone = plot(NaN*[1 1],NaN*[1 1],'Color','w');
    bRX = barh(rangePlotSingle(:,1),[mNu90Nplus1(:,1),mNu90Stack_Single(:,1)],'stacked','BarWidth',b(1).BarWidth/nSysRF);
    bBF = barh(rangePlotSingle(:,2),[mNu90Nplus1(:,1),mNu90Stack_Single(:,2)],'stacked','BarWidth',b(1).BarWidth/nSysRF);
    bRX(1).FaceAlpha = 0; bRX(1).LineStyle = 'none';
    bRX(2).LineStyle = 'none'; bRX(2).FaceAlpha =0.5; bRX(2).FaceColor = rgb('PowderBlue');
    bBF(1).FaceAlpha = 0;  bBF(1).LineStyle = 'none';
    bBF(2).LineStyle = 'none';  bBF(2).FaceAlpha =0.5; bBF(2).FaceColor = rgb('DarkCyan');
    if ~strcmp(SysBudget,'08') && ~strcmp(SysBudget,'09')
        bEL = barh(rangePlotSingle(:,3),[mNu90Nplus1(:,1),mNu90Stack_Single(:,3)],'stacked','BarWidth',b(1).BarWidth/nSysRF);
        bEL(1).FaceAlpha = 0; bEL(1).LineStyle = 'none';
        bEL(2).LineStyle = 'none'; bEL(2).FaceAlpha =0.5; bEL(2).FaceColor = rgb('Orange');
        legSingle = {'Single Contributions','Stat + Column Density and Inel. Cross Section',...
            'Stat + Magnetic Fields','Stat + Energy Loss Functions'};
        leg = legend([b,pnone,bRX(2),bBF(2),bEL(2)], legStack{:},legSingle{:});
    else
         legSingle = {'Single Contributions','Stat + Column Density and Inel. Cross Section',...
            'Stat + Magnetic Fields'};
          leg = legend([b,pnone,bRX(2),bBF(2)], legStack{:},legSingle{:});
    end
    leg.NumColumns = 2;
end
legend boxoff
yticks(rangePlot);
yticklabels(rangeLabel);
ytickangle(90)
xlabel('neutrino mass sensitivity 90% C.L. (meV)');
ylabel('Scan energy range below Endpoint');
if Nranges==3
    ylim([min(min((rangePlotSingle)))-20 max(rangePlot)+30]);
    xmin = 200;%ceil(min(min([mNu90Stack60 mNu90Stack30]))-10);
    xmax = ceil(max(max(mNu90Nplus1(:,4)))*1.01);
elseif Nranges==2
    ylim([min(min((rangePlotSingle)))-8 max(rangePlot)+15]);
    xmin = floor(min(min(mNu90Nplus1))/10)*10-10;
    xmax = ceil(max(max(mNu90Nplus1(:,4)))*1.01);
elseif Nranges==1
    ylim([min(min((rangePlotSingle)))-10 max(rangePlot)-5]);
    xmin = floor(min(min(mNu90Nplus1))/10)*10;
    xmax = ceil(max(max(mNu90Nplus1(:,4)))*1.01);
end
xlim([xmin xmax]);
xticks(xmin:10:xmax);
if contains(TD,'MTDcreator') || contains(TD,'Optimize')
    labelTime = TimeSec/(24*60*60*124/148);
else
    labelTime = TimeSec/(24*60*60);
end
title(sprintf('Response Function Uncertainty Breakdown \nNominal KATRIN %.0f days, B_a = %.0fG, B_{T/max} = %.0f %%, Background = %.0f mcps',...
    labelTime, MACE_Ba_T*1e4,WGTS_B_T/3.6*100,1e3*BKG_RateSec));
PrettyFigureFormat;
PlotFontSize = 14;
set(gca,'FontSize',PlotFontSize);
grid on;

SysUncertainties = [sprintf('Response Function: '),...
    sprintf('\\Delta \\rhod\\sigma = %.1f%%, ',100*sqrt(s.WGTS_CD_MolPerCm2_RelErr^2+s.ISXsection_RelErr^2)),...
    sprintf('\\DeltaB_a = %s, \\DeltaB_{T/max} = %.1f%%',s.MACE_Ba_T_Err_str,100*s.MACE_Bmax_T_RelErr),...
    s.ELoss];
a=annotation('textbox', [0.14 0.1 1 0.10], ...
    'String', SysUncertainties, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left');
a.FontSize=11;a.FontWeight='bold';

if strcmp(SavePlot,'ON')
  TDlabel = strrep(TDlabel1,['_',num2str(range(end)),'eV'],'');
    save_name = sprintf('%s_%0.0fd_MTD-%s_SensitivityNominal_NuMassScan_MultiBar_RFBreakdown_Bfields%.0fpercent_Ba%.0fG_SingleSys%s',SysBudget,TimeSec/86400,TDlabel,WGTS_B_T/3.6*100,MACE_Ba_T*1e4,SingleSys);
    if ~exist('../sensitivity_nominalKATRIN/plots/png/RFBreakdown/','dir')
        mkdir ../sensitivity_nominalKATRIN/plots/png/RFBreakdown/
        mkdir ../sensitivity_nominalKATRIN/plots/pdf/RFBreakdown/
        mkdir ../sensitivity_nominalKATRIN/plots/fig/RFBreakdown/
    end
    print(f5,['../sensitivity_nominalKATRIN/plots/png/RFBreakdown/',save_name,'.png'],'-dpng');
    savefig(f5,['../sensitivity_nominalKATRIN/plots/fig/RFBreakdown/',save_name,'.fig'],'compact');
    publish_figurePDF(f5,['../sensitivity_nominalKATRIN/plots/pdf/RFBreakdown/',save_name,'.pdf']);
end

end


