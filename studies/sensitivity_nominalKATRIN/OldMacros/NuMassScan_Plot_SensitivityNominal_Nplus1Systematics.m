% Plot Sensitivity for 3 ranges
% as a function of Ba (and Background)
% loop over systematics
% ----------------------------------inputs -------------------------------------------
function NuMassScan_Plot_SensitivityNominal_Nplus1Systematics(varargin)
p=inputParser;
p.addParameter('SysBudget','03',@(x)ischar(x));  % Systematics Input Information in Result folder
p.addParameter('range',[30,45,60],@(x)all(isfloat(x)));% MTD energy range: 30,45,60eV below E0 (18575eV)
p.addParameter('MACE_Ba_T',7*1e-04,@(x)isfloat(x));
p.addParameter('WGTS_B_T',0.7*3.6,@(x)isfloat(x));
p.addParameter('BKG_RateSec','',@(x)isfloat(x));%if empty, fill accordning to range and Ba ("optimal")
p.addParameter('Anchor6G','OFF',@(x)ismember(x,{'ON','OFF'})); % Background Option
p.addParameter('Anchor6GValue',335e-3,@(x)isfloat(x)); % 335=FT, [26keV, 32 keV]ROI
p.addParameter('TimeSec',3*365*24*60*60,@(x)isfloat(x));
p.addParameter('TD','',@(x)ischar(x)); %if empty, fill accordning to range and Ba ("optimal")
p.addParameter('SavePlot','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('SingleSys','ON',@(x)ismember(x,{'ON','OFF'})); % shows also single systematic contributions
p.addParameter('TDRSys','OFF',@(x)ismember(x,{'ON','OFF'})); % Show also TDR systematics budget (not for TD=DR30 (there it's shwon anyway)
p.parse(varargin{:});
SysBudget     = p.Results.SysBudget;
MACE_Ba_T     = p.Results.MACE_Ba_T;
WGTS_B_T      = p.Results.WGTS_B_T;
range         = p.Results.range;
TimeSec       = p.Results.TimeSec;
TD            = p.Results.TD;
SavePlot      = p.Results.SavePlot;
SingleSys     = p.Results.SingleSys;
TDRSys        = p.Results.TDRSys;
BKG_RateSec   = p.Results.BKG_RateSec;
Anchor6G      = p.Results.Anchor6G;
Anchor6GValue = p.Results.Anchor6GValue;

if strcmp(TD,'DR30')
    TDRSys = 'OFF'; %displayed in a different way
end
% ----------------------------------inputs end -------------------------------------------
Nranges        = numel(range);
nSys           = 4; %stat, stat+TC, stat+TC+FSD, %stat+TC+FSD+RF
nSingleSys     = 3; %stat + ...'TC','FSD','RF'
mNu90Nplus1    = zeros(Nranges,nSys);
mNu90SingleSys  = zeros(Nranges,nSingleSys);
range_SingleSys = zeros(Nranges,nSingleSys);
range_TDR       = zeros(Nranges,1);

% Data Import
for i=1:Nranges
    if isempty(TD) || contains(TD,'Sensitivity') 
        TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range(i),MACE_Ba_T*1e4);
    elseif contains(TD,'Optimize')
        TD  = sprintf('OptimizeMTD%.0f_NuMassFactor_BkgFraction_all_03',range(i));
        if strcmp(SysBudget,'08') || strcmp(SysBudget,'09')
            SysBudget = [SysBudget,'_NoELoss'];
        end
    elseif contains(TD,'MTDcreator') && contains(TD,'Ba7.0')
        TD  = sprintf('MTDcreator_E018575.0_%.0feV_B35_Ba7.0_RedF0.7_NuMF0.59_BkgF0.20',range(i));
        if strcmp(SysBudget,'08') || strcmp(SysBudget,'09')
            SysBudget = [SysBudget,'_NoELoss'];
        end
    elseif contains(TD,'eV_B35_Ba3.0_RedF1.0_NuMF0.40_BkgF0.10_B9')
        TD  = sprintf('MTDcreator_E018575.0_%0.feV_B35_Ba3.0_RedF1.0_NuMF0.40_BkgF0.10_B9',range(i));
        if strcmp(SysBudget,'08')
            SysBudget = [SysBudget,'_NoELoss'];
        end
    end
    TDlabel1 = strrep(TD,sprintf('_Ba%.0fG',MACE_Ba_T*1e4),''); %
    TDlabel = strrep(TDlabel1,'Sensitivity_','');

    save_nameNplus1 = sprintf('./results/%s_%0.0fd_SensitivityNominal_ResultsNuMassScan_MTD-%s_Ba%.0fG_Bs%.2fT_SystematicsNplus1.mat',SysBudget,round(TimeSec/(86400)),TDlabel,MACE_Ba_T*1e4,WGTS_B_T);
    %disp(save_nameNplus1);
    
    if contains(TD,'B9')
            TDlabel1 = strrep(TD,sprintf('_Ba%.0fG',MACE_Ba_T*1e4),''); %
            TDlabel = strrep(TDlabel1,'Sensitivity_','');
            save_nameNplus1 = sprintf('./results/%s_%0.0fd_SensitivityNominal_ResultsNuMassScan_MTD-MTDcreator_E018575.0_%0.0feV_B35_Ba3.0_RedF1.0_NuMF0.40_BkgF0.10_B9_Ba3G_Bs3.60T_SystematicsNplus1.mat',SysBudget,round(TimeSec/(86400)),range(i));
            %disp(save_nameNplus1);
    end
    
    if exist(save_nameNplus1,'file')==2
        dNplus1 = importdata(save_nameNplus1);
    else
        disp(save_nameNplus1); fprintf('File doesnt exist! \n');
        return
    end
    mNu90Nplus1(i,:)  = dNplus1.mNu90;
    
    if strcmp(SingleSys,'ON')
        save_nameSingleSys = sprintf('./results/%s_%0.0fd_SensitivityNominal_ResultsNuMassScan_MTD-%s_Ba%.0fG_Bs%.2fT_Systematics.mat',SysBudget,round(TimeSec/(86400)),TDlabel,MACE_Ba_T*1e4,WGTS_B_T);
        
        if contains(TD,'B9')
            save_nameSingleSys = sprintf('./results/%s_%0.0fd_SensitivityNominal_ResultsNuMassScan_MTD-MTDcreator_E018575.0_%0.0feV_B35_Ba3.0_RedF1.0_NuMF0.40_BkgF0.10_B9_Ba3G_Bs3.60T_Systematics.mat',SysBudget,round(TimeSec/(86400)),range(i));
            %disp(save_nameNplus1);
        end
        
        
        if exist(save_nameSingleSys,'file')==2
            dSingleSys = importdata(save_nameSingleSys);
        else
            fprintf('File doesnt exist! \n');
            return
        end
        mNu90SingleSys(i,:) = dSingleSys.mNu90(2:4);
        if Nranges==1
        range_SingleSys(i,:) = range(i)-((0:0.135:(nSingleSys-1)/6)+0.47);  
          range_TDR(i) = range(i)+1.5;
        elseif Nranges==2
           range_SingleSys(i,:) = range(i)-([0,2,4]+7); 
             range_TDR(i) = range(i)+3.5;
        elseif Nranges==3
        range_SingleSys(i,:) = range(i)-((0:1:nSingleSys-1)+3.5);      
          range_TDR(i) = range(i)+1.5;
        end
      
    end
    SysBudget = strrep(SysBudget,'_NoELoss','');
end

%range_SingleSys = sort(reshape(range_SingleSys,1,numel(range_SingleSys)));
%% MultiBar Plot
s = GetSysBudget('SysBudget',SysBudget);
if isempty(BKG_RateSec)
BKG_RateSec = GetBackground('MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'Anchor6G',Anchor6G,'Anchor6GValue',Anchor6GValue);
end
close;
f55 = figure('Name','MultiBarPlot','Renderer','opengl');
set(f55, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
% Convert to Stacking Format: N+1
PlotmNu90 = 1e3*mNu90Nplus1.^0.5; % convert to sensitivity on neutrino mass (meV)
t = zeros(Nranges,nSys);
rangeLabel = cell(Nranges,1);
for i=1:Nranges
    t(i,:) = Convert2Stack(PlotmNu90(i,:),PlotmNu90(i,1));
    if t(i,3)<0 % then TC+FSD is smaller than TC
        t(i,3)=0;
    end
    rangeLabel{i} = sprintf('%.0f eV',range(i));
end
% Convert to Stacking Format: N + stat
 mNuStat = sqrt(mNu90Nplus1(:,1))'*1e3;
 mNuSys = zeros(Nranges,nSingleSys);
 for i=1:nSingleSys
 mNuSys(:,i) = (sqrt(mNu90SingleSys(:,i))'*1e3-mNuStat)';
 end
 % In case only 1 range: add one invisible row (workaround barh)
if Nranges==1
    t = [t;NaN.*zeros(1,nSys)];
    range =  [range,range+1];
    rangeLabel = {rangeLabel{:};' '};
    range_SingleSys = [range_SingleSys;range_SingleSys+1];
    mNuSys = [mNuSys;NaN.*zeros(1,nSingleSys)];
    mNuStat = [mNuStat,NaN];
end

b = barh(range,t,'stacked');
b(1).LineStyle = '--';
b(2).LineStyle = 'none';b(3).LineStyle = 'none'; b(4).LineStyle = 'none';
b(1).FaceColor = rgb('White');
b(2).FaceColor = rgb('FireBrick');
b(3).FaceColor = rgb('GoldenRod');
b(4).FaceColor = rgb('CadetBlue');
leg_str = {' Statistical';'+ Theoretical Corrections';'+ Final State Distribution';'+ Response Function'};
switch SingleSys
    case 'ON'
        if Nranges~=1
        b(1).BarWidth = b(1).BarWidth/2;
        end
        hold on;
       % if Nranges==3
        SingleSysBarWidth = b(1).BarWidth/(2*nSingleSys);
        %elseif  Nranges==2
       %   SingleSysBarWidth = b(1).BarWidth/(2*nSingleSys);   
        %end
        pnone = plot(NaN*[1 1],NaN*[1 1],'Color','w');
        bTC  = barh(range_SingleSys(:,1)',[mNuStat;mNuSys(:,1)']','stacked','BarWidth',SingleSysBarWidth);
        bFSD = barh(range_SingleSys(:,2)', [mNuStat;mNuSys(:,2)']','stacked','BarWidth',SingleSysBarWidth);
        bRF  = barh(range_SingleSys(:,3)', [mNuStat;mNuSys(:,3)']','stacked','BarWidth',SingleSysBarWidth);
        bTC(2).FaceColor =rgb('FireBrick'); bTC(2).LineStyle = 'none'; bTC(2).FaceAlpha = 0.5;
        bFSD(2).FaceColor =rgb('GoldenRod'); bFSD(2).LineStyle = 'none'; bFSD(2).FaceAlpha = 0.5;
        bRF(2).FaceColor =rgb('CadetBlue'); bRF(2).LineStyle = 'none'; bRF(2).FaceAlpha = 0.5;
        bTC(1).FaceColor = 'w';bFSD(1).FaceColor = 'w'; bRF(1).FaceColor = 'w';
        bTC(1).LineStyle = 'none'; bFSD(1).LineStyle = 'none'; bRF(1).LineStyle = 'none'; bTC(1).FaceAlpha = 0;  
        leg_str = {leg_str{:},'Single Contributions','Stat + Theoretical Corrections', 'Stat + Final State Distribution','Stat + Response Function'};
       leg = legend([b,pnone,bTC(2),bFSD(2),bRF(2)],leg_str{:});
        if strcmp(TDRSys,'ON') && ~strcmp(TD,'DR30')
            hold on;
            %TDR Systemaitcs for higher range: TDR + Marco Kleesiek's value due to rhod sigma and Ba
            if Nranges==3 % 30, 45 ,60
            TDRSys90 = sqrt(sqrt(mNu90Nplus1(:,1).^2+(1.64.*[0.0213;0.0220;0.0226]).^2))*1e3-mNuStat';
            elseif Nranges==2 % 30 and 60
              TDRSys90 = sqrt(sqrt(mNu90Nplus1(:,1).^2+(1.64.*[0.0213;0.0226]).^2))*1e3-mNuStat';    
            end
            bTDR =   barh(range_TDR,[mNuStat',TDRSys90, PlotmNu90(:,4)-TDRSys90],'stacked','FaceColor',rgb('SlateGray'),'FaceAlpha',1);
            bTDR(1).BarWidth = bTC.BarWidth; bTDR(2).BarWidth = b(2).BarWidth/2;
            bTDR(1).FaceColor = 'w'; bTDR(1).LineStyle='none';bTDR(2).LineStyle='none';bTDR(1).FaceAlpha = 0;
            bTDR(3).FaceColor = 'w'; bTDR(3).LineStyle = 'none';
            leg_str = {leg_str{1:4},'TDR + Neutrino16',leg_str{5:end}};
            leg = legend([b,bTDR(2),pnone,bTC(2),bFSD(2),bRF(2)],leg_str{:});      
        end

    case 'OFF'
        leg = legend(leg_str{:});
end

leg.NumColumns = 2;
legend boxoff
yticklabels(rangeLabel);
ytickangle(90)
xlabel('neutrino mass sensitivity 90% C.L. (meV)');
ylabel('MTD range');
if strcmp(TD,'DR30')
    hold on;
    mNu90DRStat = sqrt(0.018*1.64)*1e3;
    mNu90DRSys  = sqrt(sqrt(0.018^2+0.017^2)*1.64)*1e3-mNu90DRStat;
    bTDR = barh([31;NaN],[mNu90DRStat,mNu90DRSys;NaN,NaN],'Stacked');
    bTDR(1).FaceColor = 'w'; bTDR(1).LineStyle = '--';
    bTDR(2).FaceColor = rgb('SlateGray'); bTDR(2).LineStyle = 'none';
    leg_str = {leg_str{:},'TDR04 Statistical    (\sigma = 0.018eV^2)','TDR04 Systematics (\sigma = 0.017eV^2)'}; %_{m²_\nu}
    if strcmp(SingleSys,'ON')
    leg = legend([b,pnone,bTC(2),bFSD(2),bRF(2),bTDR(1),bTDR(2)],leg_str{:});
    else
    leg = legend([b,bTDR(1),bTDR(2)],leg_str{:});   
    end
    ylim([min(range)-1.5; max(range)+1.3]);
    xlim([160 210]);
    xticks((160:5:210));
    leg.NumColumns = 3;
    % ylabel('MTD');
    yticklabels({'30 eV','30 eV TDR'});
elseif Nranges==1
    ylim([min(range)-1.5; max(range)]);
end
if ~strcmp(TD,'DR30')
    xmax = sqrt(max(mNu90Nplus1(:,end)))*1e3+10;
    xlim([200 xmax]);
    xticks((200:20:xmax));
% Thierry
if contains(TD,'MTDcreator')
%        xlim([550, 700]);
%        xticks((550:10:700));
       xmin =floor(min(mNuStat)/10)*10-10;
        xlim([xmin, xmax]);
        xticks((xmin:10:xmax));
    elseif contains(TD,'MTDcreator')
        xmin =floor(min(mNuStat)/10)*10-10;
        xlim([xmin, xmax]);
        xticks((xmin:10:xmax));
    end
end
grid on;

% Systematic Budget Annotation
SysUncertainties = [sprintf('Response Function: '),...
    sprintf('\\Delta \\rhod\\sigma = %.1f%%, ',100*sqrt(s.WGTS_CD_MolPerCm2_RelErr^2+s.ISXsection_RelErr^2)),...
    sprintf('\\DeltaB_a = %s, \\DeltaB_{T/max} = %.1f%%, ',s.MACE_Ba_T_Err_str,100*s.MACE_Bmax_T_RelErr),...
    sprintf('\nEnergy Loss Uncertainties from Eur. Phys. J. D 10, 39-52'),...
    s.ELoss,...
    sprintf('\nFinal State Distribution:'),...
    sprintf(' Normalization %.0f %%, ',s.FSDNorm_RelErr*100),...
    sprintf('Bin-to-Bin uncorrelated %.0f %% (GS), %.0f %% (ES)',100*s.FSDShapeGS_RelErr,100*s.FSDShapeES_RelErr)];
a=annotation('textbox', [0.14 0.1 1 0.12], ...
    'String', SysUncertainties, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left');
a.FontSize=11;a.FontWeight='bold';
if contains(TD,'MTDcreator') || contains(TD,'Optimize')
    labelTime = TimeSec/(24*60*60*124/148);
else
    labelTime = TimeSec/(24*60*60);
end
TDlabel2 = strrep(strrep(TDlabel1,sprintf('_%.0feV',range(end)),'Opt'),'_',' ');
title(sprintf('Nominal KATRIN %.0f days, MTD: %s \nB_a = %.0fG, B_{T/max} = %.0f %%, Background = %.0f mcps',...
    labelTime,TDlabel2, MACE_Ba_T*1e4,WGTS_B_T/3.6*100,BKG_RateSec*1e3));

PrettyFigureFormat;
set(gca,'FontSize',16);
leg.FontSize = 12;
% save plot
if strcmp(SavePlot,'ON')
    if ~exist('../sensitivity_nominalKATRIN/plots/pdf/SysBreakdown/','dir')
        mkdir ../sensitivity_nominalKATRIN/plots/pdf/SysBreakdown/
        mkdir ../sensitivity_nominalKATRIN/plots/png/SysBreakdown/
        mkdir ../sensitivity_nominalKATRIN/plots/fig/SysBreakdown/
    end
    if strcmp(TD,'DR30')
        save_name = sprintf('%s_%0.0fd_SensitivityNominalScan_SystematicBreakdownBar_%s_Ba%.0fG_B%.0fSingleSys%s',SysBudget,round(TimeSec/(86400)),strrep(TDlabel1,sprintf('_%.0feV',range(end)),'Opt'),MACE_Ba_T*1e4,WGTS_B_T/3.6*100,SingleSys);
    else
        save_name = sprintf('%s_%0.0fd_SensitivityNominalScan_SystematicBreakdownBar_%s_Ba%.0fG_B%.0fSingleSys%s_TDRSys%s',SysBudget,round(TimeSec/(86400)),strrep(TDlabel1,sprintf('_%.0feV',range(end)),'Opt'),MACE_Ba_T*1e4,WGTS_B_T/3.6*100,SingleSys,TDRSys);
    end
    print(f55,['../sensitivity_nominalKATRIN/plots/png/SysBreakdown/',save_name,'.png'],'-dpng');
    savefig(f55,['../sensitivity_nominalKATRIN/plots/fig/SysBreakdown/',save_name,'.fig'],'compact');
    publish_figurePDF(f55,['../sensitivity_nominalKATRIN/plots/pdf/SysBreakdown/',save_name,'.pdf']);

end


