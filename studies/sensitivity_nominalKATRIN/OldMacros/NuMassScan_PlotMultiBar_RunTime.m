function NuMassScan_PlotMultiBar_RunTime(varargin)
p=inputParser;
p.addParameter('SysBudget',{'08_NoELoss','07','06'},@(x)iscell(x));  % Systematics Input Information in Result folder
p.addParameter('range',60,@(x)isfloat(x));% MTD energy range: 30,45,60eV below E0 (18575eV)
p.addParameter('MACE_Ba_T',7*1e-04,@(x)isfloat(x));
p.addParameter('WGTS_B_T',0.7*3.6,@(x)isfloat(x));
p.addParameter('BKG_RateSec','',@(x)isfloat(x));%if empty, fill accordning to range and Ba ("optimal")
p.addParameter('Anchor6G','OFF',@(x)ismember(x,{'ON','OFF'})); % Background Option
p.addParameter('TimeSec',(124/148*24*60*60).*[900, 300, 42],@(x)all(isfloat(x)));
p.addParameter('TD','MTDcreator',@(x)ischar(x)); %if empty, fill accordning to range and Ba ("optimal")
p.addParameter('SavePlot','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('SingleSys','ON',@(x)ismember(x,{'ON','OFF'})); % shows also single systematic contributions

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
Anchor6G      = p.Results.Anchor6G;

% Get Background, if not specified
if isempty(BKG_RateSec)
    BKG_RateSec = GetBackground('MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'Anchor6G',Anchor6G);
end

% Get TD, if not specified
if isempty(TD) || contains(TD,'Sensitivity')
    TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,MACE_Ba_T*1e4);
elseif contains(TD,'Optimize')
    TD  = sprintf('OptimizeMTD%.0f_NuMassFactor_BkgFraction_all_03',range);
elseif contains(TD,'MTDcreator')
    TD  = sprintf('MTDcreator_E018575.0_%.0feV_B35_Ba7.0_RedF0.7_NuMF0.59_BkgF0.20',range);
end
TDlabel1 = strrep(TD,sprintf('_Ba%.0fG',MACE_Ba_T*1e4),''); %
TDlabel = strrep(TDlabel1,'Sensitivity_','');

Ntime = numel(TimeSec);
nEffects = 4; % 1- stat, 2.TC, 3.TC + FSD, 4. TC + FSD + RF
nSingleSys   = 3; % TC,  FSD, RF
% Init 
mNu90Nplus1           = zeros(Ntime,nEffects); % 1- stat, 2.tc, 3.tc and fsd, 4 tc fsd rf
mNu90Stack            = zeros(Ntime,nEffects); % conversion to stacked format
mNu90Stack_Single     = zeros(Ntime,nSingleSys);   % for single bars
timeLabel             = cell(Ntime,1);
timePlot              = zeros(Ntime,1);
timePlotSingle        = zeros(Ntime,nSingleSys);
%% Data import

for t=1:Ntime
    results_filename    = sprintf('./results/%s_%0.0fd_SensitivityNominal_ResultsNuMassScan_MTD-%s_Ba%.0fG_Bs%.2fT_Systematics.mat',SysBudget{t},TimeSec(t)/86400,TDlabel,MACE_Ba_T*1e4,WGTS_B_T);
    results_filenameNp1 = sprintf('./results/%s_%0.0fd_SensitivityNominal_ResultsNuMassScan_MTD-%s_Ba%.0fG_Bs%.2fT_SystematicsNplus1.mat',SysBudget{t},TimeSec(t)/86400,TDlabel,MACE_Ba_T*1e4,WGTS_B_T);
    if exist(results_filename,'file')==2
        d    = importdata(results_filename); %'stat','TC','FSD','RF','all','RF_EL','RF_BF','RF_RX','RF_BFRX'
        dNplus1 = importdata(results_filenameNp1);%'stat',TC, TC+FSD,TC+FSD+RF
    else
        fprintf('File doesnt exist! \n');
        return
    end
    mNu90Nplus1(t,:)  = sqrt(dNplus1.mNu90)*1e3;
    mNu90Stack(t,:) = Convert2Stack(mNu90Nplus1(t,:),mNu90Nplus1(t,1));
    
    mNu90Stack_Single(t,1) = 1e3*sqrt(d.mNu90(2))-1e3*sqrt(d.mNu90(1)); % TC
    mNu90Stack_Single(t,2) = 1e3*sqrt(d.mNu90(3))-1e3*sqrt(d.mNu90(1)); % FSD
    mNu90Stack_Single(t,3) = 1e3*sqrt(d.mNu90(4))-1e3*sqrt(d.mNu90(1)); % RF
    
    if TimeSec(t)==(124/148*24*60*60)*42
        timeLabel{t} = sprintf('42 days');
    elseif TimeSec(t)==(124/148*24*60*60)*365
        timeLabel{t} = sprintf('365 days');
    elseif TimeSec(t)==(124/148*24*60*60)*300
        timeLabel{t} = sprintf('300 days');
    elseif TimeSec(t)==(124/148*24*60*60)*900
        timeLabel{t} = sprintf('900 days');
    end
    
    timePlot(t) = t*20;
    if Ntime~=1
        SingleBarWidth = 3.2;
        nOffset = 5.5;
    elseif Ntime==1
        SingleBarWidth=4;
        nOffset = 4.1;
    end
    timePlotSingle(t,:) = timePlot(t)-((1:SingleBarWidth:3*SingleBarWidth)+nOffset);
end

%% MultiBar plot
f5 = figure('Name','MultiBar_RunTime','Renderer','opengl'); 
set(f5, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
b = barh(timePlot,mNu90Stack,'stacked','BaseValue',0);
b(1).FaceColor = rgb('White');
b(2).FaceColor = rgb('FireBrick');
b(3).FaceColor = rgb('GoldenRod');
b(4).FaceColor = rgb('CadetBlue');
b(1).LineStyle = '--'; b(2).LineStyle = 'none';b(3).LineStyle = 'none';b(4).LineStyle = 'none';
% legends
legStack = {'Statistical','+ Theoretical Corrections','+ Final States Distribution','+ Response Function'};
leg = legend(legStack{:},'Location','southeast');

if strcmp(SingleSys,'ON')
    b(1).BarWidth = b(1).BarWidth*0.6;
    hold on;
    pnone = plot(NaN*[1 1],NaN*[1 1],'Color','w');
    bTC = barh(timePlotSingle(:,1),[mNu90Nplus1(:,1),mNu90Stack_Single(:,1)],'stacked','BarWidth',b(1).BarWidth/nSingleSys);
    bFSD = barh(timePlotSingle(:,2),[mNu90Nplus1(:,1),mNu90Stack_Single(:,2)],'stacked','BarWidth',b(1).BarWidth/nSingleSys);
    bRF = barh(timePlotSingle(:,3),[mNu90Nplus1(:,1),mNu90Stack_Single(:,3)],'stacked','BarWidth',b(1).BarWidth/nSingleSys);
    bTC(1).FaceAlpha = 0; bTC(1).LineStyle = 'none';
    bTC(2).LineStyle = 'none'; bTC(2).FaceAlpha =0.5; bTC(2).FaceColor = rgb('FireBrick');
    bFSD(1).FaceAlpha = 0;  bFSD(1).LineStyle = 'none';
    bFSD(2).LineStyle = 'none';  bFSD(2).FaceAlpha =0.5; bFSD(2).FaceColor = rgb('GoldenRod');
    bRF(1).FaceAlpha = 0; bRF(1).LineStyle = 'none';
    bRF(2).LineStyle = 'none'; bRF(2).FaceAlpha =0.5; bRF(2).FaceColor = rgb('CadetBlue');
    legSingle = {'Single Contributions','Stat + Theoretical Corrections',...
        'Stat + Final States','Stat + Response Function'};
    leg = legend([b,pnone,bTC(2),bFSD(2),bRF(2)], legStack{:},legSingle{:});
    leg.NumColumns = 2;
end
legend boxoff
yticks(timePlot);
yticklabels(timeLabel);
ytickangle(90)
xlabel('neutrino mass sensitivity 90% C.L. (meV)');
ylabel('runtime');
if Ntime==3
    ylim([min(min((timePlotSingle)))-5 max(timePlot)+10]);
    if range==60
        xmin = 260;%ceil(min(min([mNu90Stack60 mNu90Stack30]))-10);
    elseif range==30 && ~contains(TD,'MTDcreator')
        xmin=310;
    elseif range==30 && contains(TD,'MTDcreator')
        xmin = 380;
    end
    xmax = ceil(max(max(mNu90Nplus1(:,4)))*1.01);
end

xlim([xmin xmax]);
if contains(TD,'MTDcreator')
    xticks(xmin:40:xmax);
else
xticks(xmin:20:xmax);
end
PrettyFigureFormat;
PlotFontSize = 14;
set(gca,'FontSize',PlotFontSize);
grid on;
if contains(TD,'Optimize') || contains(TD,'MTDcreator')
titlestr = sprintf('Samak Optimized MTD - %.0feV scan range -  Background %.0f mcps',range,BKG_RateSec*1e3);
title(titlestr);
else
    title(sprintf('MTD %s - Background %.0f mcps ',strrep(TD,'_',' '),BKG_RateSec*1e3));
end
if strcmp(SavePlot,'ON')
    if ~exist('../sensitivity_nominalKATRIN/plots/png/RunTime/','dir')
        mkdir ../sensitivity_nominalKATRIN/plots/png/RunTime/
        mkdir ../sensitivity_nominalKATRIN/plots/pdf/RunTime/
        mkdir ../sensitivity_nominalKATRIN/plots/fig/RunTime/
    end
    save_name = sprintf('RunTime_MTD-%s_BKG-%.0f',TD,BKG_RateSec*1e3);
    print(f5,['../sensitivity_nominalKATRIN/plots/png/RunTime/',save_name,'.png'],'-dpng');
    savefig(f5,['../sensitivity_nominalKATRIN/plots/fig/RunTime/',save_name,'.fig'],'compact');
    publish_figurePDF(f5,['../sensitivity_nominalKATRIN/plots/pdf/RunTime/',save_name,'.pdf']);
    
end
end