% plot script for KSN1 sterile contour plot
% ratio -> syst over stat
% Lisa, April 2020
%% settings

CL  =95;%68.27;        % 1 sigma
Npar = 1;           % confidence level condition for Npar parameter
SavePlot = 'OFF';
PlotSplines = 'ON'; % plot smooth spline instead of points -> fails for closed contours
range = [95:-5:45,41,40];%;%:-5:75);%
nGridSteps = 50;
DataType = 'Twin';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'OFF';
SysBudget = 24;
PlotStatDom = 'ON';
%% stat. only: load grid (or calculate if doesn't exist)
nContours = numel(range);

mnu4SqStat = cell(nContours,1);
sin2T4Stat = cell(nContours,1);
for i=1:nContours
    GridArg = {'range',range(i),...
        'nGridSteps',nGridSteps,...
        'chi2','chi2Stat',...
        'freePar',freePar,...
        'RunList',RunList,...
        'SmartGrid',SmartGrid,...
        'RecomputeFlag','OFF',...
        'SysBudget',SysBudget};
    [mnu4SqTmp,sin2T4Tmp,chi2Tmp,chi2_refTmp,~] = KSN1GridSearch('DataType',DataType,GridArg{:});
    [mnu4SqStat{i}, sin2T4Stat{i}] = ...
        KSN1Grid2Contour(mnu4SqTmp,sin2T4Tmp,chi2Tmp,chi2_refTmp,CL,'Mode','Old','Npar',Npar);
    
    %     if strcmp(DataType,'Twin')
    %        mnu4SqStat{i} =  cell2mat(mnu4SqStat{i});
    %         sin2T4Stat{i} = cell2mat( sin2T4Stat{i});
    %     end
end

%% stat. % syst.: load grid (or calculate if doesn't exist)
mnu4SqCM = cell(nContours,1);
sin2T4CM = cell(nContours,1);

for i=1:numel(range)
    GridArg = {'range',range(i),...
        'nGridSteps',nGridSteps,...
        'chi2','chi2CMShape',...
        'freePar',freePar,...
        'RunList',RunList,...
        'SmartGrid',SmartGrid,...
        'RecomputeFlag','OFF',...
        'SysBudget',SysBudget};
    [mnu4SqTmp,sin2T4Tmp,chi2Tmp,chi2_refTmp,~] = KSN1GridSearch('DataType',DataType,GridArg{:});
    [mnu4SqCM{i}, sin2T4CM{i}] = ...
        KSN1Grid2Contour(mnu4SqTmp,sin2T4Tmp,chi2Tmp,chi2_refTmp,CL,'Mode','Old','Npar',Npar);
    
    %     if strcmp(DataType,'Twin')
    %         mnu4SqCM{i} =  cell2mat(mnu4SqCM{i});
    %         sin2T4CM{i} = cell2mat( sin2T4CM{i});
    %     end
end
%
% mnu4SqCM = cell2mat(mnu4SqCM');
% sin2T4CM = cell2mat(sin2T4CM');
%% calculate syst only (1 sigma)
mnu4SqCommon   = cell(nContours,1);
sin2T4SystOnly = cell(nContours,1);
sin2T4StatOnly = cell(nContours,1);
for i=1:nContours
    mNuCMtmp = mnu4SqCM{i};  mNuStattmp = mnu4SqStat{i};
    sinCMtmp = sin2T4CM{i};  sinStattmp = sin2T4Stat{i};
    
    [LogicCM,LogicStat]= ismember(mNuCMtmp,mNuStattmp);
    mnu4SqCommon{i}   = mNuCMtmp(LogicCM); %equivalent to mNuStattmp(LogicStat,i);
    sin2T4StatOnly{i} = sinStattmp(LogicStat);
    sin2T4SystOnly{i} = sqrt(sinCMtmp(LogicCM).^2-sinStattmp(LogicStat).^2);
end

%% some plot arguments
Colors =  {rgb('DodgerBlue'),rgb('Orange'),rgb('DarkSlateGray'),rgb('FireBrick'),...
                rgb('Magenta'),rgb('LimeGreen'),rgb('CadetBlue'),rgb('Navy'),...
                rgb('ForestGreen'),rgb('PowderBlue'),rgb('Pink'),rgb('DarkOrange'),rgb('Black'),...
                rgb('ForestGreen'),rgb('PowderBlue'),rgb('Pink'),rgb('DarkOrange')};
LineStyles = {'-','-.',':','--','-','-.',':','--','-','-.',':','--','-','-.',':','--'};
legStr = cell(nContours+1,1);
legStr{1} = sprintf('\\sigma^2_{syst.} = \\sigma^2_{stat.}');
plotdir  = [getenv('SamakPath'),'ksn1ana/LisaSterile/plots/'];
MakeDir(plotdir);

%% interpolation
if ismember(PlotSplines,{'ON','Lin'})
    if strcmp(PlotSplines,'ON')
        intArg = 'spline';
    else
        intArg = 'nearest';
    end
    mnu4SqCommonPlot   = cell(nContours,1);
    sin2T4SystOnlyPlot = cell(nContours,1);
    sin2T4StatOnlyPlot = cell(nContours,1);
    sin2T4CMPlot       = cell(nContours,1);
    mnu4SqCMPlot       = cell(nContours,1);
    
    for i=1:nContours
        mnu4SqCommonPlot{i}   = logspace(log10(min(mnu4SqCommon{i})),log10((range(i)-5)^2),1e4);
        sin2T4SystOnlyPlot{i} = interp1(mnu4SqCommon{i},sin2T4SystOnly{i},mnu4SqCommonPlot{i},intArg);
        sin2T4StatOnlyPlot{i} = interp1(mnu4SqCommon{i},sin2T4StatOnly{i},mnu4SqCommonPlot{i},intArg);
        
         mnu4SqCMPlot{i} = logspace(log10(min(mnu4SqCM{i})),log10(max(mnu4SqCM{i})),1e3);
         sin2T4CMPlot{i}= interp1(mnu4SqCM{i},sin2T4CM{i},mnu4SqCMPlot{i},'spline');
    end
else
    mnu4SqCommonPlot   =  mnu4SqCommon;
    sin2T4SystOnlyPlot =  sin2T4SystOnly;
    sin2T4StatOnlyPlot =  sin2T4StatOnly;
    sin2T4CMPlot = sin2T4CM;
    mnu4SqCMPlot = mnu4SqCM;
end

%% get ratio syst % stat
StatDomFraction = zeros(nContours,1); % fraction of m4 that are stat. dominated;
 mnu4SqCommonLin   = cell(nContours,1);

for i=1:nContours
    mnu4SqCommonLin{i}   = linspace(min(mnu4SqCommon{i}),max(mnu4SqCommon{i}),1e5);
    sin2T4Syst_tmp = interp1(mnu4SqCommon{i},sin2T4SystOnly{i},mnu4SqCommonLin{i},'lin');
    sin2T4Stat_tmp = interp1(mnu4SqCommon{i},sin2T4StatOnly{i},mnu4SqCommonLin{i},'lin');
       
    Ratio = sin2T4Syst_tmp.^2./sin2T4Stat_tmp.^2;
    StatDomFraction(i) = sum(Ratio<=1)./numel(Ratio);
    
end
%% ratio:  syst only/ stat only
% pl = cell(nContours,1);
% GetFigure;
% pref = plot(ones(10,1),logspace(0,4,10),'k-','LineWidth',2);
% hold on;
% %Colors = (colorcube(nContours));
% for i=1:nContours  
%     pl{i}= plot(sin2T4SystOnlyPlot{i}.^2./sin2T4StatOnlyPlot{i}.^2,mnu4SqCommonPlot{i},...
%         LineStyles{i},'LineWidth',2.5,'Color',Colors{i});
%     set(gca,'YScale','log');
%     set(gca,'XScale','lin');
%     xlabel(sprintf('\\sigma^2_{syst.}/\\sigma^2_{stat.}(|U_{e4}|^2)'));
%     ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
% 
%     if strcmp(PlotStatDom,'ON')
%         legStr{i+1} = sprintf('%.0f eV range (%.0f%% stat. dominated)',...
%             range(i),100*StatDomFraction(i));
%     else
%         legStr{i+1} = sprintf('%.0f eV range',range(i));
%     end
% end
% PrettyFigureFormat;
% leg = legend([pref,pl{:}]',legStr{:},'EdgeColor',rgb('Silver'),'Location','southeast');
% if nContours>7
%     leg.NumColumns = 2;
% end
% % leg.Title.String = 'Lower fit boundary - 18574 eV';
% % leg.Title.FontWeight = 'normal';
% plotname1 = sprintf('%sksn1_StatOverSyst_%s.png',plotdir,DataType);
% print(gcf,plotname1,'-dpng','-r450');
% fprintf('save plot to %s \n',plotname1);

%% plot 2: mnu on x -axis, stat and syst seperately

Colors = (cool(nContours));
pl = cell(nContours,1);
GetFigure;
pref = plot(logspace(0,4,10),ones(10,1),'-','LineWidth',2,'Color',rgb('DimGray'));
hold on;
for i=1:nContours  
%     pl{i}= plot(mnu4SqCommonPlot{i},sin2T4StatOnlyPlot{i}.^2,...
%         LineStyles{i},'LineWidth',2.5,'Color',rgb(Colors{i}));
%     hold on;
%     pl{i}= plot(mnu4SqCommonPlot{i},sin2T4SystOnlyPlot{i}.^2,...
%         ':','LineWidth',2.5,'Color',rgb(Colors{i}));
 pl{i}= plot(mnu4SqCommonPlot{i},sin2T4SystOnlyPlot{i}.^2./sin2T4StatOnlyPlot{i}.^2,...
        LineStyles{i},'LineWidth',2.5,'Color',Colors(i,:));
   if range(i)==65
       pl{i}.LineWidth=3.5;
       pl{i}.LineStyle = '-';
       pl{i}.Color = rgb('Black');
   end
%     set(gca,'YScale','lin');
    set(gca,'XScale','log');
  ylabel(sprintf('\\sigma^2_{syst.}/\\sigma^2_{stat.}(|{\\itU}_{e4}|^2)'));
    xlabel(sprintf('{\\itm}_4^2 (eV^2)'));

    if strcmp(PlotStatDom,'ON')
        legStr{i+1} = sprintf('%.0f eV range (%.0f%% stat. dominated)',...
            range(i),100*StatDomFraction(i));
    else
        legStr{i+1} = sprintf('%.0f eV range',range(i));
    end
end
ylim([0 10])
xlim([1.3 7e3]);
PrettyFigureFormat;
leg = legend([pref,pl{:}]',legStr{:},'EdgeColor',rgb('Silver'),'Location','northwest');
% if nContours>7
%     leg.NumColumns = 2;
% end
% leg.Title.String = 'Lower fit boundary - 18574 eV';
% leg.Title.FontWeight = 'normal';
plotname3 = sprintf('%sksn1_StatOverSyst_%s.png',plotdir,DataType);
print(gcf,plotname3,'-dpng','-r450');
export_fig(gcf,strrep(plotname3,'.png','.pdf'));
fprintf('save plot to %s \n',plotname3);




%% sanityplot: syst only , stat only
SanityPlot = 'ON';
if strcmp(SanityPlot,'ON')
    myContour = 13;
    GetFigure;
    pStat= plot(mnu4SqCommonPlot{myContour},sin2T4StatOnlyPlot{myContour}.^2,'-','LineWidth',2.5);
    hold on;
    pCM = plot(mnu4SqCMPlot{myContour},sin2T4CMPlot{myContour}.^2,':','LineWidth',2.5);
    pSyst= plot(mnu4SqCommonPlot{myContour},sin2T4SystOnlyPlot{myContour}.^2,'-.','LineWidth',2.5);
    set(gca,'YScale','log');
    set(gca,'XScale','log');
    ylabel(sprintf('\\sigma^2(|U_{e4}|^2)'));
    xlabel(sprintf('{\\itm}_4^2 (eV^2)'));
    PrettyFigureFormat;
    leg = legend([pCM,pStat,pSyst],' Total',' Stat. only',' Syst. only',...
        'EdgeColor',rgb('Silver'),'Location','southwest');
    
    xlim([min(mnu4SqCommonPlot{myContour})-0.5 max(mnu4SqCommonPlot{myContour})*1.6]);
    t = title(sprintf('%.2f%% C.L. (%.0f parameter) - %.0f eV range',CL,Npar,range(myContour)),...
        'FontWeight','normal');
    plotname2 = sprintf('%sksn1_StatOverSyst_SanityPlot%.0feV_%s.png',plotdir,range(myContour),DataType);
    print(gcf,plotname2,'-dpng','-r450');
    export_fig(gcf,strrep(plotname2,'.png','.pdf'));
    fprintf('save plot to %s \n',plotname2);
end





