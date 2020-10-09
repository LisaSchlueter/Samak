DataType = 'Real';%Real';
AnaFlag = 'StackPixel';%Ring';
RingMerge = 'Full';%'None';

AllBkgPoints = 'ON';

if strcmp(AnaFlag,'Ring')
    SysBudget = 39;
    if strcmp(RingMerge,'Full')
        AnaStr = AnaFlag;
    else
        AnaStr = sprintf('Ring%s',RingMerge);
    end
else
    SysBudget = 38;
    AnaStr = AnaFlag;
end
savedir = [getenv('SamakPath'),'knm2ana/knm2_background/results/'];
savename = sprintf('%sknm2_backgroundTimeEvolution140_%s_%s.mat',...
    savedir,DataType,AnaStr);

if strcmp(AllBkgPoints,'ON')
    savename =  strrep(savename,'.mat','_allB.mat');
end

if exist(savename,'file') 
    load(savename,'FitResult','RunAnaArg','A');
else
    SigmaSq =  0.0124+0.0025;
    
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar','mNu E0 Bkg Norm',...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',38,...
        'AnaFlag',AnaFlag,...
        'chi2','chi2Stat',...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',1,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'RingMerge',RingMerge};
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(40);
    %%
    Time = A.SingleRunData.qUfrac(end,:).*A.SingleRunData.TimeSec;
    BkgPoint = A.SingleRunData.TBDIS(end,:)./Time;
    LiveTime = hours(A.SingleRunData.StartTimeStamp-A.SingleRunData.StartTimeStamp(1));
    
    Time_All   = A.SingleRunData.qUfrac.*A.SingleRunData.TimeSec;
    TBDIS_Rate =  A.SingleRunData.TBDIS./Time;
    qU = A.RunData.qU;
    
    MakeDir(savedir);
    save(savename,'Time','BkgPoint','LiveTime','Time_All','TBDIS_Rate','qU');
end

%% Extract background for 140eV subrun
figure('Units','normalized','Position',[0.1,0.1,0.6,0.4]);
if strcmp(AllBkgPoints,'ON')
    nBkg = sum(qU-18574>0)-1; % number of other bkg points
    MeanBrest  = mean(TBDIS_Rate(end-nBkg:end-1,:));
    MeanTimerest = mean(Time_All(end-nBkg:end-1,:));
%    for i=1:nBkg
% i =1;
%     eB = errorbar(LiveTime,TBDIS_Rate(end-i,:),...
%         1.112.*sqrt(TBDIS_Rate(end-i,:).*Time_All(end-i,:))./Time_All(end-i,:),...
%         '.','MarkerSize',15,'LineWidth',2,'Color',rgb('Silver'),'MarkerEdgeColor',rgb('Silver'),'CapSize',0);
    % hold on;
%     end
end

% 140 eV point
e140 = errorbar(LiveTime,BkgPoint,1.112.*sqrt(BkgPoint.*Time)./Time,...
    '.','MarkerSize',15,'LineWidth',2,'Color',rgb('Silver'),'MarkerEdgeColor',rgb('DodgerBlue'),'CapSize',0);
% linear fit to 140 eV
[par140, err140, ~,~,~] = linFit(LiveTime',BkgPoint',1.112.*(sqrt(BkgPoint.*Time)./Time)');
hold on;
% other background points
if strcmp(AllBkgPoints,'ON')
    nBkg = sum(qU-18574>0)-1; % number of other bkg points
    
    MeanBrest     = mean(TBDIS_Rate(end-nBkg:end-1,:));
    MeanTimerest  = mean(Time_All(end-nBkg:end-1,:));
    ErroMeanBrest = (1.112.*sqrt(MeanBrest.*MeanTimerest)./MeanTimerest)./sqrt(numel(LiveTime));
    eB = errorbar(LiveTime,MeanBrest,ErroMeanBrest,...
        '.','MarkerSize',15,'LineWidth',2,'Color',rgb('Orange'),'MarkerEdgeColor',rgb('Orange'),'CapSize',0);
    
    % linear fit to other point
  [parB, errB, ~,~,~] = linFit(LiveTime',MeanBrest',ErroMeanBrest');
  pFitB = plot(LiveTime,parB(1).*LiveTime+parB(2),'-','LineWidth',2,'Color',rgb('Orange'));
end

pFit = plot(LiveTime,par140(1).*LiveTime+par140(2),'-','LineWidth',2,'Color',rgb('DodgerBlue'));
PrettyFigureFormat('FontSize',20);
xlabel('Time (hours)');
ylabel('Background rate (cps)');
ylim([0.12 0.38]);
xlim([-30 max(LiveTime)+30]);
%title(sprintf('Backround subrun qU - 18574 = 140 eV'),'FontWeight','normal','FontSize',16);

if strcmp(AllBkgPoints,'ON')
 leg = legend([pFit,pFitB],sprintf('{\\itB}_{140 eV}  = %.2g \\pm %.1g   mcps/h (stat + NP)',par140(1)*1e3,err140(1)*1e3),...
     sprintf('{\\itB}_{0-45 eV} = %.2g \\pm %.1g mcps/h (stat + NP)',parB(1)*1e3,errB(1)*1e3),...
    'FontSize',16,'EdgeColor',rgb('Silver'),'Location','northwest');
else
leg = legend(pFit,sprintf('B_{140 eV} = %.2g \\pm %.1g mcps/h (stat + NP)',par140(1)*1e3,err140(1)*1e3),...
    'FontSize',16,'EdgeColor',rgb('Silver'),'Location','northwest');
end
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.7]));
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = strrep(strrep(savename,'results','plots'),'.mat','.png');
print(plotname,'-dpng','-r350');

%% plot 2: comparison with other kg points - 1 by 1
figure('Units','normalized','Position',[0.1,0.1,0.6,0.5]);
if strcmp(AllBkgPoints,'ON')
    nBkg = sum(qU-18574>0)-1; % number of other bkg points
    
    for i=1:nBkg
        subplot(nBkg,1,i)
        Bkg = TBDIS_Rate(end-i,:);
        BkgErr =   1.112.*sqrt(Bkg.*Time_All(end-i,:))./Time_All(end-i,:);
        pl = plot(linspace(0,1200,10),zeros(10,1),'-','Color',rgb('Black'),'LineWidth',1.5);
        hold on;
          plot(linspace(0,1200,10),-3.*ones(10,1),'-','Color',rgb('Silver'),'LineWidth',1);
           plot(linspace(0,1200,10),3.*ones(10,1),'-','Color',rgb('Silver'),'LineWidth',1);
         eB = plot(LiveTime,(BkgPoint-Bkg)./BkgErr,...
            '.','MarkerSize',15,'LineWidth',2,'Color',rgb('DodgerBlue'));
        PrettyFigureFormat('FontSize',22);
        ax = gca;
        mypos = ax.Position;
        ax.Position = [mypos(1) mypos(2)+0.01 mypos(3) mypos(4)+0.055];
        ylabel(sprintf('%.0f eV',qU(end-i)-18574))
        if i==nBkg
            xlabel('LiveTime')
            %ylabel(sprintf('Norm. residuals (\\sigma)'))
            t = text(-150,12,sprintf('Residuals (\\sigma)'),'Rotation',90,'FontSize',ax.XLabel.FontSize);
        else
            xticklabels('');
        end
        ylim([-5 6.5])
        yticks([-3 0 3])
        fprintf('Mean = %.1f mcps \n',1e3.*mean(Bkg-BkgPoint))
        pnone = plot(0,100,'.','Color',rgb('White'),'MarkerSize',1);
        leg = legend(pnone,sprintf('\\langle{\\itB}_{140 eV} - {\\itB}\\rangle = %.1f mcps',1e3.*mean(BkgPoint-Bkg)),...
            'Location','northeast','EdgeColor',rgb('Silver'),'FontSize',get(gca,'FontSize')-3);
        legend boxoff
        %set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.7]));
    end
    plotname = strrep(strrep(savename,'results','plots'),'.mat','_Residuals.png');
    print(plotname,'-dpng','-r350');

end


