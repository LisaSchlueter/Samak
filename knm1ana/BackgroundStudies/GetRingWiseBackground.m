% RingWise Fit Using Ring Analysis Class
MyRingList=1:12;
ReDoConstructor='ON';
ReDoFit='ON';
RunList = 'KNM1';
switch ReDoConstructor
    case 'ON'
        M = MultiRunAnalysis('RunList',RunList);
        R = RingAnalysis('RunAnaObj',M,'RingList',MyRingList);
end

%% List of Pixels per Ring
RingCounter=0;
for i=MyRingList
    RingCounter=RingCounter+1;
    PixelPerRing(RingCounter)=numel(R.MultiObj(RingCounter).PixList);
    fprintf('------------------------------------------------------------------------------------------------------\n')
    fprintf(2,'Ring %.0f: %0.f Pixels =',i,PixelPerRing)
    disp(R.MultiObj(RingCounter).PixList)
    initRingBackground(RingCounter)=R.MultiObj(RingCounter).ModelObj.BKG_RateSec_i;
end

%% Fit Individual Rings
switch ReDoFit
    case 'ON'
        R.FitRings('SaveResult','ON','RecomputeFlag','ON');
        R.PlotFits('SavePlot','ON','PlotPar',3);
end

BkgPerRingMCPS    = (R.FitResult.par(:,3)+initRingBackground(:))*1000;
BkgErrPerRingMCPS = R.FitResult.err(:,3)*1000;

%% Ring Wise Background
fig1000 = figure('Renderer','opengl');
set(fig1000,'units','normalized','pos',[0.1, 0.1,0.9,0.9]);
subplot(3,1,[1 2])
bar(R.RingList,BkgPerRingMCPS,'facecolor',rgb('IndianRed'));
hold on
errorb(R.RingList,BkgPerRingMCPS,BkgErrPerRingMCPS,...
    rgb('IndianRed'),'LineWidth',1)
hold off
xlim([0 13]);
ylabel('mcps');
title(sprintf('KNM1 %.0f runs (%0.f - %0.f) - %.0f Pixels - %.2f mcps',...
    M.nRuns,M.RunList(1),M.RunList(end),...
    numel(M.PixList),...
    sum(BkgPerRingMCPS)));
PrettyFigureFormat      
set(gca,'FontSize',20);
subplot(3,1,3)
e = bar(R.RingList,PixelPerRing);
xlim([0 13]);
ylim([2 13]);
xlabel('Ring');
ylabel('Pixels');
PrettyFigureFormat      
set(gca,'FontSize',20);

%% Ring Wise Background Renormalized to 1 pixel
BkgPerRingPerPixelMCPS=BkgPerRingMCPS./PixelPerRing';
BkgErrPerRingPerPixelMCPS=BkgErrPerRingMCPS./PixelPerRing';
fig1000 = figure('Renderer','opengl');
set(fig1000,'units','normalized','pos',[0.1, 0.1,0.9,0.9]);
bar(R.RingList,BkgPerRingPerPixelMCPS,'facecolor',rgb('IndianRed'));
hold on
errorb(R.RingList,BkgPerRingPerPixelMCPS,BkgErrPerRingPerPixelMCPS,...
    rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'))
hold off
xlim([0 13]);
ylim([2 3.1]);
xlabel('Ring');
ylabel('Pixel equivalent background (mcps)');
title(sprintf('KNM1 %.0f runs (%0.f - %0.f) - %.0f Pixels - %.2f mcps',...
    M.nRuns,M.RunList(1),M.RunList(end),...
    numel(M.PixList),...
    sum(BkgPerRingMCPS)));
PrettyFigureFormat      
set(gca,'FontSize',20);

%% Background Per Ring simulation
%% For 12 rings
%M149.ReadSingleRunData;
MyRing=1;
SubRun              = 9 ;
AnchorBackground    = 295; % mcps
AnchorPixelsPerRing = PixelPerRing';
AnchorPixels        = sum(PixelPerRing);
TimePerSubRun       = M.SingleRunData.TimeperSubRun(SubRun);
Ntrials             = 189;

SimBkgPerRing    = AnchorBackground/sum(BkgPerRingMCPS) .* ...
                   AnchorPixels/numel(M.PixList) .* ...
                   BkgPerRingPerPixelMCPS .* AnchorPixelsPerRing .* ...
                   TimePerSubRun * 1e-3;
SimBkgPerRing     = repmat(SimBkgPerRing,12,Ntrials);
Sample = poissrnd(SimBkgPerRing);

% Plot Distribution of Events in Ring MyRing
fig1000 = figure('Renderer','opengl');
set(fig1000,'units','normalized','pos',[0.1, 0.1,0.9,0.9]);
% Sim
subplot(1,2,1)
nhist(Sample(MyRing,:));
title(sprintf('Sim: Counts Per Subrun %.1f sec in Ring%.0f',TimePerSubRun,MyRing));
PrettyFigureFormat
% Data
subplot(1,2,2)
Data=sum(M.SingleRunData.TBDIS(SubRun,:,cell2mat(M.ModelObj.FPD_RingPixList(MyRing))),3);
nhist(Data);
title(sprintf('Data: Counts Per Subrun %.1f sec in Ring%.0f',TimePerSubRun,MyRing));
PrettyFigureFormat;

