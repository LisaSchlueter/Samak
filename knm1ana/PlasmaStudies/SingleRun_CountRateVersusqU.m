%
% Macro to study correlation between Rates for Sing-Run 
% and qU actual values for each subruns
% 
% Thierry Lasserre
% Last Modified, 17 May 2019
%
set(0,'DefaultFigureWindowStyle','docked')

% Read Data KNM1
ReadData = 'ON';
E0 = 18575;

switch ReadData
    case 'ON'
        option = {...
            'RunList','KNM1_m149mvRW',...
            'exclDataStart',2,...
            'fixPar','1 5 6 7 8 9 10',...
            'chi2','chi2Stat',...
            'ELossFlag','KatrinD2',...
            'StackTolerance',0.2,...
            'NonPoissonScaleFactor',1.5,...
            'Debug','ON',...
            'AnaFlag','StackPixel',...
            'PixList',[]...
            };
        SPR = MultiRunAnalysis('DataType','Real',option{:});
        SPRT= MultiRunAnalysis('DataType','Twin',option{:});
end

%% Correction Factor For Run Rates - Same For Twins & Real
rhodActivity = (SPR.SingleRunData.WGTS_MolFrac_DT/2+SPR.SingleRunData.WGTS_MolFrac_HT/2+SPR.SingleRunData.WGTS_MolFrac_TT).*SPR.SingleRunData.WGTS_CD_MolPerCm2;
RunRateCorrectionFactor = (rhodActivity./mean((SPR.SingleRunData.WGTS_MolFrac_DT/2+SPR.SingleRunData.WGTS_MolFrac_HT/2+SPR.SingleRunData.WGTS_MolFrac_TT))./mean(SPR.SingleRunData.WGTS_CD_MolPerCm2)).^-1;
            
%% Extract Rates for Single Runs
rate  = zeros(numel(SPR.RunData.qU),numel(SPR.RunList));
ratet = zeros(numel(SPRT.RunData.qU),numel(SPRT.RunList));

for i=1:1:numel(SPR.RunData.qU)
        rate(i,:) = SPR.SingleRunData.TBDIS(i,:)./SPR.SingleRunData.TimeperSubRun(i,:);
        ratet(i,:)= SPRT.SingleRunData.TBDIS(i,:)./SPRT.SingleRunData.TimeperSubRun(i,:);
end

%% Loop on qU bins
Nbins = 1:2:40;
for qUbin=Nbins
    
    % REAL: Build Reference Rates w & w/o correction
    Rateref  = rate(qUbin,:);
    RaterefC = rate(qUbin,:).*SPR.R200RateErosCorrection;
    qUref    = SPR.SingleRunData.qU(qUbin,:);
    % Fit Corrected Rates Versus qU
    [pC,SC] = polyfit(qUref-E0,RaterefC,1);
    [y_fitC,deltaC] = polyval(pC,qUref-E0,SC);
    steC = sqrt(diag(inv(SC.R)*inv(SC.R')).*SC.normr.^2./SC.df);
    % Fit Not Corrected Rates Versus qU
    [p,S] = polyfit(qUref-E0,Rateref,1);
    [y_fit,delta] = polyval(p,qUref-E0,S);
    ste = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df);
    
    % TWINS: Build Reference Rates w & w/o correction
    RateTref  = ratet(qUbin,:);
    RateTrefC = ratet(qUbin,:).*SPRT.R200RateErosCorrection;

    % Fit Corrected Rates Versus qU
    [pTC,STC] = polyfit(qUref-E0,RateTrefC,1);
    [y_fitTC,deltaTC] = polyval(pTC,qUref-E0,STC);
    steTC = sqrt(diag(inv(STC.R)*inv(STC.R')).*STC.normr.^2./STC.df);
    % Fit Not Corrected Rates Versus qU
    [pT,ST] = polyfit(qUref-E0,RateTref,1);
    [y_fitT,deltaT] = polyval(pT,qUref-E0,ST);
    steT = sqrt(diag(inv(ST.R)*inv(ST.R')).*ST.normr.^2./ST.df);
    
    
   
    % Plot Rate Verus qU
    figure
    
    % REAL
    subplot(2,4,[1 2 3])
    plot(qUref-E0,RaterefC,'s',...
        'LineWidth',2,'Color',rgb('DarkGreen'),'MarkerFaceColor',rgb('DarkGreen'),'MarkerSize',8);
    ylabel('rate (cps)');
    xlabel(sprintf('qUref-%.1f',E0));
    PrettyFigureFormat
    set(gca,'FontSize',14);
    hold on;
    plot(qUref-E0,Rateref,'s',...
        'LineWidth',2,'Color',rgb('Gray'),'MarkerFaceColor',rgb('Gray'),'MarkerSize',8);
    plot(qUref-E0,y_fit,'r-.','LineWidth',5,'Color',rgb('Gray'));
    plot(qUref-E0,y_fitC,'r-.','LineWidth',5,'Color',rgb('DarkGreen'));
    hold off
    title(sprintf('Data: C: slope=%.2f+-%.2f cps/V - NC: slope=%.2f+-%.2f cps/V',pC(1),steC(1),p(1),ste(1)));
    subplot(2,4,4)
    hnc = histogram(Rateref,'FaceColor',rgb('Gray'));
    hold on
    hc = histogram(RaterefC,'FaceColor',rgb('DarkGreen'));
    legend([hnc hc],sprintf('Raw - \\sigma=%.1f cps',std(Rateref)),sprintf('Corrected - \\sigma=%.1f cps',std(RaterefC)),'FontSize',10); legend boxoff;
    xlabel('rate (cps)');
    hold off
    PrettyFigureFormat
    set(gca,'FontSize',14);
    
    % TWIN
    subplot(2,4,[5 6 7])
    plot(qUref-E0,RateTrefC,'s',...
        'LineWidth',2,'Color',rgb('DarkOrange'),'MarkerFaceColor',rgb('DarkOrange'),'MarkerSize',8);
    ylabel('rate (cps)');
    xlabel(sprintf('qUref-%.1f',E0));
    PrettyFigureFormat
    set(gca,'FontSize',14);
    hold on;
    plot(qUref-E0,RateTref,'s',...
        'LineWidth',2,'Color',rgb('Gray'),'MarkerFaceColor',rgb('Gray'),'MarkerSize',8);
    plot(qUref-E0,y_fitT,'r-.','LineWidth',5,'Color',rgb('Gray'));
    plot(qUref-E0,y_fitTC,'r-.','LineWidth',5,'Color',rgb('DarkOrange'));
    hold off
    title(sprintf('Twin MC: C: slope=%.2f+-%.2f cps/V - NC: slope=%.2f+-%.2f cps/V',pTC(1),steTC(1),pT(1),steT(1)));
    subplot(2,4,8)
    hnc = histogram(RateTref,'FaceColor',rgb('Gray'));
    hold on
    hc = histogram(RateTrefC,'FaceColor',rgb('DarkOrange'));
    legend([hnc hc],sprintf('Raw - \\sigma=%.1f cps',std(RateTref)),sprintf('Corrected - \\sigma=%.1f cps',std(RateTrefC)),'FontSize',10); legend boxoff;
    xlabel('rate (cps)');
    hold off
    PrettyFigureFormat
    set(gca,'FontSize',14);
end

set(0,'DefaultFigureWindowStyle','normal')
