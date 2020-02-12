%
% Pile-up dependent efficiency
% Data by Sanshiro / Léo
%
% ROIqU: Retarding Potential
% pileupEffCorr: Pile-up Efficiency Correction
% 
% Last Updated: 26/07/2018
%

% Two runs as input data - 40806 and 40769
fileRun=load('Run40806_rate_HV_pileup.dat');
ROIqU          = fileRun(:,1);
RateFPDcps     = fileRun(:,2);
RatePixelcps   = fileRun(:,3);
pileupEffCorr  = fileRun(:,4);

% Interpolation FPD rate
pileupEffCorrFPDInterp = @(ratefpd) interp1(RateFPDcps,pileupEffCorr,ratefpd,'spline');
rateFPD = linspace(min(RateFPDcps),max(RateFPDcps),1000);

% Interpolation Pixel rate
pileupEffCorrPixelInterp = @(ratepixel) interp1(RatePixelcps,pileupEffCorr,ratepixel,'spline');
ratePixel = linspace(min(RatePixelcps),max(RatePixelcps),1000);

% Plot / FPD rate
p=Plot(RateFPDcps,pileupEffCorr,RateFPDcps,pileupEffCorrFPDInterp(RateFPDcps));
p.XLabel = 'FPD rate (cps)';
p.YLabel = 'Correction Factor';
p.Title  = 'Pile-up Efficiency Correction - FPD-wise';

% Plot / Pixel rate
p=Plot(RatePixelcps,pileupEffCorr,RatePixelcps,pileupEffCorrPixelInterp(RatePixelcps));
p.XLabel = 'Pixel rate (cps)';
p.YLabel = 'Correction Factor';
p.Title  = 'Pile-up Efficiency Correction - Pixel-wise';

% FIT - Pixel-wise
[xData, yData] = prepareCurveData( RatePixelcps, pileupEffCorr );

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'Bisquare';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData);
legend( h, 'pileupEffCorr vs. RatePixelcps', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel RatePixelcps
ylabel pileupEffCorr
grid on
PrettyFigureFormat

% 
% Plot fit with data.
figure(1)
errorbar(RatePixelcps,pileupEffCorr,pileupEffCorr*0,'-s','MarkerSize',5,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
plot(RatePixelcps,fitresult(RatePixelcps),'LineWidth',2);
hold off
xlabel('Pixel Rate (cps)');
ylabel('Correction Factor');
legend('data','1-order polynomial fit','Location','NorthEast')
title('Pile-Up - Rate Dependent Efficiency - Fit');
grid on
PrettyFigureFormat