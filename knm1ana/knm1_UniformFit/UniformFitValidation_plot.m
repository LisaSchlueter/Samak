% load file from "UniformFitValidation" and plot results
% calculate pixel-wise response function
% fit (average) MC twins with pixel-wise response function

nPix = 117;
RunList = 'KNM1';
savedir = [getenv('SamakPath'),'knm1ana/knm1_UniformFit/results/'];
savenameRF = [savedir,sprintf('ResponseFunction_Pixelwise%.0f_%s.mat',nPix,RunList)];
savenameFit = [savedir,sprintf('FitWithResponseFunction_Pixelwise%.0f_%s.mat',nPix,RunList)];

dRF = importdata(savenameRF); % pixel-wise response function
dFit = importdata(savenameFit); % pixel-wise response function

%%  plot Fit result histogram
GetFigure;
Parameter = 1;
hpixel = histogram(dFit.FitPar(:,Parameter)-dFit.FitParStack(Parameter),'FaceColor',rgb('SkyBlue'),'EdgeColor','none','FaceAlpha',1); 
hpixel.BinWidth = 0.0005;
hold on;
pmean = plot(mean(dFit.FitPar(:,Parameter)-dFit.FitParStack(Parameter)).*ones(2,1),[0,30],'-','Color',rgb('Orange'),'LineWidth',3);
PrettyFigureFormat('FontSize',24);
xlabel(sprintf('\\Delta{\\itm}_\\nu^2 (eV^{ 2})'));
ylabel('Occurence');
leg = legend([hpixel, pmean],sprintf('Fits with pixel-wise response function'),...
    sprintf('\\langleFits with pixel-wise response function\\rangle'));  
leg.Location = 'northwest';
PrettyLegendFormat(leg);
hold off
ylim([0 max(hpixel.BinCounts)+2.2]);

plotdir = strrep(savedir,'results','plots');
plotname = [plotdir,sprintf('FitResult_RFPixelwise%.0f_%s.pdf',nPix,RunList)];
export_fig(plotname);

fprintf('Standard deviation = %.1e eV^2 \n',std(dFit.FitPar(:,Parameter)));

%% plot response function

GetFigure;
Residuals = (mean(dRF.plotRF));%-RFStack(:,27)')./(RFStack(:,27)');
Residuals(isnan(Residuals)) = 0;  Residuals(isinf(Residuals)) = 0;
Residuals(Residuals>1e5) = 0;
[l,a] = boundedline(dRF.Te-dRF.qU(27),Residuals,std(dRF.plotRF));
l.Color = rgb('Orange'); l.LineStyle = '-'; l.LineWidth = 2;
a.FaceAlpha = 1; a.FaceColor = rgb('SkyBlue');
hold on;
pstack2 = plot(dRF.Te-dRF.qU(27),dRF.RFStack(:,27),':','Color',rgb('Red'),'LineWidth',2);
xlabel(sprintf('{\\itE}_{kin} - \\langle{\\itqU}\\rangle (eV)'));
ylabel('Transmission probability');
PrettyFigureFormat('FontSize',24);
leg2 = legend([a,l,pstack2],sprintf('Pixel-wise: 1\\sigma band'),....
    sprintf('\\langlePixel-wise \\rangle'),...
    sprintf('Average pixel')); 
leg2.Title.String = 'Response function'; leg2.Title.FontWeight = 'normal'; 
PrettyLegendFormat(leg2);
leg2.FontSize = get(gca,'FontSize');
leg2.Location='northwest';

xlim([-0.15,3.3]);
ylim([0 0.82]);

plotname2 = [plotdir,sprintf('ResponseFunction_RFPixelwise%.0f_%s.pdf',nPix,RunList)];
export_fig(plotname2);