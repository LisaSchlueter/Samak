%
% Investigate KNM2 Non-Poissonian Component of Background
% T. Lasserre
% Last Modified: 14/10/2019
%

%%
% Read REAL Run-wise Fit & Plot/Study Background Distribution
%
if ~exist('Real','var')
    RunList = 'KNM2_Prompt';
    Real = MultiRunAnalysis('RunList',RunList,...
        'DataType','Real',...
        'minuitOpt','min;minos','fitter','minuit',...
        'exclDataStart',12,...
        'DopplerEffectFlag','OFF',...
        'RadiativeFlag','ON',...
        'RingMerge','Full');
    qUmin = round(Real.ModelObj.qU(Real.exclDataStart));
    qUmax = round(Real.ModelObj.qU(end));
    range = round(Real.ModelObj.qU(Real.exclDataStart)-Real.ModelObj.Q_i);
end
%% Load SingleRunObj
if isempty(Real.SingleRunObj)
    Real.LoadSingleRunObj
end
% Load Single Fit Results
%Real.FitRunList('Recompute','ON');

%%
% For each Run:
% Select all points (subruns) above E0 (5 points per scan)
% Build Main Distributions
% - Counts
% - Count Errors
% - Time
% -  Rate
%
qUrange = Real.RunData.qU>18573.7;
BkgSubRun      = squeeze(Real.SingleRunData.TBDIS(qUrange,:)); BkgSubRun=reshape(BkgSubRun,1,numel(BkgSubRun));
BkgSubRunError = squeeze(Real.SingleRunData.TBDISE(qUrange,:)); BkgSubRunError=reshape(BkgSubRunError,1,numel(BkgSubRunError));
BkgSubRunTime  = squeeze(Real.SingleRunData.TimeSec.*Real.SingleRunData.qUfrac(qUrange,:)); BkgSubRunTime=reshape(BkgSubRunTime,1,numel(BkgSubRunTime));

Threshold      = 1e4;
BkgSubRun      = BkgSubRun(BkgSubRun<Threshold);
BkgSubRunError = BkgSubRunError(BkgSubRun<Threshold);
BkgSubRunTime  = BkgSubRunTime(BkgSubRun<Threshold);

BkgSubRunRate  = BkgSubRun./BkgSubRunTime*1e3;

% Plot SubRun Background Distributions
myMainTitle = sprintf('KATRIN - KNM2 Subrun-wise Background - %d Runs - [%.0f - %0.f] eV',...
    numel(Real.RunList),qUmin,qUmax);
maintitle   = myMainTitle;
savefile    = sprintf('plots/BackgroundSubRunWise_%d_%s_%.0feVbelowE0_1.png',...
    numel(Real.RunList),Real.chi2,abs(range));
fig1      = figure('Name','KATRIN - KNM2 Subrun-wise Background','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=18;a.FontWeight='bold';
subplot(2,2,1)
histogram(BkgSubRun,'FaceColor',Real.PlotColor,'Normalization','probability','LineWidth',2);
xlabel('Counts Per Subrun (mcps)'); ylabel('Frequency');
PrettyFigureFormat
subplot(2,2,2)
histogram(BkgSubRunError,'FaceColor',Real.PlotColor,'Normalization','probability','LineWidth',2);
xlabel('Count errors Per Subrun (mcps)'); ylabel('Frequency');
PrettyFigureFormat
subplot(2,2,3)
histogram(BkgSubRunTime,'FaceColor',Real.PlotColor,'Normalization','probability','LineWidth',2);
xlabel('Time Per Subrun (second)'); ylabel('Frequency');
PrettyFigureFormat
subplot(2,2,4)
histogram(BkgSubRunRate,'FaceColor',Real.PlotColor,'Normalization','probability','LineWidth',2);
xlabel('Rate Per Subrun (mcps)'); ylabel('Frequency');
PrettyFigureFormat
export_fig(gcf,savefile,'-m3');

%% Display
t = PrintTable(sprintf('%s RunList Subrun-wise Background Measurements',Real.StackFileName));
t.addRow('minimum value',min(BkgSubRun),'counts');
t.addRow('average value',mean(BkgSubRun),'counts');
t.addRow('maximum value',max(BkgSubRun),'counts');
t.addRow('','','');
t.addRow('minimum uncertainty on counts',min(BkgSubRunError),'counts');
t.addRow('average uncertainty on counts',mean(BkgSubRunError),'counts');
t.addRow('maximum uncertainty on counts',max(BkgSubRunError),'counts');
t.addRow('','','');
t.addRow('minimum subrun time',min(BkgSubRunTime),'seconds');
t.addRow('average subrun time',mean(BkgSubRunTime),'seconds');
t.addRow('maximum subrun time',max(BkgSubRunTime),'seconds');
t.addRow('','','');
t.addRow('','','');
t.addRow('minimum rate  ',min(BkgSubRunRate),'mcps');
t.addRow('average rate  ',mean(BkgSubRunRate),'mcps');
t.addRow('maximum rate  ',max(BkgSubRunRate),'mcps');
t.addRow('','','');
t.display;
%t.HasHeader = true;
t.Format = 'tex';
t.Caption =sprintf('%s RunList Subrun-wise Background Measurements',Real.StackFileName);
t.print;

%%
% Build a table of Counts
% All Counts are renormalized to an averaged subrun time
myMainTitle = sprintf('KATRIN - KNM2 Subrun-wise Background - %d Runs - [%.0f - %0.f] eV',...
    numel(Real.RunList),qUmin,qUmax);
maintitle   = myMainTitle;
savefile    = sprintf('plots/BackgroundSubRunWise_%d_%s_%.0feVbelowE0_2.png',...
    numel(Real.RunList),Real.chi2,abs(range));
fig1      = figure('Name','KATRIN - KNM2 Subrun-wise Background','NumberTitle','off','rend','painters','pos',[10 10 1000 1000]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=18;a.FontWeight='bold';
BkgSubRunR = BkgSubRun./BkgSubRunTime.*mean(BkgSubRunTime);
% Real
subplot(2,1,1)
h = histogram(BkgSubRunR,40,'FaceColor',Real.PlotColor,'Normalization','probability','LineWidth',2);
pdG = fitdist(BkgSubRunR','Normal');
histfit(BkgSubRunR,40,'Normal');
PrettyFigureFormat; % pretty display
ylabel('Subruns');
%ylim([0 250]);
subplot(2,1,2)
h = histogram(BkgSubRunR,40,'FaceColor',Real.PlotColor,'Normalization','probability','LineWidth',2);
pdN = fitdist(BkgSubRunR','poisson');
histfit(BkgSubRunR,40,'Poisson');
PrettyFigureFormat; % pretty display
xlabel(sprintf('Subrun-wise Background Counts (mcps) in %.1f sec',mean(BkgSubRunTime)));
ylabel('Subruns');
%ylim([0 250]);
export_fig(gcf,savefile,'-m3');

%% Nice Plot to illustrate the background Rate Systematics
myMainTitle = sprintf('KATRIN - KNM2 Subrun-wise Background - %d subruns - [%.0f - %0.f] eV',...
    numel(BkgSubRunR),qUmin,qUmax);
maintitle   = myMainTitle;
savefile    = sprintf('plots/BackgroundSubRunWise_%d_%s_%.0feVbelowE0_3.png',...
    numel(Real.RunList),Real.chi2,abs(range));
fig1      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
    'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

h = histogram(BkgSubRunR,40,'Normalization','pdf',...
    'FaceColor',rgb('DodgerBlue'),'LineWidth',2,'FaceAlpha',0.7);
xlabel(sprintf('Background counts in %.0f sec',mean(BkgSubRunTime)));
ylabel('Frequency');
PrettyFigureFormat
% Gaussian PDF
pdfG = @(b) 1/sqrt(2*pi*pdG.sigma.^2) .* exp(-(b-pdG.mu).^2/(2*pdG.sigma.^2));
% Poisson PDF
pdfP = @(b) 1/sqrt(2*pi*pdN.lambda) .* exp(-(b-pdN.lambda).^2/(2*pdN.lambda));
hold on
%b = (h.BinEdges(1:end-1)+h.BinWidth/2);
b = linspace(min(BkgSubRunR),max(BkgSubRunR),100);
g=plot(b,pdfG(b),'Color',rgb('IndianRed'),'LineWidth',4);
p=plot(b,pdfP(b),'Color',rgb('GoldenRod'),'LineWidth',4);
leg=legend([h g p],...
    sprintf('Background subscans'),...
    sprintf('Gaussian \\sigma=%.2f counts',...
    (pdG.sigma)),sprintf('Poisson \\sigma=%.2f counts',sqrt(pdN.lambda)),...
    'location','northwest');
leg.Color = 'none'; legend boxoff;
hold off
PrettyFigureFormat
                        set(gca,'FontSize',24);
export_fig(gcf,savefile,'-m3');
