%
% Investigate KNM1 Background Rate Systematics
% T. Lasserre
% Last Modified: 19/6/2019
%

%%
% Read REAL Run-wise Fit & Plot/Study Background Distribution
%
if ~exist('Real','var')
    RunList = 'KNM1';
    Real = MultiRunAnalysis('RunList',RunList,...
        'DataType','Real',...
        'minuitOpt','min;minos','fitter','minuit',...
        'exclDataStart',14);
    qUmin = round(Real.ModelObj.qU(Real.exclDataStart));
    qUmax = round(Real.ModelObj.qU(end));
    range = round(Real.ModelObj.qU(Real.exclDataStart)-Real.ModelObj.Q_i);
    Real.GetPlotColor;
end
% Load SingleRunObj
% if isempty(Real.SingleRunObj)
%     Real.LoadSingleRunObj
% end
% Load Single Fit Results
Real.FitRunList('Recompute','ON');

%%
% Read TWIN Run-wise Fit & Plot/Study Background Distribution
%
if ~exist('Twin','var')
    RunList = 'KNM1';
    Twin = MultiRunAnalysis('RunList',RunList,...
        'DataType','Twin',...
        'minuitOpt','min;minos','fitter','minuit',...
        'exclDataStart',14);
    qUmin = round(Twin.ModelObj.qU(Twin.exclDataStart));
    qUmax = round(Twin.ModelObj.qU(end));
    range = round(Twin.ModelObj.qU(Twin.exclDataStart)-Twin.ModelObj.Q_i);
    Twin.GetPlotColor;
end
% Load SingleRunObj
% if isempty(Twin.SingleRunObj)
%     Twin.LoadSingleRunObj
% end
% Load Single Fit Results
Twin.FitRunList('Recompute','ON');

%% Run-wise Background Distribution TWIN/REAL
myMainTitle=[sprintf('KATRIN - KNM1 Background - %d Runs - [%.0f - %0.f] eV',...
    numel(Real.RunList),qUmin,qUmax)];
maintitle=myMainTitle;
savefile=sprintf('plots/Background_RateSystematics_%d_%s_%.0feVbelowE0_1.png',...
    numel(Real.RunList),Real.chi2,abs(range));
fig1 = figure('Name','KATRIN - KNM1 Background','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=18;a.FontWeight='bold';
% Real
subplot(2,1,1)
h = histogram(Real.SingleRun_FitResults.chi2Stat.B*1e3,9,'FaceColor',Real.PlotColor);
xc=h.BinEdges(1:end-1)+h.BinWidth/2;
yc=h.BinCounts;
bar(xc,yc,'FaceColor',Real.PlotColor);
PrettyFigureFormat; % pretty display
fw = ezfit(xc,yc,'pdf(b) = n*exp(-(b-m)^2/(2*s^2));m=290');
showfit(fw,'fitcolor',Real.PlotColor,'fitlinewidth',5,'dispfitlegend','off');
ylabel('Real Runs');
grid on;
subplot(2,1,2)
h = histogram(Twin.SingleRun_FitResults.chi2Stat.B*1e3,9,'FaceColor',Twin.PlotColor);
xc=h.BinEdges(1:end-1)+h.BinWidth/2;
yc=h.BinCounts;
bar(xc,yc,'FaceColor',Twin.PlotColor);
PrettyFigureFormat; % pretty display
fw = ezfit(xc,yc,'pdf(b) = n*exp(-(b-m)^2/(2*s^2));m=290');
showfit(fw,'fitcolor',Twin.PlotColor,'fitlinewidth',5,'dispfitlegend','off');
xlabel('Run-wise Background Rate (mcps)');
ylabel('Twin Runs');
export_fig(gcf,savefile,'-m3');

%% Display KNM1 background
fprintf('RunList: %s - Number of runs: %0.f\n',Real.StackFileName, numel(Real.RunList));
fprintf('Minimum Background: %0.1f mcps\n',min(Real.SingleRun_FitResults.chi2Stat.B*1e3));
fprintf('Maximum Background: %0.1f mcps\n',max(Real.SingleRun_FitResults.chi2Stat.B*1e3));
fprintf('Mean Background: %0.1f mcps\n',mean(Real.SingleRun_FitResults.chi2Stat.B*1e3));
fprintf('Background Standard Deviation: %0.1f mcps\n',std(Real.SingleRun_FitResults.chi2Stat.B*1e3));
fprintf('Mean Background from Gaussian Fit: %0.1f mcps\n',fw.m(1));
fprintf('Standard Deviation Background from Gaussian Fit: %0.1f mcps\n',abs(fw.m(3)));

%% Stacked Background Fit: Asimov
savefile=sprintf('plots/Background_RateSystematics_%d_%s_%.0feVbelowE0_2.png',...
    numel(Real.RunList),Real.chi2,abs(range));
Twin.fixPar='5 6 7 8 9 10 11';
Twin.Fit; Twin.PlotFit;
fprintf('Background from Stacked Fit: %0.2f mcps\n',1e3*(Twin.FitResult.par(3)+Twin.ModelObj.BKG_RateSec_i));
fprintf('Background Error from Stacked Fit: %0.2f mcps\n',1e3*Twin.FitResult.err(3));
export_fig(gcf,savefile,'-m3');

%% Stacked Background Fit: Random Stat Fluct
savefile=sprintf('plots/Background_RateSystematics_%d_%s_%.0feVbelowE0_3.png',...
    numel(Real.RunList),Real.chi2,abs(range));
Twin.fixPar='5 6 7 8 9 10 11';
Twin_TBDIS_Asimov    = Twin.RunData.TBDIS;
Twin_TBDISE_Asimov   = Twin.RunData.TBDISE;
Twin.RunData.TBDIS   = Twin_TBDIS_Asimov + randn(numel(Twin.RunData.qU),1).*Twin_TBDISE_Asimov;
Twin.RunData.TBDISE  = sqrt(Twin.RunData.TBDIS);
Twin.Fit; Twin.PlotFit;
Twin.RunData.TBDIS   = Twin_TBDIS_Asimov;
Twin.RunData.TBDISE  = Twin_TBDISE_Asimov;
fprintf('Background from Stacked Fit: %0.2f mcps\n',1e3*(Twin.FitResult.par(3)+Twin.ModelObj.BKG_RateSec_i));
fprintf('Background Error from Stacked Fit: %0.2f mcps\n',1e3*Twin.FitResult.err(3));
export_fig(gcf,savefile,'-m3');
