load('TwinStackFitMonteCarloMinos500_40eVrange_statsysfluctated.mat');
%% Format Results
mass2    = cell2mat(arrayfun(@(x) x.par(1),StackFitResults,'UniformOutput',false));
e0       = cell2mat(arrayfun(@(x) x.par(2) + MRL.ModelObj.Q_i,StackFitResults,'UniformOutput',false));
B        = cell2mat(arrayfun(@(x) x.par(3) + MRL.ModelObj.BKG_RateAllFPDSec,StackFitResults,'UniformOutput',false));
N        = cell2mat(arrayfun(@(x) x.par(4),StackFitResults,'UniformOutput',false));
mass2Err = cell2mat(arrayfun(@(x) x.err(1),StackFitResults,'UniformOutput',false));
e0Err    = cell2mat(arrayfun(@(x) x.err(2) ,StackFitResults,'UniformOutput',false));
BErr     = cell2mat(arrayfun(@(x) x.err(3) + MRL.ModelObj.BKG_RateAllFPDSec,StackFitResults,'UniformOutput',false));
NErr     = cell2mat(arrayfun(@(x) x.err(4),StackFitResults,'UniformOutput',false));
chi2     = cell2mat(arrayfun(@(x) x.chi2min,StackFitResults,'UniformOutput',false));
ErrMat   = (arrayfun(@(x) x.errmat(1:4,1:4),StackFitResults,'UniformOutput',false));

Nsim = 1000;
range    = round(MRL.ModelObj.qU(MRL.exclDataStart)-MRL.ModelObj.Q_i);
qUmin = round(MRL.ModelObj.qU(MRL.exclDataStart));
qUmax = round(MRL.ModelObj.qU(end));
    
%% Plot Results: Neutrino Mass Squared
myMainTitle=[sprintf('KATRIN - KNM1 Uniform Stacked Toy Monte Carlo - %d Runs - [%.0f - %0.f] eV - Stat+Sys',...
    numel(MRL.RunList),qUmin,qUmax)];
maintitle=myMainTitle;
savefile=sprintf('plots/KNM1_StackedToyMC_%d_%s_%.0feVbelowE0-1.png',...
    numel(MRL.RunList),MRL.chi2,abs(range));
fig1 = figure('Name','KATRIN - KNM1 Uniform Stacked Fit','NumberTitle','off','rend','painters','pos',[10 10 1400 500]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=18;a.FontWeight='bold';
% Neutrino Mass Squared
subplot(1,2,1)
h1=histogram(mass2(mass2Err>0 & mass2Err<2.5 & chi2<50),'LineWidth',2,'normalization', 'cdf','FaceAlpha',0.5,'FaceColor',MRL.PlotColor);h1.NumBins=15;
ylabel(sprintf('Probability (%0.f trials)',Nsim));
xlabel('Neutrino Mass Squared (eV^2)');
legstr = sprintf('mean = %0.2f eV^2 \n \\sigma = %0.2f eV^2',...
    mean(mass2(mass2Err>0. & mass2Err<2.5 & chi2<50)),std(mass2(mass2Err>0. & mass2Err<2.5 & chi2<50)));
leg = legend([h1],legstr,'Location','northwest');
leg.Color = 'none'; legend boxoff;
PrettyFigureFormat;set(gca,'FontSize',18);
% Neutrino MAss Squared Error
subplot(1,2,2)
h1=histogram(mass2Err(mass2Err>0.3 & mass2Err<2.5 & chi2<50),'LineWidth',2,'normalization', 'cdf','FaceAlpha',0.5,'FaceColor',MRL.PlotColor);h1.NumBins=15;
ylabel(sprintf('Probability (%0.f trials)',Nsim));
xlabel('Neutrino Mass Squared 1\sigma Error (eV^2)');
legstr = sprintf('mean = %0.2f eV^2 \n \\sigma = %0.2f eV^2',...
    mean(mass2Err(mass2Err>0. & mass2Err<2.5 & chi2<50)),std(mass2Err(mass2Err>0. & mass2Err<2.5 & chi2<50)));
leg = legend([h1],legstr,'Location','northwest');
leg.Color = 'none'; legend boxoff;
PrettyFigureFormat;set(gca,'FontSize',18);
export_fig(gcf,savefile,'-m3');
