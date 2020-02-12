%
% Fit of the Neutrino Mass Squared
% For N x Twin Monte Carlo Simulations
%

Nsim = 1000;

if ~exist('Twin','var')
    RunList = 'KNM1';
    Twin = MultiRunAnalysis('RunList',RunList,'DataType','Twin','minuitOpt','min;minos','fitter','minuit','chi2','chi2CMShape','TwinBias_mnuSq',3);
    Twin.fixPar='5 6 7 8 9 10 11';
    Twin_TBDIS_Asimov   = Twin.RunData.TBDIS;
    Twin_TBDISE_Asimov  = Twin.RunData.TBDISE;
    Twin.exclDataStart  = 14;
    dof = numel(Twin.RunData.qU(Twin.exclDataStart:end)) - 4;
    range    = round(Twin.ModelObj.qU(Twin.exclDataStart)-Twin.ModelObj.Q_i);
end

%%
for countfit = 1:Nsim
    
    % Initiazlization
    Twin.ModelObj.SetFitBias(0);
    %Twin.exclDataStart  = 17;
    range    = round(Twin.ModelObj.qU(Twin.exclDataStart)-Twin.ModelObj.Q_i);
    qUmin = round(Twin.ModelObj.qU(Twin.exclDataStart));
    qUmax = round(Twin.ModelObj.qU(end));
    Twin.fixPar         = '5 6 7 8 9 10 11';
    
    %Stat Fluctuations
    %Twin.RunData.TBDIS   = Twin_TBDIS_Asimov + randn(numel(Twin.RunData.qU),1).*Twin_TBDISE_Asimov;
    
    %Stat+Syst Fluctuations
    Twin.RunData.TBDIS = mvnrnd(Twin_TBDIS_Asimov,Twin.FitCMShape,1)';
    Twin.RunData.TBDISE  = sqrt(abs(Twin.RunData.TBDIS));
    % Fit
    tic
    Twin.Fit;
    toc
    
    % Store Results
    StackFitResults(countfit) = Twin.FitResult;
    
    % Display
    fprintf('Fit Stacked %s %s:  %.0f %% completed \n',...
        Twin.DataType,Twin.StackFileName,countfit/Nsim*100);
end

%% Format Results
mass2    = cell2mat(arrayfun(@(x) x.par(1),StackFitResults,'UniformOutput',false));
e0       = cell2mat(arrayfun(@(x) x.par(2) + Twin.ModelObj.Q_i,StackFitResults,'UniformOutput',false));
B        = cell2mat(arrayfun(@(x) x.par(3) + Twin.ModelObj.BKG_RateAllFPDSec,StackFitResults,'UniformOutput',false));
N        = cell2mat(arrayfun(@(x) x.par(4),StackFitResults,'UniformOutput',false));
mass2Err = cell2mat(arrayfun(@(x) x.err(1),StackFitResults,'UniformOutput',false));
e0Err    = cell2mat(arrayfun(@(x) x.err(2) ,StackFitResults,'UniformOutput',false));
BErr     = cell2mat(arrayfun(@(x) x.err(3) + Twin.ModelObj.BKG_RateAllFPDSec,StackFitResults,'UniformOutput',false));
NErr     = cell2mat(arrayfun(@(x) x.err(4),StackFitResults,'UniformOutput',false));
chi2     = cell2mat(arrayfun(@(x) x.chi2min,StackFitResults,'UniformOutput',false));
ErrMat   = (arrayfun(@(x) x.errmat(1:4,1:4),StackFitResults,'UniformOutput',false));

% Save Results Server
%save('/home/iwsatlas1/lasserre/Samak2.0/knm1ana/knm1_ToyMCsensitivity/TwinStackFitMonteCarloMinos500_40eVrange_statsysfluctated.mat','StackFitResults');
%save('TwinStackFitMonteCarloMinos500_40eVrange_statsysfluctated.mat','StackFitResults');
%save('TwinStackFitMonteCarloMinos500_40eVrange_statsysfluctated_3ev2.mat','StackFitResults');

%% Plot Results: Neutrino Mass Squared
myMainTitle=[sprintf('KATRIN - KNM1 Uniform Stacked Toy Monte Carlo - %d Runs - [%.0f - %0.f] eV - Stat+Sys',...
    numel(Twin.RunList),qUmin,qUmax)];
maintitle=myMainTitle;
savefile=sprintf('plots/KNM1_StackedToyMC_%d_%s_%.0feVbelowE0-1.png',...
    numel(Twin.RunList),Twin.chi2,abs(range));
fig1 = figure('Name','KATRIN - KNM1 Uniform Stacked Fit','NumberTitle','off','rend','painters','pos',[10 10 1400 500]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=18;a.FontWeight='bold';
% Neutrino Mass Squared
subplot(1,2,1)
h1=histogram(mass2(mass2Err>0 & mass2Err<2.5 & chi2<50),'LineWidth',2,'normalization', 'probability','FaceAlpha',0.5,'FaceColor',Twin.PlotColor);h1.NumBins=15;
ylabel(sprintf('Probability (%0.f trials)',Nsim));
xlabel('Neutrino Mass Squared (eV^2)');
legstr = sprintf('mean = %0.2f eV^2 \n \\sigma = %0.2f eV^2',...
    mean(mass2(mass2Err>0. & mass2Err<2.5 & chi2<50)),std(mass2(mass2Err>0. & mass2Err<2.5 & chi2<50)));
leg = legend([h1],legstr,'Location','northwest');
leg.Color = 'none'; legend boxoff;
PrettyFigureFormat;set(gca,'FontSize',18);
% Neutrino MAss Squared Error
subplot(1,2,2)
h1=histogram(mass2Err(mass2Err>0.3 & mass2Err<2.5 & chi2<50),'LineWidth',2,'normalization', 'probability','FaceAlpha',0.5,'FaceColor',Twin.PlotColor);h1.NumBins=15;
ylabel(sprintf('Probability (%0.f trials)',Nsim));
xlabel('Neutrino Mass Squared 1\sigma Error (eV^2)');
legstr = sprintf('mean = %0.2f eV^2 \n \\sigma = %0.2f eV^2',...
    mean(mass2Err(mass2Err>0. & mass2Err<2.5 & chi2<50)),std(mass2Err(mass2Err>0. & mass2Err<2.5 & chi2<50)));
leg = legend([h1],legstr,'Location','northwest');
leg.Color = 'none'; legend boxoff;
PrettyFigureFormat;set(gca,'FontSize',18);
export_fig(gcf,savefile,'-m3');

%% Plot Results: Endpoint
myMainTitle=[sprintf('KATRIN - KNM1 Uniform Stacked Toy Monte Carlo - %d Runs - [%.0f - %0.f] eV - Stat+Sys',...
    numel(Twin.RunList),qUmin,qUmax)];
maintitle=myMainTitle;
savefile=sprintf('plots/KNM1_StackedToyMC_%d_%s_%.0feVbelowE0-2.png',...
    numel(Twin.RunList),Twin.chi2,abs(range));
fig1 = figure('Name','KATRIN - KNM1 Uniform Stacked Fit','NumberTitle','off','rend','painters','pos',[10 10 1400 500]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=18;a.FontWeight='bold';
% Endpoint
subplot(1,2,1)
h1=histogram(e0(mass2Err>0. & mass2Err<2.5 & chi2<50),'LineWidth',2,'normalization', 'probability','FaceAlpha',0.5,'FaceColor',Twin.PlotColor);h1.NumBins=15;
ylabel(sprintf('Probability (%0.f trials)',Nsim));
xlabel('Effective Endpoint (eV)');
legstr = sprintf('mean = %0.2f eV \n \\sigma = %0.2f eV',...
    mean(e0(mass2Err>0. & mass2Err<2.5 & chi2<50)),std(e0(mass2Err>0. & mass2Err<2.5 & chi2<50)));
leg = legend([h1],legstr,'Location','northwest');
leg.Color = 'none'; legend boxoff;
PrettyFigureFormat;set(gca,'FontSize',18);
% Endpoint Error
subplot(1,2,2)
h1=histogram(e0Err(mass2Err>0.2 & mass2Err<2.5 & chi2<50),'LineWidth',2,'normalization', 'probability','FaceAlpha',0.5,'FaceColor',Twin.PlotColor);h1.NumBins=15;
ylabel(sprintf('Probability (%0.f trials)',Nsim));
xlabel('Effective Endpoint (eV)');
legstr = sprintf('mean = %0.2f eV \n \\sigma = %0.2f eV',...
    mean(e0Err(mass2Err>0. & mass2Err<2.5 & chi2<50)),std(e0Err(mass2Err>0. & mass2Err<2.5 & chi2<50)));
leg = legend([h1],legstr,'Location','northwest');
leg.Color = 'none'; legend boxoff;
PrettyFigureFormat;set(gca,'FontSize',18);
export_fig(gcf,savefile,'-m3');

%% Chi2 & Correlation Matrix
myMainTitle=[sprintf('KATRIN - KNM1 Uniform Stacked Toy Monte Carlo - %d Runs - [%.0f - %0.f] eV - Stat+Sys',...
    numel(Twin.RunList),qUmin,qUmax)];
maintitle=myMainTitle;
savefile=sprintf('plots/KNM1_StackedToyMC_%d_%s_%.0feVbelowE0-3.png',...
    numel(Twin.RunList),Twin.chi2,abs(range));
fig1 = figure('Name','KATRIN - KNM1 Uniform Stacked Fit','NumberTitle','off','rend','painters','pos',[10 10 1400 500]);
a=annotation('textbox', [0 0.91 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=18;a.FontWeight='bold';
% chi2
subplot(1,2,1)
h1=histogram(chi2(mass2Err>0. & mass2Err<2.5 & chi2<50),'LineWidth',2,'normalization', 'probability','FaceAlpha',0.5,'FaceColor',Twin.PlotColor);h1.NumBins=15;
hold on
h2=histogram(chi2rnd(dof-1,1e6,1),'LineWidth',2,'normalization', 'probability','FaceAlpha',0.,'FaceColor','White','EdgeColor','Red');
h2.BinEdges=h1.BinEdges;
h1=histogram(chi2(mass2Err>0. & mass2Err<2.5 & chi2<50),'LineWidth',2,'normalization', 'probability','FaceAlpha',0.5,'FaceColor',Twin.PlotColor);h1.NumBins=15;
hold off
ylabel(sprintf('Probability (%0.f trials)',Nsim));
xlabel('\chi ^2');
legstr = sprintf('mean = %0.2f \n \\sigma = %0.2f',...
    mean(chi2(mass2Err>0. & mass2Err<2.5 & chi2<50)),std(chi2(mass2Err>0. & mass2Err<2.5 & chi2<50)));
leg = legend([h1 h2],legstr,sprintf('\\chi ^2 statistics \n %.0f dof',dof),'Location','northeast');
leg.Color = 'none'; legend boxoff;
PrettyFigureFormat;set(gca,'FontSize',18);
% correlation
subplot(1,2,2)
h1=scatter(mass2(mass2Err>0. & mass2Err<2.5 & chi2<50 ),e0(mass2Err>0. & mass2Err<2.5 & chi2<50),'s','MarkerFaceColor',Twin.PlotColor,'MarkerEdgeColor',Twin.PlotColor);
ylabel('Endpoint (eV)');
xlabel('Neutrino Mass Squared (eV^2)');
PrettyFigureFormat;set(gca,'FontSize',18);
export_fig(gcf,savefile,'-m3');

%%
% A = MultiRunAnalysis('RunList','KNM1','chi2','chi2Stat','exclDataStart',17,'fixPar','5 6 7 8 9 10 11'); %take 17 for covariance matrix normalization
% A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
% A.RunData.TBDIS = A.ModelObj.TBDIS;
% A.RunData.TBDISE = sqrt(A.ModelObj.TBDIS);
% A.Fit; 
% mNuSqStat    = A.FitResult.par(1);
% mNuSqErrStat = A.FitResult.err(1);
% sA = sqrt(mNuSqErrStat*1.64);
% fprintf('%.0feV range \n',range);
% fprintf('sensitivity mNu 90%%C.L. = %.3f (stat) \n',sA);

%% Asimov Sensitivity
Twin.RunData.TBDIS   = Twin_TBDIS_Asimov;
Twin.RunData.TBDISE  = Twin_TBDISE_Asimov;
Twin.Fit;
mNuSqStat    = Twin.FitResult.par(1);
mNuSqErrStat = Twin.FitResult.err(1);
sA = sqrt(mNuSqErrStat*1.64);
fprintf('%.0feV range \n',range);
fprintf('sensitivity mNu 90%%C.L. = %.3f (stat) \n',sA);

%% Toy MC Verus Asimov Sensitivity
myMainTitle=[sprintf('KATRIN - KNM1 Uniform Stacked Toy Monte Carlo - %d Runs - [%.0f - %0.f] eV - Stat+Sys',...
    numel(Twin.RunList),qUmin,qUmax)];
maintitle=myMainTitle;
savefile=sprintf('plots/KNM1_StackedToyMC_%d_%s_%.0feVbelowE0-4.png',...
    numel(Twin.RunList),Twin.chi2,abs(range));
fig1 = figure('Name','KATRIN - KNM1 Uniform Stacked Fit','NumberTitle','off','rend','painters','pos',[10 10 1000 500]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=18;a.FontWeight='bold';
s90=sqrt(mass2Err(mass2Err>0.15 & mass2Err<3.5  & chi2<50)*1.64);
h=histogram(s90,'LineWidth',2,'normalization', 'probability','FaceAlpha',0.5,'FaceColor',Twin.PlotColor);h.NumBins=15;
hold on
sAsimov=line([sA sA],[0 max(h(1).Values)],'Color',rgb('Red'),'LineWidth',2);
s2sigma=line([mean(s90)-2*std(s90) mean(s90)-2*std(s90)],[0 max(h(1).Values)/2],'Color',rgb('Orange'),'LineWidth',2);
s2sigma=line([mean(s90)+2*std(s90) mean(s90)+2*std(s90)],[0 max(h(1).Values)/2],'Color',rgb('Orange'),'LineWidth',2);
hold off
ylabel(sprintf('Probability (%0.f trials)',Nsim));
xlabel('m_{\nu} - 90%CL sensitivity (eV)');
legstr = sprintf('mean = %0.2f eV\n \\sigma = %0.2f eV',...
    mean(s90),std(s90));
leg = legend([h s2sigma sAsimov ],legstr,sprintf('2\\sigma interval \n [%.1f %.1f] eV',...
    mean(s90)-2*std(s90),mean(s90)+2*std(s90)),sprintf('Asimov\n %0.2f eV',sA),'Location','northwest');
leg.Color = 'none'; legend boxoff;
PrettyFigureFormat;set(gca,'FontSize',18);
export_fig(gcf,savefile,'-m3');
