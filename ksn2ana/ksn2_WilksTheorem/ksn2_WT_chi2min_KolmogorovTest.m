% Test of Wilk's theorem (coverage)
% chi2min distribution (best fit chi2): Kolmogorov-Smirnov Test
Hypothesis = 'H2';
SavePlt = 'ON';
MergeNew = 'ON';
RmDuplicates = 'ON';
InterpMode = 'lin';

switch Hypothesis
    case 'H0'
        NrandMC = 419;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        randMC_new  = 1:1250;
        MergeStr = sprintf('_MergeNew%.0f',numel(randMC_new));
    case 'H1'
        randMC = [1:1500];
  Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
        chi2 = 'chi2CMShape';
        MergeNew = 'OFF'; % nothing new
        NrandMC = numel(randMC);
    case 'H2'
        randMC      = 1:1500;
        Twin_sin2T4 = 0.07;
        Twin_mNu4Sq = 20;
        chi2        = 'chi2CMShape';
        MergeNew    = 'OFF'; % nothing new
end

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_Interp%s_%.0fsamples%s_RmDouble%s.mat',...
        savedir,InterpMode,NrandMC,MergeStr,RmDuplicates);
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_Interp%s_%.0fsamples.mat',...
        savedir,Twin_mNu4Sq,Twin_sin2T4,InterpMode,NrandMC);
end


if exist(savefile,'file')
    load(savefile);
    fprintf('load file from %s \n',savefile);
else
     fprintf('file does not exist: %s \n',savefile);
     return
end


%% calculate KS-Test for best fit
dof = 23;
chi2min    = sort(chi2_bf);
Chi2CDFEmp = arrayfun(@(x) sum(chi2min<=x)./numel(chi2min),chi2min); % empirical cdf


Chi2CDFTheo = chi2cdf(chi2min,dof);                                  % theoretical cdf
[h,p,ksstat,cv] = kstest(chi2min,'CDF',[chi2min,Chi2CDFTheo]);

% plot cdf
GetFigure;
pEmp =plot(chi2min,Chi2CDFEmp,'-.','LineWidth',2);
hold on;
x = linspace(0,max(chi2_bf),1e3);
pTheo = plot(x,chi2cdf(x,dof),'-','LineWidth',2);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\chi^2_{min} (%.0f dof)',dof));
ylabel(sprintf('Cumulative probability'));
resultsStr = sprintf('KS test: p-value = %.2g',p);
title(resultsStr,'FontWeight','normal','FontSize',get(gca,'FontSize'))
legend([pTheo,pEmp],sprintf(' \\chi^2 distribution with %.0f dof',dof),...
    sprintf(' Empirical distribution (%.0f samples)',numel(chi2min)),...
    'EdgeColor',rgb('Silver'),'Location','southeast');
ylim([-0.05 1.05])
xlim([0 max(x)])
%% save
if strcmp(SavePlt,'ON')
    plotnameChi2KS = strrep(strrep(savefile,'results','plots'),'.mat','_KStest_Chi2BfCDF.png');
    print(gcf,plotnameChi2KS,'-dpng','-r450');
    fprintf('save plot to %s \n',plotnameChi2KS);
end