% Test of Wilk's theorem (coverage)
% chi2null distribution (best fit chi2): Kolmogorov-Smirnov Test
Hypothesis = 'H0';
switch Hypothesis
    case 'H0'
        NrandMC = 1e3;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
    case 'H1'
        NrandMC = 439;
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
end
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_%.0fsamples.mat',savedir,NrandMC);
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_%.0fsamples.mat',savedir,Twin_mNu4Sq,Twin_sin2T4,NrandMC);
end

if exist(savefile,'file')
    load(savefile);
    fprintf('load file from %s \n',savefile);
else
     fprintf('file does not exist: %s \n',savefile);
     return
end


%% calculate KS-Test for best fit
dof = 25;
chi2_null    = sort(chi2_null);
Chi2CDFEmp = arrayfun(@(x) sum(chi2_null<=x)./numel(chi2_null),chi2_null); % empirical cdf


Chi2CDFTheo = chi2cdf(chi2_null,dof);                                  % theoretical cdf
[h,p,ksstat,cv] = kstest(chi2_null,'CDF',[chi2_null,Chi2CDFTheo]);

% plot cdf
GetFigure;
pEmp =plot(chi2_null,Chi2CDFEmp,'-.','LineWidth',2);
hold on;
x = linspace(0,max(chi2_null),1e3);
pTheo = plot(x,chi2cdf(x,dof),'-','LineWidth',2);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\chi^2_{null} (%.0f dof)',dof));
ylabel(sprintf('Cumulative probability'));
resultsStr = sprintf('KS test: p-value = %.2g',p);
title(resultsStr,'FontWeight','normal','FontSize',get(gca,'FontSize'))
legend([pTheo,pEmp],sprintf(' \\chi^2 distribution with %.0f dof',dof),...
    sprintf(' Empirical distribution (%.0f samples)',numel(chi2_null)),...
    'EdgeColor',rgb('Silver'),'Location','southeast');
ylim([-0.05 1.05])

%% save
plotnameChi2KS = strrep(strrep(savefile,'results','plots'),'.mat','_KStest_Chi2NullCDF.png');
print(gcf,plotnameChi2KS,'-dpng','-r450');
fprintf('save plot to %s \n',plotnameChi2KS);