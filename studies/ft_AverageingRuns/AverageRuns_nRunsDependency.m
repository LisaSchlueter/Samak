%-----------------------------------------------------------------------------------------------------
% Get standard error of mean as a function of combined runs
%-----------------------------------------------------------------------------------------------------
addpath(genpath('../../../Samak2.0'));
close all; clear;
CorrCoeff = 0.9;
exclDataStart = 7;
TD ='StackCD100all';%'StackCD100_3hours';
sort = 'down';
MRA = MultiRunAnalysis('RunList',TD,...
                    'fixPar','1 5 6','exclDataStart',exclDataStart,...
                     'DataEffCorr','ROI+PileUp','ringCutFlag','ex2','chi2','chi2CM');
RunListall = MRA.StackedRuns;                 
nRuns   = zeros(numel(RunListall),1);
SysCM   = cell(numel(RunListall),1);
Mean    =  zeros(numel(RunListall),1);
MeanErr = zeros(numel(RunListall),1);
chi2min = zeros(numel(RunListall),1);

MRA.FitRunList_AveragePar; close;
E0all     =  MRA.RunList_FitResults.chi2CMall.('E0');
E0Errall  =  MRA.RunList_FitResults.chi2CMall.('E0Err');
E0corr    =  MRA.RunList_FitResults.chi2CMcorr.('E0');
E0Errcorr =  MRA.RunList_FitResults.chi2CMcorr.('E0Err');
E0stat    =  MRA.RunList_FitResults.chi2Stat.('E0');
E0Errstat =  MRA.RunList_FitResults.chi2Stat.('E0Err');

for i=1:numel(RunListall)
    if strcmp(sort,'up')
        nRuns(i) = numel(RunListall)-i+1;
        MRA.RunList_FitResults.chi2CMall.('E0') = E0all(1:(end-i+1));
        MRA.RunList_FitResults.chi2CMall.('E0Err') = E0Errall(1:(end-i+1));
        MRA.RunList_FitResults.chi2CMcorr.('E0') = E0corr(1:(end-i+1));
        MRA.RunList_FitResults.chi2CMcorr.('E0Err') = E0Errcorr(1:(end-i+1));
        MRA.RunList_FitResults.chi2Stat.('E0') = E0stat(1:(end-i+1));
        MRA.RunList_FitResults.chi2Stat.('E0Err') = E0Errstat(1:(end-i+1));
    elseif strcmp(sort,'down')
        nRuns(i) = i;
        MRA.RunList_FitResults.chi2CMall.('E0') = E0all((end-i+1):end);
        MRA.RunList_FitResults.chi2CMall.('E0Err') = E0Errall((end-i+1):end);
        MRA.RunList_FitResults.chi2CMcorr.('E0') = E0corr((end-i+1):end);
        MRA.RunList_FitResults.chi2CMcorr.('E0Err') = E0Errcorr((end-i+1):end);
        MRA.RunList_FitResults.chi2Stat.('E0') = E0stat((end-i+1):end);
        MRA.RunList_FitResults.chi2Stat.('E0Err') = E0Errstat((end-i+1):end);
    end
SysCM{i} = MRA.ComputeSysCMRunList('CorrCoeff',CorrCoeff,'plotDist','OFF');


fprintf('----------BEGIN FIT MINUIT-------------- \n');
% Init Fit Parameter
if strcmp(sort,'up')
parMean = mean(E0all(1:(end-i+1)));
Data = E0all(1:(end-i+1));
elseif strcmp(sort,'down')
parMean = mean(E0all((end-i+1):end));
Data = E0all((end-i+1):end);
end
% Minuit Arguments
tmparg = sprintf(['fix %s;set pri -10;'...
    'migrad minos'],'');
% Minuit Input: Init Fit Parameter, Data, Covariance Matrix
Args = {parMean, {Data, SysCM{i}}, '-c', tmparg};
[par, err, chi2, errmat] = fminuit('chi2meanErr',Args{:});
dof = numel(E0all(1:(end-i+1)))-1;
Mean(i) = par+MRA.ModelObj.Q_i;
MeanErr(i) = err;
chi2min(i) = chi2;
end

%%
fig111 = figure(111);
set(fig111, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.7, 0.8]);
plot(nRuns,MeanErr,...
     'o--','MarkerFaceColor',rgb('CadetBlue'),'MarkerEdgeColor',rgb('SteelBlue'),'MarkerSize',10,'LineWidth',3);
ylabel('\sigma_{E0eff} (stat + sys) (eV)');
if strcmp(sort,'up')
xlabel(sprintf('number of runs (run %.0f upwards)',MRA.StackedRuns(1)));
elseif strcmp(sort,'down')
 xlabel(sprintf('number of runs (run %.0f downwards)',MRA.StackedRuns(end)));   
end
legend(['standard error of mean E0_{eff} for n runs']);
PrettyFigureFormat;
set(gca,'FontSize',18);
xlim([0.5 max(nRuns)+0.5]);
ylim([min(MeanErr) max(MeanErr)]);
title(sprintf('Standard Error of Mean for Correlated E0 Measurements (\\rho=%.1f) \n Runs: %s, %.0feV below E0',CorrCoeff,strrep(TD,'_','-'),MRA.ModelObj.Q_i-(MRA.ModelObj.qU(exclDataStart))));

savename = sprintf('AverageRuns_nRunsDependency_%s_%.0fbelowE0_%.1fCorrCoeff_%s',TD,MRA.ModelObj.Q_i-(MRA.ModelObj.qU(exclDataStart)),CorrCoeff,sort);
publish_figurePDF(fig111,['./plots/pdf/',savename,'.pdf']);
print(fig111,['./plots/png/',savename,'.png'],'-dpng');
savefig(fig111,['./plots/fig/',savename,'.fig']);
%%
%save(['./results/Result',savename,'.mat'],'Mean','MeanErr','chi2min','dof','TD','nRuns','RunListall');