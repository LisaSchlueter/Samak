% Test of Wilk's theorem (coverage)
%  DeltaChi2Crit convergence test
% cumulative pdf with critical delta chi2 for 95%CL.
Hypothesis = 'H0';
InterpMode = 'lin';
SavePlt = 'ON';
MergeNew = 'ON';
RmDuplicates = 'ON';

switch Hypothesis
    case 'H0'
        randMC =[1001:1260,1294:1300,1349:1500];%11:1e3;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        chi2 = 'chi2CMShape';
        randMC_new  = 1:1251;
    case 'H1'
        randMC = 1:1e3;
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
        chi2 = 'chi2Stat';
        MergeNew = 'OFF'; % nothing new
end

if strcmp(MergeNew,'ON')
    MergeStr = sprintf('_MergeNew%.0f',numel(randMC_new));
    NrandMC = numel(randMC)+numel(randMC_new);
else
    MergeStr = '';
    NrandMC = numel(randMC);
end

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];

if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_Interp%s_%.0fsamples%s_RmDouble%s.mat',...
        savedir,InterpMode,numel(randMC),MergeStr,RmDuplicates);
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_Interp%s_%.0fsamples%s_RmDouble%s.mat',...
        savedir,Twin_mNu4Sq,Twin_sin2T4,InterpMode,numel(randMC),MergeStr,RmDuplicates);
end

if exist(savefile,'file')
    d = importdata(savefile);
    fprintf('load file from %s \n',savefile);
else
    fprintf('file does not exist: %s \n',savefile);
    return
end

savefileConv = strrep(savefile,'.mat','_Convergence.mat');
if exist(savefileConv,'file')
    dC = importdata(savefileConv);
    fprintf('load convergence file from %s \n',savefileConv);
else
    % calculate sample evolution
    nSamplesConv = 2:1:numel(d.chi2_delta);
    Chi2Crit = zeros(numel(nSamplesConv),1);
    
    for i=1:numel(nSamplesConv)
        DeltaChi2_tmp = sort(d.chi2_delta(1:nSamplesConv(i)));
        DeltaChi2CDF_tmp = arrayfun(@(x) sum(DeltaChi2_tmp<=x)./numel(DeltaChi2_tmp),DeltaChi2_tmp);
        [DeltaChi2CDFquantile,ia] = unique(DeltaChi2CDF_tmp);
        Chi2Crit(i) = interp1(DeltaChi2CDFquantile,DeltaChi2_tmp(ia),0.95,'lin');
    end
    save(savefileConv,'Chi2Crit','nSamplesConv');
    dC = importdata(savefileConv);
end

%% plot
GetFigure;
plot(dC.nSamplesConv,dC.Chi2Crit,'-','LineWidth',2,'MarkerSize',10);
PrettyFigureFormat;
xlabel('Number of samples');
ylabel(sprintf('\\Delta\\chi^2_{crit.}'));
xlim([0 max(nSamplesConv)])
plotname = strrep(strrep(savefile,'results','plots'),'.mat','_DeltaChi2Crit_Convergence.png');
print(gcf,plotname,'-dpng','-r450');
fprintf('save plot to %s \n',plotname);
 