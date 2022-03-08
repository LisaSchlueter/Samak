% twins for the three periods seperatetly

range = 40;

NP = 1.112;
savedir = [getenv('SamakPath'),'knm2ana/knm2_TimeEvolution/results/'];
savefile = sprintf('%sknm2_RunwiseFits_%.0feV_NP%.3g.mat',savedir,range,NP);

if exist(savefile,'file')
    load(savefile)
    fprintf('load file %s \n',savefile);
else
    MakeDir(savedir);
    BKG_PtSlope = 3*1e-06;
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'fixPar','E0 Bkg Norm',...
        'DataType','Real',...
        'FSDFlag','KNM2_0p1eV',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'BKG_PtSlope',3*1e-06,...
        'chi2','chi2Stat',...
        'DopplerEffectFlag','FSD',...
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',NP,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'BKG_PtSlope',BKG_PtSlope};
    
    % read data and set up model
    M = MultiRunAnalysis(RunAnaArg{:});
    M.exclDataStart = M.GetexclDataStart(range);
    %% fit all runs one by one (or load fit results from file)
    FitResults = M.FitRunList;
    save(savefile,'M','FitResults');
    fprintf('save file %s \n',savefile);
end

%%  plots
PltStyle = {'DispHist','ON','HideGaps','OFF','TimeStyle','date','saveplot','pdf'};
%M.PlotFitRunList('Parameter','B',PltStyle{:},'LinFit','ON');
%M.PlotFitRunList('Parameter','E0','ShowRWPeriods','ON',PltStyle{:},'YLim',18573+[-0.05 1.7]);
%M.PlotFitRunList('Parameter','N','ShowRWPeriods','OFF',PltStyle{:},'DisplayStyle','Rel','YLim',[-0.07 0.095]);
M.PlotFitRunList('Parameter','pVal','ShowRWPeriods','OFF',PltStyle{:},'YLim',[-0.05 1.05]);

%%
M.PlotFitRunList('DisplayStyle','Abs','Parameter','RhoD','HideGaps','OFF','YLim',[4.20,4.25].*1e17)
M.PlotFitRunList('DisplayStyle','Abs','Parameter','T2','HideGaps','OFF','ShowRWPeriods','ON','Saveplot','pdf');
M.PlotFitRunList('saveplot','pdf','DisplayStyle','Rel','Parameter','E0','YLim',[-0.65,0.8],'DispHist','ON');
M.PlotFitRunList('saveplot','pdf','DisplayStyle','Abs','Parameter','B','YLim','','DispHist','ON');
M.PlotFitRunList('saveplot','pdf','DisplayStyle','Abs','Parameter','N','YLim',[0.89,1.05],'DispHist','ON');

%% evaluate E0 stability
ndof = numel(FitResults.E0)-1;
chi2E0 = (wmean(FitResults.E0,1./FitResults.E0Err.^2)-FitResults.E0).^2./FitResults.E0Err.^2;
pE0 = 1-chi2cdf(sum(chi2E0),ndof);
[hE0, pE0r] =runstest(FitResults.E0); % runs test test for randomness, h = 0 -> accept null hypothesis, that x is in random order
% E0 rw period-wise
Idx1 = 1:171;   % RW1: 1-171
Idx2 = 172:268; % RW2: 172-268
Idx3 = 269:361;% RW3: 267:361

chi2E01 = (wmean(FitResults.E0(Idx1),1./FitResults.E0Err(Idx1).^2)-FitResults.E0(Idx1)).^2./FitResults.E0Err(Idx1).^2;
pE01 = 1-chi2cdf(sum(chi2E01),numel(Idx1)-1);

chi2E02 = (wmean(FitResults.E0(Idx2),1./FitResults.E0Err(Idx2).^2)-FitResults.E0(Idx2)).^2./FitResults.E0Err(Idx2).^2;
pE02 = 1-chi2cdf(sum(chi2E02),numel(Idx3)-1);

chi2E03 = (wmean(FitResults.E0(Idx3),1./FitResults.E0Err(Idx3).^2)-FitResults.E0(Idx3)).^2./FitResults.E0Err(Idx3).^2;
pE03 = 1-chi2cdf(sum(chi2E03),numel(Idx3)-1);


[hE01, pE0r1] =runstest(FitResults.E0(Idx1)); % runs test test for randomness, h = 0 -> accept null hypothesis, that x is in random order
[hE02, pE0r2] =runstest(FitResults.E0(Idx2)); % runs test test for randomness, h = 0 -> accept null hypothesis, that x is in random order
[hE03, pE0r3] =runstest(FitResults.E0(Idx3)); % runs test test for randomness, h = 0 -> accept null hypothesis, that x is in random order


fprintf('E0 stability: \n');
fprintf('All scans: chi2/cdof = %.1f/%.0f, p =%.4g, runs test p = %.2f \n',sum(chi2E0),ndof,pE0,pE0r);
fprintf('Rw1 scans: chi2/cdof = %.1f/%.0f, p =%.4g,  runs test p = %.2f \n',sum(chi2E0),numel(Idx1)-1,pE01,pE0r1);
fprintf('Rw2 scans: chi2/cdof = %.1f/%.0f,  p =%.4g, runs test p = %.2f \n',sum(chi2E0),numel(Idx2)-1,pE02,pE0r2);
fprintf('Rw3 scans: chi2/cdof = %.1f/%.0f,  p =%.4g, runs test p = %.2f \n',sum(chi2E0),numel(Idx3)-1,pE03,pE0r3);


%% normalization
chi2N = (wmean(FitResults.N,1./FitResults.NErr.^2)-FitResults.N).^2./FitResults.NErr.^2;
pN = 1-chi2cdf(sum(chi2N),ndof);
[hN, pNr] =runstest(FitResults.N); % runs test test for randomness

%% background stability
chi2B = (wmean(FitResults.B,1./FitResults.BErr.^2)-FitResults.B).^2./FitResults.BErr.^2;
pB = 1-chi2cdf(sum(chi2B),ndof);
[hB, pBr] =runstest(FitResults.B); % runs test test for randomness

LiveTime = hours(M.SingleRunData.StartTimeStamp-M.SingleRunData.StartTimeStamp(1));        
[coeff,coeffErr,chi2,dof]= linFit(LiveTime',1e3.*FitResults.B,1e3.*FitResults.BErr);
pBlin = 1-chi2cdf(chi2,dof);

%% runstest with time-slope corrected counts
RateBcorr = 1e3.*FitResults.B'-(coeff(1).*LiveTime+coeff(2)-1e3.*wmean(FitResults.B,1./FitResults.BErr.^2));
[hBcorr, pBrcorr] =runstest(RateBcorr); % runs test test for randomness


%% pvalue: ks test
cdf_u = makedist('Uniform');
[a,b]=kstest(FitResults.pValue,'CDF',cdf_u);
[hpVal, ppVal] =runstest(FitResults.pValue); % runs test test for randomness


