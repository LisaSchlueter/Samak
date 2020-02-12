% combining periods with drifting RW voltage
% based on Fake MC
% January 2020, Lisa

% reference file defines most KATRIN parameters
InitFile = @ref_FakeRun_KNM2_CD84_2hours; %84 column density, 2 hours

% linear drift of rear wall voltage
nRuns    = sum([121,95,92]);
<<<<<<< HEAD
maxRWdrift = 0.4;
%maxRWdrift = 1;
=======
%maxRWdrift = 0.4;
maxRWdrift = 1;
>>>>>>> 2680d1afe0066666933d7782d14d252046b5598f
RWpotential  = maxRWdrift/nRuns.*(1:nRuns);
TwinBias_Q = 18573.60+RWpotential;
meanE0  = mean(TwinBias_Q);
RecomputeFlag = 'OFF';
Plots = 'OFF';
%% load or calculate
savedir  = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];
savename = [savedir,sprintf('knm2_FakeRunRWcombi_Drift_%.0fRuns_MaxDrift_%.0fmeV.mat',...
    nRuns,1e3*maxRWdrift)];

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
else
    CommonArg = {'RunNr',1,...% has no meaning
        'DataType','Fake',...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2',...
        'exclDataStart',11,... % 11==40eV range (28 subruns)
        'chi2','chi2Stat',...
        'RingMerge','Full',...
        'minuitOpt','min;migrad',...
        'NonPoissonScaleFactor',1,...
        'AnaFlag','StackPixel',...
        'fixPar','mNu E0 Bkg Norm',...
        'RingList',1:12};
    
    %% Calculate Fake Data
    TBDIS_runwise = zeros(38,nRuns);
    for i=1:nRuns
        progressbar(i/nRuns);
       TBDIS_runwise(:,i) = CalculateFakeRun(CommonArg,InitFile,TwinBias_Q(i));  
    end
  
    % stacked
    TBDIS =sum(TBDIS_runwise,2);
    
    %% Stacked Model: same settings, but 308 times the measurement time
    InitFileCombi =  @ref_FakeRun_KNM2_CD84_308x2hours; %84 column density, 308 runs a 2 hours
    M = RunAnalysis(CommonArg{:},'FakeInitFile',InitFileCombi);
    M.RunData.TBDIS = TBDIS;
    M.Fit;
    FitResult_uncorr = M.FitResult;
    %% Stacked with correction:
  
    RecWidth = max(RWpotential)-min(RWpotential);
    FSDArg = {'Dist','Rect','Sigma',RecWidth};
    M.ModelObj.LoadFSD(FSDArg{:},'SanityPlot','OFF');
    M.ModelObj.ComputeTBDDS;
    M.ModelObj.ComputeTBDIS;
    M.Fit;
    FitResult_corr = M.FitResult;
    
    save(savename,'FitResult_uncorr','FitResult_corr','FSDArg','M','InitFileCombi','InitFile',...
        'CommonArg','TwinBias_Q','nRuns','TBDIS','TBDIS_runwise');
end
%% result: 
fprintf('- Linear drift of %.0f meV ofer %.0f runs--------------- \n',maxRWdrift*1e3,nRuns);
fprintf('- Result without FSD correction------------------- \n');
<<<<<<< HEAD
fprintf('Delta mNuSq = %.2g eV2 \n',FitResult_uncorr.par(1));
=======
fprintf('Delta mNuSq = %.3f eV2 \n',FitResult_uncorr.par(1));
>>>>>>> 2680d1afe0066666933d7782d14d252046b5598f
fprintf('Delta E0    = %.2g eV \n',FitResult_uncorr.par(2)+M.ModelObj.Q_i-meanE0);
fprintf('------------------------------------------------- \n');

fprintf('- Result with FSD correction--------------------- \n');
fprintf('Delta mNuSq = %.2g eV2 \n',FitResult_corr.par(1));
fprintf('Delta E0    = %.2g eV \n',FitResult_corr.par(2)+M.ModelObj.Q_i-meanE0);
fprintf('------------------------------------------------- \n');

%% plots
if strcmp(Plots,'ON')
    % plot modified FSD
    M.ModelObj.LoadFSD(FSDArg{:},'SanityPlot','ON','ZoomPlot','ON');
    
    % plot RW voltage evolution over time
    f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
    h1 = plot(1:nRuns,TwinBias_Q-meanE0,'Color',rgb('DodgerBlue'),'LineWidth',3);
    PrettyFigureFormat('FontSize',25);
    ylabel('Scan-wise rear wall voltage (V)');
    xlabel('Scan number');
    xlim([1,nRuns]);
    ylim([min(TwinBias_Q-meanE0), max(TwinBias_Q-meanE0)]);
    export_fig(f1,strrep(strrep(savename,'results','plots'),'.mat','.pdf'));
    
    
    % plot RW voltage distribution
    f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
    h1 = histogram(TwinBias_Q-meanE0,'FaceColor',rgb('DodgerBlue'),'LineWidth',1);
    PrettyFigureFormat('FontSize',25);
    ylabel('Occurrence');
    xlabel('Plasma potential (V)');
    xlim([min(TwinBias_Q-meanE0)-0.03, max(TwinBias_Q-meanE0)+0.03]);
    export_fig(f1,strrep(strrep(savename,'results','plots'),'.mat','_hist.pdf'));
    
end

%% auxillary function
function TBDIS = CalculateFakeRun(CommonArg,InitFile,E0)
D  = RunAnalysis(CommonArg{:},'FakeInitFile',InitFile);
D.ModelObj.Q_i = E0;
D.ModelObj.ComputeTBDDS;
D.ModelObj.ComputeTBDIS;
TBDIS = D.ModelObj.TBDIS;
end
