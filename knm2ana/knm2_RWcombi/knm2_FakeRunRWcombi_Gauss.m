% combining periods with Gaussian distributed RW voltage
% based on Fake MC
% January 2020, Lisa

% reference file defines most KATRIN parameters
InitFile = @ref_FakeRun_KNM2_CD84_2hours; %84 column density, 2 hours

% linear drift of rear wall voltage
nRuns    = sum([121,95,92]);
SigmaRW = 0.25; % flucutuation of RW voltage (std) in eV
TwinBias_Q = 18573.7+randn(nRuns,1).*SigmaRW;

RecomputeFlag = 'OFF';
Plots = 'ON';

%% load or calculate
savedir  = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];
savename = [savedir,sprintf('knm2_FakeRunRWcombi_Gaus_%.0fRuns_SigmaRW_%.0fmeV.mat',...
    nRuns,1e3*SigmaRW)];

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
    meanE0  = mean(TwinBias_Q);
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
    FSDArg = {'Dist','Gauss','Sigma',std(TwinBias_Q)};
    M.ModelObj.LoadFSD(FSDArg{:},'SanityPlot','OFF');
    M.ModelObj.ComputeTBDDS;
    M.ModelObj.ComputeTBDIS;
    M.Fit;
    FitResult_corr = M.FitResult;
     meanE0  = mean(TwinBias_Q);
     
    save(savename,'FitResult_uncorr','FitResult_corr','FSDArg','M','InitFileCombi','InitFile',...
        'CommonArg','SigmaRW','TwinBias_Q','nRuns','TBDIS','TBDIS_runwise');
end
%% result: 
fprintf('------------------------------------------------- \n');
fprintf('- Gaussian RW fluc. of %.0f meV ofer %.0f runs----- \n',SigmaRW*1e3,nRuns);
fprintf('- Result without FSD correction------------------ \n');
fprintf('Delta mNuSq = %.2g eV2 \n',FitResult_uncorr.par(1));
fprintf('Delta E0    = %.2g eV \n',FitResult_uncorr.par(2)+M.ModelObj.Q_i-meanE0);
fprintf('------------------------------------------------- \n');

fprintf('- Result with FSD correction--------------------- \n');
fprintf('Delta mNuSq = %.2g eV2 \n',FitResult_corr.par(1));
fprintf('Delta E0    = %.2g eV \n',FitResult_corr.par(2)+M.ModelObj.Q_i-meanE0);
fprintf('------------------------------------------------- \n');

%% sanity plots
if strcmp(Plots,'ON')
    % plot modified FSD
    M.ModelObj.LoadFSD(FSDArg{:},'SanityPlot','ON','ZoomPlot','ON');
    
    % plot RW voltage distribution
    f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
    h1 = histogram(TwinBias_Q-meanE0,'FaceColor',rgb('DodgerBlue'),'FaceAlpha',1);
    PrettyFigureFormat('FontSize',25);
    xlabel('Scan-wise rear wall voltage (V)');
    ylabel('Occurrence');
    export_fig(f1,strrep(strrep(savename,'results','plots'),'.mat','.pdf'));
end
%% auxillary function
function TBDIS = CalculateFakeRun(CommonArg,InitFile,E0)
D  = RunAnalysis(CommonArg{:},'FakeInitFile',InitFile);
D.ModelObj.Q_i = E0;
D.ModelObj.ComputeTBDDS;
D.ModelObj.ComputeTBDIS;
TBDIS = D.ModelObj.TBDIS;
end
