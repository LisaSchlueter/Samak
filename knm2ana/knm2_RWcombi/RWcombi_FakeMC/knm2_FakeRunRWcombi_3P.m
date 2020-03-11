% combining periods with 3 different RW settings
% based on Fake MC
% January 2020, Lisa

% reference file defines most KATRIN parameters
InitFile = @ref_FakeRun_KNM2_CD84_2hours; %84 column density, 2 hours

% 3 different rear wall settings
nRuns    = [121,95,92];
FakeE0   = 18573.7+[0.02,0.07,-0.1];
%FakeE0   = 18573.7+[0.4,0.6,-1];
%FakeE0   = 18573.7+[1,2,-3];
meanE0  = wmean(FakeE0,nRuns);
RecomputeFlag = 'OFF';
Plots = 'OFF';

%% load or calculate
savedir  = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];
savename = [savedir,sprintf('knm2_FakeRunRWcombi_3P_%.0fRuns_E0%.0fmeV_%.0fmeV_%.0fmeV.mat',...
    sum(nRuns),1e3*(FakeE0(1)-18573.7),1e3*(FakeE0(2)-18573.7),1e3*(FakeE0(3)-18573.7))];

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
    TBDIS1 = CalculateFakeRun(CommonArg,InitFile,FakeE0(1),nRuns(1));
    TBDIS2 = CalculateFakeRun(CommonArg,InitFile,FakeE0(2),nRuns(2));
    TBDIS3 = CalculateFakeRun(CommonArg,InitFile,FakeE0(3),nRuns(3));
    % stacked
    TBDIS = TBDIS1 + TBDIS2+ TBDIS3;
    
    %% Stacked Model: same settings, but 308 times the measurement time
    InitFileCombi =  @ref_FakeRun_KNM2_CD84_308x2hours; %84 column density, 308 runs a 2 hours
    M = RunAnalysis(CommonArg{:},'FakeInitFile',InitFileCombi,'TwinBias_Q',18573.70);
    M.RunData.TBDIS = TBDIS;
    M.Fit;
    FitResult_uncorr = M.FitResult;
    %% Stacked with correction:
    MultiPos  = FakeE0-meanE0;
    MultiWeights = [nRuns(1),nRuns(2),nRuns(3)]./sum(nRuns);
    FSDArg = {'Sigma',0.001,'MultiPos',MultiPos,'MultiWeights',MultiWeights,'Dist','Gauss'};
    M.ModelObj.LoadFSD(FSDArg{:});
    M.ModelObj.ComputeTBDDS;
    M.ModelObj.ComputeTBDIS;
    M.Fit;
    FitResult_corr = M.FitResult;
    
    save(savename,'FitResult_uncorr','FitResult_corr','FSDArg','M','InitFileCombi','InitFile',...
        'CommonArg','FakeE0','nRuns','TBDIS','TBDIS1','TBDIS2','TBDIS3');
end
%% result: 
fprintf('- Result without FSD correction------------------- \n');
fprintf('Delta mNuSq = %.2g eV2 \n',FitResult_uncorr.par(1));
fprintf('Delta E0    = %.2g eV \n',FitResult_uncorr.par(2)+M.ModelObj.Q_i-meanE0);
fprintf('------------------------------------------------- \n');
fprintf('- Result with FSD correction--------------------- \n');
fprintf('Delta mNuSq = %.2g eV2 \n',FitResult_corr.par(1));
fprintf('Delta E0    = %.2g eV \n',FitResult_corr.par(2)+M.ModelObj.Q_i-meanE0);
fprintf('------------------------------------------------- \n');

%% plot 
if strcmp(Plots,'ON')
    %  plotmodified FSD
    M.ModelObj.LoadFSD(FSDArg{:},'SanityPlot','ON','ZoomPlot','ON');
    
    % plot RW voltage distribution
    RWvoltages = [repmat(FakeE0(1)-meanE0,[1,nRuns(1)]),...
        repmat(FakeE0(2)-meanE0,[1,nRuns(2)]),...
        repmat(FakeE0(3)-meanE0,[1,nRuns(3)])];
    f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
    h1 = histogram(RWvoltages,'FaceColor',rgb('DodgerBlue'),'FaceAlpha',1);
    PrettyFigureFormat('FontSize',25);
    xlabel('Scan-wise rear wall voltage (V)');
    ylabel('Occurrence');
    export_fig(f1,strrep(strrep(savename,'results','plots'),'.mat','.pdf'));
end

%% auxillary function
function TBDIS = CalculateFakeRun(CommonArg,InitFile,E0,nRuns)
D  = RunAnalysis(CommonArg{:},'FakeInitFile',InitFile,'TwinBias_Q',E0);
TBDIS = nRuns*D.RunData.TBDIS;
end

