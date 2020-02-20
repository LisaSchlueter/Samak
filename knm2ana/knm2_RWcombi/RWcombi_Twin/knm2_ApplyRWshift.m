% test strategy to account for different endpoints using FSD
% fist step: known shift from MC

%% settings
E0Mode ='3Period'; % runwise, 3Period
RunList = 'KNM2_Prompt';
FSDFlag = 'BlindingKNM2';
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
chi2 = 'chi2Stat';
DataType = 'Twin';
range = 40;
Twin_SameqUFlag = 'OFF';

CommonArg = {'FSDFlag',FSDFlag,...
    'ELossFlag',ELossFlag,...
    'AnaFlag',AnaFlag,...
    'chi2',chi2};

RecomputeFlag = 'ON';
%% load if possible
savedir = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];
savename = [savedir,sprintf('knm2_ApplyRWshift_%s_%s_%.0feV.mat',E0Mode,RunList,range)];

if strcmp(Twin_SameqUFlag,'ON')
    savename = strrep(savename,'.mat','SameqU.mat');
end

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
else
    %% define twin E0 bias
    switch E0Mode
        case 'Runwise'
            TwinBias_Q = 'Fit';
            % get run-wise endpoint
            D = MultiRunAnalysis(CommonArg{:},'RunList',RunList,'DataType','Real','fixPar','E0 Bkg Norm');
            D.exclDataStart = D.GetexclDataStart(range);
            D.PlotFitRunList('DispHist','Separate','DisplayStyle','Rel','SavePlot','ON')
            CommonArg = {CommonArg{:},'exclDataStart',D.exclDataStart};
            D.FitRunList;
            E0 = D.SingleRun_FitResults.chi2Stat.E0;
            E0Err = D.SingleRun_FitResults.chi2Stat.E0Err;
            % broaden FSDs by standard deviation of endpoint distribution
            FSDArg = {'Sigma',std(E0)};
            % set up model and calculate twins
            M = MultiRunAnalysis(CommonArg{:},'RunList',RunList,...
                'DataType','Twin','TwinBias_Q',TwinBias_Q,'Twin_SameqUFlag',Twin_SameqUFlag,...
                'fixPar','mNu E0 Bkg Norm');
            
        case '3Period' %3 rear wall periods
            % read data and set up model
            savename3P = [savedir,'knm2_FitE0_3Periods.mat'];
            if exist(savename3P,'file') && strcmp(RecomputeFlag,'OFF')
                load(savename3P)
            else
                D1 = MultiRunAnalysis(CommonArg{:},'RunList','KNM2_RW1','DataType','Real','fixPar','E0 Bkg Norm');
                D1.exclDataStart = D1.GetexclDataStart(range);
                CommonArg = {CommonArg{:},'exclDataStart',D1.exclDataStart};
                D1.Fit;
                D2 = MultiRunAnalysis(CommonArg{:},'RunList','KNM2_RW2','DataType','Real','fixPar','E0 Bkg Norm');
                D2.Fit;
                D3 = MultiRunAnalysis(CommonArg{:},'RunList','KNM2_RW3','DataType','Real','fixPar','E0 Bkg Norm');
                D3.Fit;
                
                TwinBias_Q = [(D1.FitResult.par(2)+D1.ModelObj.Q_i).*ones(1,D1.nRuns),...
                    (D2.FitResult.par(2)+D2.ModelObj.Q_i).*ones(1,D2.nRuns),...
                    (D3.FitResult.par(2)+D3.ModelObj.Q_i).*ones(1,D3.nRuns,1)];
                E0    = [D1.FitResult.par(2),D2.FitResult.par(2),D3.FitResult.par(2)]+D1.ModelObj.Q_i;
                E0Err =[D1.FitResult.err(2),D2.FitResult.err(2),D3.FitResult.err(2)];
                meanE0 = wmean(E0,1./E0Err.^2);
                
                MultiWeights = [D1.nRuns,D2.nRuns,D3.nRuns]./sum([D1.nRuns,D2.nRuns,D3.nRuns]);
                MultiPos  = E0-meanE0; %shift
                
                save(savename3P,'TwinBias_Q','E0','E0Err','meanE0','MultiWeights','MultiPos','D1','D2','D3','CommonArg');
            end
            % use 3 gaussians shifted by difference of real wall potential
            % weighted with number of runs in real wall setting
            % gaussian width from subrun-wise delta qU
            %set up model and calculate twins
            M = MultiRunAnalysis(CommonArg{:},'RunList',RunList,'DataType','Twin','TwinBias_Q',TwinBias_Q,'fixPar','mNu E0 Bkg Norm');
            
            FSDArg = {'Sigma',mean(std(M.SingleRunData.qU')),'MultiPos',MultiPos,'MultiWeights',MultiWeights};
    end
    
    
    %% unbroadened FSD
    M.Fit
    FitResult_ref = M.FitResult;
    %% boraden FSD
    M.ModelObj.LoadFSD(FSDArg{:});
    M.ModelObj.ComputeTBDDS;
    M.ModelObj.ComputeTBDIS;
    %% fit
    M.Fit;
    FitResult = M.FitResult;
    %% save
    save(savename,'FitResult','FitResult_ref','FSDArg','M','TwinBias_Q','CommonArg','RunList');
    
    if strcmp(E0Mode,'Runwise')
          save(savename,'E0','E0Err','D','-append');
    else
    end
end

%% result
switch E0Mode
    case '3Period'
        meanQ = mean(TwinBias_Q);
    case 'Runwise'
        meanQ = mean(E0);
end
fprintf('-----------Regular FSD:-------------------------\n')
fprintf('Neutrino mass squared shift %.2f meV^2 \n',1e3*FitResult_ref.par(1));
fprintf('Endpoint shift from MC truth %.2f meV \n',1e3*(meanQ-M.ModelObj.Q_i-FitResult_ref.par(2)))

fprintf('-----------Broadened FSD:-------------------------\n')
fprintf('Neutrino mass squared shift %.2f meV^2 \n',1e3*FitResult.par(1)); 
fprintf('Endpoint shift from MC truth %.2f meV \n',1e3*(meanQ-M.ModelObj.Q_i-FitResult.par(2)))
%% plots
%M.ModelObj.LoadFSD(FSDArg{:},'BinningFactor',1,'SanityPlot','ON')
%M.ModelObj.ComputeTBDDS;
%M.ModelObj.ComputeTBDIS;
%M.Fit;

