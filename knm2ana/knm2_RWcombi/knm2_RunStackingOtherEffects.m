RunList = 'KNM2_Prompt';    % all runs
fixPar = 'mNu E0 Bkg Norm'; % free parameter
DataType = 'Twin';
FSDFlag = 'BlindingKNM2';
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel';     % uniform FPD
range = 40;                 % in eV below E0
chi2 = 'chi2Stat';

CommonArg = {'RunList',RunList,...                           % argument for RunAnalysis class, used later
    'fixPar',fixPar,...
    'DataType',DataType,...
    'FSDFlag',FSDFlag,...
    'ELossFlag',ELossFlag,...
    'AnaFlag',AnaFlag,'chi2',chi2,...
    'TwinBias_Q','Fit'};

%% twins with common endpoint, rest is runwise
A1 = MultiRunAnalysis(CommonArg{:},'TwinBias_Q',18573.7);     % constant endpoint
A1.exclDataStart = A1.GetexclDataStart(range); % get correct range
CommonArg = {CommonArg{:},'exclDataStart',A1.exclDataStart};  % add to default arguments
A1.ModelObj.LoadFSD('Sigma',0.01)
A1.Fit;

%% twins with common endpoint and common qU (subrun-wise), rest is runwise
A2 = MultiRunAnalysis(CommonArg{:},'TwinBias_Q',18573.7,'Twin_SameqUFlag','ON');     % constant endpoint
A2.Fit;

%% twins with common endpoint and common qUfrac (subrun-wise), rest is runwise
A3 = MultiRunAnalysis(CommonArg{:},'TwinBias_Q',18573.7,'Twin_SameqUfracFlag','ON');   % constant endpoint
A3.Fit;

%% twins with common endpoint and common column density (subrun-wise), rest is runwise
A4 = MultiRunAnalysis(CommonArg{:},'TwinBias_Q',18573.7,'Twin_SameCDFlag','ON');     % constant endpoint
A4.Fit;
%% result
fprintf('-----------Same endpoint:-------------------------\n')
fprintf('Neutrino mass squared shift %.2f meV^2 \n',1e3*A1.FitResult.par(1));
fprintf('Endpoint shift from MC truth %.2f meV \n',1e3*(A1.FitResult.par(2)));

fprintf('-----------Same endpoint + same qU:-------------------------\n')
fprintf('Neutrino mass squared shift %.2f meV^2 \n',1e3*A2.FitResult.par(1)); 
fprintf('Endpoint shift from MC truth %.2f meV \n',1e3*(A2.FitResult.par(2)));








