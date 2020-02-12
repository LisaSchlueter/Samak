R = MultiRunAnalysis('RunList','StackCD100_3hours','chi2','chi2CMShape',...
    'DataEffCor','RunSummary','exclDataStart',7);
SysEffects = struct(...
                'RF_EL','ON',...  % Response Function(RF) EnergyLoss
                'RF_BF','ON',...  % RF B-Fields
                'RF_RX','ON',...  % Column Density, inel cross ection
                'FSD','ON',...
                'TASR','ON',...
                'TCoff_RAD','ON',...
                'TCoff_OTHER','ON',...
                'DOPoff','OFF');
R.ComputeCM('DataDriven','OFF','SysEffects',SysEffects);
R.Fit;

FitResults = cell(4,1);
chi2min = zeros(4,1);
E0 = zeros(4,1);
E0err = zeros(4,1);

FitResults{1} = R.FitResult;
chi2min(1) = R.FitResult.chi2min;
E0(1)      = R.FitResult.par(2);
E0err(1)   = R.FitResult.err(2);

SysEffects = struct(...
                'RF_EL','OFF',...  % Response Function(RF) EnergyLoss
                'RF_BF','ON',...  % RF B-Fields
                'RF_RX','ON',...  % Column Density, inel cross ection
                'FSD','ON',...
                'TASR','ON',...
                'TCoff_RAD','ON',...
                'TCoff_OTHER','ON',...
                'DOPoff','OFF');
R.ComputeCM('DataDriven','OFF','SysEffects',SysEffects);
R.Fit;
FitResults{2} = R.FitResult;
chi2min(2) = R.FitResult.chi2min;
E0(2)      = R.FitResult.par(2);
E0err(2)   = R.FitResult.err(2);
%%
R = MultiRunAnalysis('RunList','StackCD100_3hours','chi2','chi2CMShape',...
    'DataEffCor','RunSummary','exclDataStart',7);
R.ComputeCM('DataDriven','OFF')
R.Fit;
FitResults{3} = R.FitResult;
chi2min(3) = R.FitResult.chi2min;
E0(3)      = R.FitResult.par(2);
E0err(3)   = R.FitResult.err(2);

%%
R = MultiRunAnalysis('RunList','StackCD100_3hours','chi2','chi2CMShape',...
    'DataEffCor','RunSummary','exclDataStart',7);
SysEffects = struct(...
                'RF_EL','OFF',...  % Response Function(RF) EnergyLoss
                'RF_BF','ON',...  % RF B-Fields
                'RF_RX','ON',...  % Column Density, inel cross ection
                'FSD','ON',...
                'TASR','ON',...
                'TCoff_RAD','ON',...
                'TCoff_OTHER','ON',...
                'DOPoff','OFF');
R.ComputeCM('DataDriven','OFF','SysEffects',SysEffects);
R.Fit;
FitResults{4} = R.FitResult;
chi2min(4) = R.FitResult.chi2min;
E0(4)      = R.FitResult.par(2);
E0err(4)   = R.FitResult.err(2);
E0 = E0+R.ModelObj.Q_i;