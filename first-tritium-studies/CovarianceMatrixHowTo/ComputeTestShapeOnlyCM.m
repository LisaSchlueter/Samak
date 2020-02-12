RunList = 40667;
MRA = MultiRunAnalysis('RunList',RunList,'exclDataStart',1);

% Define Covariance Matrix properties
myEffects = struct('FSD','ON',...
                   'TASR','OFF',...
                   'RF_EL','OFF',...  % Response Function(RF) EnergyLoss
                   'RF_BF','OFF',...  % RF B-Fields
                   'RF_RX','OFF'); % choose 1 Effect only!  
MRA.ComputeCM('SysEffect',myEffects,'nTrials',500,'Recompute','ON','DataDriven','ON');
MRA.FitCM_Obj.ComputeCM_FSD;
%MRA.FitCM_Obj.ComputeCM_TASR;
%MRA.FitCM_Obj.ComputeCM_RF;
MRA.FitCM_Obj.PlotCM('ConvergenceTest','OFF','saveplot','OFF');

%% Generate 10000 IBD based on MRA.FitCM_Obj.CovMat
tbdisXXXs = MRA.FitCM_Obj.GenerateTBDISmultisims('nMultisims',10000);

% Plot Distribution of Normalizations of each TBDIS
v  = sum(tbdisXXXs');       % normalization of each simulation
vm = mean(sum(tbdisXXXs')); % mean normalization factor 
figure(555)
subplot(2,2,1)
histogram(v./vm);
title('Multisim normalizations');
xlabel('Normalization');
ylabel('multisims');
PrettyFigureFormat

% Normalize each new simulation to the normalization factor 
% All spectra have now the same number of events
tbdisXXXsNorm = tbdisXXXs./v'.*vm;

% Recompute New Covariance Matrix - Shape ONLY via MC
 CMsoMC = cov(tbdisXXXsNorm);
% figure(666)
% imagesc(CMsoMC);
% colorbar;
% PrettyFigureFormat

% Compute Fractional Shape Only MC Covariance Matrix
CMsoMCFrac = bsxfun(@rdivide,CMsoMC,mean(tbdisXXXsNorm)');  %divides 1st row of CovMat by TBDIS(1), etc...
CMsoMCFrac = bsxfun(@rdivide,CMsoMCFrac,mean(tbdisXXXsNorm));   %divides columnwise
subplot(2,2,3)
imagesc(MRA.FitCM_Obj.CovMatFrac);
title('Original Frac');
xlabel('qU bin'); ylabel('qU bin');
colorbar;
PrettyFigureFormat
pbaspect([1 1 1])
subplot(2,2,4)
imagesc(CMsoMCFrac);
title('New Shape Only Frac');
xlabel('qU bin'); ylabel('qU bin');
colorbar;
PrettyFigureFormat
pbaspect([1 1 1])

% Plot Difference CovMat - CovMatsoMC
subplot(2,2,2)
imagesc(MRA.FitCM_Obj.CovMatFrac-CMsoMCFrac);
title(' CovMatFrac - MC CovMatShapeOnlyFrac');
xlabel('qU bin'); ylabel('qU bin');
colorbar;
pbaspect([1 1 1])
PrettyFigureFormat

%% Systematics Envelop (diagonal)
figure(222)
[SysLine SysArea] = boundedline(MRA.ModelObj.qU-MRA.ModelObj.Q_i, ones(MRA.ModelObj.nqU,1),...
    sqrt(diag(MRA.FitCM_Obj.CovMat))./(MRA.nRuns.*MRA.ModelObj.TBDIS),'-b');
hold on
[SysLine SysArea] = boundedline(MRA.ModelObj.qU-MRA.ModelObj.Q_i, ones(MRA.ModelObj.nqU,1),...
    sqrt(diag(CMsoMC))./(MRA.nRuns.*MRA.ModelObj.TBDIS),'-r');
hold off