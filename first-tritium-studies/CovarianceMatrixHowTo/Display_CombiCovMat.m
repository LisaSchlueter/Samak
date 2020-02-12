% Display a Single Covariance Matrix
%RunList = [40531,40538:40543,40603,40604,40610:40613,40667:40693];
%exclDataStart=1;
%MRA = MultiRunAnalysis('RunList',RunList,'exclDataStart',exclDataStart);
MRA = RunAnalysis('RunNr',40668);
%Define Covariance Matrix properties
myEffects = struct(...
               'RF_EL','ON',...  % Response Function(RF) EnergyLoss
               'RF_BF','ON',...  % RF B-Fields
               'RF_RX','ON',...  % Column Density, inel cross ection
               'FSD','ON',...
               'TASR','ON',...
               'TCoff_RAD','ON',...
               'TCoff_OTHER','ON',...
               'DOPoff','OFF');
           
MRA.ComputeCM('SysEffect',myEffects,'nTrials',1000',...
                  'Recompute','OFF','DataDriven','OFF');
MRA.FitCM_Obj.PlotCM('ConvergenceTest','OFF','saveplot','ON');
%MRA.FitCM_Obj.PlotStack('plotStat','OFF','saveplot','OFF');

%strfile = sprintf('CombiCM-NormShape-Excl%g',MRA.exclDataStart);
%saveas(110,['plots/' strfile '.png']);

%% Decomposition
figqU=figure('Renderer','opengl');
set(figqU, 'Units', 'normalized', 'Position', [0.1, 0.2, 0.7, 0.7]);
subplot(2,2,1)
imagesc(MRA.FitCM_Obj.CovMatFracNorm);pbaspect([1 1 1])
xlabel('qU bin');ylabel('qU bin');
colorbar
title('Combi decomposition: normalization term')
%PrettyFigureFormat
subplot(2,2,2)
imagesc(log(abs(MRA.FitCM_Obj.CovMatFracNorm)));pbaspect([1 1 1])
xlabel('qU bin');ylabel('qU bin');
colorbar
title('log(abs(normalization term))')
%PrettyFigureFormat
subplot(2,2,3)
imagesc(MRA.FitCM_Obj.CovMatFracShape);pbaspect([1 1 1])
xlabel('qU bin');ylabel('qU bin');
colorbar
title('Combi decomposition: shape + mixed terms')
%PrettyFigureFormat
subplot(2,2,4)
imagesc(log(abs(MRA.FitCM_Obj.CovMatFracShape)));pbaspect([1 1 1])
xlabel('qU bin');ylabel('qU bin');
colorbar
title('log(abs(shape + mixed terms))')
%PrettyFigureFormat
saveas(42,'plots/CMdecompose_Combi.png');

%% Stackplot Save
strfile = sprintf('Stack-CombiCM-NormShape-Excl%g',200);
saveas(999,['plots/' strfile '.png']);