% Display a Single Covariance Matrix
%RunList = [40538:40543,40603,40604,40610:40613,40667:40693];
RunList = [40673:40680];
MRA = MultiRunAnalysis('RunList',RunList,'exclDataStart',1);

%Define Covariance Matrix properties
myEffects = struct('TCoff_RAD','ON','TCoff_OTHER','ON'); % choose 1 Effect only!  
MRA.InitializeCM('SysEffect',myEffects,'nTrials',1,...
                  'Recompute','ON','DataDriven','ON');
              
MRA.FitCM_Obj.ComputeCM_TCoff;
MRA.FitCM_Obj.PlotCM('ConvergenceTest','OFF','saveplot','OFF');

saveas(110,'plots/CMplot_TCoff.png');

%% Decomposition
figqU=figure(42)
set(figqU, 'Units', 'normalized', 'Position', [0.1, 0.2, 0.7, 0.7]);
subplot(2,2,1)
imagesc(MRA.FitCM_Obj.CovMatFracNorm);pbaspect([1 1 1])
xlabel('qU bin');ylabel('qU bin');
colorbar
title('TCoff decomposition: normalization term')
PrettyFigureFormat
subplot(2,2,2)
imagesc(log(abs(MRA.FitCM_Obj.CovMatFracNorm)));pbaspect([1 1 1])
xlabel('qU bin');ylabel('qU bin');
colorbar
title('log(abs(normalization term))')
PrettyFigureFormat
subplot(2,2,3)
imagesc(MRA.FitCM_Obj.CovMatFracShape);pbaspect([1 1 1])
xlabel('qU bin');ylabel('qU bin');
colorbar
title('RF decomposition: shape + mixed terms')
PrettyFigureFormat
subplot(2,2,4)
imagesc(log(abs(MRA.FitCM_Obj.CovMatFracShape)));pbaspect([1 1 1])
xlabel('qU bin');ylabel('qU bin');
colorbar
title('log(abs(shape + mixed terms))')
PrettyFigureFormat
saveas(42,'plots/CMdecompose_TCoff.png');

%% Decomposition Correlation Matrix
figqU=figure(422)
set(figqU, 'Units', 'normalized', 'Position', [0.1, 0.2, 0.9, 0.5]);
subplot(1,2,1)
corplot(MRA.FitCM_Obj.CovMatFracNorm);pbaspect([1 1 1])
xlabel('qU bin');ylabel('qU bin');
colorbar
title('TCoff decomposition: normalization term')
PrettyFigureFormat
subplot(1,2,2)
corplot(MRA.FitCM_Obj.CovMatFracShape);pbaspect([1 1 1])
xlabel('qU bin');ylabel('qU bin');
colorbar
title('log(abs(normalization term))')
PrettyFigureFormat
saveas(422,'plots/CorMdecompose_TCoff.png');

