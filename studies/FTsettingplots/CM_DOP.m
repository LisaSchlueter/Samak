addpath(genpath('../../../Samak2.0'));
close all

% Build FSD Covariance Matrix with a given TD
% Must use Flat Time Distribution
A = ref_FTMTD('TD','Flat200');
B = ref_FTMTD('TD','Flat200','DopplerEffectFlag','matConv');

% Fix the Grid for interpolation
[X,Y] = meshgrid(A.qU);

% Define your CM Class with your number of trials, choose SysEffect,...
myEffects = struct(...
    'DOPoff','ON'); 
C = CovarianceMatrix('StudyObject',A,...
    'AltObject',B,...   
    'nTrials',0,...
    'RunFlag','OFF', 'nRuns',1,...
    'SysEffect',myEffects,...
    'RecomputeFlag','ON');
C.ComputeCM_DOPoff;
C.PlotCM
V = C.CovMatFrac;

% Create a TBDIS with another TD 
C = ref_FTMTD('TD','FT-TL4');

% Interpolate Doppler CM
[Xq,Yq] = meshgrid(C.qU);
VI = interp2(X,Y,V,Xq,Yq,'makima');

%% Plot Results
figure(1)
subplot(1,2,1)
imagesc(V);
colormap(flipud(gray));
colorbar
xlabel('qU bin'); ylabel('qU bin');
title('original DOPoff Flat200')
PrettyFigureFormat
pbaspect([1 1 1])
subplot(1,2,2)
imagesc(VI);
colormap(flipud(gray)); 
colorbar
xlabel('qU bin'); ylabel('qU bin');
title('interpolated DOPoff - FT-TL4')
PrettyFigureFormat
pbaspect([1 1 1])
set(gcf, 'Position', [100, 100, 1000, 500])
saveas(gcf,'CM_DOPoff.png')

% Test with correlation Matrix 
figure(2)
subplot(1,2,1)
corplot(V);
xlabel('qU bin'); ylabel('qU bin');
title('original DOPoff Flat200')
pbaspect([1 1 1])
subplot(1,2,2)
corplot(VI);
xlabel('qU bin'); ylabel('qU bin');
title('interpolated DOPoff - FT-TL4')
pbaspect([1 1 1])
set(gcf, 'Position', [100, 100, 1000, 500])
saveas(gcf,'CMCor_DOPoff.png')
