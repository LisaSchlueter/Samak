%
% ?Resizing of the FSD+TASR covariance matrix 
%  from the FT-TL3 MTD to the FT-TL4 MTD 
%
%  Thierry Lasserre
%  May 16 2018
%

addpath(genpath('../../../Samak2.0'));
close all

% Build Covariance Matrix with a given TD
A = ref_FTMTD('TD','FT-TL3');
[X,Y] = meshgrid(A.qU);

% Define your CM Class with your number of trials, choose SysEffect,...
myEffects = struct(...
    'RF_EL','OFF',...        % Response Function(RF) EnergyLoss
    'RF_BF','OFF',...        % RF B-Fields
    'RF_RX','OFF',...        % RF Column Density, Cross Section
    'TASR','ON',...         % Tritium Activity Stack Runs);
    'FSD','ON',...           % FSD Noralization
    'DOPoff','OFF',...       % Doppler (if switched off)
    'TCoff','ON',...        % thepretical corrections (if switched off)
    'TCoff_RAD','OFF',...    % theretical corrections (if switched off)
    'TCoff_OTHER','OFF',...  % theretical corrections (if switched off)
    'BM1S','OFF');           % Background
C = CovarianceMatrix(...
    'StudyObject',A,...
    'nTrials',1000,...
    'RunFlag','OFF','nRuns',1,...
    'SysEffect',myEffects,...
    'WGTS_TASR_RelErr',5e-3,...
    'RecomputeFlag','ON');

% Compute / Plot CM
C.ComputeCM;
V = C.CovMatFrac;

% Create a TBDIS with another TD 
B = ref_FTMTD('TD','FT-TL4');

% Interpolate CM
[Xq,Yq] = meshgrid(B.qU);
VI = interp2(X,Y,V,Xq,Yq,'makima');

%% Plot Results
figure(1)
subplot(1,2,1)
imagesc(V);
colormap(flipud(gray));
xlabel('qU bin'); ylabel('qU bin');
title('original CM FT-TL3')
PrettyFigureFormat
pbaspect([1 1 1])
subplot(1,2,2)
imagesc(VI);
colormap(flipud(gray)); 
xlabel('qU bin'); ylabel('qU bin');
title('interpolated CM - FT-TL4')
PrettyFigureFormat
pbaspect([1 1 1])
set(gcf, 'Position', [100, 100, 1200, 500])
saveas(gcf,'CMresizing_originalFT3.png')

% Test with correlation Matrix 
figure(2)
subplot(1,2,1)
corplot(V);
xlabel('qU bin'); ylabel('qU bin');
title('original CM - FT-TL3')
pbaspect([1 1 1])
subplot(1,2,2)
corplot(VI);
xlabel('qU bin'); ylabel('qU bin');
title('interpolated CM - FT-TL4')
pbaspect([1 1 1])
set(gcf, 'Position', [100, 100, 1200, 500])
saveas(gcf,'CMresizing_interplatedFT4.png')

