addpath(genpath('../../../Samak2.0'));

A =InitKatrin_ft_RF('TD','FT-TL4');
A.TimeSec = 24*60*60;
myEffects = struct(...
                 'RF_EL','OFF',...  % Response Function(RF) EnergyLoss
                 'RF_BF','OFF',...  % RF B-Fields
                 'RF_RX','OFF',...  % RF Column Density, Cross Section
                 'TASR','ON',...  % Tritium Activity Stack Runs);
                 'FSD','OFF',...   % FSD                
                 'DOPoff','OFF',...%,...% Doppler (ON/OFF only)
                 'TCoff_RAD','OFF',...
                 'TCoff_OTHER','OFF');   % Theoretical corrections (ON/OFF only)
             
C = CovarianceMatrix('StudyObject',A, 'nTrials',100,'RunFlag','OFF', 'nRuns',1, 'SysEffect',myEffects);
C.RecomputeFlag = 'ON';
C.ComputeCM_TASR;
%C.ComputeCM;
%C.PlotStack('plotStat','ON','saveplot','ON');