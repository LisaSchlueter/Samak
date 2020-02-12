addpath(genpath('../../../Samak2.0'));
set(0,'DefaultFigureWindowStyle','normal');

%% 1. Define your Study Obect (TBD Object)
A =ref_RunSummaries_StackPix(40666,'ex2');
A.TimeSec = 24*60*60;
%% 2. Choose systematic effect(s)
myEffects = struct(...
                 'RF_EL','OFF',...  % Response Function(RF) EnergyLoss
                 'RF_BF','OFF',...  % RF B-Fields
                 'RF_RX','OFF',...  % RF Column Density, Cross Section
                 'TASR','OFF',...  % Tritium Activity Stack Runs);
                 'FSD','ON',...   % FSD
                 'BM1S','OFF',...  % Background
                 'DOPoff','OFF',...%,...% Doppler (ON/OFF only)
                 'TCoff','OFF');   % Theoretical corrections (ON/OFF only)
             
%% 3. Define your CM Class with your number of trials, choose SysEffect,...
C = CovarianceMatrix('StudyObject',A, 'nTrials',100,'RunFlag','OFF', 'nRuns',1, 'SysEffect',myEffects);
%% 4. Chose your uncertainties
NormFlag  = 'OFF';
ShapeFlag = 'ON';
            
%% 5.(optional) Enforces Recomputing all Lookup Tables
% -If unsure, if computed correctly before -> RecomputeFlag = 'ON'
% -If CM is already computed previously with a different RunTime than your
% Study Object: Choose RecomputeFlag = 'OFF' 
%-> CM will be normalized automatically with the fractional CM
C.RecomputeFlag = 'ON';

%% 6.Compute Covariance Matrix for your SysEffect
%C.ComputeCM_FSD('NormFlag',NormFlag,'ShapeFlag',ShapeFlag);   %Compute CM for Final States Distribution Effect
C.ComputeCM_FSD();   %Compute CM for Final States Distribution Effect

%% 8. Plots
C.PlotCM;       %Option: 'saveplot','ON' -->saves plot in directory "plots" (you need to create this in your study)
C.DecomposeCM;  % Plot Decomposed CM   ; Option: 'saveplot','ON

%% 9.(optional) Tests: Rank and Convergence
C.ComputeRank;    %Compute the Rank of the CM as a first sanity check
C.ConvergenceTest('filename',C.CovMatFile,'Nmax',C.nTrials); 
