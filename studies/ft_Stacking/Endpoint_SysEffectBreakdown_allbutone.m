% ------------------------------------------------------------------------------------------------%
% Switch on systematic effects to see effect on fitted parameter
% All One - but 1 Systematic Effect
% ------------------------------------------------------------------------------------------------%

RunList = [40538:40543,40603,40604,40610:40613,40667:40693];
MRA = MultiRunAnalysis('RunList',RunList,'fixPar','1 5 6','chi2','chi2CM',...
    'exclDataStart',9);

% Select Systematic Effects used in Fit
myEffects = struct(...
    'RF_EL','ON',...  % Response Function(RF) EnergyLoss
    'RF_BF','ON',...  % RF B-Fields
    'RF_RX','ON',...  % Column Density, inel cross ection
    'FSD','ON',...
    'TASR','ON',...
    'TCoff_RAD','ON',...
    'TCoff_OTHER','ON');
MRA.InitializeCM('SysEffects',myEffects,...
    'DataDriven','ON',...
    'WGTS_CD_MolPerCm2_RelErr',0.05);
E0    = zeros(7,1);
E0Err = zeros(7,1);
for i=1:7
    MRA.SimulateStackRuns;  %reset to init ModelObject
    if i==1
        MRA.FitCM_Obj.SysEffect.RF_EL = 'OFF';         % Response Function OFF
        MRA.FitCM_Obj.SysEffect.RF_BF = 'OFF';
        MRA.FitCM_Obj.SysEffect.RF_RX = 'OFF';
    elseif i==2
        MRA.FitCM_Obj.SysEffect.RF_EL = 'ON';          % FSD OFF
        MRA.FitCM_Obj.SysEffect.RF_BF = 'ON';
        MRA.FitCM_Obj.SysEffect.RF_RX = 'ON';
        MRA.FitCM_Obj.SysEffect.FSD = 'OFF';
    elseif i==3
        MRA.FitCM_Obj.SysEffect.FSD = 'ON';            % TASR OFF
        MRA.FitCM_Obj.SysEffect.TASR = 'OFF';
    elseif i==4
        MRA.FitCM_Obj.SysEffect.TASR = 'ON';           % Radiatice Corrections OFF
        MRA.FitCM_Obj.SysEffect.TCoff_RAD = 'OFF';
    elseif i==5
        MRA.FitCM_Obj.SysEffect.TCoff_RAD = 'ON';      % Other Theo. Correction OFF
        MRA.FitCM_Obj.SysEffect.TCoff_OTHER = 'OFF';
    elseif i==6
        MRA.FitCM_Obj.SysEffect.TCoff_OTHER = 'ON';     % all ON
        MRA.FitCM_Obj.SysEffect.TCoff_RAD = 'ON';
        MRA.FitCM_Obj.SysEffect.TASR = 'ON';
        MRA.FitCM_Obj.SysEffect.RF_EL = 'ON';
        MRA.FitCM_Obj.SysEffect.RF_BF = 'ON';
        MRA.FitCM_Obj.SysEffect.RF_RX = 'ON';
        MRA.FitCM_Obj.SysEffect.FSD = 'ON';
    elseif i==7
         MRA.chi2 ='chi2Stat';                         % Statitcs only (all OFF)
    end
    MRA.ComputeCM('InitNormFit','ON');
    MRA.Fit;
    E0(i)    = MRA.ModelObj.Q_i+MRA.FitResult.par(2);
    E0Err(i) = MRA.FitResult.err(2);
    %close all;
    %MRA.PlotFit('ResidualsFlag','ON','saveplot','ON','Mode','Count');
end
fig10 = figure(10);
errorbar([1:7]',E0,E0Err,'Color',rgb('CadetBlue'),'LineWidth',3,'LineStyle','--');
ylabel('E0 (eV)');
xticks([1 2 3 4 5 6 7]);
xticklabels({'Response Function','Final States','Activity Fluctuation','Radiative Corrections','Other Theo. Corrections','Combi CM','Statistical'});
%xticklabels({'RF','FSD','TASR','TCoffRad','TCoffOTHER','Combi CM','Stat Only'})
xlim([0.5 7.5])
xlabel('Covariance Matrix for ALL effects, but...');
PrettyFigureFormat
set(gca,'FontSize',16)
Runtitle1 = sprintf('KATRIN First Tritium %.0feV below E0 eV',MRA.ModelObj.qU(MRA.exclDataStart)-18575);
Runtitle2 = sprintf('\n %uRuns Stacked %.0f - %.0f',numel(MRA.StackedRuns),MRA.RunList(1),MRA.RunList(end));
title([Runtitle1,Runtitle2])
