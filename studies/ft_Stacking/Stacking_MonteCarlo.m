% ------------------------------------------------------------------------%
% Fit to Simulated Data: Stacking Runs Test
% L.Schlueter (June 2018)
% ------------------------------------------------------------------------%

RunList = [40538:40543, 40603,40604,40610:40613,40667:40693]; %
A = MultiRunAnalysis('RunList',RunList,'AnaFlag','StackPixel',...
    'chi2','chi2Stat','ringCutFlag','ex2','fix','');
%FitCM = A.FitCM; %Stack + normal CM

% Load Sinle Runs
if isempty(A.SingleRunObj)
    A.LoadSingleRunObj;
end

% Compute Stacking 
TBDIS_Stack = zeros(A.ModelObj.nqU,1);
for i=1:length(A.StackedRuns)
 TBDIS_Stack = TBDIS_Stack + A.SingleRunObj{i,1}.TBDIS;
end

% Fit Simulation; not Data
A.RunData.qU = A.ModelObj.qU;
A.RunData.TBDIS = TBDIS_Stack;

%%
chi2 = 'Stat';

switch chi2
     case 'StackCMOnly'
        A.chi2 = 'chi2CM';
        A.FitCM = (A.ModelObj.TBDIS(:).*(A.StackCM_Obj.CovMatFrac).*A.ModelObj.TBDIS(:)') + diag(A.ModelObj.TBDIS(:));
    case 'CM' % CM + Stack 
        A.chi2 = 'chi2CM';
        A.ComputeCM;
        A.ComputeCM_Stack;
    case 'Stat'
        A.chi2 = 'chi2Stat';       
end

close all;
A.Fit;
A.PlotFit('ResidualsFlag','ON','saveplot','OFF');


% TBDISmultisimsys = mvnrnd(A.ModelObj.TBDIS,A.ModelObj.TBDIS(:).*(A.FitStackCMFrac).*A.ModelObj.TBDIS(:)',100);  
% TBDISmultisimstat = mvnrnd(A.ModelObj.TBDIS,A.FitCMStat,100);
% figure(8);
% stairs((TBDISmultisimsys'./A.ModelObj.TBDIS));
% hold on;
% stairs(TBDISmultisimstat'./A.ModelObj.TBDIS);


