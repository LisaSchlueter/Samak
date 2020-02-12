% Display a Single Covariance Matrix
%RunList = [40531,40538:40543,40603,40604,40610:40613,40667:40693];
MRA = MultiRunAnalysis('RunList','StackCD100all','exclDataStart',1);

%Define Covariance Matrix properties
myEffects = struct('FSD','ON'); % choose 1 Effect only!  
MRA.ComputeCM('SysEffect',myEffects,'nTrials',100',...
                  'Recompute','OFF','DataDriven','OFF','StackCM','OFF');
              
MRA.FitCM_Obj.ComputeCM_FSD;
MRA.FitCM_Obj.PlotCM('ConvergenceTest','ON','saveplot','OFF');