% compute stacking covariance matrix for different qUErr
if ~exist('M','var')
    M = MultiRunAnalysis('RunList','KNM1','chi2','chi2Stat',...
        'fixPar','5 6 7 8 9 10 11',...
        'DataType','Twin',...
        'FSDFlag','Sibille0p5eV',...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'exclDataStart',14);
end

M.ComputeCM('SysEffects',struct('Stack','ON'),'BkgCM','OFF');
Stack_qUErr_i = M.FitCM_Obj.Stack_qUErr;
Stack_qUfracRelErr_i = M.FitCM_Obj.Stack_qUfracRelErr;
%%
Stack_qUErr = [0.1,0.5,1];
mNuSqErr = zeros(numel(Stack_qUErr)+1,1);
M.chi2 = 'chi2Stat';
M.Fit;
mNuSqErr(1) = M.FitResult.err(1);

M.chi2 = 'chi2CMShape';
CM = M.FitCM_Obj;
CM.nTrials= 1000;

for i=1:numel(Stack_qUErr)
M.ComputeCM('SysEffects',struct('Stack','ON'),'BkgCM','OFF',...
  'Stack_qUErr',Stack_qUErr(i).*Stack_qUErr_i,'DataDriven','OFF',...
  'Stack_qUfracRelErr',Stack_qUfracRelErr_i)
M.Fit;
mNuSqErr(i+1) = M.FitResult.err(1);
end

%%
plot(1:numel(mNuSqErr),sqrt(mNuSqErr.^2-mNuSqErr(1)^2),'x');


