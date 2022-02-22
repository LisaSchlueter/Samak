% mnusq uncertainty as a function of slope constraints

% settings
DataType = 'Twin';
range = 40;

%MaxSlopeCpsPereV_v = [99,(30:-5:5).*1e-06,1e-06];
MaxSlopeCpsPereV_v = [1e9,[5,15].*1e-06];

savedir = [getenv('SamakPath'),'knm1ana/knm1_FitBkgSlope/results/'];

savename = [savedir,sprintf('knm1_BkgSlopeErr_%s_%.0feVrange_BkgConstraints_%.0fto%.0fmuCpsPerS_%.0ffits.mat',...
    DataType,range,1e6.*min(MaxSlopeCpsPereV_v),1e6.*max(MaxSlopeCpsPereV_v),numel(MaxSlopeCpsPereV_v))];

if exist(savename,'file')
    load(savename);
else
    MakeDir(savedir);

    % Init Model Object and covariance matrix object
     A = MultiRunAnalysis('RunList','KNM1',...
        'chi2','chi2Stat',...
        'DataType',DataType,...
        'fixPar','mNu E0 Norm Bkg',...
        'NonPoissonScaleFactor',1,...
        'SysBudget',22,...
        'minuitOpt','min ; minos',...
        'FSDFlag','Sibille0p5eV',...
        'ELossFlag','KatrinT2',...
        'AngularTFFlag','OFF',...
        'SynchrotronFlag','ON',...
        'RadiativeFlag','ON',...
        'DopplerEffectFlag','OFF');
    
    A.exclDataStart =    A.GetexclDataStart(range);
    
      FitResults = cell(numel(MaxSlopeCpsPereV_v)+1,1);
      FitCMFracShape = cell(numel(MaxSlopeCpsPereV_v)+1,1);
      mNuSqErr = zeros(numel(MaxSlopeCpsPereV_v)+1,1);
    %% Fit: stat only + background slope fixed for reference
    A.chi2 = 'chi2Stat';
    A.InitModelObj_Norm_BKG;
    A.Fit;
    FitResults{1} = A.FitResult;
    mNuSqErr(1) = A.FitResult.err(1);
    FitCMFracShape{1} = A.FitCMFracShape;
    %% Fit: stat only + bkg cm
    A.chi2 = 'chi2CMShape';
    
    for i=1:numel(MaxSlopeCpsPereV_v)
        A.InitModelObj_Norm_BKG;
        A.ComputeCM('SysEffects',struct('FSD','OFF'),...
            'BkgCM','ON','BkgPtCM','OFF',...
            'BkgMode','SlopeFit',...
            'MaxSlopeCpsPereV',MaxSlopeCpsPereV_v(i));
         FitCMFracShape{i+1} = A.FitCMFracShape;
        A.InitModelObj_Norm_BKG;
        A.Fit;
        FitResults{i+1} = A.FitResult;
        mNuSqErr(i+1) = A.FitResult.err(1);
        A.ModelObj.SetFitBias(0);
    end
   
  
   
    save(savename,'FitResults','mNuSqErr','MaxSlopeCpsPereV_v','FitCMFracShape');
end

%% display
if any(mNuSqErr(:)<mNuSqErr(1))
    Idx = find(mNuSqErr(:)<mNuSqErr(1));
    mNuSqErr(Idx) = mNuSqErr(1);
end
mNuSq_sys = sqrt(mNuSqErr(2:end).^2-mNuSqErr(1).^2);
GetFigure;
plot(MaxSlopeCpsPereV_v,mNuSq_sys,'x')
