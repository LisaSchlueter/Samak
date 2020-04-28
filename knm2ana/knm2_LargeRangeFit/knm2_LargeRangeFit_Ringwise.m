function [Diffs_E0, E0, E0Err,FitResults,Rate_m,Rate_d,qU_d,Time] = knm2_LargeRangeFit_Ringwise
%% large range fit for single rings

%% Common settings
savedir = [getenv('SamakPath'),'tritium-data/fit/Knm2/SingleRingFit/LargeRange/'];
MakeDir(savedir);

ROIstr = 'Default';
Knm2AnaArg = {'fixPar','E0 Bkg Norm',...
    'DataType','Real',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'RingMerge','Full',...
    'NonPoissonScaleFactor',1.112,...
    'ROIFlag',ROIstr,...
    'AngularTFFlag','ON',...
    'SynchrotronFlag','ON'};

%% loop over all rings and periods
nRings = 4;
FitResults = cell(3,1);
E0         = zeros(3,nRings);
E0Err      = zeros(3,nRings);
Rate_d     = zeros(3,nRings,11);
Rate_m     = zeros(3,nRings,11);
qU_d       = zeros(3,nRings,11);
Time       = zeros(3,nRings,11);
for j = 1:3
    RunList = ['KNM2_RW' num2str(j)];
    start = 90;               % fit range in eV below endpoint
    stop  = 45;               % fit range in eV below endpoint
    savename = sprintf('%sSingleRingFit_%s_%.0feV_%.0feV.mat',savedir,RunList,start,stop);
    
    if exist(savename,'file') 
        file = importdata(savename);
        E0(j,:)       = file.FitResult.par(:,2);
        E0Err(j,:)    = file.FitResult.err(:,2);
        FitResults{j} = file.FitResult; % save all fit results
        qU_d(j,:,:)   = file.qU(file.includeDataIndex,:)';
        Rate_m(j,:,:) = file.Rate_model(file.includeDataIndex,:)';
        Rate_d(j,:,:) = file.Rate_data(file.includeDataIndex,:)';
        Time(j,:,:)   = file.TimeSubRunSec(file.includeDataIndex,:)';
        fprintf('load fit result from %s \n',savename);
    else
        D   = MultiRunAnalysis('RunList',RunList,Knm2AnaArg{:});
        % find good fit range
        
        BkgSubRun = numel(D.ModelObj.qU); % background subrun (last one)
        exclDataStart_tmp = find((D.ModelObj.qU)>=18574-start,1);
        exclDataStop_tmp  = find((D.ModelObj.qU)>=18574-stop,1);
        includeDataIndex = [exclDataStart_tmp:exclDataStop_tmp,BkgSubRun];
        D.exclDataStart = includeDataIndex;
        
        % init ring object
        R   = RingAnalysis('RunAnaObj',D,'RingList',1:nRings);
        R.FitRings('SaveResult','OFF',...
            'RecomputeFlag','ON',...  % load from storage or recalculate
            'AsymErr','OFF');
        FitResult = R.FitResult;
        Q_i   = D.ModelObj.Q_i;
        
        qU            = zeros(D.ModelObj.nqU,nRings);
        Rate_data     = zeros(D.ModelObj.nqU,nRings);
        Rate_model    = zeros(D.ModelObj.nqU,nRings);
        TimeSubRunSec = zeros(D.ModelObj.nqU,nRings);
       
        for r=1:nRings
            qU(:,r) = R.MultiObj(1).RunData.qU;
            Rate_data(:,r) = R.MultiObj(r).RunData.TBDIS./(R.MultiObj(r).RunData.qUfrac.*R.MultiObj(r).RunData.TimeSec);
            Rate_model(:,r) = R.MultiObj(r).ModelObj.TBDIS./(R.MultiObj(r).ModelObj.qUfrac.*R.MultiObj(r).ModelObj.TimeSec);
            TimeSubRunSec(:,r) = (R.MultiObj(r).ModelObj.qUfrac.*R.MultiObj(r).ModelObj.TimeSec);
        end
        
        % save
        save(savename,'FitResult','Q_i','Knm2AnaArg','RunList','includeDataIndex',...
                      'qU','Rate_data','Rate_model','TimeSubRunSec');
        fprintf('save fit result to %s \n',savename)
        E0(j,:)       = R.FitResult.par(:,2);
        E0Err(j,:)    = R.FitResult.err(:,2);
        FitResults{j} = R.FitResult; % save all fit results
        Rate_m(j,:,:) = Rate_model(includeDataIndex,:)';
        Rate_d(j,:,:) = Rate_data(includeDataIndex,:)';
        qU_d(j,:,:)   = qU(includeDataIndex,:)';
        Time(j,:,:)   = TimeSubRunSec(includeDataIndex,:)';
    end
end

% result with respect to rear period 2, ring 1
Diffs_E0 = E0 - E0(2,1);

%% sanity
chi2min = cell2mat(cellfun(@(x) x.chi2min,FitResults,'UniformOutput',false)');
dof = cell2mat(cellfun(@(x) x.dof,FitResults,'UniformOutput',false)');
pVal = 1-chi2cdf(chi2min,dof);
end
