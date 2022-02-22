% Test fit of background slope
% compare background slope systematics:
% - with a covariance matrix
% - nuissance parameter + pull Term
mNuSqErr = zeros(3,1);
%% set up model: 
% settings
RunList               = 'KNM1';
exclDataStart         = 13; % 40eV range = 27 subruns
BkgConstrain = 'OFF';
BkgRange = -5; %eV with respect to endpoint, background points used for cov mat

savedir = [getenv('SamakPath'),'knm1ana/knm1_FitBkgSlope/results/'];
DataType = 'Real';
savename = [savedir,sprintf('knm1_FitBkgSlope_%s_%s_%.0f_BkgConstrain%s_%.0feVBkgRange.mat',...
    RunList,DataType,exclDataStart,BkgConstrain,BkgRange)];

if exist(savename,'file')
    load(savename);
else
    MakeDir(savedir);
    % Init Model Object and covariance matrix object
    A = MultiRunAnalysis('RunList',RunList,...
        'chi2','chi2Stat','DataType',DataType,...
        'exclDataStart',exclDataStart,...
        'fixPar','5 6 7 8 9 10 11',...
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',1.064,...
        'minuitOpt','min ; migrad',...
        'FSDFlag','Sibille',...
        'ELossFlag','KatrinT2',...
        'SysBudget',22);
    %% Fit: stat only + background slope fixed for reference
    A.fixPar = '5 6 7 8 9 10 11 12';
    A.chi2 = 'chi2Stat';
    A.Fit;
    A.ModelObj.SetFitBias(0);
    mNuSqErr(1) = A.FitResult.err(1);
    %% Fit: stat only + with free background slope
    A.fixPar = '5 6 7 8 9 10 11';
    switch BkgConstrain
        case 'OFF'
            A.pullFlag = 3; % no pull
        case 'ON'
            A.pullFlag = 6;
    end
    A.chi2 = 'chi2Stat';
    A.Fit('CATS','OFF');
    mNuSqErr(2) = A.FitResult.err(1);
    A.ModelObj.SetFitBias(0);
    %% Fit: stat + background slope systematics, background slope fixed
    A.fixPar = '5 6 7 8 9 10 11 12';
    A.chi2 = 'chi2CMShape';
    switch BkgConstrain
        case 'ON'
            MaxSlopeCpsPereV=15*1e-06;
        case 'OFF'
            MaxSlopeCpsPereV = 1e9;
    end
    A.ComputeCM('BkgCM','ON','SysEffects',struct('FSD','OFF'),...
        'MaxSlopeCpsPereV',MaxSlopeCpsPereV,'BkgRange',BkgRange);
    A.Fit;
    mNuSqErr(3) = A.FitResult.err(1);
    A.ModelObj.SetFitBias(0);
    %% save
    save(savename,'mNuSqErr','MaxSlopeCpsPereV');
end
%%
fprintf('--------------------------------------------------------------\n')
fprintf('Nu-mass squared sensitivity increase with respect to stat + bkg slope fixed \n');
fprintf('Data Type: %s\n',DataType);
fprintf('Background constrain: %s (%.3g mcps/eV) \n',BkgConstrain,MaxSlopeCpsPereV);
fprintf('Background slope as free fit parameter: %.3f eV^2 \n',sqrt(mNuSqErr(2)^2-mNuSqErr(1)^2))
fprintf('Background slope fixed + cov mat: %.3f eV^2 \n',sqrt(mNuSqErr(3)^2-mNuSqErr(1)^2))
fprintf('--------------------------------------------------------------\n')

