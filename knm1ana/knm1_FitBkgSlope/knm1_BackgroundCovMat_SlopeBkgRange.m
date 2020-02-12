% Background slope
% check systematic uncertainty as a function of spectrum points, used for
% background slope covariance matrix calculation

%% set up model: 
% settings
RunList               = 'KNM1';
exclDataStart         = 14; % 40eV range = 27 subruns
BkgConstrain = 'OFF';

savedir = [getenv('SamakPath'),'knm1ana/knm1_FitBkgSlope/results/'];
DataType = 'Twin';
savename = [savedir,sprintf('knm1_BackgroundCovMat_SlopeBkgRange_%s_%s_%.0f_BkgConstrain%s.mat',...
    RunList,DataType,exclDataStart,BkgConstrain)];

if exist(savename,'file')
    load(savename);
else
    MakeDir(savedir); % create results dir, if it doesnt exist yet
    
    % qU points in KNM1 MTD close to endpoint (eV): -4.6, -2.6, 2.4, 7.4, 17.4, 27.4, 47.4
    BkgRange = [-5,-3,0,3]; %eV with respect to endpoint, background points used for cov mat
    mNuSqCMErr = zeros(numel(BkgRange),1);
    
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
    mNuSqStatErr(1) = A.FitResult.err(1);
    
    %% Fit: stat + background slope systematics, background slope fixed
    A.fixPar = '5 6 7 8 9 10 11 12';
    A.chi2 = 'chi2CMShape';
    switch BkgConstrain
        case 'ON'
            MaxSlopeCpsPereV=15*1e-06;
        case 'OFF'
            MaxSlopeCpsPereV = 1e9;
    end
    
    for i=1:numel(BkgRange)
        A.ComputeCM('BkgCM','ON','SysEffects',struct('FSD','OFF'),...
            'MaxSlopeCpsPereV',MaxSlopeCpsPereV,'BkgRange',BkgRange(i));
        A.Fit;
        mNuSqCMErr(i) = A.FitResult.err(1);
        A.ModelObj.SetFitBias(0);
    end
    save(savename,'mNuSqStatErr','mNuSqCMErr','MaxSlopeCpsPereV','BkgRange');
end

%%
fig1 = figure('Units','normalized','pos',[0.1,0.1,0.5,0.5]);
plot(BkgRange,sqrt(mNuSqCMErr.^2-mNuSqStatErr^2),'x-','LineWidth',3,'MarkerSize',9);
PrettyFigureFormat('FontSize',24);
xlabel('Lower boundary below the endpoint (eV)');
switch DataType
    case 'Real'
        ylabel(sprintf('1\\sigma uncertainty on m_\\nu^2 (eV^2)'));
        
    case 'Twin'
        ylabel(sprintf('1\\sigma sensitivity on m_\\nu^2 (eV^2)'));
end
title('background range - slope systematics')
grid on;



