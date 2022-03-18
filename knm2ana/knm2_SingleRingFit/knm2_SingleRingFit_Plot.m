range     = 40;
freePar   = 'mNu E0 Bkg Norm';
chi2      = 'chi2Stat+';
DataType  = 'Real';
RingMerge = 'None';
BKG_PtSlope = 3*1e-06;
FSDFlag   = 'KNM2';
RecomputeFlag = 'OFF';

savedir = [getenv('SamakPath'),'knm2ana/knm2_SingleRingFit/results/'];
savename = sprintf('%sknm2_SingleRingFit_%s_%s_%s_BkgPT%.2g_%.0feV_%s_Ring%s.mat',...
    savedir,DataType,strrep(freePar,' ',''),FSDFlag,BKG_PtSlope.*1e6,range,chi2,RingMerge);

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
    fprintf('load %s \n',savename);
else
    %% settings
    SigmaSq =  0.0124+0.0025;
    if strcmp(chi2,'chi2Stat')
        NP = 1;
    elseif strcmp(chi2,'chi2Stat+')
        NP = 1.112;
        chi2 = 'chi2Stat';
    else
        NP = 1.112;
    end
    RunAnaArg = {...
        'RunList','KNM2_Prompt',...  % define run number -> see GetRunList
        'fixPar',freePar,...         % free Parameter !!
        'DataType',DataType,...              % Real, Twin or Fake
        'FSDFlag',FSDFlag,...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'RingMerge',RingMerge,...             % 'Full' == 4 Pseudo rings
        'fitter','minuit',...
        'minuitOpt','min;minos',...
        'chi2',chi2,...
        'NonPoissonScaleFactor',NP,...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','ON',...
        'SysBudget',38,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'DopplerEffectFlag','FSD',...
        'BKG_PtSlope',BKG_PtSlope};
    
    
    %% read data and set up model: MultiRunAnalysis
    A = MultiRunAnalysis(RunAnaArg{:}); % object of class MultiRunAnalysis
    A.exclDataStart = A.GetexclDataStart(range);
    
    %% start ringwise analysis
    
    R = RingAnalysis('RunAnaObj',A,'RingList',A.RingList); % object of class RingAnalysis
    
    %% fit every ring - one after the other
    R.FitRings('SaveResult','ON',...
        'RecomputeFlag',RecomputeFlag,...  % load from storage or recalculate
        'AsymErr','OFF');                 % asymmetric from scan and more correct uncertainties -> only for mNuSq
    
    save(savename,'R');
end
%% display
CommonArg = {'SavePlot','OFF',...
    'Blind','OFF',...       % show relative or absolute results
    'PlotMode','Abs'};
%%
Q_i = R.MultiObj(1).ModelObj.Q_i ;
if strcmp(RingMerge,'Full')
    YLim = [-2.5,3;Q_i-0.13 Q_i+0.13;1.65 2.3;0.98 ,0.994];
else
     YLim = [-5.3,7.5;Q_i-0.27, Q_i+0.42;1.55 2.45;0.9652 ,1.008];
end
%%
R.PlotFits(CommonArg{:},...
    'PlotPar',1,...        % 1 == neutrino mass, 2 == E0
    'YLim',YLim(1,:),... % force y-axis to good limits
    'linFit','ON');
%%
    R.PlotFits(CommonArg{:},...
        'PlotPar',2,...        % 1 == neutrino mass, 2 == E0
        'YLim',YLim(2,:),... % force y-axis to good limits
        'linFit','ON');        % show linear fit
  %%  
  R.PlotFits(CommonArg{:},...      % show relative or absolute values
        'PlotPar',3,...        % 1 == neutrino mass, 2 == E0
        'YLim',YLim(3,:),... % force y-axis to good limits
        'linFit','ON');        % show linear fit
%%
 R.PlotFits(CommonArg{:},...      % show relative or absolute values
        'PlotPar',4,...        % 1 == neutrino mass, 2 == E0
        'YLim',YLim(4,:),... % force y-axis to good limits
        'linFit','ON');        % show linear fit
%% some statistics
Par = 4;
x = R.FitResult.par(:,Par);
if Par==1
    xerr = 0.5.*(R.FitResult.errPos(:,Par)-R.FitResult.errNeg(:,Par));
else
    xerr = R.FitResult.err(:,Par);
end
chi2const = sum((x-wmean(x,1./xerr.^2)).^2./xerr.^2);
pconst = 1-chi2cdf(chi2const,numel(x)-1);
fprintf('Compatible width no radial dependence at p=%.2f (chi2/dof = %.1f/%.0f) \n',pconst,chi2const,numel(x)-1);


