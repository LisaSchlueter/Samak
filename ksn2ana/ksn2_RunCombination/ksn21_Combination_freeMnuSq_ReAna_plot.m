%% plot with combined ksn1+2 free nu-mass
chi2          = 'chi2CMShape';
DataType      = 'Twin';
nGridSteps    = 30;
freePar       = 'mNu E0 Bkg Norm';
range         = 40;
savedir = [getenv('SamakPath'),'SterileAnalysis/GridSearchFiles/Combi/',DataType,'/'];
MakeDir(savedir)
savename = sprintf('%sKSN12Combi_ReAna_GridSearch_%s_%s_Uniform_%s_%.0fnGrid.mat',savedir,DataType,strrep(freePar,' ',''),chi2,nGridSteps);

% load file
if exist(savename,'file')
    d = importdata(savename);
    fprintf('load results from file %s \n',savename)
else
    fprintf('cannot find file %s \n',savename)
    return
end

%%
savename2 = sprintf('%sKSN12Combi_ReAna_Interp1Obj.mat',savedir);
if exist(savename2,'file')
    A = importdata(savename2);
else
    % int obejct for interpolatio and plot
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar','mNu E0 Norm Bkg',...%free par
        'SysBudget',40,...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'RadiativeFlag','ON',...
        'FSDFlag','KNM2_0p1eV',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'chi2',chi2,...
        'NonPoissonScaleFactor',1.112,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'TwinBias_Q',18573.7,...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',3*1e-06,...
        'TwinBias_BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
    A = MultiRunAnalysis(RunAnaArg{:});
    save(savename2,'A');
end
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};

S = SterileAnalysis(SterileArg{:});

%%
S.sin2T4 = d.sin2T4;
S.mNu4Sq = d.mnu4Sq;
S.chi2 = d.chi2;
S.chi2_ref = d.chi2_ref;
S.chi2_Null = d.FitResult_Null.chi2min;
S.mNuSq = cell2mat(cellfun(@(x) x.par(1),d.FitResults,'UniformOutput',0));
S.E0 = cell2mat(cellfun(@(x) x.par(2),d.FitResults,'UniformOutput',0));
if strcmp(DataType,'Twin') % 2 fits didn't converge....
    S.chi2(21,11) = NaN;
    S.chi2(25,11) = NaN;
end

if strcmp(DataType,'Twin')
    BF = 'OFF';
else
    BF = 'ON';
end
    S.InterpMode = 'spline';
    S.Interp1Grid('Maxm4Sq',40^2);
S.ContourPlot('BestFit',BF,'SavePlot','OFF');

if S.chi2_Null<S.chi2_bf
    chi2min_bf       = S.chi2_Null;
    sin2T4_bf       = 0;
    mNu4Sq_bf       = 0;
else
    sin2T4_bf       = S.sin2T4_bf;
    mNu4Sq_bf       = S.mNu4Sq_bf;
    chi2min_bf       = S.chi2_bf;
end
sin2T4_contour  = S.sin2T4_contour;
mNu4Sq_contour  = S.mNu4Sq_contour;

% append contour
save(savename,'sin2T4_contour','mNu4Sq_contour','sin2T4_bf','mNu4Sq_bf','-append');

