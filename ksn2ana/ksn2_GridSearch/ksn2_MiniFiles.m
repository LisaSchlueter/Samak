% save contours analysis into mini files


chi2 = 'chi2CMShape';
DataType = 'Real';
range = 40;
freePar = 'mNu E0 Norm Bkg';

savedir = [getenv('SamakPath'),'SterileAnalysis/MiniFiles/KSN2/'];
MakeDir(savedir);
savefile = sprintf('%sKSN2contour_%s_%s_%s_%.0feV.mat',savedir,DataType,strrep(freePar,' ',''),chi2,range);

if exist(savefile,'file')
else
    %% configure RunAnalysis object
    if strcmp(chi2,'chi2Stat')
        NonPoissonScaleFactor = 1;
    elseif  strcmp(chi2,'chi2CMShape')
        NonPoissonScaleFactor = 1.112;
    end
    
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar',freePar,...%free par
        'SysBudget',40,...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'RadiativeFlag','ON',...
        'FSDFlag','KNM2_0p5eV',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'chi2',chi2,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'TwinBias_Q',18573.7,...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',3*1e-06,...
        'TwinBias_BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
    A = MultiRunAnalysis(RunAnaArg{:});
    %% configure Sterile analysis object
    SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',...
        'range',range,...
        'LoadGridArg',{'ExtmNu4Sq','ON','mNu4SqTestGrid',5}};
    
    %%
    S = SterileAnalysis(SterileArg{:});
    S.InterpMode = 'spline';
    %% load grid file + interpolate + cotour and best fit
    
    if contains(freePar,'mNu') && strcmp(DataType,'Real')
        S.nGridSteps = 50;
        S.LoadGridFile('mNu4SqTestGrid',5);
    elseif contains(freePar,'mNu') && strcmp(DataType,'Twin')
        S.nGridSteps = 40;
        S.LoadGridFile('mNu4SqTestGrid',5);
    elseif ~contains(freePar,'mNu') && strcmp(DataType,'Twin')
        S.nGridSteps = 30;
        S.LoadGridFile('mNu4SqTestGrid',5);
    elseif ~contains(freePar,'mNu') && strcmp(DataType,'Real')
        S.nGridSteps = 30;
        S.LoadGridFile('ExtmNu4Sq','ON','mNu4SqTestGrid',5,'Extsin2T4','ON');
        S.InterpMode = 'Mix';
    end
    
    S.Interp1Grid;
    S.ContourPlot('BestFit','ON','SavePlot','OFF');
    sin2T4_contour  =  S.sin2T4_contour;
    mNu4Sq_contour  =  S.mNu4Sq_contour;
    FitResults_Null = S.FitResults_Null;
    
    if strcmp(DataType,'Real')
        [DeltaChi2, SignificanceCL,SignificanceSigma] = S.CompareBestFitNull;
        sin2T4_bf      = S.sin2T4_bf;
        mNu4Sq_bf      =  S.mNu4Sq_bf;
        mNuSq_bf       =  S.mNuSq_bf;
        chi2_bf        = S.chi2_bf;
    end
    
    save(savefile,'sin2T4_contour','mNu4Sq_contour','FitResults_Null');
    if strcmp(DataType,'Real')
        save(savefile,'sin2T4_bf','mNu4Sq_bf','mNuSq_bf','chi2_bf',....
            'SignificanceCL','SignificanceSigma','DeltaChi2','-append')
    end
    fprintf('save to file %s \n',savefile);
end
