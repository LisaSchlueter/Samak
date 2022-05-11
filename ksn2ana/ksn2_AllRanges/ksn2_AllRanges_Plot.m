% ksn2 calculate chi2 grid search for all ranges
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_AllRanges/results/'];
freePar = 'mNu E0 Norm Bkg';
savename = sprintf('%sksn2_Object_%s.mat',savedir,strrep(freePar,' ',''));
if exist(savename,'file')
    load(savename)
else
    %% settings that might change
    chi2 = 'chi2CMShape';
    DataType = 'Real';
    range = 40;
    
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
    Real = MultiRunAnalysis(RunAnaArg{:});
    %% configure Sterile analysis object
    
    
    %% define fit ranges
    RangeStandard = Real.GetexclDataStart(40);
    qU = round(Real.RunData.qU-18574);
    ranges = sort(round(-qU(1:RangeStandard)));
    
    if strcmp(DataType,'Real')
        if contains(freePar,'mNu')
            nGridSteps = 45;
            LoadGridArg = {'ExtmNu4Sq','OFF','mNu4SqTestGrid',5} ;
            ranges =[ranges(1), ranges(end)];
        else
            nGridSteps = 30;
            LoadGridArg = {'ExtmNu4Sq','ON','mNu4SqTestGrid',5} ;
        end
    else
        if contains(freePar,'mNu')
            nGridSteps = 40;
            ranges = ranges(end);
        else
            nGridSteps = 30;
        end
        LoadGridArg = {'ExtmNu4Sq','OFF','mNu4SqTestGrid',5} ;
    end
    
    SterileArg = {'RunAnaObj',Real,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',...
        'LoadGridArg',LoadGridArg,...
        'nGridSteps',nGridSteps };
    S = SterileAnalysis(SterileArg{:});
    MakeDir(savedir);
    save(savename,'Real','S','RunAnaArg','SterileArg','qU','ranges');
end

HoldOn = 'OFF';
pHandles = cell(numel(ranges),1);


if contains(freePar,'mNu')
    
    for i=1:numel(ranges)
        S.range = ranges(i);
        if ranges(i) == 90
            
            %% 90 eV range only
            S.LoadGridFile(S.LoadGridArg{:});
            S.Interp1Grid('Maxm4Sq',50^2);
            BestFit = S.FindBestFit;
            
            S.LoadGridFile(S.LoadGridArg{:});
            S.Interp1Grid;
            S.chi2_ref = BestFit.chi2min;
            S.ContourPlot('BestFit','ON','ReCalcBF','OFF','HoldOn','OFF','Color',S.PlotColors{1},'LineStyle',S.PlotLines{1});
            mNu4Sq_contour = S.mNu4Sq_contour;
            sin2T4_contour = S.sin2T4_contour;
            FitResults_Null = S.FitResults_Null;
            
            %
            S.LoadGridFile(S.LoadGridArg{:});
            S.Interp1Grid('Maxm4Sq',60^2,'Minm4Sq',150);
            BestFit_2ndMin = S.FindBestFit; % second (lower) minima
            
            savename_tmp = sprintf('%sksn2_AllRanges_chi2CMShape_%s_%.0feV.mat',savedir,strrep(freePar,' ',''),ranges(i));
            save(savename_tmp,'mNu4Sq_contour','sin2T4_contour','BestFit','BestFit_2ndMin','FitResults_Null');
            fprintf('save %s \n',savename_tmp);
        else
            
            S.LoadGridFile(S.LoadGridArg{:});
            S.Interp1Grid;
            BestFit = S.FindBestFit;
            S.ContourPlot('BestFit','ON','ReCalcBF','OFF','HoldOn','OFF','Color',S.PlotColors{1},'LineStyle',S.PlotLines{1});
            mNu4Sq_contour = S.mNu4Sq_contour;
            sin2T4_contour = S.sin2T4_contour;
            FitResults_Null = S.FitResults_Null;
            
            savename_tmp = sprintf('%sksn2_AllRanges_chi2CMShape_%s_%.0feV.mat',savedir,strrep(freePar,' ',''),ranges(i));
            save(savename_tmp,'mNu4Sq_contour','sin2T4_contour','BestFit','FitResults_Null');
            fprintf('save %s \n',savename_tmp);
        end
        
    end
else
    
    for i=1:numel(ranges)
        Real.exclDataStart = Real.GetexclDataStart(ranges(i));
        % configure Sterile analysis object
        
        S.range = ranges(i);
        S.LoadGridFile(S.LoadGridArg{:});
        S.Interp1Grid;
        
        S.ContourPlot('BestFit','ON','HoldOn',HoldOn,'Color',S.PlotColors{i},'LineStyle',S.PlotLines{i});
        HoldOn = 'ON';
        
        mNu4Sq_contour = S.mNu4Sq_contour;
        sin2T4_contour = S.sin2T4_contour;
        BestFit = S.FindBestFit;
        FitResults_Null = S.FitResults_Null;
        
        if ranges(i)==40 && ~contains(freePar,'mNu')
            S.LoadGridFile(S.LoadGridArg{:},'Extsin2T4','ON');
            S.Interp1Grid('Maxm4Sq',34^2);
            BestFit = S.FindBestFit;
        elseif ranges(i)==90 && ~contains(freePar,'mNu')
            S.LoadGridFile(S.LoadGridArg{:});
            S.Interp1Grid('Maxm4Sq',5000);
            S.ContourPlot('BestFit','ON','HoldOn','OFF','Color',S.PlotColors{2},'LineStyle',S.PlotLines{i});
            BestFit = S.FindBestFit;
            S.LoadGridFile(S.LoadGridArg{:});
            S.Interp1Grid;
            S.chi2_ref = BestFit.chi2min;
            S.ContourPlot('BestFit','ON','ReCalcBF','OFF','HoldOn',HoldOn,'Color',S.PlotColors{3},'LineStyle',S.PlotLines{i});
            mNu4Sq_contour = S.mNu4Sq_contour;
            sin2T4_contour = S.sin2T4_contour;
            
            
        end
        
        savename_tmp = sprintf('%sksn2_AllRanges_chi2CMShape_%s_%.0feV.mat',savedir,strrep(freePar,' ',''),ranges(i));
        save(savename_tmp,'mNu4Sq_contour','sin2T4_contour','BestFit','FitResults_Null');
        fprintf('save %s \n',savename_tmp);
        
    end
    
end

%%
% S.RunAnaObj.fixPar = ConvertFixPar('freePar','mNu E0 Bkg Norm','nPar',S.RunAnaObj.nPar);
% S.LoadGridFile(S.LoadGridArg{:});
%     S.Interp1Grid;
%     S.ContourPlot('BestFit','ON','HoldOn',HoldOn,'Color',S.PlotColors{i},'LineStyle',S.PlotLines{i});
