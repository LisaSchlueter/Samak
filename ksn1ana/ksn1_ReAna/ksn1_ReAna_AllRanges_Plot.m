% ksn2 calculate chi2 grid search
%% settings that might change
DataType              = 'Real';
chi2                  = 'chi2CMShape';
freePar               = 'mNu E0 Norm Bkg';

savedir = [getenv('SamakPath'),'ksn1ana/ksn1_ReAna/results/'];
savename = sprintf('%sksn1_Object_%s.mat',savedir,strrep(freePar,' ',''));
if exist(savename,'file')
    load(savename)
else
    
    if contains(freePar,'mNu')
        nGridSteps = 40;
    else
        nGridSteps = 30;
    end
    
    if strcmp(chi2,'chi2Stat')
        NonPoissonScaleFactor = 1;
    elseif  strcmp(chi2,'chi2CMShape')
        NonPoissonScaleFactor = 1.064;
    end
    
    RunAnaArg = {'RunList','KNM1',...
        'chi2',chi2,...
        'DataType',DataType,...
        'fixPar',freePar,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'SysBudget',200,...
        'minuitOpt','min ; minos',...
        'FSDFlag','KNM2_0p1eV',...
        'ELossFlag','KatrinT2A20',...
        'AngularTFFlag','ON',...
        'SynchrotronFlag','ON',...
        'RadiativeFlag','ON',...
        'DopplerEffectFlag','FSD',...
        'BKG_PtSlope',-2.2*1e-06};
    
    
    %% configure RunAnalysis object
    Real = MultiRunAnalysis(RunAnaArg{:});
    
    % define fit ranges
    RangeStandard = Real.GetexclDataStart(40);
    qU = round(Real.RunData.qU-18575);
    ranges = sort(round(-qU(1:RangeStandard)));
    if contains(freePar,'mNu')
        ranges = ranges(end);
    end
    
    SterileArg = {'RunAnaObj',Real,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'nGridSteps',nGridSteps,...
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',...
        'range',ranges(1)};
    
    S = SterileAnalysis(SterileArg{:});
    S.LoadGridArg = {'ExtmNu4Sq','OFF','mNu4SqTestGrid',2};
    save(savename,'Real','S','RunAnaArg','SterileArg','qU','ranges');
    
end
HoldOn = 'OFF';
pHandles = cell(numel(ranges),1);


if contains(freePar,'mNu')
    
    S.range = 40;
    S.LoadGridFile(S.LoadGridArg{:});
    S.Interp1Grid('Maxm4Sq',34^2);
    BestFit = S.FindBestFit;
    S.LoadGridFile(S.LoadGridArg{:});
    S.Interp1Grid('Maxm4Sq',36^2);
    S.chi2_ref = BestFit.chi2min;
    S.ContourPlot('BestFit','ON','ReCalcBF','OFF','HoldOn','OFF','Color',S.PlotColors{1},'LineStyle',S.PlotLines{1});
    
    mNu4Sq_contour = S.mNu4Sq_contour;
    sin2T4_contour = S.sin2T4_contour;
    FitResults_Null = S.FitResults_Null;
    savename_tmp = sprintf('%sksn1_AllRanges_chi2CMShape_%s_%.0feV.mat',savedir,strrep(freePar,' ',''),S.range);
    save(savename_tmp,'mNu4Sq_contour','sin2T4_contour','BestFit','FitResults_Null');
    fprintf('save %s \n',savename_tmp);
    
    S.range = 94;
    %% 90 eV range only
    %     S.LoadGridFile(S.LoadGridArg{:});
    %     S.Interp1Grid('Maxm4Sq',50^2);
    %     BestFit = S.FindBestFit;
    %% find first minimum and contour
    S.LoadGridFile(S.LoadGridArg{:});
    S.Interp1Grid;
    BestFit = S.FindBestFit;
    S.ContourPlot('BestFit','ON','HoldOn','OFF','Color',S.PlotColors{1},'LineStyle',S.PlotLines{1});
    mNu4Sq_contour = S.mNu4Sq_contour;
    sin2T4_contour = S.sin2T4_contour;
    FitResults_Null = S.FitResults_Null;
    
    % get 99% C.L. 
    S.ConfLevel = 99;
     S.ContourPlot('BestFit','ON','HoldOn','ON','Color',S.PlotColors{2},'LineStyle',S.PlotLines{2});
     mNu4Sq_contour99 = S.mNu4Sq_contour;
    sin2T4_contour99 = S.sin2T4_contour;
    S.ConfLevel = 95;
  
    % get also position of 2nd minimum
    S.LoadGridFile(S.LoadGridArg{:});
    S.Interp1Grid('Maxm4Sq',2000);
   % S.ContourPlot('BestFit','ON','HoldOn','ON','Color',S.PlotColors{2},'LineStyle',S.PlotLines{2});
    BestFit_2 = S.FindBestFit; % second (lower) minimum
    plot(BestFit_2.sin2T4,BestFit_2.mNu4Sq,'p','MarkerSize',15,'MarkerFaceColor',S.PlotColors{1},'Color',S.PlotColors{1});
    
    savename_tmp = sprintf('%sksn1_AllRanges_chi2CMShape_%s_%.0feV.mat',savedir,strrep(freePar,' ',''),ranges(1));
    save(savename_tmp,'mNu4Sq_contour','sin2T4_contour','BestFit','FitResults_Null',...
       'BestFit_2','mNu4Sq_contour99','sin2T4_contour99');
    fprintf('save %s \n',savename_tmp);
    
    
else
    
    for i=1:numel(ranges)
         if ranges(i)==94  
             S.nGridSteps = 45;
         end
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
        
        
        savename_tmp = sprintf('%sksn1_AllRanges_chi2CMShape_%s_%.0feV.mat',savedir,strrep(freePar,' ',''),ranges(i));
        save(savename_tmp,'mNu4Sq_contour','sin2T4_contour','BestFit','FitResults_Null');
        fprintf('save %s \n',savename_tmp);
        
        if ranges(i)==94  
            %% 94 eV range: plot also 99% C.L., because of closed contour 
            S.ConfLevel = 99;
            S.ContourPlot('BestFit','ON','HoldOn','OFF','Color',S.PlotColors{2},'LineStyle',S.PlotLines{2});
            mNu4Sq_contour99 = S.mNu4Sq_contour;
            sin2T4_contour99 = S.sin2T4_contour;
           
            
            % get second minimum
            S.LoadGridFile(S.LoadGridArg{:});
            S.Interp1Grid('Maxm4Sq',1e3);
            S.ContourPlot('BestFit','ON','HoldOn','ON','Color',S.PlotColors{t},'LineStyle',S.PlotLines{t});
            BestFit_2 = S.FindBestFit;
            
            
            save(savename_tmp,'mNu4Sq_contour99','sin2T4_contour99','BestFit_2','-append');
            fprintf('append 99%% exclusion to  %s \n',savename_tmp);
            
            S.nGridSteps = 30; % reset
             S.ConfLevel = 95;
        end
    end
    
end


