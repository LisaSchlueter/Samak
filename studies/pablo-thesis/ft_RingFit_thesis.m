%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ring fit analysis
% ---------------------------------------------------------------------- %
% This study shows the results of the fitted parameters of the KATRIN
% First Tritium Data ring-wise, for stacked runs
%
% Pablo I. Morales Guzm�n  TUM/MPP
% Last update : 06/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../../../Samak2.0'));

    % Select runs to analyze (from 40667 to 40693 are 3h FT runs
    RunList = 'StackCD100all';
    % Select the rings you want to analyze, ex. [1,5,8], [4:10], etc.
    ringlist = [1:10];

    % Apply any pulls you wish (recommended for FT: 2 on mnu, none on the others)
    pulls    = [Inf,Inf,Inf,Inf,Inf,Inf];
    
    % Choose data range to analyze
    ranges = [9,7]; % 9 = -200 eV; 1 = -1600 eV
    Nranges = length(ranges);
    
    % Choose fixed parameters (optional)
    fixPar   = '1 5 6';
    
    % Choose fitter
    fitter = 'minuit';
%load('data/StackedRingsStackedRunsResultsThesis2');

for uncertainty = 1:2
    % Choose type of fit (chi2)
    switch uncertainty
        case 1
            chi2name = 'chi2Stat';
        case 2
            chi2name = 'chi2CM';
    end
    
    for range = 1:Nranges
        
        dataStart = ranges(range);
        
        % Call Analysis class
        MRA = MultiRunAnalysis('AnaFlag','Ring','RingList',ringlist,...
            'chi2',chi2name,'RunList',RunList,...
            'fitter',fitter,'exclDataStart',dataStart,...
            'fixPar',fixPar,'pulls',pulls,...
            'DataEffcorr','OFF',...
            'EffCorrFlag','ON'); %'ROI+PileUp'
        
        % Do the fits
        MRA.FitAllRings();
        
        % Save the results for the plotting class
        FitResult.par     = MRA.AllFitResults{1};
        FitResult.err     = MRA.AllFitResults{2};
        FitResult.chi2min = MRA.AllFitResults{3};
        FitResult.dof     = MRA.AllFitResults{5};
        
        GeneralTable = [FitResult.par,FitResult.err,FitResult.chi2min,FitResult.dof];
        
        switch uncertainty
            case 1 % Statistical
                switch range
                    case 1
                        TableShortStat = GeneralTable;
                        
                    case 2
                        TableMedStat = GeneralTable;
                        
                    case 3
                        TableLongStat = GeneralTable;
                        
                end
            case 2 % systematics with CovMat
                switch range
                    case 1
                        TableShortSys = GeneralTable;
                        
                    case 2
                        TableMedSys = GeneralTable;
                        
                    case 3
                        TableLongSys = GenerlaTable;
                end
        end
        
    end
    
end

save('data/StackedRingsStackedRunsResultsThesis2');









