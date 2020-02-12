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



% Select runs to analyze (from 40667 to 40693 are 3h FT runs)
RunList = 40667;
% Select the rings you want to analyze, ex. [1,5,8], [4:10], etc.
%PixelList = [1:124];
PixelList = [1:100,102:112,114:123];


% Choose data range to analyze
ranges = [9,7,1]; % 9 = -200 eV; 1 = -1600 eV
Nranges = length(ranges);

% Choose fixed parameters (optional)
fixPar   = '1 5 6';

% Choose fitter
fitter = 'minuit';

for uncertainty = 1:1
    % Choose type of fit (chi2)
    switch uncertainty
        case 1
            chi2name = 'chi2P';
        case 2
            chi2name = 'chi2CM';
    end
    
    for range = 1:Nranges
        
        dataStart = ranges(range);
        
        % Call Analysis class
        MRA = MultiRunAnalysis('AnaFlag','SinglePixel','PixList',PixelList,...
            'chi2',chi2name,'RunList',RunList,...
            'fitter',fitter,'exclDataStart',dataStart,...
            'fixPar',fixPar,...
            'DataEffcorr','OFF'); %'ROI+PileUp'
        
        % Do the fits
        MRA.FitAllSinglePixels();
        
        % Save the results for the plotting class
        FitResult.par     = MRA.AllFitResults{1};
        FitResult.err     = MRA.AllFitResults{2};
        FitResult.chi2min = MRA.AllFitResults{3};
        FitResult.dof     = MRA.AllFitResults{5};
        
        GeneralTable = [FitResult.par,FitResult.err,FitResult.chi2min,FitResult.dof,PixelList'];
        
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

save('data/R40667SinglePix');




