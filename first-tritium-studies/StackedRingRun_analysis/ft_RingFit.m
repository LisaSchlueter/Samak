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
runlist = 'StackCD100all';

% Select the rings you want to analyze, ex. [1,5,8], [4:10], etc.
ringlist = [1:2];

% Apply any pulls you wish (recommended for FT: 2 on mnu, none on the others)
pulls    = [Inf,Inf,Inf,Inf,Inf,Inf];

% Choose type of fit (chi2)
chi2name = 'chi2Stat'; 

% Choose data range to analyze
dataStart = 1; % 9 = -200 eV; 1 = -1600 eV

% Choose fixed parameters (optional)
fixPar   = '1 5 6'; 

% Choose fitter
fitter = 'minuit';

% Call Analysis class
MRA = MultiRunAnalysis('AnaFlag','Ring','RingList',ringlist,...
    'chi2',chi2name,'RunList',runlist,...
    'fitter',fitter,'exclDataStart',dataStart,...
    'fixPar',fixPar,'pulls',pulls,...
    'DataEffcorr','OFF'); %'ROI+PileUp'

% Do the fits
MRA.FitAllRings();

% Save the results for the plotting class
FitResult.par     = MRA.AllFitResults{1};
FitResult.err     = MRA.AllFitResults{2};
FitResult.chi2min = MRA.AllFitResults{3};
FitResult.dof     = MRA.AllFitResults{5};

% Call the plotting class and do the plots
P = PLOTC('Xdata',MRA.RunData.qU,'Ydata',MRA.RunData.TBDIS,...
    'ModelObj',MRA.ModelObj,'saveplot','OFF',...
    'titleFlag','Stack','RunList',MRA.RunList,'startqU',dataStart,...
    'FitResult',FitResult,'RingList',ringlist);
P.StackFileName = MRA.StackFileName;
P.plotRingIterative();









