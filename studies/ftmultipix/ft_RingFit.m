%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ring fit analysis
% ---------------------------------------------------------------------- %
% This study shows the results of the fitted parameters of the KATRIN
% First Tritium Data ring-wise, for stacked runs
%
% Pablo I. Morales Guzmán  TUM/MPP
% Last update : 06/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../../../Samak2.0'));

% Select runs to analyze (from 40667 to 40693 are 3h FT runs)
runlist = [40667,40668,40669,40670,40671,40672,40673,40674,40675,40676,...
    40677,40678,40679,40680,40681,40682,40683,40684,40685,40686,40687,...
    40688,40689,40690,40691,40692,40693];

% Select the rings you want to analyze, ex. [1,5,8], [4:10], etc.
ringlist = [5:5];

% Apply any pulls you wish (recommended for FT: 2 on mnu, none on the others)
pulls = [2^2,Inf,Inf,Inf];

% Choose type of fit (chi2)
chi2name = 'chi2Stat'; 

% Choose fixed parameters (optional)
fixPar = '';

% Choose fitter
fitter = 'minuit';

% (fixed parameters only work with minuit)
if ~isempty(fixPar); fitter = 'minuit'; end

% Call Analysis class 
MRA = MultiRunAnalysis('AnaFlag','Ring','RingList',ringlist,...
    'chi2',chi2name,'RunList',runlist,...
    'fitter',fitter,'exclDataStart',9,'fixPar',fixPar,'pulls',pulls);

% Do the fits
 MRA.FitAllRings();

% Save the results for the plotting class
FitResult.par = MRA.AllFitResults{1};
FitResult.err = MRA.AllFitResults{2};
FitResult.chi2min = MRA.AllFitResults{3};
FitResult.dof = MRA.AllFitResults{5};

% Call the plotting class and do the plots
P = PLOTC('Xdata',MRA.RunData.qU,'Ydata',MRA.RunData.TBDIS,...
    'ModelObj',MRA.ModelObj,'CovMat',MRA.FitCM,'saveplot','export',...
    'titleFlag','Stack','RunList',MRA.RunList,'startqU',9,...
    'FitResult',FitResult,'RingList',ringlist);

P.plotRingIterative();















