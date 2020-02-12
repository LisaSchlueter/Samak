%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ring fits analysis
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

% Select the rings you want to analyze
ringlist = [1:13];

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
MRA = MultiRunAnalysis('AnaFlag','Ring','ringlist',ringlist,...
    'chi2',chi2name,'RunList',runlist,...
    'fitter',fitter,'exclDataStart',9,'fixPar',fixPar,'pulls',pulls);

% Do the fits
MRA.FitAllRings();

% Do the plots
