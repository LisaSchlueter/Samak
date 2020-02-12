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
ringlist = [1:4];

% Apply any pulls you wish (recommended for FT: 2 on mnu, none on the others)
pulls = [2^2,Inf,Inf,Inf,inf,inf];

% Choose data range to analyze
dataStart = 1; % 9 = -200 eV

% Choose type of fit (chi2)
chi2name = 'chi2Stat'; 

% Choose fixed parameters (optional)
fixPar = '1 5 6';

% Choose fitter
fitter = 'minuit';

% (fixed parameters only work with minuit)
if ~isempty(fixPar); fitter = 'minuit'; end

% Memory allocation for variables in for loop
saverhoDmin = zeros(length(ringlist),1);
saverhoDUnc_low = zeros(length(ringlist),1);
saverhoDUnc_up = zeros(length(ringlist),1);
plotname = cell(length(ringlist),1);

% Call plot class for the saving function 
P = PLOTC('saveplot','export');

% Do the rho_d scans, ring by ring
for ri = 1:length(ringlist)
    % Call the class with only the selected ring
    MRA = MultiRunAnalysis('AnaFlag','Ring','ring',ringlist(ri),...
    'chi2',chi2name,'RunList',runlist,...
    'fitter',fitter,'exclDataStart',dataStart,'fixPar',fixPar,'pulls',pulls);

    
    % Do the rho_d scan and save the relevant results
    MRA.RhoDScan();
    saverhoDmin(ri) = MRA.CDmin;
    saverhoDUnc_low(ri) = MRA.rhoDlowerUnc;
    saverhoDUnc_up(ri) = MRA.rhoDupperUnc;
    
    % Save the plot
    P.savename = sprintf('RhoDScan_ring%.0d.pdf',ringlist(ri));
    plotname{ri} = ['plots/',P.savename];
    P.savePlotFun();
    close all
end

% Call the plotting class and do the plots
P = PLOTC('saveplot','export','RunList',MRA.RunList,...
    'titleFlag','RhoDScanRing',...
    'RingList',ringlist);

P.plotRingRhoDScan(saverhoDmin,saverhoDUnc_low,saverhoDUnc_up,plotname);














