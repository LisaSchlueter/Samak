addpath(genpath('../../../Samak2.0'));
close all;
clear;

% Select runs to analyze
%   For MultiRunAnalysis:
%   StackCD100all: All 100 % CD scans up/down
%   StackCD100up: All 100 % CD scans up
%   StackCD100down: All 100 % CD scans down
%   StackCD100random: All 100 % CD scans random
RunList = 'StackCD100all';

% Select Segmentation of the detector
% Options: StackPixel, SinglePixel, MultiPixel, Ring
Segmentation = 'StackPixel';

% Select ranges:
%   9 Short down to -200 eV from E0
%   7 Medium down to -400 eV from E0
%   1 Long down to -1600 eV from E0
ranges = [9,7,1];
Nranges = length(ranges);

% Select pixels (just neccessary for pixel segmentation)
pixellist = 1;

% Apply any pulls you wish t the fitted parameters
% The order is: mnu, E0, B, N, FSD Gr, FSD Ex
pulls    = [Inf,Inf,Inf,Inf,Inf,Inf];

% Choose fixed parameters
% '1 5 6' (usual choice): mnu, FSD Gr, and FSD E fixed
fixPar   = '1 5 6';

% Choose fitter
% Options: minuit, matlab
fitter = 'minuit';

% Choose type of fit
% Options: chi2Stat, chi2CM
chi2name = 'chi2CM';

% Allocate variables
FitResults = zeros(Nranges,10);


for range = 1:Nranges
    % Choose data range to analyze
    dataStart = ranges(range);
    
    MRA = MultiRunAnalysis('AnaFlag',Segmentation,...
        'chi2',chi2name,'RunList',RunList,...
        'fitter',fitter,'exclDataStart',dataStart,...
        'fixPar',fixPar,'pulls',pulls,...
        'DataEffcorr','OFF',...
        'EffCorrFlag','ON'); %'ROI+PileUp'
    
    MRA.Fit();
    
    mnu = MRA.FitResult.par(1);
    e0 = MRA.FitResult.par(2);
    bck = MRA.FitResult.par(3:3+length(pixellist)-1);
    norm = MRA.FitResult.par(3+length(pixellist):end-2);
    mnu_e = MRA.FitResult.err(1);
    e0_e = MRA.FitResult.err(2);
    bck_e = MRA.FitResult.err(3:3+length(pixellist)-1);
    norm_e = MRA.FitResult.err(3+length(pixellist):end-2);
    chi2 = MRA.FitResult.chi2min;
    dof = MRA.FitResult.dof;
    FitResults(range,:) = [mnu,e0,bck,norm,mnu_e,e0_e,bck_e,norm_e,chi2,dof];
    
end






