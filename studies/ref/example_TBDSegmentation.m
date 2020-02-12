%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic examples on how to use the different segmentations of FPD in Samak
% Notes 14/05/2018: 
% Only simpsons integration work
%
% TO DO 
% - Include all theoretical corrections in CumFrac filename
%
%
% Pablo I. Morales G. TUM/MPP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../../../Samak2.0'));
clear;

% REMEMBER: only simspons integration works for now !

%% No segmentation

% Calculates average magnetic field in analyzing plane
% Full column density
% Full background


% Just choose option of FPD_Segmentation to OFF
FPD_Segmentation_No = 'OFF';

% Create object
A_NoSeg = TBD('FPD_Segmentation',FPD_Segmentation_No,'FPD_Pixel',100,'IStype','SIMPFAST');
% (It does not matter what you write in FPD_Pixel, samak will change it to 1
% and give a warning)

% Calculate the differential spectrum
A_NoSeg.ComputeTBDDS();

% Calculate integral spectrum
A_NoSeg.ComputeTBDIS();

%% Multi Pixel

% Calculates magnetic field in analyzing plane for the specific pixel
% Pixel-wise column density
% Pixel-wise background

% Create object
% Just choose option of FPD_Segmentation to MULTIPIXEL
FPD_Segmentation_MP = 'MULTIPIXEL';

% Give the pixels you want to analyze as a list
% Ex: FPD_Pixel = 1:148 (all pixels)
% Ex: FPD_Pixel = 23:56
% Ex: FPD_Pixel = 5
% Ex: FPD_Pixel = [1,5,8,23,45,67,90,111,147]

FPD_Pixel_MP = [1,5,8,23,24,25,26,27,28,45,67,90,111,147];

% Background can be included as XmasData or FLAT
% BKG_Flag = 'XmasData' reads a file with the info 
% BKG_Flag = 'ON' and 'BKG_Type' = FLAT assigns the same amount of background to each pixel
% Optionally give the background for whole detector and it will be
% distributed within the pixels

opt_bkg_MP = {...
    'BKG_Flag','ON',...
    'BKG_Type','FLAT',...
    'BKG_RateAllFPDSec',300e-3}; % 300 mcps in this case, so 300/148 for each pixel

A_MP = TBD('FPD_Segmentation',FPD_Segmentation_MP,'FPD_Pixel',FPD_Pixel_MP,opt_bkg_MP{:},'IStype','SIMPFAST');


% Calculate the differential spectrum (done to all pixels selected, ComputeTBDDSallPixels no longer necessary)
A_MP.ComputeTBDDS();

% Calculate integral spectrum (done to all pixels selected, ComputeTBDDSallPixels no longer necessary)
A_MP.ComputeTBDIS();

%% Single Pixel
% Calculates magnetic field in analyzing plane for the specific pixel
% Pixel-wise column density
% Pixel-wise background

% Single pixel is basically the same as using multipixel and giving just
% one pixel in FPD_Pixel. The difference is that it makes sure just a
% single pixel is used
% (It does not matter what you write in FPD_Pixel, samak will change it to 1
% and give a warning)
% If you are doing single pixel analysis, best practice would be to change
% this option to SINGLEPIXEL

FPD_Segmentation_SP = 'SINGLEPIXEL';

FPD_Pixel = 139;

opt_bkg_SP = {...
    'BKG_Flag','ON',...
    'BKG_Type','FLAT',...
    'BKG_RateAllFPDSec',300e-3}; % 300 mcps in this case, so 300/148 for the pixel selected

A_SP = TBD('FPD_Segmentation',FPD_Segmentation_SP,'FPD_Pixel',FPD_Pixel,opt_bkg_SP{:},'IStype','SIMPFAST');


%% (Optional) for fit:
wannafit = true;
if wannafit
% First add fluctuations to the integral spectrum
A_MP.AddStatFluctTBDIS();

% Assign data
Data = [A_MP.qU(:), A_MP.TBDIS(:), A_MP.TBDISE(:)];


% Initialize parameters (VERY IMPORTANT! for multipixel fit, a bad initialization could cause the fitter not to converge)
i_mnu      = 0;
i_Q        = 0;
i_B        = A_MP.TBDIS(end,:)/(A_MP.qUfrac(end)*A_MP.TimeSec);
i_N        = zeros(1,A_MP.nPixels);

DoF = A_MP.nqU*A_MP.nPixels - length(i_B) - length(i_N) - 2;

ParIni = [i_mnu i_Q i_B i_N];
        
         
% Use a different class without the fluctuation to fit (put background to 0)
A_MPfit = TBD('FPD_Segmentation',FPD_Segmentation_MP,'FPD_Pixel',FPD_Pixel_MP,'IStype','SIMPFAST',opt_bkg_MP{:},'BKG_RateAllFPDSec',0);
A_MPfit.ComputeTBDDS();
A_MPfit.ComputeTBDIS();

% Fit with Matlab
fprintf('----------BEGIN FIT-------------- \n');
options = optimoptions('fminunc','Algorithm','quasi-newton',...
    'OptimalityTolerance',1e-7,'StepTolerance',1e-8,...
    'FiniteDifferenceType','central',...
    'MaxFunctionEvaluations',1e6,'UseParallel',true,...
    'Display','iter');
DataTBD = {Data,A_MPfit};
TBDfun = @(xx) Chi2example(xx,DataTBD);
[par,chi2min,exitflag,output,grad,Hessian] = fminunc(TBDfun,ParIni,options);
errmat = 0.5*Hessian;
varcov = inv(errmat);
err = sqrt(diag(varcov))';
norms_fit = par(3+A_MPfit.nPixels:3+2*A_MPfit.nPixels-1);
bcks_fit = par(3:2+A_MPfit.nPixels);

% Display fir results

mnuSq_report = A_MPfit.mnuSq_i+par(1);
mnu_report = sqrt(abs(mnuSq_report));
err_mnu = sqrt(abs((mnuSq_report + err(1)))) - mnu_report;
BCK_fit = sum(bcks_fit);
BCK_err = sqrt(sum(bcks_fit.^2));
N_ave = mean(abs(norms_fit));
N_err = sqrt(sum(norms_fit.^2))/A_MPfit.nPixels;
fprintf(2,'----------------------------------------------\n');
fprintf('===============================================\n');
fprintf('  m^2       = %g +/- %g eV^2\n', mnuSq_report,err(1));
fprintf('  m         = %g +/- %g eV\n', mnu_report,err_mnu);
fprintf('  dQ        = %g +/- %g eV\n', par(2),err(2));
fprintf('  B total   = %g +/- %g cps\n', BCK_fit,BCK_err);
fprintf('  dN average= %g +/- %g\n', N_ave,N_err);
fprintf('  Chi2/dof  = %g/%g\n', chi2min, DoF);
fprintf('===============================================\n');
 
% Plot data from fit

figure(1)
pixel_norms = NaN(148,1);
pixel_norms(A_MPfit.FPD_Pixel) = norms_fit;
FPDViewer(pixel_norms,'ReDrawSkeleton','ON');
title(['Normalization bias TEST']);

figure(2)
pixel_bcks = NaN(148,1);
pixel_bcks(A_MPfit.FPD_Pixel) = bcks_fit;
FPDViewer(pixel_bcks,'ReDrawSkeleton','ON');
title(['Background TEST']);

end


function chi2 = Chi2example(p,DATA)
DAT = DATA{1};
Afit = DATA{2};

y = DAT(:,2);    % bin content of non-null bins
z = DAT(:,3);   % uncertainties bin content of non-null bins

Afit.ComputeTBDDS(...
    'mSq_bias',p(1),...
    'E0_bias',p(2),...
    'N_bias',p(Afit.nPixels+3:Afit.nPixels+3+Afit.nPixels-1),...
    'B_bias',p(3:Afit.nPixels+2));
% p(3:A.nPixels+2)p(A.nPixels+3:A.nPixels+3+A.nPixels-1)
Afit.ComputeTBDIS();

m = Afit.TBDIS(:);

% Poisson deviance chi square.
chi2 = sum(( y - m ).^2 ./ z.^2);
end




