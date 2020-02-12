%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Endpoint sensitivity study
% Detector as 1 pixel equivalent
% Configuration for First Tritium May

% P. Morales 2018
% Last update 16/03/2018


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initialization
clear;
addpath(genpath('../../../Samak2.0'));
%Parameters
nfit            = 1;
runtimes        = (2)*60*60*24;
%Ntimes          = length(runtimes);
Ntimes = 5;

nTeBinFac       = 5;
time_dist       = 'Flat30';
Fitter          = 'Matlab';
nPixels         = 30;

FSD             = 'DOSS';
FSDadd = ''; FSDplot = 'only DT';

BCK             = 'XmasData';

Mode            = 'Read';

display         = 'OFF';
FitFlag         = 'ON';
savespectrum    = 'OFF';
plotresults     = 'OFF';
plotdist        = 'OFF';
saveresults     = 'OFF';
savedist        = 'OFF';


% Parametrization: True Value
mnuSq_t = (0.00)^2;

% Init
opt_bin = {...
    'nTeBinningFactor',nTeBinFac};

opt_wgtsmace = {...
    'KTFFlag','Compute'}; %MACE+WGTSIS %SSCW_DCOfficial %SSCW_DCperqU

opt_bkg = {...
    'BKG_Flag',BCK};

opt_fsd= {'DTFSD',FSD};

%Initialize matrices
fit_p = zeros(nfit,10);
muhat = zeros(Ntimes,1);
sigmahat = zeros(Ntimes,1);
fit_xi2min = zeros(nfit,1);

opt_katrin = {...
    'TD',time_dist,...
    'TimeSec',runtimes,...
    'Mode',Mode};

A = InitKATRINE0_allPixels(...
    opt_katrin{:},...
    opt_wgtsmace{:},...
    opt_fsd{:},...
    opt_bkg{:},...
    opt_bin{:},...
    'mnuSq_i',mnuSq_t,...
    'nPixels',nPixels);

    A.ComputeTBDDSallPixels();
    A.ComputeTBDISallPixels();
    model = A.TBDISallPixels;
    
    A.AddStatFluctTBDISallPixels();
    Data = [...
        reshape(repmat(A.qU,1,A.nPixels),1,A.nqU*A.nPixels)' , ...
        reshape((A.TBDISallPixels),1,A.nqU*A.nPixels)', ...
        reshape(sqrt((A.TBDISallPixels)),1,A.nqU*A.nPixels)'];
    
    figure(3)
    title(sprintf('TBD Integral Spectra and Fit -  %g Pixels',A.nPixels));
    hold on
    Rfit  = reshape(model(1:A.nqU*A.nPixels),A.nqU,A.nPixels);
    Rdata    = reshape(Data(1:A.nqU*A.nPixels,2),A.nqU,A.nPixels);
RdataErr = reshape(Data(1:A.nqU*A.nPixels,3),A.nqU,A.nPixels);
    ribbon(Rfit);colormap('summer')
    hold on
stem3(Rdata,'Marker','s',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'LineStyle','none');
hold off
ylabel('qU')
xlabel('Pixel')
zlabel('Integral Spectrum Fit')
hold off
view([-100 45])
%publish_figure(3,'figures/fit_tbdallpixelssim3.eps');

OptionZ.Duration = 5;
OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'WellMadeVid',OptionZ)