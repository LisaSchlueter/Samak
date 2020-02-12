function [N , B] = fit_tbdallpixels(varargin)
%
% Fit TBD all Pixels Simulation
% 148 pixels  
% 
% Th. Lasserre - CEA Saclay
% Last Updated: Feb. 12 2018
%

close all

% Parser
p = inputParser;
p.addParameter('Run',110,@(x)isfloat(x));
p.addParameter('MyMode','Simulation',@(x)ismember(x,{'Bootcamp','Simulation'}));
p.addParameter('Chi2Type','G',@(x)ismember(x,{'G','P'}));
p.addParameter('pub','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('displayPlot','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('StatFluct','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('nPixels',148,@(x)isfloat(x) && x<150);
p.addParameter('mnuSq_t',(0)^2,@(x)isfloat(x));

p.parse(varargin{:});
Run         =    p.Results.Run;
MyMode      =    p.Results.MyMode;
Chi2Type    =    p.Results.Chi2Type;
pub         =    p.Results.pub;
display     =    p.Results.display;
displayPlot =    p.Results.displayPlot;
StatFluct   =    p.Results.StatFluct;
nPixels     =    p.Results.nPixels;
mnuSq_t     =    p.Results.mnuSq_t;

ParnPixels = zeros(nPixels,15);
%progressbar('...Samak TBD Fit Pixel-Wise Backgrounds ...');
for p = 1:1:nPixels
    %progressbar(p/nPixels);
    fprintf(2,'--------------------------------------------------------------\n')
    fprintf(2,'Run %g - Fitting Pixel %0.2f \n',Run,p)
    [parA, err, chi2min, ndof ] =  fit_tbdpixel(...
        'Run',Run,'MyMode',MyMode,'mnuSq_t',mnuSq_t,...
        'Fitter','Minuit','Pixel',p,'StatFluct',StatFluct,...
        'display',display,'displayPlot',displayPlot,...
        'Chi2Type',Chi2Type);
    ParnPixels(p,:) = [p parA err chi2min ndof];
end

% Output for further fit
B = ParnPixels(:,5)';
N = ParnPixels(:,4)';

% Save Results
strtmp = sprintf('./results/bc_fitallpixels_run%g_%s_fix12.mat',Run,Chi2Type);
save(strtmp,'ParnPixels');

%% Build pdf file with all fits
switch pub
    case 'ON'
        PATH = getenv('PATH');
        setenv('PATH', [PATH ':/usr/local/bin/']);
        cd figures
        ls bc_fit_run*.eps;
        command = 'gs -sDEVICE=pdfwrite -sOutputFile="fittbdallpixels.pdf" -dNOPAUSE -dEPSCrop -c "<</Orientation 2>> setpagedevice" -f bc_fit_run*.eps -c quit';
        unix(command);
        unix('rm bc_fit_run*.eps');
        mycommand1 = sprintf('mv fittbdallpixels.pdf bc_fitallpixels_run%g.pdf',Run);
        unix(mycommand1);
        mycommand2 = sprintf('open bc_fitallpixels_run%g_%s.pdf',Run, Chi2Type);
        unix(mycommand2); 
        cd ..
end

end


