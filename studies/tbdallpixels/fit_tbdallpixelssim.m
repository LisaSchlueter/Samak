function [par,err,chi2min,ndof] = fit_tbdallpixelssim(varargin)
%
%            FIT INTEGRAL SPECTRUM
%                   Tritium
%          All Pixels Simultaneously
%
%          Backgrounds are fixed parameters
%
%  - Background fixed with pixel-wise fit : 148 par. fixed
%  - A common fit for: E0 and NuMassSq
%  - nPixels free normalizations
%
% Th. Lasserre - CEA Saclay
% Last Updated: Feb. 16 2018               
%

% Initialization
clear par ; clear parnomix;
tic;
% Parser
p = inputParser;
p.addParameter('Run',110,@(x)isfloat(x));
p.addParameter('MyMode','Simulation',@(x)ismember(x,{'Bootcamp','Simulation'}));
p.addParameter('pub','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('StatFluct','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('nPixels',148,@(x)isfloat(x) && x<150);
p.addParameter('FixNorms','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('Backgrounds','READ',@(x)ismember(x,{'FIT','READ'}));
p.addParameter('mnuSq_t',0.,@(x)isfloat(x));
p.addParameter('FPDView','OFF',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});
Run         =    p.Results.Run;
MyMode      =    p.Results.MyMode;
pub         =    p.Results.pub;
display     =    p.Results.display;
StatFluct   =    p.Results.StatFluct;
nPixels     =    p.Results.nPixels;
FixNorms    =    p.Results.FixNorms;
Backgrounds =    p.Results.Backgrounds;
mnuSq_t     =    p.Results.mnuSq_t;
FPDView     =    p.Results.FPDView;

% Parametrization: True Value
switch MyMode
    case 'Bootcamp'
        switch Run
            case 0
                A = init_bootcampR0('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat0.mat');
                figprefix = 'bc_fitallsim_run0';
            case 1
                A = init_bootcampR1('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat1.mat');
                figprefix = 'bc_fitallsim_run1';
            case 2
                A = init_bootcampR2('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat2.mat');
                figprefix = 'bc_fitallsim_run2';
            case 3
                A = init_bootcampR3('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat3.mat');
                figprefix = 'bc_fitallsim_run3';
            case 4
                A = init_bootcampR4('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat4.mat');
                figprefix = 'bc_fitallsim_run4';
            case 5
                A = init_bootcampR5('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat5.mat');
                figprefix = 'bc_fitallsim_run5';
            case 6
                A = init_bootcampR6('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat6.mat');
                figprefix = 'bc_fitallsim_run6';
            case 7
                A = init_bootcampR7('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat7.mat');
                figprefix = 'bc_fitallsim_run7';
            case 8
                A = init_bootcampR8('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat8.mat');
                figprefix = 'bc_fitallsim_run8';
            case 9
                A = init_bootcampR9('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat9.mat');
                figprefix = 'bc_fitallsim_run9';
            case 10
                A = init_bootcampR10('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat10.mat');
                figprefix = 'bc_fitallsim_run10';
            case 110
                A = init_bootcampR110('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat1_10.mat');
                figprefix = 'bc_fitallsim_run1_10';
            case 35410
                A = init_bootcampR35410('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat35410.mat');
                figprefix = 'bc_fitallsim_run35410';
            case 35411
                A = init_bootcampR35411('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat35411.mat');
                figprefix = 'bc_fitallsim_run35411';
            case 35412
                A = init_bootcampR35412('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat35412.mat');
                figprefix = 'bc_fitallsim_run35412';
            case 35413
                A = init_bootcampR35413('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat35413.mat');
                figprefix = 'bc_fitallsim_run35413';
            case 35420
                A = init_bootcampR35420('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat35420.mat');
                figprefix = 'bc_fitallsim_run35420';
            case 35422
                A = init_bootcampR35422('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat35422.mat');
                figprefix = 'bc_fitallsim_run35422';
            case 3541035422
                A = init_bootcampR3541035422('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat3541035422.mat');
                figprefix = 'bc_fitallsim_run35410_35422';
        end
        A.TimeSec = D.TimeSec; A.ComputeNormFactorTBDDS;
    case 'Simulation'
        A = init_tbdallpixels('Pixel',0,'nPixels',nPixels);
end
A.mnuSq_i = mnuSq_t;

% Init
switch MyMode
    case 'Bootcamp'
        quIndex        = 1:3:A.nPixels*3;
        countIndex     = quIndex + 1;
        coutErrIndex   = quIndex + 2;
%         for i=1:1:A.nPixels
%             A.qUpixel(i,:) = D.TBDISallPixels(:,quIndex(i));
%         end
        A.qUpixel(1:A.nPixels,:) = D.TBDISallPixels(:,quIndex(1:A.nPixels))'; 
        switch display
            case 'ON'
                fprintf(2,'----------------------------------------\n');
                fprintf(2,'qU Matrix - %g Pixels\n',A.nPixels);
                disp(num2str(A.qUpixel(1:A.nPixels,:)));
                fprintf(2,'----------------------------------------\n');
        end
    case 'Simulation'
        A.qUpixel(1:A.nPixels,:) = repmat(A.qU',A.nPixels,1);
        
end

% Fit Backgrounds
switch Backgrounds
    case 'FIT'
        [~ , tmp] = fit_tbdallpixels(...
            'Run',Run,'display','ON','displayPlot','OFF',...
            'nPixels',nPixels,'StatFluct','OFF','mnuSq_t',mnuSq_t);
        switch MyMode
            case 'Simulation'
                A.BKG_RateSecallPixels_i = A.BKG_RateAllFPDSec./A.nPixels+tmp;
            otherwise
                A.BKG_RateSecallPixels_i = tmp;
        end
    case 'READ'
        fitr = load(sprintf('results/bc_fitallpixels_run%g_G.mat',Run));
        switch MyMode
            case 'Simulation'
                A.BKG_RateSecallPixels_i = A.BKG_RateAllFPDSec./A.nPixels+fitr.ParnPixels(:,4)';
            otherwise
                A.BKG_RateSecallPixels_i = fitr.ParnPixels(:,4)';
        end
end

switch display
    case 'ON'
        fprintf(2,'----------------------------------------\n');
        fprintf(2,'Background Vector - %g Pixels\n',A.nPixels);
        disp(num2str(A.BKG_RateSecallPixels_i(1:A.nPixels)));
        fprintf(2,'----------------------------------------\n');
end
 
% Init Normalization
fitr = load(sprintf('results/bc_fitallpixels_run%g_G.mat',Run));
A.BKG_RateSecallPixels_i = fitr.ParnPixels(1:A.nPixels,4)';

% Init
i_mnuSq    = A.mnuSq_i;%(0.2*randn)^2;
i_E0       = 0;%0.05*randn;
i_N        = fitr.ParnPixels(1:A.nPixels,5)';

% Data
switch MyMode
    case 'Bootcamp'
        Data = [...
            reshape(D.TBDISallPixels(:,quIndex(1:A.nPixels)),1,A.nqU*A.nPixels)',...
            reshape(D.TBDISallPixels(:,countIndex(1:A.nPixels)),1,A.nqU*A.nPixels)',...
            reshape(D.TBDISallPixels(:,coutErrIndex(1:A.nPixels)),1,A.nqU*A.nPixels)'];
        DataL = {Data,A};
    case 'Simulation'
        A.ComputeTBDDSallPixels(); A.ComputeTBDISallPixels();
        switch StatFluct
            case 'ON'
                AddStatFluctTBDISallPixels(A);
        end
        Data = [...
            reshape(repmat(A.qU,1,A.nPixels),1,A.nqU*A.nPixels)' , ...
            reshape((A.TBDISallPixels),1,A.nqU*A.nPixels)', ...
            reshape(sqrt((A.TBDISallPixels)),1,A.nqU*A.nPixels)'];
        DataL = {Data,A};
end
switch display
    case 'ON'
        fprintf(2,'----------------------------------------\n');
        fprintf(2,'%s Data For Fit\n',MyMode);
        disp(num2str(Data));
        fprintf(2,'----------------------------------------\n');
end

% Initializing Fit
ParIni    = [i_mnuSq i_E0 i_N];


% Build Name Parameter List
parnames = ['mSq E0'];
parnamescor = [sprintf('''mSq'',''E0''')];
for k=1:1:A.nPixels
tmp         = sprintf(' N%g',k); 
parnames    = [parnames tmp];
tmpcor      = sprintf(',''N%g''',k);
parnamescor = [parnamescor tmpcor];
end

% Build Fixed Parameters Mist
switch FixNorms
    case 'ON'
        FixParStr = [];
        FixPar    = [3:A.nPixels+2]; % Normalization
        for l=3:1:A.nPixels+2
            FixParStr = [FixParStr sprintf('; fix %g',l)];
        end
        npar       = 2;
    case 'OFF'
        npar       = 2+A.nPixels;
        FixParStr = [];
end

tmparg = ['set pri 1;' FixParStr ' ; min'];
Args = {ParIni, DataL,'-n', parnames,'-c',tmparg};

switch display
    case 'ON'
        fprintf(2,'----------------------------------------\n');
        fprintf(2,'Minuit Fitting Options \n');
        disp(Args);
        fprintf(2,'----------------------------------------\n');
end

[par, err, chi2min, errmat] = pminuit('chi2_tbd_allpixelssim',Args);

switch display
    case 'ON'
        fprintf(2,'----------------------------------------\n');
        fprintf(2,'Data / Fit: Spectra \n');
        disp([ Data(:,2) model_tbdallpixelssim(par,A)]);
        fprintf(2,'----------------------------------------\n');
end

switch display
    case 'ON'
        fprintf(2,'--------------------------------------------------------------\n');
        fprintf('  Processing ....\n');
        fprintf(2,'--------------------------------------------------------------\n');
        fprintf(2,'  mnuSq \t= %.3f \t \t %g \t eV^2 \n', (i_mnuSq+par(1)),err(1));
        fprintf(2,'  E0 \t= %.3f \t \t %g \t eV \n',  (i_E0+par(2)),err(2));
        for i=1:1:A.nPixels
        fprintf(2,'  N%g \t= %.3f \t \t %g \t  \n',  i,  par(2+i),err(2+i));
        end
        ndof = A.nqU*A.nPixels-npar;
        fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,ndof);
        fprintf(2,'--------------------------------------------------------------\n');
end

%% Save Results
strtmp = sprintf('./results/bc_fitallpixels_run%g_sim_%gPix.mat',Run,nPixels);
save(strtmp,'par','err','chi2min','errmat');

%% Plot Results
figure(1)
subplot(2,1,1)
hdata = errorbar(Data(1:A.nqU*A.nPixels,1),Data(1:A.nqU*A.nPixels,2),Data(1:A.nqU*A.nPixels,3),...
    'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
hold on
model = model_tbdallpixelssim(par,A);

for (i=0:1:A.nPixels-1)
    hfit1 = plot(Data(i*A.nqU+1:i*A.nqU+A.nqU,1),model(i*A.nqU+1:i*A.nqU+A.nqU),...
        'Color','Red','LineWidth',1,'LineStyle','-');
end

hold off
grid on
xlabel('qU (eV)','FontSize',10);
ylabel('Counts','FontSize',10);
title(sprintf('KATRIN TBD - %s - Pixel 1 to %.0f - Doppler %s',A.TD,A.nPixels,A.DopplerEffectFlag));
set(gca,'FontSize',12);
set(gca,'yscale','lin');
mydata = sprintf('Data: NuMassSq=%.4f eV^2 - E0=%.4f meV \n',A.mnuSq_i,A.Q);
myfit = sprintf('Fit: mSq=%.4f \\pm %.4f eV - E0=%.4f \\pm %.4f eV',(i_mnuSq+par(1)),err(1),i_E0+par(2),err(2));
mychi2 = sprintf('\\chi2 / dof=%.1f/%.0f\n',chi2min,A.nqU*A.nPixels-2);
legend([hdata  hfit1],mychi2,myfit,'Location','NorthEast') ;% legend(a,'boxoff');
axis([min(A.qU) max(A.qU)+1 0.7*min(Data(:,2)) max(Data(:,2))*1.2])

subplot(2,1,2)
hdata = errorbar(Data(1:A.nqU*A.nPixels,1),Data(1:A.nqU*A.nPixels,2)-model(1:A.nqU*A.nPixels),Data(1:A.nqU*A.nPixels,3),...
    'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
hold on
line([Data(1,1) Data(end,1)],[0 0 ],'LineStyle','--','Color','Red');
hold off
grid on
xlabel('qU (eV)','FontSize',10);
ylabel('Residuals','FontSize',10);
set(gca,'FontSize',12);
set(gca,'yscale','lin');
%axis([min(A.qU) max(A.qU)+1 min(Data(1:A.nqU*A.nPixels,2)-model(1:A.nqU*A.nPixels))*2 max(Data(1:A.nqU*A.nPixels,2)-model(1:A.nqU*A.nPixels))*2])
%publish_figure(1,'figures/fit_tbdallpixelssim1.eps');

figure(2)
title(sprintf('Correlation Matrix of %g Fitted Parameters',npar));
        [h cor] = corplot(errmat(1:numel(ParIni),1:numel(ParIni)));
%         xticklabels({'mSq','E0','N1','N2','N3'})
%         yticklabels({'mSq','E0','N1','N2','N3'})
%        disp(parnamescor);
%publish_figure(2,'figures/fit_tbdallpixelssim2.eps');

figure(3)
title(sprintf('TBD Integral Spectra and Fit -  %g Pixels',A.nPixels));
hold on
[x, y] = meshgrid(A.qU,1:1:A.nPixels);
Rfit  = reshape(model(1:A.nqU*A.nPixels),A.nqU,A.nPixels);
Rdata = reshape(Data(1:A.nqU*A.nPixels,2),A.nqU,A.nPixels);
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

figure(4)
title(sprintf('TBD Integral Spectra and Fit Residuals -  %g Pixels',A.nPixels));
hold on
[x, y] = meshgrid(A.qU,1:1:A.nPixels);
Rfit     = reshape(model(1:A.nqU*A.nPixels),A.nqU,A.nPixels);
Rdata    = reshape(Data(1:A.nqU*A.nPixels,2),A.nqU,A.nPixels);
RdataErr = reshape(Data(1:A.nqU*A.nPixels,3),A.nqU,A.nPixels);
ribbon(Rfit-Rfit);colormap('summer')
hold on
stem3((Rdata-Rfit)./RdataErr,'Marker','s',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'LineStyle','none');
hold off
ylabel('qU')
xlabel('Pixel')
zlabel('Integral Spectrum Fit Residuals')
hold off
view([180 0])
%publish_figure(4,'figures/fit_tbdallpixelssim4.eps');

switch pub
    case 'ON'
        %% Build pdf file with all fits
        PATH = getenv('PATH');
        setenv('PATH', [PATH ':/usr/local/bin/']);
        cd figures
        command = 'gs -sDEVICE=pdfwrite -sOutputFile="fit_tbdallpixelssim.pdf" -dNOPAUSE -dEPSCrop -c "<</Orientation 0>> setpagedevice" -f fit_tbdallpixelssim*.eps -c quit';
        unix(command);
        unix('rm fit_tbdallpixelssim*.eps');
        cd ..
        unix('open figures/fit_tbdallpixelssim.pdf');
        close all
end
toc;
end


