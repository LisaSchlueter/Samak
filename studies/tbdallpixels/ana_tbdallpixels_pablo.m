function [ListGoodPixels] = ana_tbdallpixels_pablo(varargin)
%
% Read / Analyse Tritium Fit
% All Pixels
%
% T. Lasserre, 
% Last Updated: Feb. 12 2018
%

% Initialization
clear par ; clear parnomix;

% Parser
p = inputParser;
p.addParameter('Run',110,@(x)isfloat(x));
p.addParameter('fign',1,@(x)isfloat(x) && x>0);
p.addParameter('pub',@(x)ismember(x,{'ON','OFF','YES'}));
p.addParameter('display','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('nPixels',148,@(x)isfloat(x) && x<150);
p.addParameter('FPD','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('ReDrawSkeleton','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('file','bc_fitallpixels_run110_sim_148Pix.mat',@(x)ismember(x,{'tbd_allpixels_v1.mat','bc_fitallpixels_run0'}));
p.parse(varargin{:});

Run                      =    p.Results.Run;
fign                     =    p.Results.fign;
pub                      =    p.Results.pub;
display                  =    p.Results.display;

FPD                      =    p.Results.FPD;
ReDrawSkeleton           =    p.Results.ReDrawSkeleton;
file                     =    p.Results.file;

% Open
%strfile = sprintf('./results/bc_fitallpixels_run%g_G.mat',Run);
strfile = sprintf('./results/%s',file);
f = importdata(strfile);
% Bad Pixel List
%BadPixels  =  [10 77 100 106 107 130];
BadPixels  =  [100 130];
GoodPixelBinary = true(148);
switch display
    case 'ON'
        fprintf(2,'BAD Pixels - Begin ...\n');
        disp(num2str(BadPixels));
        fprintf(2,'BAD Pixels - End ...\n');
end

switch display
    case 'ON'
        fprintf(2,'Good Pixels Binary - Begin ...\n');
        disp(num2str(GoodPixelBinary'));
        fprintf(2,'Good Pixels Binary - End ...\n');
end
% Build Reduced List of Good Pixels - 148 Entries OR LESS
ListGoodPixels = ([1:1:148].*GoodPixelBinary'); ListGoodPixels(ListGoodPixels==0) = [];
switch display
    case 'ON'
        fprintf(2,'Good Pixels - Begin ...\n');
        disp(num2str(ListGoodPixels));
        fprintf(2,'Good Pixels - End ...\n');
end

%% Chi2
str_title = sprintf('\\chi^2 ');
MeanChi2 = mean((f(ListGoodPixels,14))');

figure(1)
% Histo
title(str_title);
nhist(f(ListGoodPixels,14),'binfactor',3,'color','qualitative','box','box'); grid on;
xlabel('\chi^2'); ylabel('Entries');
publish_figure(1,'figures/ana_tbdallpixels01.eps');
% Map
figure(2)
set(gcf, 'Position', [100, 100, 600, 500])
switch FPD
    case 'ON'
        if numel(BadPixels)>0
            FPDViewer(f(:,14)','ReDrawSkeleton',ReDrawSkeleton,'outliers',BadPixels);
        else
            FPDViewer(f(:,14)','ReDrawSkeleton',ReDrawSkeleton);
            title(str_title);
        end
        title(str_title);
end
publish_figure(2,'figures/ana_tbdallpixels02.eps');

%% NuMassSq
[MeanNuMassSq , MeanNuMassSqError ]   = swmean(f(ListGoodPixels,2),(1./f(ListGoodPixels,8).^2)); % .*escalefactor(ListGoodPixels)');
fprintf(2,'Weighted Mean  NuMassSq = %.3f +- %.3f eV ...\n',MeanNuMassSq,MeanNuMassSqError);

% Histo
figure(50)
subplot(2,2,1)
str_title = sprintf('NuMassSq (eV^2) ');
title(str_title);
nhist(f(ListGoodPixels,2),'binfactor',3,'color','qualitative','box','box'); grid on;
xlabel('NuMassSq (eV^2)'); ylabel('Entries');

% Map
figure(11)
set(gcf, 'Position', [100, 100, 600, 500])
switch FPD
    case 'ON'
        str_title = sprintf('NuMassSq (meV) ');
        if numel(BadPixels)>0
            FPDViewer(f(:,2)','ReDrawSkeleton',ReDrawSkeleton,'outliers',BadPixels);
        else
            FPDViewer(f(:,2)','ReDrawSkeleton',ReDrawSkeleton);
            title(str_title);
        end
end
publish_figure(11,'figures/ana_tbdallpixels04.eps');

% Map - Weigthed Mean
figure(12)
set(gcf, 'Position', [100, 100, 600, 500])
switch FPD
    case 'ON'
        str_title = sprintf('NuMassSq - %.2f eV^2 ',MeanNuMassSq);
        if numel(BadPixels)>0
            FPDViewer(f(:,2)' - MeanNuMassSq,'ReDrawSkeleton',ReDrawSkeleton,'outliers',BadPixels);
        else
            FPDViewer(f(:,2)' - MeanNuMassSq,'ReDrawSkeleton',ReDrawSkeleton);
            title(str_title);
        end
        title(str_title);
end
publish_figure(12,'figures/ana_tbdallpixels05.eps');

%% Endpoint
[MeanWidth , MeanWidthError ]       = swmean(f(ListGoodPixels,3),(1./f(ListGoodPixels,9).^2)); % .*escalefactor(ListGoodPixels)');
fprintf(2,'Weighted Mean Endpoint = %.3f +- %.3f eV ...\n',MeanWidth,MeanWidthError);

% Histo
figure(50)
subplot(2,2,2)
str_title = sprintf('Endpoint-18575 (meV) ');
title(str_title);
nhist(f(ListGoodPixels,3),'binfactor',3,'color','qualitative','box'); grid on;
xlabel('Endpoint (eV)'); ylabel('Entries');

% Map
figure(21)
set(gcf, 'Position', [100, 100, 600, 500])
switch FPD
    case 'ON'
        str_title = sprintf('Endpoint-18575 (meV) ');
        if numel(BadPixels)>0
            FPDViewer(f(:,3)','ReDrawSkeleton',ReDrawSkeleton,'outliers',BadPixels);
        else
            FPDViewer(f(:,3)','ReDrawSkeleton',ReDrawSkeleton);
            title(str_title);
        end
        title(str_title);
end
publish_figure(21,'figures/ana_tbdallpixels06.eps');


% Map - Weigthed Mean
figure(22)
set(gcf, 'Position', [100, 100, 600, 500])
switch FPD
    case 'ON'
        str_title = sprintf('Endpoint - %.2f eV ',MeanWidth);
        if numel(BadPixels)>0
            FPDViewer(f(:,3)' - MeanWidth,'ReDrawSkeleton',ReDrawSkeleton,'outliers',BadPixels);
        else
            FPDViewer(f(:,3)' - MeanWidth,'ReDrawSkeleton',ReDrawSkeleton);
        end
        title(str_title);
end
publish_figure(22,'figures/ana_tbdallpixels07.eps');

%% Background
[MeanBackground , MeanBackgroundError]         = swmean(f(ListGoodPixels,4),(1./f(ListGoodPixels,10).^2)); % .*escalefactor(ListGoodPixels)');
fprintf(2,'Weighted Mean Background = %.3f +- %.3f cps ...\n',MeanBackground,MeanBackgroundError);
%disp(num2str(f(ListGoodPixels,5)'));

% Histo
figure(50)
subplot(2,2,3)
str_title = sprintf('Mean Background (mcps) ');
title(str_title);
nhist(f(ListGoodPixels,4),'binfactor',3,'color','qualitative','box'); grid on;
xlabel('Background (mcps)'); ylabel('Entries');

% Map
figure(31)
set(gcf, 'Position', [100, 100, 600, 500])
switch FPD
    case 'ON'
        str_title = sprintf('Background (mcps) ');
        if numel(BadPixels)>0
            FPDViewer(f(:,4)','ReDrawSkeleton',ReDrawSkeleton,'outliers',BadPixels);
        else
            FPDViewer(f(:,4)','ReDrawSkeleton',ReDrawSkeleton);
            title(str_title);
        end
title(str_title);
end
publish_figure(31,'figures/ana_tbdallpixels08.eps');

% Map - Weigthed Mean
figure(32)
set(gcf, 'Position', [100, 100, 600, 500])
switch FPD
    case 'ON'
        str_title = sprintf('Background - %.2f mcps ',MeanBackground);
        if numel(BadPixels)>0
            FPDViewer(f(:,4)','ReDrawSkeleton',ReDrawSkeleton,'outliers',BadPixels);
        else
            FPDViewer(f(:,4)'- MeanBackground,'ReDrawSkeleton',ReDrawSkeleton);
            title(str_title);
        end
        title(str_title);
end
publish_figure(32,'figures/ana_tbdallpixels09.eps');

%% Normalization
[MeanNormalization , MeanNormalizationError]        = swmean(f(ListGoodPixels,5),(1./f(ListGoodPixels,11).^2)); % .*escalefactor(ListGoodPixels)');
fprintf(2,'Weighted Mean Normalization = %.2f +- %.2f cps ...\n',MeanNormalization,MeanNormalizationError);
%disp(num2str(f(ListGoodPixels,4)'));

% Histo
figure(50)
subplot(2,2,4)
str_title = sprintf('Mean Normalization ');
title(str_title);
nhist(f(ListGoodPixels,5),'binfactor',3,'color','qualitative','box'); grid on;
xlabel('Normalization (cps)'); ylabel('Entries');
publish_figure(50,'figures/ana_tbdallpixels03.eps');

% Map
figure(41)
set(gcf, 'Position', [100, 100, 600, 500])
switch FPD
    case 'ON'
        str_title = sprintf('Normalization ');
if numel(BadPixels)>0
            FPDViewer(f(:,5)','ReDrawSkeleton',ReDrawSkeleton,'outliers',BadPixels);
        else
            FPDViewer(f(:,5)'- MeanBackground,'ReDrawSkeleton',ReDrawSkeleton);
            title(str_title);
end
title(str_title);
end
publish_figure(41,'figures/ana_tbdallpixels10.eps');

% Map - Weigthed Mean
figure(42)
set(gcf, 'Position', [100, 100, 600, 500])
switch FPD
    case 'ON'
        str_title = sprintf('Normalization - %.2f cps ',MeanNormalization);
        if numel(BadPixels)>0
            FPDViewer(f(:,5)','ReDrawSkeleton',ReDrawSkeleton,'outliers',BadPixels);
        else
            FPDViewer(f(:,5)'- MeanNormalization,'ReDrawSkeleton',ReDrawSkeleton);
            title(str_title);
        end
        title(str_title);
end
publish_figure(42,'figures/ana_tbdallpixels11.eps');

%% Build pdf file with all fits
PATH = getenv('PATH');
setenv('PATH', [PATH ':/usr/local/bin/']);
cd figures
command = 'gs -sDEVICE=pdfwrite -sOutputFile="tbdfit-allpixels-fpd.pdf" -dNOPAUSE -dEPSCrop -c "<</Orientation 0>> setpagedevice" -f ana_tbdallpixels*.eps -c quit';
unix(command);
unix('rm ana_tbdallpixels*.eps');
mycommand1 = sprintf('mv tbdfit-allpixels-fpd.pdf bc_fitallpixels_run%g-fpdviewer.pdf',Run);
unix(mycommand1);
mycommand2 = sprintf('open bc_fitallpixels_run%g-fpdviewer.pdf',Run);
unix(mycommand2);
cd ..

% disp([ListGoodPixels' f(ListGoodPixels,2)-MeanNuMassSq f(ListGoodPixels,3)-MeanWidth f(ListGoodPixels,4)-MeanBackground])
% [maxP maxPI] = max(f(ListGoodPixels,2)-MeanNuMassSq);
% [minP minPI] = min(f(ListGoodPixels,2)-MeanNuMassSq);
% [maxW maxWI] = max(f(ListGoodPixels,3)-MeanWidth) ;
% [minW minWI] = min(f(ListGoodPixels,3)-MeanWidth) ;
% [maxA maxAI] = max(f(ListGoodPixels,4));
% [minA minAI] = min(f(ListGoodPixels,4));
% disp([ListGoodPixels(minPI) minP ListGoodPixels(maxPI) maxP] );
% disp([ListGoodPixels(minWI) minW ListGoodPixels(maxWI) maxW]);
% disp([ListGoodPixels(minAI) minA ListGoodPixels(maxAI) maxA]);
end
