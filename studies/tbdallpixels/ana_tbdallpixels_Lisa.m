function [ListGoodPixels] = ana_tbdallpixels_Lisa(varargin)
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
p.addParameter('Run',3541035422,@(x)isfloat(x));
p.addParameter('fign',1,@(x)isfloat(x) && x>0);
p.addParameter('pub',@(x)ismember(x,{'ON','OFF','YES'}));
p.addParameter('display','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('nPixels',148,@(x)isfloat(x) && x<150);
p.addParameter('FPD','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('ReDrawSkeleton','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('file', 'bc_fitallpixels_run3541035422_P_fix12.mat',@(x)ismember(x,{'tbd_allpixels_v1.mat','bc_fitallpixels_run0'}));
p.parse(varargin{:});

Run                      =    p.Results.Run;
fign                     =    p.Results.fign;
pub                      =    p.Results.pub;
display                  =    p.Results.display;
display                  =    p.Results.display;
FPD                      =    p.Results.FPD;
ReDrawSkeleton           =    p.Results.ReDrawSkeleton;
file                     =    p.Results.file;

% Open
strfile = sprintf('./results/%s',file);
f = importdata(strfile);
chi2tmp = f(:,14); 
%Exclude bad Pixels
for i = 148:-1:1
    if abs(chi2tmp(i))>1000
    fprintf('Excluding Pixel %u \n', i);
    f(i,:)= NaN;
    end
end
chi2 = f(:,14);
goodPixels = f(:,1); %pixel number
%par = f(:,2:5); %parameter values ( 2 numass, 3 Endpoint, 4 N, 5 B)
par = f(:,4:5); %parameter values (4 N, 5 B)
parerr = f(:,10:11); %uncertainties 
ndof = (f(1,15)+4)*6+2;
%ndof =f(1,15)+2;
%% Chi2
str_title = sprintf('\\chi^2 (%u ndof)',ndof);
MeanChi2 = nanmean(chi2);

figure(1)
% Histo
title(str_title);
nhist(chi2,'binfactor',3,'color','qualitative','box','box'); grid on;
xlabel('\chi^2'); ylabel('Entries');
publish_figure(1,'figures/ana_tbdallpixels01.eps');
% Map
figure(2)
set(gcf, 'Position', [100, 100, 600, 500])
switch FPD
    case 'ON'
      FPDViewer(chi2,'ReDrawSkeleton',ReDrawSkeleton);
      title(str_title);
end
publish_figure(2,'figures/ana_tbdallpixels02.eps');

%% Background
%load('./data/background_run35410.mat');
d0 = importdata('./data/Background/BackgroundRate_35410.txt');
BKGData = d0.data(:,2);
d1 = importdata('./data/Background/BackgroundRate_35410.txt');
BKGData = BKGData + d1.data(:,2);
d2 = importdata('./data/Background/BackgroundRate_35410.txt');
BKGData = BKGData + d2.data(:,2);
d3 = importdata('./data/Background/BackgroundRate_35410.txt');
BKGData = BKGData + d3.data(:,2);
d4 = importdata('./data/Background/BackgroundRate_35410.txt');
BKGData = BKGData + d4.data(:,2);
d5 = importdata('./data/Background/BackgroundRate_35410.txt');
BKGData = BKGData + d5.data(:,2);
BKGData = BKGData./6;

%exclude NaN to calculate weighted mean
Background_tmp =par(:,2); 
BackgroundErr_tmp = parerr(:,2);
Background_tmp(isnan(Background_tmp)) = []; 
BackgroundErr_tmp(isnan(BackgroundErr_tmp)) = [];
[MeanBackground , MeanBackgroundError]= swmean(Background_tmp, (1./BackgroundErr_tmp.^2)); % .*escalefactor(ListGoodPixels)');
fprintf(2,'Weighted Mean Background = %.3f +- %.3f cps ...\n',MeanBackground,MeanBackgroundError);

% Histo
figure(50)
subplot(2,1,1)
str_title = sprintf('Mean Background (mcps) ');
title(str_title);
nhist(par(:,2),'binfactor',3,'color','qualitative','box'); grid on;
xlabel('Background (mcps)'); ylabel('Entries');

% Map
figure(31)
set(gcf, 'Position', [100, 100, 600, 500])
switch FPD
    case 'ON'
        str_title = sprintf('Background (mcps) ');
        FPDViewer(par(:,2),'ReDrawSkeleton',ReDrawSkeleton);
        title(str_title);
end
publish_figure(31,'figures/ana_tbdallpixels08.eps');

figure(33)
set(gcf, 'Position', [100, 100, 600, 500])
switch FPD
    case 'ON'
        str_title = sprintf('Background Δ(Fit-Data) (mcps)  ');
        FPDViewer(par(:,2)-BKGData*1e3,'ReDrawSkeleton',ReDrawSkeleton);
        title(str_title);
end
publish_figure(33,'figures/ana_tbdallpixels12.eps');



% Map - Weigthed Mean
figure(32)
set(gcf, 'Position', [100, 100, 600, 500])
switch FPD
    case 'ON'
        str_title = sprintf('Background - %.3f (mcps) ',MeanBackground);
        FPDViewer(par(:,2)-MeanBackground,'ReDrawSkeleton',ReDrawSkeleton);
        title(str_title);
end
publish_figure(32,'figures/ana_tbdallpixels09.eps');

%% Normalization

%exclude NaN to calculate weighted mean
Norm_tmp =par(:,1); 
NormErr_tmp = parerr(:,1);
Norm_tmp(isnan(Norm_tmp)) = []; 
NormErr_tmp(isnan(NormErr_tmp)) = [];
[MeanNormalization , MeanNormalizationError]= ...
    swmean(Norm_tmp, (1./NormErr_tmp.^2));
fprintf(2,'Weighted Mean Normalization = %.2f +- %.2f cps ...\n',MeanNormalization,MeanNormalizationError);
%disp(num2str(f(ListGoodPixels,4)'));

% Histo
figure(50)
subplot(2,1,2)
str_title = sprintf('Mean Normalization ');
title(str_title);
nhist(par(:,1),'binfactor',4,'color','qualitative','box'); grid on;
xlabel('Normalization (cps)'); ylabel('Entries');
publish_figure(50,'figures/ana_tbdallpixels03.eps');

% Map
figure(41)
set(gcf, 'Position', [100, 100, 600, 500])
switch FPD
    case 'ON'
        str_title = sprintf('Normalization ');
        FPDViewer(par(:,1),'ReDrawSkeleton',ReDrawSkeleton);
        title(str_title);
end
publish_figure(41,'figures/ana_tbdallpixels10.eps'); % -MeanBackground

% Map - Weigthed Mean
figure(42)
set(gcf, 'Position', [100, 100, 600, 500])
switch FPD
    case 'ON'
        str_title = sprintf('Normalization - %.3f cps ',MeanNormalization);
        FPDViewer(par(:,1)- MeanNormalization,'ReDrawSkeleton',ReDrawSkeleton);
        title(str_title);
        
end
publish_figure(42,'figures/ana_tbdallpixels11.eps');

%% NuMassSq
% [MeanNuMassSq , MeanNuMassSqError ]   = swmean(f(ListGoodPixels,2),(1./f(ListGoodPixels,8).^2)); % .*escalefactor(ListGoodPixels)');
% fprintf(2,'Weighted Mean  NuMassSq = %.3f +- %.3f eV ...\n',MeanNuMassSq,MeanNuMassSqError);
% 
% % Histo
% figure(50)
% subplot(2,2,1)
% str_title = sprintf('NuMassSq (eV^2) ');
% title(str_title);
% nhist(f(ListGoodPixels,2),'binfactor',3,'color','qualitative','box','box'); grid on;
% xlabel('NuMassSq (eV^2)'); ylabel('Entries');
% 
% % Map
% figure(11)
% set(gcf, 'Position', [100, 100, 600, 500])
% switch FPD
%     case 'ON'
%         str_title = sprintf('NuMassSq (meV) ');
%         if numel(BadPixels)>0
%             FPDViewer(f(:,2)','ReDrawSkeleton',ReDrawSkeleton,'outliers',BadPixels);
%         else
%             FPDViewer(f(:,2)','ReDrawSkeleton',ReDrawSkeleton);
%             title(str_title);
%         end
% end
% publish_figure(11,'figures/ana_tbdallpixels04.eps');
% 
% % Map - Weigthed Mean
% figure(12)
% set(gcf, 'Position', [100, 100, 600, 500])
% switch FPD
%     case 'ON'
%         str_title = sprintf('NuMassSq - %.2f eV^2 ',MeanNuMassSq);
%         if numel(BadPixels)>0
%             FPDViewer(f(:,2)' - MeanNuMassSq,'ReDrawSkeleton',ReDrawSkeleton,'outliers',BadPixels);
%         else
%             FPDViewer(f(:,2)' - MeanNuMassSq,'ReDrawSkeleton',ReDrawSkeleton);
%             title(str_title);
%         end
%         title(str_title);
% end
% publish_figure(12,'figures/ana_tbdallpixels05.eps');

%% Endpoint
% [MeanWidth , MeanWidthError ]       = swmean(f(ListGoodPixels,3),(1./f(ListGoodPixels,9).^2)); % .*escalefactor(ListGoodPixels)');
% fprintf(2,'Weighted Mean Endpoint = %.3f +- %.3f eV ...\n',MeanWidth,MeanWidthError);
% 
% % Histo
% figure(50)
% subplot(2,2,2)
% str_title = sprintf('Endpoint-18575 (meV) ');
% title(str_title);
% nhist(f(ListGoodPixels,3),'binfactor',3,'color','qualitative','box'); grid on;
% xlabel('Endpoint (eV)'); ylabel('Entries');
% 
% % Map
% figure(21)
% set(gcf, 'Position', [100, 100, 600, 500])
% switch FPD
%     case 'ON'
%         str_title = sprintf('Endpoint-18575 (meV) ');
%         if numel(BadPixels)>0
%             FPDViewer(f(:,3)','ReDrawSkeleton',ReDrawSkeleton,'outliers',BadPixels);
%         else
%             FPDViewer(f(:,3)','ReDrawSkeleton',ReDrawSkeleton);
%             title(str_title);
%         end
%         title(str_title);
% end
% publish_figure(21,'figures/ana_tbdallpixels06.eps');
% 
% 
% % Map - Weigthed Mean
% figure(22)
% set(gcf, 'Position', [100, 100, 600, 500])
% switch FPD
%     case 'ON'
%         str_title = sprintf('Endpoint - %.2f eV ',MeanWidth);
%         if numel(BadPixels)>0
%             FPDViewer(f(:,3)' - MeanWidth,'ReDrawSkeleton',ReDrawSkeleton,'outliers',BadPixels);
%         else
%             FPDViewer(f(:,3)' - MeanWidth,'ReDrawSkeleton',ReDrawSkeleton);
%         end
%         title(str_title);
% end
% publish_figure(22,'figures/ana_tbdallpixels07.eps');

%% Build pdf file with all fits
PATH = getenv('PATH');
setenv('PATH', [PATH ':/usr/local/bin/']);
cd figures
command = 'gs -sDEVICE=pdfwrite -sOutputFile="tbdfit-allpixels-fpd.pdf" -dNOPAUSE -dEPSCrop -c "<</Orientation 0>> setpagedevice" -f ana_tbdallpixels*.eps -c quit';
unix(command);
unix('rm ana_tbdallpixels*.eps');
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
