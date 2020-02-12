function [ListGoodPixels] = KrLineReadAllPixels(varargin)
%
% Read / Analyse KrLine Fit
% All Pixels
%
% T. Lasserre, August 2017
%

% Initialization
clear par ; clear parnomix;
addpath('fminuit'); addpath('krypton-data'); addpath('tools'); %digits(6);

% Parser
p = inputParser;
p.addParameter('fign',1,@(x)isfloat(x) && x>0);
p.addParameter('pub',@(x)ismember(x,{'ON','OFF','YES'}));
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('FPD','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('file','Data_KrL332HSAllPixelsDopplerONRipplesON_v2.mat',...
@(x)ismember(x,{'Data_KrL332AllPixelsDopplerONRipplesON_v2.mat',...
'Data_KrL332AllPixelsDopplerONRipplesON.mat','Data_KrL332HSAllPixelsDopplerONRipplesON.mat'}));
p.parse(varargin{:});

fign       =    p.Results.fign;
pub        =    p.Results.pub;
display    =    p.Results.display;
FPD        =    p.Results.FPD;
file       =    p.Results.file;

% Variables
Pixels     = ones(1,148);

% Open
f = importdata(file);


% Find the Line Type
tmp = strfind(file,'KrL332');
if tmp == 6
    %BadPixels  = [123 124 125 126 128 129 134 135 136 137 138 139 140 141 147]+1; % L3-32 Line
    %BadPixels  =  [20 28 37 74 99 100 111 112 113 114 123 124 125 126 127 128 129 130 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148]; % L3-32 Line
    %HS
    BadPixels  =  [23 100 112 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148];
    TD         = 'KrL3_32'; strline = 'L3-32';
    fprintf(2,'Krypton L3-32 Line...\n');
end
tmp = strfind(file,'KrK32');
if tmp  >= 5
    % BadPixels  = [6 14 26 27 34 40 56 68 81 95 97 99 100 103 109 111 112 113 114 123 124 125 126 127 128 129 130 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148];  % K32 Line
    BadPixels  = [6 26 34 36 38 40 41 56 64 68 69 70 71 80 81 84 85 91 92 93 94 95 97 99 100 101 103 109 110 111 112 113 114 117 122 123 124 125 126 127 128 129 130 131 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148];
    TD         = 'KrK32'; strline = 'K32';
    fprintf(2,'Krypton K-32 Line ...\n');
end
if tmp == []
    fprintf(2,'Krypton Line Not Found - Exit...\n');
    return;
end
switch display
    case 'ON'
        fprintf(2,'BAD Pixels - Begin ...\n');
        disp(num2str(BadPixels));
        fprintf(2,'BAD Pixels - End ...\n');
end

% Doppler
tmp = strfind(file,'DopplerON');  if tmp>0  fprintf(2,'Doppler Effect ON ...\n'); end
tmp = strfind(file,'DopplerOFF'); if tmp>0  fprintf(2,'Doppler Effect OFF ...\n'); end

% HV Ripples
tmp = strfind(file,'RipplesON');  if tmp>0  fprintf(2,'HV Ripples ON  ...\n'); end
tmp = strfind(file,'RipplesOFF'); if tmp>0  fprintf(2,'HV Ripples OFF ...\n'); end

% Correct Energy Scale
% A=InitKrKATRIN('TD',TD);  escalefactor = A.Pixel_MACE_Ba_VCorr;

% Build List of Good Pixels - 148 Entries
GoodPixelBinary = Pixels; GoodPixelBinary(BadPixels)=0; GoodPixelBinary = GoodPixelBinary';
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
str_title = sprintf('\\chi^2 - Data - Krypton %s Line - Good Pixels',strline);

MeanChi2 = mean((f(ListGoodPixels,11))');

% Histo
figure(1)
title(str_title);
nhist(f(ListGoodPixels,11),'binfactor',3,'color','qualitative'); grid on;
xlabel('\chi^2'); ylabel('Entries');

% Map
switch FPD
    case 'ON'
FPDViewer(f(:,11)','outliers',BadPixels);
title(str_title);
end

%% Line Position
f(:,1) = f(:,1) * (1972.4531/1972.5);
[MeanPosition  MeanPositionError ]   = swmean(f(ListGoodPixels,1),(1./f(ListGoodPixels,6).^2)); % .*escalefactor(ListGoodPixels)');
%MeanPositionError = 0; % .*escalefactor(ListGoodPixels)');
fprintf(2,'%s Weighted Mean Line Position = %.3f +- %.3f eV ...\n',strline,MeanPosition,MeanPositionError);

% Histo
figure(50)
str_title = sprintf('Line Position - Data - Krypton %s Line - Good Pixels',strline);
title(str_title);
nhist(f(ListGoodPixels,1),'binfactor',3,'color','qualitative'); grid on;
xlabel('Line Position (eV)'); ylabel('Entries');

% Map
switch FPD
    case 'ON'
        str_title = sprintf('Line Position (eV) - Data - Krypton %s Line - Good Pixels',strline);
FPDViewer(f(:,1)','outliers',BadPixels);
title(str_title);
end

% Map - Weigthed Mean
switch FPD
    case 'ON'
        str_title = sprintf('Line Position - %.2f eV - Data - Krypton %s Line - Good Pixels',MeanPosition,strline);
FPDViewer(f(:,1)' - MeanPosition,'outliers',BadPixels);
title(str_title);
end

%% Line Width
%MeanWidth         = mean(f(ListGoodPixels,2)); % .*escalefactor(ListGoodPixels)');
[MeanWidth MeanWidthError ]       = swmean(f(ListGoodPixels,2),(1./f(ListGoodPixels,7).^2)); % .*escalefactor(ListGoodPixels)');
%MeanWidthError    = 0; % .*escalefactor(ListGoodPixels)');
fprintf(2,'%s Weighted Mean Line Width = %.3f +- %.3f eV ...\n',strline,MeanWidth,MeanWidthError);

% Histo
figure(100)
str_title = sprintf('Line Width (eV) - Data - Krypton %s Line - Good Pixels',strline);
title(str_title);
nhist(f(ListGoodPixels,2),'binfactor',3,'color','qualitative'); grid on;
xlabel('Line Width (eV)'); ylabel('Entries');

% Map
switch FPD
    case 'ON'
        str_title = sprintf('Line Width (eV) - Data - Krypton %s Line - Good Pixels',strline);
FPDViewer(f(:,2)','outliers',BadPixels);
title(str_title);
end


% Map - Weigthed Mean
switch FPD
    case 'ON'
        str_title = sprintf('Line Width - %.2f eV - Data - Krypton %s Line - Good Pixels',MeanWidth,strline);
FPDViewer(f(:,2)' - MeanWidth,'outliers',BadPixels);
title(str_title);
end

%% Line Amplitude
[MeanAmplitude MeanAmplitudeError]         = swmean(f(ListGoodPixels,3),(1./f(ListGoodPixels,8).^2)); % .*escalefactor(ListGoodPixels)');
%MeanAmplitudeError     = 0; % .*escalefactor(ListGoodPixels)');
fprintf(2,'%s Weighted Mean Amplitude = %.3f +- %.3f cps ...\n',strline,MeanAmplitude,MeanAmplitudeError);

% Histo
figure(150)
str_title = sprintf('Mean Amplitude (ps) - Data - Krypton %s Line - Good Pixels',strline);
title(str_title);
nhist(f(ListGoodPixels,3),'binfactor',3,'color','qualitative'); grid on;
xlabel('Amplitude (cps)'); ylabel('Entries');

% Map
switch FPD
    case 'ON'
        str_title = sprintf('Amplitude (cps) - Data - Krypton %s Line - Good Pixels',strline);
FPDViewer(f(:,3)','outliers',BadPixels);
title(str_title);
end

% Map - Weigthed Mean
switch FPD
    case 'ON'
        str_title = sprintf('Amplitude - %.2f cps - Data - Krypton %s Line - Good Pixels',MeanAmplitude,strline);
FPDViewer(f(:,3)' - MeanAmplitude,'outliers',BadPixels);
title(str_title);
end

%% Background
[MeanBackground MeanBackgroundError]        = swmean(f(ListGoodPixels,4),(1./f(ListGoodPixels,9).^2)); % .*escalefactor(ListGoodPixels)');
%MeanBackgroundError    = 0; % .*escalefactor(ListGoodPixels)');
fprintf(2,'%s Weighted Mean Background = %.2f +- %.2f cps ...\n',strline,MeanBackground,MeanBackgroundError);

% Histo
figure(200)
str_title = sprintf('Mean Background (cps) - Data - Krypton %s Line - Good Pixels',strline);
title(str_title);
nhist(f(ListGoodPixels,4),'binfactor',3,'color','qualitative'); grid on;
xlabel('Background (cps)'); ylabel('Entries');

% Map
switch FPD
    case 'ON'
        str_title = sprintf('Background (cps) - Data - Krypton %s Line - Good Pixels',strline);
FPDViewer(f(:,4)','outliers',BadPixels);
title(str_title);
end

% Map - Weigthed Mean
switch FPD
    case 'ON'
        str_title = sprintf('Background - %.2f cps - Data - Krypton %s Line - Good Pixels',MeanBackground,strline);
FPDViewer(f(:,4)' - MeanBackground,'outliers',BadPixels);
title(str_title);
end

disp([ListGoodPixels' f(ListGoodPixels,1)-MeanPosition f(ListGoodPixels,2)-MeanWidth f(ListGoodPixels,3)-MeanAmplitude])
[maxP maxPI] = max(f(ListGoodPixels,1)-MeanPosition);
[minP minPI] = min(f(ListGoodPixels,1)-MeanPosition);
[maxW maxWI] = max(f(ListGoodPixels,2)-MeanWidth) ;
[minW minWI] = min(f(ListGoodPixels,2)-MeanWidth) ;
[maxA maxAI] = max(f(ListGoodPixels,3));
[minA minAI] = min(f(ListGoodPixels,3));
disp([ListGoodPixels(minPI) minP ListGoodPixels(maxPI) maxP] );
disp([ListGoodPixels(minWI) minW ListGoodPixels(maxWI) maxW]);
disp([ListGoodPixels(minAI) minA ListGoodPixels(maxAI) maxA]);
end