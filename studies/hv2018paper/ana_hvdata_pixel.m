function [ListGoodPixels] = ana_hvdata_pixel(varargin)
%
% Read / Analyse KrLine Fit
% All Pixels
%
% T. Lasserre, August 2017
%

% Initialization
clear par ; clear parnomix;

% Parser
p = inputParser;
p.addParameter('fign',1,@(x)isfloat(x) && x>0);
p.addParameter('pub',@(x)ismember(x,{'ON','OFF','YES'}));
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('FPD','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('file','hvdata_l332_40pixels_v3.mat',...
@(x)ismember(x,{'hvdata_l332_40pixels_v1.mat','hvdata_l332_40pixels_v2.mat','hvdata_l332_40pixels_v3.mat','hvdata_l332_40pixels_v4.mat'}));
p.parse(varargin{:});

fign       =    p.Results.fign;
pub        =    p.Results.pub;
display    =    p.Results.display;
FPD        =    p.Results.FPD;
file       =    p.Results.file;

% Variables
Pixels     = ones(1,148);

% Open
strfile = sprintf('./results/%s',file);
f = importdata(strfile);f(41:148,1:13)=0;
% Bad Pixel List
BadPixels  =  [40:1:148];

switch display
    case 'ON'
        fprintf(2,'BAD Pixels - Begin ...\n');
        disp(num2str(BadPixels));
        fprintf(2,'BAD Pixels - End ...\n');
end

% Build List of Good Pixels - 148 Entries
strline = 'L3-32';
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

MeanChi2 = mean((f(ListGoodPixels,12))');

% Histo
figure(1)
title(str_title);
nhist(f(ListGoodPixels,12),'binfactor',3,'color','qualitative','box','box'); grid on;
xlabel('\chi^2'); ylabel('Entries');

% Map
switch FPD
    case 'ON'
fprintf('##############################\n');
f(:,12)'        
fprintf('##############################\n');
FPDViewer(f(:,12)','outliers',BadPixels);
title(str_title);
end

%% Line Position
%f(:,1) = f(:,1) * (1972.4531/1972.5);
[MeanPosition  MeanPositionError ]   = swmean(f(ListGoodPixels,2),(1./f(ListGoodPixels,7).^2)); % .*escalefactor(ListGoodPixels)');
%MeanPositionError = 0; % .*escalefactor(ListGoodPixels)');
fprintf(2,'%s Weighted Mean Line Position = %.3f +- %.3f eV ...\n',strline,MeanPosition,MeanPositionError);

% Histo
figure(50)
str_title = sprintf('Line Position - Data - Krypton %s Line - Good Pixels',strline);
title(str_title);
nhist(f(ListGoodPixels,2),'binfactor',3,'color','qualitative','box','box'); grid on;
xlabel('Line Position (eV)'); ylabel('Entries');

% Map
switch FPD
    case 'ON'
        str_title = sprintf('Line Position (eV) - Data - Krypton %s Line - Good Pixels',strline);
FPDViewer(f(:,2)','outliers',BadPixels);
title(str_title);
end

% Map - Weigthed Mean
switch FPD
    case 'ON'
        str_title = sprintf('Line Position - %.2f eV - Data - Krypton %s Line - Good Pixels',MeanPosition,strline);
FPDViewer(f(:,2)' - MeanPosition,'outliers',BadPixels);
title(str_title);
end

%% Line Width
%MeanWidth         = mean(f(ListGoodPixels,2)); % .*escalefactor(ListGoodPixels)');
[MeanWidth MeanWidthError ]       = swmean(f(ListGoodPixels,3),(1./f(ListGoodPixels,8).^2)); % .*escalefactor(ListGoodPixels)');
%MeanWidthError    = 0; % .*escalefactor(ListGoodPixels)');
fprintf(2,'%s Weighted Mean Line Width = %.3f +- %.3f eV ...\n',strline,MeanWidth,MeanWidthError);

% Histo
figure(100)
str_title = sprintf('Line Width (eV) - Data - Krypton %s Line - Good Pixels',strline);
title(str_title);
nhist(f(ListGoodPixels,3),'binfactor',3,'color','qualitative','box'); grid on;
xlabel('Line Width (eV)'); ylabel('Entries');

% Map
switch FPD
    case 'ON'
        str_title = sprintf('Line Width (eV) - Data - Krypton %s Line - Good Pixels',strline);
FPDViewer(f(:,3)','outliers',BadPixels);
title(str_title);
end


% Map - Weigthed Mean
switch FPD
    case 'ON'
        str_title = sprintf('Line Width - %.2f eV - Data - Krypton %s Line - Good Pixels',MeanWidth,strline);
FPDViewer(f(:,3)' - MeanWidth,'outliers',BadPixels);
title(str_title);
end

%% Line Amplitude
[MeanAmplitude MeanAmplitudeError]         = swmean(f(ListGoodPixels,4),(1./f(ListGoodPixels,9).^2)); % .*escalefactor(ListGoodPixels)');
%MeanAmplitudeError     = 0; % .*escalefactor(ListGoodPixels)');
fprintf(2,'%s Weighted Mean Amplitude = %.3f +- %.3f cps ...\n',strline,MeanAmplitude,MeanAmplitudeError);
disp(num2str(f(ListGoodPixels,4)'));

% Histo
figure(150)
str_title = sprintf('Mean Amplitude (ps) - Data - Krypton %s Line - Good Pixels',strline);
title(str_title);
nhist(f(ListGoodPixels,4),'binfactor',3,'color','qualitative','box'); grid on;
xlabel('Amplitude (cps)'); ylabel('Entries');

% Map
switch FPD
    case 'ON'
        str_title = sprintf('Amplitude (cps) - Data - Krypton %s Line - Good Pixels',strline);
FPDViewer(f(:,4)','outliers',BadPixels);
title(str_title);
end

% Map - Weigthed Mean
switch FPD
    case 'ON'
        str_title = sprintf('Amplitude - %.2f cps - Data - Krypton %s Line - Good Pixels',MeanAmplitude,strline);
FPDViewer(f(:,4)' - MeanAmplitude,'outliers',BadPixels);
title(str_title);
end

%% Offset
[MeanOffset MeanOffsetError]        = swmean(f(ListGoodPixels,5),(1./f(ListGoodPixels,10).^2)); % .*escalefactor(ListGoodPixels)');
%MeanOffsetError    = 0; % .*escalefactor(ListGoodPixels)');
fprintf(2,'%s Weighted Mean Offset = %.2f +- %.2f cps ...\n',strline,MeanOffset,MeanOffsetError);
disp(num2str(f(ListGoodPixels,5)'));

% Histo
figure(200)
str_title = sprintf('Mean Offset (cps) - Data - Krypton %s Line - Good Pixels',strline);
title(str_title);
nhist(f(ListGoodPixels,5),'binfactor',3,'color','qualitative','box'); grid on;
xlabel('Offset (cps)'); ylabel('Entries');

% Map
switch FPD
    case 'ON'
        str_title = sprintf('Offset (cps) - Data - Krypton %s Line - Good Pixels',strline);
FPDViewer(f(:,5)','outliers',BadPixels);
title(str_title);
end

% Map - Weigthed Mean
switch FPD
    case 'ON'
        str_title = sprintf('Offset - %.2f cps - Data - Krypton %s Line - Good Pixels',MeanOffset,strline);
FPDViewer(f(:,5)' - MeanOffset,'outliers',BadPixels);
title(str_title);
end

% disp([ListGoodPixels' f(ListGoodPixels,2)-MeanPosition f(ListGoodPixels,3)-MeanWidth f(ListGoodPixels,4)-MeanAmplitude])
% [maxP maxPI] = max(f(ListGoodPixels,2)-MeanPosition);
% [minP minPI] = min(f(ListGoodPixels,2)-MeanPosition);
% [maxW maxWI] = max(f(ListGoodPixels,3)-MeanWidth) ;
% [minW minWI] = min(f(ListGoodPixels,3)-MeanWidth) ;
% [maxA maxAI] = max(f(ListGoodPixels,4));
% [minA minAI] = min(f(ListGoodPixels,4));
% disp([ListGoodPixels(minPI) minP ListGoodPixels(maxPI) maxP] );
% disp([ListGoodPixels(minWI) minW ListGoodPixels(maxWI) maxW]);
% disp([ListGoodPixels(minAI) minA ListGoodPixels(maxAI) maxA]);
end