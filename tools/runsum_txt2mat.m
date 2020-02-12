function srm = runsum_txt2mat(run)
%%
%
% P. I. Morales Guzmán, December 2017
%

%#ok<*NASGU>
%run = 1;

if nargin == 0
   run = 0;
end


%% Generate path to find files
addpath(genpath('../../../Samak2.0'));
addpath(genpath('../../Samak2.0'));
addpath(genpath('../samak'));

%% Import Data into workspace

% Sub runs times
td = importdata(['tritium-bootcamp/subruntimes.dat']);
qUfrac = td/sum(td); 

% 1st col: retarding potential; 2nd col: counts; 3rd col: error
TBDISallPixels = importdata(['tritium-bootcamp/TBDIS.dat']);

if TBDISallPixels(1,1) > TBDISallPixels(end,1)
    TBDISallPixels = flipud(TBDISallPixels);
    td = flipud(td);
    qUfrac = flipud(qUfrac);
end
    

[nqU,npixel] = size(TBDISallPixels);

npixel = npixel/3;

% Content order:
%   1 Measuring time
%   2 Temperature
%   3 Column Density
%   4 Magnetic Field source
%   5 Magnetic Field Pinch
%   6 Tritium Purity
%   7 Detector Efficiency
%   8 Final States Model
%   9 Doppler Mode
%   10 Theoretical Corretions
%   (Potential Offset at pixel 3)
%   (Magnetic Field at pixel 3)

filename_runconfig = 'tritium-bootcamp/RunConfig.txt';
formatSpec_runconfig = '%s%[^\n\r]';
fileID_runconfig = fopen(filename_runconfig,'r');
runconfig = textscan(fileID_runconfig, formatSpec_runconfig,...
    'Delimiter', {''}, 'TextType', 'string','ReturnOnError', false);
fclose(fileID_runconfig);

if length(runconfig{1}) < 10
   runconfig = [runconfig{1};" "];
else
   runconfig = runconfig{1};
end

runconfig_num = str2double(runconfig);

TimeSec         = sum(td); TimeYear = TimeSec/(365.25*24*60*60);
WGTS_Temp       = runconfig_num(2);
WGTS_CD_MolPerCm2= runconfig_num(3)*1e-4;
WGTS_B_T        = runconfig_num(4);
MACE_Bmax_T     = runconfig_num(5);
WGTS_Tp         = runconfig_num(6);
FPD_MeanEff     = runconfig_num(7);

% Maybe these below changes
switch runconfig(8)
    case 'saenz'
        TTFSD = 'SAENZ';
        DTFSD = 'DOSS';
        HTFSD = 'SAENZ';
    case 'doss01'
        TTFSD = 'DOSS';
        DTFSD = 'DOSS';
        HTFSD = 'SAENZ';
    case 'doss'
        TTFSD = 'DOSS';
        DTFSD = 'OFF';
        HTFSD = 'OFF';
end

% Maybe these below changes
switch runconfig(9)
    case 'thermal'
        DopplerEffectFlag = 'ON';
    case 'bulk'
        DopplerEffectFlag = 'ON';
    case 'combined'
        DopplerEffectFlag = 'ON';
    case 'disabled'
        DopplerEffectFlag = 'OFF';
end

ScreeningFlag       = 'OFF';
FiniteExtChargeFlag = 'OFF';
WintFiniteSizeFlag  = 'OFF';
EEexchangeFlag      = 'OFF';
RecoilCoulombFlag   = 'OFF';
RadiativeFlag       = 'OFF';
RecoilWmVmAFlag     = 'OFF';

% Check
switch runconfig(10)
    case 'G'
        RadiativeFlag = 'ON';
end

% 1st column: potential offset; 2nd column: magnetic field analysis plane
Offsets = importdata(['tritium-bootcamp/Offsets.dat']);

PotentialOffset_all = Offsets(:,1);
MACE_Ba_TallPixels = Offsets(:,2);


clear runconfig runconfig_num filename_runconfig formatSpec_runconfig
clear fileID_runconfig Offsets ans
currentdirectory = pwd;
idcs = strfind(currentdirectory,'samak');

 save([currentdirectory(1:idcs-1),'samak/tritium-bootcamp/runmat',num2str(run),'.mat'],'-v7.3','-nocompression',...
     '-regexp','^(?!(idcs|currentdirectory)$).')


%% Create structure
srm = load(['tritium-bootcamp/runmat',num2str(run),'.mat']);

end
