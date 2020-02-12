function [qU krdatamap] = readKrData(varargin)
%
%             Read Krypton Data
%            To fit 83mKr Lines
%
% Input   : Line
% Output  :
%     qU
%       * retarding energy qU in eV
%
%     krdatamap
%       * count rate in cps
%       * uncertainty on count rate in cps
%
%          Th. Lasserre - CEA Saclay
%                August 2017
%

% Initialization
addpath(genpath('../../../samak'));

% Parser
p = inputParser;
p.addParameter('fign',1,@(x)isfloat(x) && x>0);
p.addParameter('pub','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('TD','DummyKrdc1',@(x)ismember(x,{'KrK32','KrL3_32','KrL3_32_HS','MoSFitterJuly2017','KrL3_32_Satellites','DummyKrdc1'}));
p.parse(varargin{:});
fign       =    p.Results.fign;
pub        =    p.Results.pub;
display    =    p.Results.display;
TD         =    p.Results.TD;

switch TD
    case 'KrL3_32'
        folder = '33051_L3-32';
        %folder = '33051_L3-32_Simulation'; %When Fitting Simulation
        %folder = '33051_L3-32_Simulation_nofluct'; %no stat. fluctuations
        prefix = '33051';
    case 'KrL3_32_HS'
        folder = '33196_L3-32_high_stat_satellite';
        prefix = '33196';
    case 'KrK32'
        folder = '33054_K-32';
        prefix = '33054';
    case 'MoSFitterJuly2017'
        folder = 'MoSFitterJuly2017';
        prefix = 'output_int_L3-32';
 case  'KrL3_32_Satellites'
        folder = '33196_L3-32_high_stat_satellite';
        prefix = '33196';
    case 'DummyKrdc1'
        folder  = 'DataChallenge18V1';
        prefix  = 'kdc_v1_';
end

switch TD
    case {'KrL3_32','KrK32'}
        krdatamap = zeros(148,2,31);
    case 'KrL3_32_HS'
        krdatamap = zeros(148,2,150);
    case 'KrL3_32_Satellites'
        krdatamap = zeros(148,2,150);
    case 'DummyKrdc1'
        krdatamap = zeros(1,2,43);
end

switch TD
    case {'KrL3_32','KrL3_32_HS','KrK32','KrL3_32_Satellites'}   
        for(i=0:1:147)
            Pixel = num2str(i);
            krfile = sprintf('../../krypton-data/%s/%s_channel%s.txt',folder,prefix,Pixel);
            krplot = sprintf('../../krypton-data/figs/%s_channel%s.eps',folder,prefix,Pixel);
            krtitle = sprintf('%s Line - Pixel %s',TD,Pixel);
            %fprintf('File %s - %s Line - Pixel %s\n',krfile,TD,Pixel);
            data = importdata(krfile);
            
            % read qU
            if (i==0) qU      = flip(data(:,1)); end   
            % read pixel
            krdatamap(i+1,1,:)    = flip(data(:,2));
            krdatamap(i+1,2,:)    = flip(data(:,3));
            
            
        end
    case 'MoSFitterJuly2017'
        krdatamap = zeros(1,2,31);
        krfile = sprintf('../../krypton-data/%s/%s.txt',folder,prefix);
        krplot = sprintf('../../krypton-data/figs/%s.eps',folder,prefix);
        krtitle = sprintf('%s',TD);
        data = importdata(krfile);
        
        % read qU
        qU                    = data(:,1);
        
        % read pixel
        krdatamap(1,1,:)      = (data(:,2));
        krdatamap(1,2,:)      = (data(:,3));
    case 'DummyKrdc1'
        %modelname = 'ssc';
        %modelname = 'fitness';
        %modelname = 'samak_relativistic';   
        %modelname= 'ssc_with_revised_amp_no_fluctuation'; % 1day subruns
        %modelname = 'fitness_1day_asimov'; 
        modelname= 'samak_1day_asimov';
        krfile = sprintf('../../../krypton-data/%s/%s%s.txt',folder,prefix,modelname);
        data = importdata(krfile); 
        qU                    = data(:,1);
        krdatamap(1,1,:)      = (data(:,2));
        krdatamap(1,2,:)      = (data(:,3));
end
