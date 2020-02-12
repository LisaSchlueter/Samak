% 
% Read Analysis Plane Data
% 
% Simulation done 
%   - for 2.7G --> in the analysis plane
%   - for 6G --> in the analysis plane
% 
% Pixel
% E-field values, Volt
% B-field value, Tesla
% 

clear all; close all;

addpath(genpath('../../../Samak2.0'));
addpath(genpath('../../tools/masumhabib-PlotPub-23bcfed/lib'));

% Switch Gauss Flag
%Ba = 2.7;  % G
Ba = 6.0;   % https://neutrino.ikp.kit.edu:8080/FirstTritium/180507_170821/GlobalSimulation-FirstTritium-PeriodSummary_May2018-WQ_18600V_6.0G.ktf

% data files - from Martin Slezak, July 2017
if Ba == 2.7 
dataEtmp   = importdata('simulation_2.7G_elPotential_pixels.txt');
dataBtmp   = importdata('simulation_2.7G_magField_pixels.txt');
% Method 1
dataE      = -(3e4-abs(dataEtmp));     % V
dataB      = (dataBtmp-(Ba*1e-4))*1e4; % Gauss
% Method 2
% dataE = -(3e4-abs(dataEtmp))/3e4;
% dataB = -(2.7e-4-abs(dataBtmp))/2.7e-4;
elseif Ba == 6
Etmp        = load('MACE_Ea_TCorrWQ6G.mat');
dataEtmp    = Etmp.Ea;
Btmp   = load('MACE_Ba_TCorrWQ6G.mat');
dataBtmp    = Btmp.Ba;
dataE = -(18600-abs(dataEtmp));  %dataE=dataE-mean(dataE);   % mV
dataB = (dataBtmp)*1e4-6.0;  %dataB=dataB-mean(dataB);           % mGauss
else
    return;
end

% Inputs
fign = 5; pub  = 0;

% FDP Rings <--> Pixel Association
pixel    = linspace(1,148,148)';
ring{1}  = [0:3]';
ring{2}  = [4:15]';
ring{3}  = [16:27]';
ring{4}  = [28:39]';
ring{5}  = [40:51]';
ring{6}  = [52:63]';
ring{7}  = [64:75]';
ring{8}  = [76:87]';
ring{9}  = [88:99]';
ring{10} = [100:111]';
ring{11} = [112:123]';
ring{12} = [124:135]';
ring{13} = [136:147]';

AplanePEB      = [pixel dataE dataB];
AplaneE_mean   = mean(dataE);
AplaneB_mean   = mean(dataB);

%% E-field values, Volt
figure(fign)
plt1 = Plot(pixel,dataE);
strptl1 = sprintf('Electric field corrections (V) - %g Gauss',Ba);
plt1.Title  = strptl1; % plot title
plt1.YLabel = '\DeltaE (Volt)'; % xlabel
plt1.XLabel = 'Pixel'; %ylabel
plt1.ShowBox = 'on';
plt1.FontSize = 16;
plt1.Legend   = {'With Alignement corrections'};
plt1.export('MACEEinhomogeneity-pixels.png');

%%
figure(fign+1)
FPDViewer(dataE,'ReDrawSkeleton','ON');
title(strptl1);

% B-field values, Tesla
figure(fign+2)
plt2 = Plot(pixel,dataB);
strptl2 = sprintf('Magnetic field corrections (G) - %g Gauss',Ba);
plt2.Title  = strptl2; % plot title
plt2.YLabel = '\DeltaB (Gauss)'; % xlabel
plt2.XLabel = 'Pixel'; %ylabel
plt2.ShowBox = 'on';
plt2.FontSize = 16;
plt2.Legend   = {'With Alignement corrections'};
plt2.export('MACEBinhomogeneity-pixels.png');

%%
figure(fign+3)
FPDViewer(dataB,'ReDrawSkeleton','ON');
title(strptl2);

return;

%% ringwise data handling
for i=1:1:13
    AplaneEBring{i} = [ring{i}+1 dataE(ring{i}+1) dataB(ring{i}+1)];
    AplaneEBring_mean{i} = [i mean(dataE(ring{i}+1)) mean(dataB(ring{i}+1))];
end
AplaneEBring_mean = reshape(cell2mat(AplaneEBring_mean),[3,13])';


%% E-field - ring-wise
figure(fign+2)
%h = plot(AplaneEBring_mean(:,1,:),AplaneEBring_mean(:,2,:),'LineWidth',2,'Color','Black')
h =bar(-(AplaneEBring_mean(:,2,:)),'red');
axis([0 14 (1-1e-6)*min(-AplaneEBring_mean(:,2,:)) (1+1e-6)*max(-AplaneEBring_mean(:,2,:))]);
grid on
xlabel('ring','FontSize',14);
ylabel('E-field in Volt','FontSize',14);
title('E-field in the MAC-E spectrometer analysis plane','FontSize',14)
set(gca,'FontSize',12);
strl = sprintf('<E> = %.2s Volt',AplaneE_mean);
a = legend(h,strl); % legend(a,'boxoff');
PrettyFigureFormat;
if pub>0
    publish_figure(fign,'efieldvalues-ring.eps')
end

%% B-field - ring-wise
figure(fign+3)
%h = plot(AplaneEBring_mean(:,1,:),AplaneEBring_mean(:,3,:),'LineWidth',2,'Color','Black');
h =bar(AplaneEBring_mean(:,3,:));
axis([0 14 0.995*min(AplaneEBring_mean(:,3,:)) 1.005*max(AplaneEBring_mean(:,3,:))]);
grid on
xlabel('ring','FontSize',14);
ylabel('B-field in Tesla','FontSize',14);
title('B-field in the MAC-E spectrometer analysis plane','FontSize',14)
set(gca,'FontSize',12);
strl = sprintf('<B> = %.2s Tesla',AplaneB_mean);
a = legend(h,strl); % legend(a,'boxoff');
PrettyFigureFormat;
if pub>0
    publish_figure(fign,'bfieldvalues-ring.eps')
end