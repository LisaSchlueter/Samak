clear; close all;

% Add path to samak folder
addpath(genpath('../../../Samak2.0'));

% load the christmas background as template
%BCK = load('C:\Users\navar\Documents\Samak2.0\inputs\BCK\BackgroundRate_35023-35110.mat');
MAG = load('C:\Users\navar\Documents\Samak2.0\inputs\WGTSMACE\MACE_Ba_TCorr6G_FT.mat');
%ELE = load('C:\Users\navar\Documents\Samak2.0\tritium-data\mat\40667mpix.mat');

%FPD = FPDViewer(BCK.bck(:,2));
%FPD = FPDViewer(nan(1,13));
FPD = FPDViewer(MAG.Ba);
%FPD = FPDViewer(ELE.qU(23,:));

screensize = (get(0, 'Screensize'));
%set(FPD, 'Position', [1,1,1000/1.5,800/1.5]);
%colorbar('off')
title('Magnetic Field variations in the analyzing plane at 6 G');
%title('Electrostatic Field variations in the analyzing plane at -18573 V');


%export_fig('plots/DetectorPixelOFF.pdf','-pdf');
publish_figurePDF(gcf,'plots/DetectorMAG.pdf');


