%
% Read VFT Data
% Taken on May 19 2018
%
% Run 40257 - down
% Run 40258 - up
%
% Run 40259 - down
% Run 40260 - up
%
% Run 40263 - down
% Run 40264 - up
%
% Run 40265 - down
% Run 40266 - up  - 16th qU value at 15 kV - Run excluded
%

clear all;close all

%data
runcount = 0;
run40257 = importdata('data/spectrum_40257_outerRingsExcl.txt'); runcount=runcount+1;
model40257 = ref_run40257(); model40257.ComputeTBDDS; model40257.ComputeTBDIS; 
run40258 = importdata('data/spectrum_40258_outerRingsExcl.txt'); runcount=runcount+1;
model40258 = ref_run40258(); model40258.ComputeTBDDS; model40258.ComputeTBDIS;
run40259 = importdata('data/spectrum_40259_outerRingsExcl.txt'); runcount=runcount+1;
model40259 = ref_run40259(); model40259.ComputeTBDDS; model40259.ComputeTBDIS;
run40260 = importdata('data/spectrum_40260_outerRingsExcl.txt'); runcount=runcount+1;
model40260 = ref_run40260(); model40260.ComputeTBDDS; model40260.ComputeTBDIS;
run40263 = importdata('data/spectrum_40263_outerRingsExcl.txt'); runcount=runcount+1;
model40263 = ref_run40263(); model40263.ComputeTBDDS; model40263.ComputeTBDIS;
run40264 = importdata('data/spectrum_40264_outerRingsExcl.txt'); runcount=runcount+1;
model40264 = ref_run40264(); model40264.ComputeTBDDS; model40264.ComputeTBDIS; 
run40265 = importdata('data/spectrum_40265_outerRingsExcl.txt'); runcount=runcount+1;
model40265 = ref_run40265(); model40265.ComputeTBDDS; model40265.ComputeTBDIS;
%run40266 = importdata('data/spectrum_40266.txt'); runcount=runcount+1;

%% Read and Order qU 
qUmatrix      = zeros(runcount,26);
qUmatrix(1,:) = flip(run40257(1:end-2,1));
qUmatrix(2,:) = (run40258(3:end,1));
qUmatrix(3,:) = flip(run40259(:,1));
qUmatrix(4,:) = (run40260(:,1));
qUmatrix(5,:) = flip(run40263(:,1));
qUmatrix(6,:) = (run40264(:,1));
qUmatrix(7,:) = flip(run40265(:,1));
%qUmatrix(8,:) = (run40266(:,1));

TBDISmatrix      = zeros(runcount,26);
%% Data
% TBDISmatrix(1,:) = flip(run40257(1:end-2,2));
% TBDISmatrix(2,:) = (run40258(3:end,2));
% TBDISmatrix(3,:) = flip(run40259(:,2));
% TBDISmatrix(4,:) = (run40260(:,2));
% TBDISmatrix(5,:) = flip(run40263(:,2));
% TBDISmatrix(6,:) = (run40264(:,2));
% TBDISmatrix(7,:) = flip(run40265(:,2));
%TBDISmatrix(8,:) = (run40266(:,2));
%% Simulation
TBDISmatrix(1,:) = (model40257.TBDIS./(model40257.qUfrac*model40257.TimeSec)-flip(run40257(1:end-2,2)))./flip(run40257(1:end-2,2)); 
TBDISmatrix(2,:) = (model40258.TBDIS./(model40257.qUfrac*model40257.TimeSec)-run40258(3:end,2))./run40258(3:end,2);
TBDISmatrix(3,:) = (model40259.TBDIS./(model40257.qUfrac*model40257.TimeSec)-flip(run40259(:,2)))./flip(run40259(:,2));
TBDISmatrix(4,:) = (model40260.TBDIS./(model40257.qUfrac*model40257.TimeSec)-run40260(:,2))./run40260(:,2);
TBDISmatrix(5,:) = (model40263.TBDIS./(model40257.qUfrac*model40257.TimeSec)-flip(run40263(:,2)))./flip(run40263(:,2));
TBDISmatrix(6,:) = (model40264.TBDIS./(model40257.qUfrac*model40257.TimeSec)-run40264(:,2))./run40260(:,2);
TBDISmatrix(7,:) = (model40265.TBDIS./(model40257.qUfrac*model40257.TimeSec)-flip(run40265(:,2)))./flip(run40265(:,2));
% Mean
MeanTBDIS = mean(TBDISmatrix,1);

%% Plot Spread qU Voltage over Runs
figure(9999)
hold on
for i = 1:1:7
    %disp(TBDISmatrix(i,:)-MeanTBDIS(1,:))
    %h(i) = errorbar(qUmatrix(i,:),(TBDISmatrix(i,:)-MeanTBDIS(1,:))./MeanTBDIS(1,:),sqrt(TBDISmatrix(i,:))./MeanTBDIS(1,:));
    h(i) = errorbar(qUmatrix(i,:)-18575,(TBDISmatrix(i,:)-MeanTBDIS(1,:)),sqrt(TBDISmatrix(i,:)));
    %h(i) = errorbar(qUmatrix(i,:)-18575,abs(TBDISmatrix(i,:)),sqrt(TBDISmatrix(i,:))*0);

end

hold off
clear plt2; plt2 = Plot();
plt2str     = sprintf('VFT 19-05-2018');
plt2.LineStyle  = 'none';
plt2.Title  = plt2str; % plot title
plt2.XLabel = 'qU index'; % xlabel
%plt2.YLabel = 'counts - mean(counts)'; %ylabel
plt2.YLabel = 'counts '; %ylabel
plt2.YScale = 'lin'; % 'linear' or 'log'
plt2.XScale = 'lin'; % 'linear' or 'log'
plt2.FontSize = 16;
%plt2.XLim = [-40 40];
%plt2.Legend = {'4257','4259','4263','4265'};
plt2.Legend = {'4257','4258','4259','4260','4263','4264','4265'};
%plt2.Legend = {'4257','4259','4263','4265'};
%plt2.Legend = {'4258','4260','4264'};
plt2.LegendLoc = 'SouthEast';
plt2title   = sprintf('VFTstack-1.png');
plt2.export(plt2title);