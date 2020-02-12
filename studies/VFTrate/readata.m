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

runcount = 0;
run40257 = importdata('data/spectrum_40257_outerRingsExcl.txt'); runcount=runcount+1;
run40258 = importdata('data/spectrum_40258_outerRingsExcl.txt'); runcount=runcount+1;
run40259 = importdata('data/spectrum_40259_outerRingsExcl.txt'); runcount=runcount+1;
run40260 = importdata('data/spectrum_40260_outerRingsExcl.txt'); runcount=runcount+1;
run40263 = importdata('data/spectrum_40263_outerRingsExcl.txt'); runcount=runcount+1;
run40264 = importdata('data/spectrum_40264_outerRingsExcl.txt'); runcount=runcount+1;
run40265 = importdata('data/spectrum_40265_outerRingsExcl.txt'); runcount=runcount+1;
run40266 = importdata('data/spectrum_40266_outerRingsExcl.txt'); runcount=runcount+1;

%% Read and Order qU 
qUmatrix      = zeros(runcount,26);
qUmatrix(1,:) = flip(run40257(1:end-2,1));
qUmatrix(2,:) = (run40258(3:end,1));
qUmatrix(3,:) = flip(run40259(:,1));
qUmatrix(4,:) = (run40260(:,1));
qUmatrix(5,:) = flip(run40263(:,1));
qUmatrix(6,:) = (run40264(:,1));
qUmatrix(7,:) = flip(run40265(:,1));
qUmatrix(8,:) = (run40266(:,1));

%% Get Mean & std
qUmean   = mean(qUmatrix);
qUmedian = median(qUmatrix);
qUstd    = std(qUmatrix);
fprintf(2, '8 Runs - Mean qU\n');
disp([num2str(qUmean'), num2str(qUstd')]);
for i=1:1:26
    fprintf(2,'Averaged qU-E_0: %.2f V \t Median qU-E_0: %.2f V \t  std qU: %.2f V \n',...
        qUmean(i)-18575,qUmedian(i)-18575,qUstd(i));
end
    
% Plot All subruns together
qUdisp = squeeze(squeeze(qUmatrix - repmat(qUmean,runcount,1)));
qUdisp = reshape(qUdisp,runcount*26,1);
nhist(qUdisp,100);
xlabel('qU setpoint - qU (V)');
ylabel('subruns');
PrettyFigureFormat

%% Plot Spread qU Voltage over Runs
figure(9999)
hold on
for i=1:1:runcount
    h(i) = errorbar(qUmean-18575,qUmean-qUmatrix(i,:),qUmean.*0);
    %errorbar(qUmean,qUmatrix(i,:),qUmean.*0);
end
hold off
clear plt2; plt2 = Plot();
plt2str     = sprintf('VFT 19-05-2018 - Stacking of %g runs',runcount);
plt2.LineStyle  = 'none';
plt2.Title  = plt2str; % plot title
plt2.XLabel = 'qU index'; % xlabel
plt2.YLabel = 'qU_{mean} - qU_{subruns}'; %ylabel
plt2.YScale = 'lin'; % 'linear' or 'log'
plt2.XScale = 'lin'; % 'linear' or 'log'
plt2.FontSize = 16;
plt2.XLim = [-500, 100];
plt2.Legend = {'4257','4258','4259','4260','4263','4264','4265','4266'};
%plt2.Legend = {'4257','4259','4263','4265'};
%plt2.Legend = {'4258','4260','4264'};
plt2.LegendLoc = 'NorthWest';
plt2title   = sprintf('VFTstack-1.png');
plt2.export(plt2title);

%% Read qUfrac
tmp=load('VFT2.mat');
qUfrac = tmp.qUfrac;
RunTime = 1860;
% Build Time Distributions
TD = 'Run40257';
qU = flip(run40257(1:end-2,1));
save('Run40257.mat','TD','qU','qUfrac','RunTime');
% Build Time Distributions
TD = 'Run40258';
qUmatrix(2,:) = (run40258(3:end,1));
save('Run40258.mat','TD','qU','qUfrac','RunTime');
% Build Time Distributions
TD = 'Run40259';
qUmatrix(3,:) = flip(run40259(:,1));
save('Run40259.mat','TD','qU','qUfrac','RunTime');
% Build Time Distributions
TD = 'Run40260';
qUmatrix(4,:) = (run40260(:,1));
save('Run40260.mat','TD','qU','qUfrac','RunTime');
% Build Time Distributions
TD = 'Run40263';
qUmatrix(5,:) = flip(run40263(:,1));
save('Run40263.mat','TD','qU','qUfrac','RunTime');
% Build Time Distributions
TD = 'Run40264';
qUmatrix(6,:) = (run40264(:,1));
save('Run40264.mat','TD','qU','qUfrac','RunTime');
% Build Time Distributions
TD = 'Run40265';
qUmatrix(7,:) = flip(run40265(:,1));
save('Run40265.mat','TD','qU','qUfrac','RunTime');
% Build Time Distributions
TD = 'Run40266';
qUmatrix(8,:) = (run40266(:,1));
save('Run40266.mat','TD','qU','qUfrac','RunTime');


