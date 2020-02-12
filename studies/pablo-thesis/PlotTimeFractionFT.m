clear; close all;

% Add path to samak folder
addpath(genpath('../../../Samak2.0'));

RunList = load('RunListFTFullCD.mat');
RunList = RunList.RunListFTFullCD;
RunList((RunList == 40610) | (RunList == 40612) | (RunList == 40539)) = [];

mpix = ''; % mpix
ringCutFlag = 'ex2'; %ex2

for r = 1:0;length(RunList)
    run_name = [num2str(RunList(r)),'ex2.mat'];
    R = load(run_name);
    if any(diff(R.qU) < 0)
       R.qU  = sort(R.qU);
       R.qUfrac = sort(R.qUfrac); 
    end
    
    hold on
    plot(R.qU,R.qUfrac)
    %pause
    hold off
    
end


RR = ref_RunAnalysis(40667,mpix,ringCutFlag);


figure(1)
RR.PlotTD;

figure(2)
RR.TD = 'FT-TL4';
RR.PlotTD;
