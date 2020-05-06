% contour with randomized twins
%% settings
CL = 0.82;
range = 95;%
nGridSteps = 25;
chi2Str = 'chi2CMShape';
DataType = 'Twin';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'OFF';
RandMC = 152:499;%[1:151,500:643];
nContours = numel(RandMC);
mnu4Sq   = cell(nContours,1);
sin2T4   = cell(nContours,1);
chi2     = cell(nContours,1);
chi2_ref = cell(nContours,1);
savefile = cell(nContours,1);
DeltaChi2 = zeros(nContours,1);
%% load grid (or calculate if doesn't exist)
for i=RandMC
    progressbar(i/nContours)
    [mnu4Sq{i},sin2T4{i},chi2{i},chi2_ref{i},savefile{i}] = KSN1GridSearch('range',range,...
        'nGridSteps',nGridSteps,...
        'chi2',chi2Str,...
        'DataType',DataType,...
        'freePar','E0 Bkg Norm',...
        'RunList',RunList,...
        'SmartGrid',SmartGrid,...
        'RecomputeFlag','OFF',...
        'RandMC',i);
    
    chi2tmp = chi2{i};
    DeltaChi2(i) = chi2tmp(1,1)-min(min(chi2tmp));
end
%%
ClosedLog95 =  DeltaChi2>= GetDeltaChi2(0.95,2);
ClosedFrac95 = sum(ClosedLog95)/nContours;
fprintf('%.0f%% C.L. : fraction of significant best fits = %.1f   (%.0f out of %.0f)\n',...
    95,ClosedFrac95*100,sum(ClosedLog95),nContours);
ClosedLog82 =  DeltaChi2>= GetDeltaChi2(0.84,2);
ClosedFrac82 = sum(ClosedLog82)/nContours;
fprintf('%.0f%% C.L. : fraction of significant best fits = %.1f  (%.0f out of %.0f)\n',...
    84,ClosedFrac82*100,sum(ClosedLog82),nContours);
%%
if strcmp(chi2Str,'chi2Stat')
    chi2Label = 'stat. only';
else
    chi2Label = 'stat. and syst.';
end
%% plot some
PlotFlag = 'OFF';
if strcmp(PlotFlag,'ON')
    for i=RandMC
        progressbar(i/nContours)
        titleStr = sprintf('Randomized MC data %.0f (%s) %.0f eV range',i,chi2Label,range);
        SaveAs = sprintf('RandomizedMC/KSN1_GridSearch_KNM1_RandomizedTwin%.0f_%s_%.0feVrange_%s_%.0fnGrid_Grid.png',...
            i,chi2Str,range,strrep(freePar,' ',''),nGridSteps);
        PlotArg ={'mnu4Sq',mnu4Sq{i},...
            'sin2T4',sin2T4{i},...
            'chi2',chi2{i},'chi2_ref',chi2_ref{i},...
            'CL',CL,...
            'titleStr',titleStr,...
            'SaveAs',SaveAs};
        KSN1GridPlot(PlotArg{:},'nInter',1e3);
        close all
    end
end
%%

 %KSN1ContourPlot(PlotArg{:},'LineStyle','-','PlotSplines','ON');
 
 