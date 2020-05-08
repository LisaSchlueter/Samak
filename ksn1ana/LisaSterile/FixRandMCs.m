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
RandMC = [1:151,500:643];
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
    
    extraStr = sprintf('_RandMC%.0f',i);
    savedir = [getenv('SamakPath'),'ksn1ana/LisaSterile/results/RandomizedMC/'];
    savefile = sprintf('%sKSN1_GridSearch_%s_%s_%s_%.0feVrange_%s_%.0fnGrid%s.mat',...
        savedir,RunList,DataType,strrep(freePar,' ',''),range,chi2Str,nGridSteps,extraStr);
    
    d = importdata(savefile);
    RunAnaArg = {'RunList',RunList,...
        'fixPar',freePar,...
        'DataType',DataType,...
        'FSDFlag','SibilleFull',...
        'ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel',...
        'chi2',chi2Str,...
        'ROIFlag','Default',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF',...
        'ISCSFlag','Edep',...
        'NonPoissonScaleFactor',1.064,...
        'TwinBias_Q',18573.70,...
        'SysBudget',22};
    
    T = MultiRunAnalysis(RunAnaArg{:});
    T.exclDataStart = T.GetexclDataStart(range);
    T.RunData.TBDIS = d.TBDIS_mc;
    T.RunData.TBDISE = sqrt(d.TBDIS_mc);
    
    T.Fit;
    FitResults_Null = T.FitResults;
    save(savefile,'FitResults_Null','-append');
end
