RunList = 'KNM1';
% Init Model Object and covariance matrix object
if ~exist('TwinFSD','var')
   TwinFSD = MultiRunAnalysis(...
       'RunList',RunList,...
       'chi2','chi2CMShape',...
       'DataType','Real',...
       'exclDataStart',14,...
       'NonPoissonScaleFactor',1.064);
end

%% Systematics With FSD Offset 
TwinFSD.fixPar = '5 6 7 8 11'; % fix only qUOffset!
SysEffects            = struct('TASR','ON','FSD','ON','RF_RX','ON','RF_EL','ON','RF_BF','ON','BkgShape','ON','TCoff_RAD','ON','TCoff_OTHER','ON','Stack','ON');
BkgCM                 = 'ON';
TwinFSD.ComputeCM('SysEffects',SysEffects,'BkgCM','ON','FSDNorm_RelErr',0);

%% Fit
TwinFSD.Fit;
TwinFSD.PlotFit;

%% Scan of Fit of FSD Offset as a function of exclDataStart
TwinFSD.qUScanFSD_TT('firstPoint',2,'lastPoint',20);

%% RhoD Scan, with Pgs as a free parameter
%[parRD, errRD, chi2RD, WGTS_CD_MolPerCm2_local, CD_bestfit] = TwinFSD.RhoDScan;