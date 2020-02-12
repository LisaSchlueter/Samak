function  [mNu_Fit, parScan, errScan, chi2minScan, dof ,mNu90,mNumin] = NuMassScan_SensitivityNominal_Systematics_Loop(varargin)
p = inputParser;
p.addParameter('SysBudget','03',@(x)ischar(x));
p.addParameter('MACE_Ba_T',7*1e-04,@(x)isfloat(x));
p.addParameter('WGTS_B_T',0.7*3.6,@(x)isfloat(x));
p.addParameter('range',60,@(x)isfloat(x));
p.addParameter('TimeSec',3*365*24*60*60,@(x)isfloat(x));
p.addParameter('RecomputeFlag','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('TD','',@(x)ischar(x)); %if empty, fill accordning to range and Ba ("optimal")
p.addParameter('ScanPrcsn',0.02,@(x)isfloat(x));
p.addParameter('mNuStop',0.6,@(x)isfloat(x));
p.addParameter('Scan','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('BKG_RateSec','',@(x)isfloat(x)); %if empty, filled according to Ba
p.addParameter('Q_i',18575,@(x)isfloat(x));
p.addParameter('PlotFit','OFF',@(x)ischar(x));
p.addParameter('SaveResults','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('FPD_MeanEff',0.9,@(x)isfloat(x) && x>=0);
p.addParameter('Anchor6G','OFF',@(x)ismember(x,{'ON','OFF'})); % Background Option
p.addParameter('Anchor6GValue',335e-3,@(x)isfloat(x)); % 335=FT, [26keV, 32 keV]ROI

p.parse(varargin{:});

SysBudget     = p.Results.SysBudget;
MACE_Ba_T     = p.Results.MACE_Ba_T;
WGTS_B_T      = p.Results.WGTS_B_T;
range         = p.Results.range;
TimeSec       = p.Results.TimeSec;
RecomputeFlag = p.Results.RecomputeFlag;
TD            = p.Results.TD;
ScanPrcsn     = p.Results.ScanPrcsn;
mNuStop       = p.Results.mNuStop;
Scan          = p.Results.Scan;
BKG_RateSec   = p.Results.BKG_RateSec;
Q_i           = p.Results.Q_i;
PlotFit       = p.Results.PlotFit;
SaveResults   = p.Results.SaveResults;
FPD_MeanEff  = p.Results.FPD_MeanEff; 
Anchor6G     = p.Results.Anchor6G;
Anchor6GValue = p.Results.Anchor6GValue;

% Compute Numass sensitivity using Scan Method
% Inside this loop:  loop over Sys Effects and save results for all sys effects to file
% Output: numel(MACE_Ba_T) files
addpath(genpath('../../../Samak2.0'));
mySysEffects  = {'TC','FSD','RF','all','RF_EL','RF_BF','RF_RX','RF_BFRX'};
nFitMax = 20;
if isempty(TD)
TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,MACE_Ba_T*1e4);
end
TDlabel = strrep(TD,sprintf('_Ba%.0fG',MACE_Ba_T*1e4),''); %
TDlabel = strrep(TDlabel,'Sensitivity_','');
save_name = sprintf('./results/%s_%0.0fd_SensitivityNominal_ResultsNuMassScan_MTD-%s_Ba%.0fG_Bs%.2fT_Systematics.mat',SysBudget,round(TimeSec/(86400)),TDlabel,MACE_Ba_T*1e4,WGTS_B_T);

%Init: gather results
chi2minScan  = zeros(numel(mySysEffects)+1,nFitMax); % chi2min distribution
mNu_Fit      = zeros(numel(mySysEffects)+1,nFitMax);
mNu90        = zeros(numel(mySysEffects)+1,1);                  % sensitivity on mnu^2 (90% C.L.)
mNumin       = zeros(numel(mySysEffects)+1,1);                  % mass with minimal chi2
parScan      = zeros(numel(mySysEffects)+1,6,nFitMax);
errScan      = zeros(numel(mySysEffects)+1,6,nFitMax);

%Load File if possible
if exist(save_name,'file')==2 && strcmp(RecomputeFlag,'OFF')
    sprintf('File already exist. Do you want to recompute? \n')
    d = importdata(save_name);
    parScan     = d.parScan;
    errScan     = d.errScan;
    chi2minScan = d.chi2minScan;
    mNu90       = d.mNu90;
    mNumin      = d.mNumin;
    dof         = d.dof;
    try
        mNu_Fit     = d.mNu_Fit;
    catch
        
    end
    return
else
    %% Do Scan
    if isempty(BKG_RateSec)
        BKG_RateSec = GetBackground('MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T);
    end
    
    Scan_Arg = {'TimeSec',TimeSec,...
        'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'Scan',Scan,'BKG_RateSec',BKG_RateSec,...
        'range',range,'ScanPrcsn',ScanPrcsn,'mNuStop',mNuStop,'plotFit',PlotFit,'TD',TD,'Q_i',Q_i,...
        'FPD_MeanEff',FPD_MeanEff};
    %stat
    [mNu_Fit(1,:), parScan(1,:,:), errScan(1,:,:), chi2minScan(1,:,:),~,mNu90(1),mNumin(1)] = ...
        NuMassScan_SensitivityNominal(Scan_Arg{:},'chi2','chi2Stat');
    
    % Systematics
    for i=1:numel(mySysEffects)
        [mNu_Fit(i+1,:), parScan(i+1,:,:), errScan(i+1,:,:), chi2minScan(i+1,:,:), dof ,mNu90(i+1),mNumin(i+1)] = ...
            NuMassScan_SensitivityNominal(Scan_Arg{:},'chi2','chi2CM',...
            'SysEffect',mySysEffects{i},'SysBudget',SysBudget,'PlotCM','OFF');
    end
    % save
    if strcmp(SaveResults,'ON')
        save(save_name,'parScan','errScan','chi2minScan','mNu_Fit','ScanPrcsn','mNu90', 'mNumin','TD','mySysEffects','dof','WGTS_B_T','MACE_Ba_T','BKG_RateSec','TimeSec');
    end
end
end
