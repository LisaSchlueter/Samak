% Compute Numass sensitivity using Scan Method
% Loop over Systematic Effects, Ba, SysBudgets,...
% Inside this loop:  loop over Sys Effects and save results for all sys effects to file
% Output: results in /results folder files
addpath(genpath('../../../Samak2.0'));
clear;
SysBudget_all = {'03'};
MACE_Ba_T_all = (3:12).*1e-04;    % B-Field analyzing plane. Background is scaled automatically in SensitivityStudy_NominalKATRIN
WGTS_B_T = 0.7*3.6;                        % MTD energy range: 30,45,60eV below E0 (18575eV)
TimeSec = 3*365*24*60*60;
RecomputeFlag = 'ON';
TDFlag = 'Opt'; % 'Opt'
ScanPrcsn = 0.02; % Scan Precision (Acceptable Delta chi2)
mNuStop = 0.6; % Stop mNu eV for neutrino mass scan
nFitMax = 20;
Q_i = 18575;
Scan = 'ON';
Anchor6G = 'OFF';
FPD_MeanEff = 0.9;
PlotFit = 'OFF';
switch TDFlag
    case 'Opt'
        range_all = [30,45,60];
    case 'IsoStat'
        range_all =  [30, 60, 90,120];
end
for sys=1:numel(SysBudget_all)
    SysBudget = SysBudget_all{sys};
    for r=1:numel(range_all)
        range = range_all(r);     
        for b=1:numel(MACE_Ba_T_all)
            MACE_Ba_T = MACE_Ba_T_all(b);
            % Do Systematics Loop
            switch TDFlag
                case 'IsoStat'
                    TD = sprintf('IsoStat%.0f',range);
                case 'Opt'
                    TD  = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,MACE_Ba_T*1e4);
            end
            BKG_RateSec = GetBackground('MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'Anchor6G',Anchor6G);
            LoopArg = {'TimeSec',TimeSec,...
                'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'Q_i',Q_i,...
                'TD',TD,'range',range,'RecomputeFlag',RecomputeFlag,...
                'BKG_RateSec',BKG_RateSec,'SysBudget',SysBudget,...
                'Scan',Scan,'PlotFit',PlotFit,'saveResults','ON',...
                'Anchor6G',Anchor6G,'FPD_MeanEff',FPD_MeanEff};
            [mNu_Fit, parScan, errScan, chi2minScan, dof ,mNu90,mNumin] = ...
                NuMassScan_SensitivityNominal_Systematics_Loop(LoopArg{:});
        end
    end
end

