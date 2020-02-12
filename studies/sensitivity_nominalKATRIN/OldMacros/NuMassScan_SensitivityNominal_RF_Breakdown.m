% Compute Numass sensitivity using Scan Method
% Loop over Ba
% Inside this loop:  loop over Sys Effects and save results for all sys effects to file
% Output: numel(MACE_Ba_T) files
clear;
addpath(genpath('../../../Samak2.0'));
SysBudget = '02';
WGTS_B_T = 3.6*0.7;
MACE_Ba_T_all = (3:12).*1e-04;
range_all = [30,45,60]; % MTD energy range: 30,45,60eV below E0 (18575eV)
TDFlag = 'Opt'; %'IsoStat'
RecomputeFlag = 'ON';
ScanPrcsn = 0.02; % Scan Precision (Acceptable Delta chi2)
mNuStop = 0.6; % Stop mNu eV for neutrino mass scan

for r=1:numel(range_all)
    range = range_all(r);
    %for b=1:numel(MACE_Ba_T_all)
        MACE_Ba_T = 7*1e-04;%MACE_Ba_T_all(b);
        
        switch TDFlag
            case 'IsoStat'
                TD = sprintf('IsoStat%.0f',range);
                save_name = sprintf('./results/%s_SensitivityNominal_ResultsNuMassScan_IsoStat%.0feV_Ba%.0fG_Bs%.2fT_Systematics_RFcomponents.mat',SysBudget,range,MACE_Ba_T*1e4,WGTS_B_T);
            case 'Opt'
                TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,MACE_Ba_T*1e4);
                save_name = sprintf('./results/%s_SensitivityNominal_ResultsNuMassScan_MTD%.0feV_Ba%.0fG_Bs%.2fT_Systematics_RFcomponents.mat',SysBudget,range,MACE_Ba_T*1e4,WGTS_B_T);
        end
        
        
        %Load File if possible
        if exist(save_name,'file')==2 && strcmp(RecomputeFlag,'OFF')
            sprintf('File already exist. Do you want to recompute? \n')
            load(save_name);
        else
            %% Do Scan
            %Init: gather results
            mySysEffects  = {'RF_EL','RF_BF','RF_RX','RF_BFRX','RF'};
            nFitMax = 20;
            TimeSec      = 3*365*24*60*60;
            chi2minScan  = zeros(numel(mySysEffects)+1,nFitMax); % chi2min distribution
            mNu90        = zeros(numel(mySysEffects)+1,1);                  % sensitivity on mnu^2 (90% C.L.)
            mNumin       = zeros(numel(mySysEffects)+1,1);                  % mass with minimal chi2
            parScan      = zeros(numel(mySysEffects)+1,6,nFitMax);
            errScan      = zeros(numel(mySysEffects)+1,6,nFitMax);
            %stat
            [~, parScan(1,:,:), errScan(1,:,:), chi2minScan(1,:,:),~,mNu90(1),mNumin(1)] = ...
                NuMassScan_SensitivityNominal('chi2','chi2Stat','TimeSec',TimeSec,...
                'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,...
                'range',range,'ScanPrcsn',ScanPrcsn,'mNuStop',mNuStop,...
                'plotFit','OFF');
            
            % Systematics
            for i=1:numel(mySysEffects)
                [~, parScan(i+1,:,:), errScan(i+1,:,:), chi2minScan(i+1,:,:), dof ,mNu90(i+1),mNumin(i+1)] = ...
                    NuMassScan_SensitivityNominal('chi2','chi2CM','SysEffect',mySysEffects{i},...
                    'TimeSec',TimeSec,'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,...
                    'range',range,'ScanPrcsn',ScanPrcsn,'mNuStop',mNuStop,...
                       'plotFit','OFF','plotCM','OFF', 'SysBudget',SysBudget);
            end
            BKG_RateSec = GetBackground('MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T);
            save(save_name,'parScan','errScan','chi2minScan','mNu90','ScanPrcsn','mNumin','TD','mySysEffects','dof','WGTS_B_T','MACE_Ba_T','BKG_RateSec','TimeSec');
        end
   % end
end
