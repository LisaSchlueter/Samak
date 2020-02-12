clear;
RecomputeFlag = 'ON';
%SysEffect     = 'stat';
SysEffect     = 'all';
SysBudget     = '08'; %defines systematic budget in ref_CovarianceMatrix_...
chi2          = 'chi2CMShape';
% Input:
MACE_Ba_T = 7*1e-04;
BsBmRed   = 0.7;
WGTS_B_T  = 3.6*BsBmRed;
TimeSec   = 900*86400;
range = 60;                                     % MTD range (eV below 18575)
NuMassFactorV = [1e-6 0.2 0.4 0.6 0.8];
BkgFractionV  = [1e-6 0.1 0.2 0.3 0.4];           
Anchor6GValue = 60e-3;

% Init Gather results
chi2min = zeros(numel(NuMassFactorV),numel(BkgFractionV));
mNu90   = zeros(numel(NuMassFactorV),numel(BkgFractionV)); %sensitivity of neutrino mass squared 90% C.L.
mNumin  = zeros(numel(NuMassFactorV),numel(BkgFractionV));

%% do stat only --> loop over time
save_file = sprintf('./results/OptimizeMTD%.0f_NuMassFactor_BkgFraction_%s_%s_%0.fmcps.mat',range,SysEffect,SysBudget,Anchor6GValue*1e3);
if exist(save_file,'file')==2 && strcmp(RecomputeFlag,'OFF')
    load(save_file);
else
    for i=1:numel(NuMassFactorV)
        for j=1:numel(BkgFractionV)
            [~,~, TD]  = MTDcreator(...
                'NuMassFactor',NuMassFactorV(i),'BkgFraction',BkgFractionV(j),...
                'BqU',[5 10 20],'Range',range,...
                'MACE_Ba_T',MACE_Ba_T,'BsBmRed',BsBmRed,...
                'MTD_Plot','OFF','NuMassSignal_Plot','OFF',...
                'Anchor6G','ON','Anchor6GValue',Anchor6GValue);
            %        [~,~,~,~,dof,tmp,~] = NuMassScan_SensitivityNominal(...
            [mnuSq_i_Fit, par, err, tmp2, dof,tmp,mNumin] =  NuMassScan_SensitivityNominal(...
                'SysBudget',SysBudget,'TD',TD,'chi2',chi2,'TimeSec',TimeSec,...
                'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,...
                'Scan','OFF','ScanPrcsn',0.02,'Anchor6G','ON','Anchor6GValue',Anchor6GValue);
            mNu90(i,j)=tmp;
            save(save_file,'TimeSec','mNu90','MACE_Ba_T','WGTS_B_T','','TD','SysEffect','range','chi2','dof');
        end
    end
    save(save_file,'TimeSec','mNu90','MACE_Ba_T','WGTS_B_T','','TD','SysEffect','range','chi2','dof');
end

%% Plot Results
Results=load(save_file);
figure(1)
strtitle = sprintf('KATRIN %g d - Ba=%.1f G- Bs=%.1f T- Bm=%.1f T - %s (%s) - IsoStatBoost %0.f eV',...
    round(TimeSec./86500),MACE_Ba_T*1e4,WGTS_B_T,WGTS_B_T/3.6*6,SysEffect, SysBudget, range);
imagesc(sqrt(Results.mNu90));
xticks([1:numel(BkgFractionV)])
xticklabels(BkgFractionV);
xlabel('Background Time Fraction');
yticks([1:numel(NuMassFactorV)])
yticklabels(NuMassFactorV);
ylabel('\nu-mass Signal Boost');
colormap(flipud(jet));
PrettyFigureFormat
title(strtitle,'FontSize',12);
hcb = colorbar;
ylabel(hcb,'90% CL Upper Limit')