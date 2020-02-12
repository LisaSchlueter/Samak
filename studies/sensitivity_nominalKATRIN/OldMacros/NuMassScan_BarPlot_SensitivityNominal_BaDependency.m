% Plot Sensitivity for 3 ranges
% as a function of Ba (and Background)
% loop over systematics
% doesnt give the same as combi CM, becasue some effects partly cancel?
clear;
SysBudget = '03';  % Systematics Input Information in Result folder
Anchor6G  = 'OFF';
s = GetSysBudet(SysBudget); %information structure
switch SysBudget
    case '00'
        nSys = 6;
    case {'01','03'}
        nSys = 5; %no TASR
end
MACE_Ba_T_all = (3:12).*1e-04;         % B-Field analyzing plane. Background is scaled automatically in SensitivityStudy_NominalKATRIN
range_all = [30,45,60];    % MTD energy range: 30,45,60eV below E0 (18575eV)
WGTS_B_T = 0.7*3.6;
%Init: gather results
mNu90all     = zeros(numel(range_all),numel(MACE_Ba_T_all),nSys);                  % sensitivity on mnu^2 (90% C.L.)
TimeSec = 3*365*24*60*60;  

for i=1:numel(range_all)
    range = range_all(i);
    for b=1:numel(MACE_Ba_T_all)
        MACE_Ba_T = MACE_Ba_T_all(b);
        TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,MACE_Ba_T*1e4);
        save_name = sprintf('./results/%s_SensitivityNominal_ResultsNuMassScan_MTD%s_Bs%.2fT_Systematics.mat',SysBudget,strrep(TD,'Sensitivity_',''),WGTS_B_T);
        if exist(save_name,'file')==2
            d = importdata(save_name);
        else
            fprintf('File doesnt exist! \n');
            return
        end  
        mNu90all(i,b,:)  = d.mNu90; %first is statistics, last one is all systematics combined
    end
end
SysEffects = d.mySysEffects;
%% Plot
MACE_Ba_T = 7*1e-04;
BKG = GetBackground('MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'Anchor6G',Anchor6G);

f66 = figure(66);
set(f66, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);

mNu90  = (squeeze(mNu90all(:,find(MACE_Ba_T_all==MACE_Ba_T),:))); %sensitivites for given Ba
SysErr90 = sqrt(mNu90(:,:).^2-mNu90(:,1).^2); %SysErr on mNu^2 90% C.L.
SysErr90(:,1) = mNu90(:,1); %stat error
SysErr90(:,end) = 0; 
t30 = SysErr90(1,:).^2;
t45 = SysErr90(2,:).^2;
t60 = SysErr90(3,:).^2;
b = barh(range_all,[t30;t45;t60],'stacked');
b(1).LineStyle = '--';
b(2).LineStyle = 'none';b(3).LineStyle = 'none'; b(4).LineStyle = 'none';b(5).LineStyle = 'none';
b(1).FaceColor = rgb('White');
b(2).FaceColor = rgb('FireBrick');
b(3).FaceColor = rgb('GoldenRod');
b(4).FaceColor = rgb('CadetBlue');
b(5).FaceColor = rgb('Navy');
legend({'Stat','TC','FSD','RF'})
legend boxoff
yticklabels({sprintf('30 eV');sprintf('45 eV');sprintf('60 eV')});
ytickangle(90)
xlabel(sprintf('\\sigma^2_{sys}(m_{\\nu^2}) 90%% C.L. (eV^2)'));
ylabel('Scan energy range below Endpoint');
PrettyFigureFormat;
set(gca,'FontSize',14);
SysUncertainties = [sprintf('Response Function: '),...
    sprintf('\\Delta \\rhod\\sigma = %.1f%%, ',100*sqrt(s.WGTS_CD_MolPerCm2_RelErr^2+s.ISXsection_RelErr^2)),...
    sprintf('\\DeltaB_a = %.0f \\muT, \\DeltaB_{T/max} = %.1f%%, ',1e6*s.MACE_Ba_T_Err,100*s.MACE_Bmax_T_RelErr),...
    sprintf('Energy Loss Uncertainties from Eur. Phys. J. D 10, 39–52'),...
    sprintf('\nFinal State Distribution:'),...
    sprintf(' Normalization %.0f %%, ',s.FSDNorm_RelErr*100),...
    sprintf('Bin-to-Bin uncorrelated %.0f %% (GS), %.0f %% (ES)',100*s.FSDShapeGS_RelErr,100*s.FSDShapeES_RelErr)];
a=annotation('textbox', [0.14 0.1 1 0.1], ...
    'String', SysUncertainties, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left');
a.FontSize=11;a.FontWeight='bold';
title(sprintf('Nominal KATRIN %.0f years \nB_a = %.0fG, B_{T/max} = %.0f %%, Background = %.0f mcps',...
    TimeSec/(365*24*60*60), MACE_Ba_T*1e4,WGTS_B_T/3.6*100,BKG*1e3));

save_name = sprintf('%s_SensitivityNominal_NuMassScan_SystematicBreakdown_SysmNuSq',SysBudget,WGTS_B_T/3.6*100);
print(f66,['./plots/png/',save_name,'.png'],'-dpng');