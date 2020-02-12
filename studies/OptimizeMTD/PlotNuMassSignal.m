addpath(genpath('../../../Samak2.0'));
E0=18575;

%Simulation Model
range        = '30';
TD           = 'Flat30_fineBin';%['Flat',range];
MACE_Ba_T    = 7e-04;
WGTS_B_T     = 0.7*3.6;
MACE_Bmax_T  = (WGTS_B_T/3.6)*6;
BkgReduction = 1;
FPD_ROIlow = 14;
% 10mcps
S10mcps = ref_TBD_NominalKATRIN('BKG_RateAllFPDSec',0.01,'MACE_Ba_T',MACE_Ba_T,...
    'WGTS_B_T',WGTS_B_T,'MACE_Bmax_T',MACE_Bmax_T,...
    'TimeSec',5*180*86400,'TD',TD,'mnuSq_i',(0)^2);
S10mcps.ComputeTBDDS; S10mcps.ComputeTBDIS;

% Nominal
BKG_RateSec = GetBackground('MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'FPD_ROIlow',FPD_ROIlow)*1e-03;
S = ref_TBD_NominalKATRIN('BKG_RateAllFPDSec',BKG_RateSec./BkgReduction,'MACE_Ba_T',MACE_Ba_T,...
    'WGTS_B_T',WGTS_B_T,'MACE_Bmax_T',MACE_Bmax_T,...
    'TimeSec',5*180*86400,'TD',TD,'mnuSq_i',0);
S.ComputeTBDDS; S.ComputeTBDIS;

%%
f1 = figure('Renderer','opengl')
set(f1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]); %'OuterPosition', [0, 0.04, 1, 0.96]);
S10mcpsM = ref_TBD_NominalKATRIN('BKG_RateAllFPDSec',0.01,'MACE_Ba_T',MACE_Ba_T,...
    'WGTS_B_T',WGTS_B_T,'MACE_Bmax_T',MACE_Bmax_T,...
    'TimeSec',5*180*86400,'TD',TD,'mnuSq_i',(0.5)^2);
S10mcpsM.ComputeTBDDS; S10mcpsM.ComputeTBDIS;
qUinter = linspace(min(S10mcpsM.qU-S10mcpsM.Q),max(S10mcpsM.qU-S10mcpsM.Q),1000);
Ratiointer = interp1(S10mcpsM.qU-S10mcpsM.Q,(S10mcpsM.TBDIS./S10mcps.TBDIS),qUinter,'spline');
pNom = plot(qUinter,Ratiointer,'LineWidth',5,'Color','Blue');
%plot(S10mcpsM.qU-S10mcpsM.Q,(S10mcpsM.TBDIS./S10mcps.TBDIS),'LineWidth',5,'Color','Blue');
hold on
pBKG = cell(5,1);
BKG = zeros(5,1);
legEntry = cell(5,1);
for i=3:2:12
    BKG((i-1)/2) = GetBackground('MACE_Ba_T',i*1e-4,'WGTS_B_T',WGTS_B_T,'FPD_ROIlow',FPD_ROIlow)*1e-03;
    Sref{i} = ref_TBD_NominalKATRIN('BKG_RateAllFPDSec',BKG((i-1)/2)./BkgReduction,'MACE_Ba_T',i*1e-4,...
        'WGTS_B_T',WGTS_B_T,'MACE_Bmax_T',MACE_Bmax_T,...
        'TimeSec',5*180*86400,'TD',TD,'mnuSq_i',(0.)^2);
    Sm{i} = ref_TBD_NominalKATRIN('BKG_RateAllFPDSec',BKG((i-1)/2)./BkgReduction,'MACE_Ba_T',i*1e-4,...
        'WGTS_B_T',WGTS_B_T,'MACE_Bmax_T',MACE_Bmax_T,...
        'TimeSec',5*180*86400,'TD',TD,'mnuSq_i',(0.5)^2);
    Sm{i}.ComputeTBDDS; Sm{i}.ComputeTBDIS;
    Sref{i}.ComputeTBDDS; Sref{i}.ComputeTBDIS;
    pBKG{(i-1)/2} = plot(Sm{i}.qU-S.Q,(Sm{i}.TBDIS./Sref{i}.TBDIS),'LineWidth',5,'Color',[0.1*i/2 0 0]);
    %    plot(Sm{i}.qU-S.Q,-(Sm{i}.TBDIS./S.TBDIS-1),'LineWidth',5,'Color','Black');
    legEntry{(i-1)/2} = sprintf('BKG = %.0fmcps, Ba = %.0fG',BKG((i-1)/2)*1e3,i);
end

%% design settings
grid on;
hold off
xlabel('retarding potential qU - 18575 (V)'); % xlabel
ylabel('signal amplitude'); %ylabel
xlim([-30 5]);
ylim([0.967 1]);
yticks(0.97:0.01:1)
PrettyFigureFormat
set(gca,'FontSize',18);
l1 = sprintf('BKG = 10 mcps,  Ba = 3G');
leg = legend([pBKG{:},pNom],{legEntry{1:end},l1},...
    'Location','southwest');
legend boxoff
%% save
print(f1,'./plots/NuMassSignal.png','-dpng');
publish_figurePDF(f1,'./plots/NuMassSignal.pdf')