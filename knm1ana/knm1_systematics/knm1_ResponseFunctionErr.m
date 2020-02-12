% Diagnosis plots for sys uncertainty on response function (eloss, b-fields, cd, cross section)
% since uncertainty on product: oNly 1 uncertainty given -> has same effect (in formula always product)

%% settings
RunList = 'KNM1';
nTrials = 1000;
RecomputeFlag = 'OFF';

% Init Model Object and covariance matrix object
if ~exist('A','var')
    A = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','DataType','Twin','SysBudget',2);
end

SysEffects = struct('RF_RX','ON','RF_EL','OFF','RF_BF','OFF');
A.ComputeCM('SysEffects',SysEffects,'nTrials',nTrials,'BkgCM','OFF')
CM = A.FitCM_Obj;
CM.ComputeCM_RF;
%% load response function info :
if strcmp(CM.SysEffect.RF_EL,'ON')
    maxE_el = 500;
else
    maxE_el = 9288;
end
ELossBinStep = 0.2;
RFBinStep = 0.04;


if strcmp(CM.SysEffect.RF_BF,'ON') && strcmp(CM.SysEffect.RF_RX,'ON') && strcmp(CM.SysEffect.RF_EL,'ON')
    rf_path = [getenv('SamakPath'),sprintf('/inputs/CovMat/RF/LookupTables/RF/')];
    rf_filename = sprintf('RFInfo_%u-Trials_%u-NIS_%s_%g-molPercm2_%.2gerr_xsection%gerr_%.2fT-Bmax_%.2fT-Bs_BT%gerr_%.0fG-Ba_%.2gerr_%s_elossBinning-%.geV-%.2geVStep_RFbinStep-%.2f.mat',...
        CM.nTrials, CM.StudyObject.NIS, CM.StudyObject.TD, CM.StudyObject.WGTS_CD_MolPerCm2,...
        CM.WGTS_CD_MolPerCm2_RelErr,CM.ISXsection_RelErr,...
        CM.StudyObject.MACE_Bmax_T,CM.StudyObject.WGTS_B_T,CM.WGTS_B_T_RelErr,...
        floor(CM.StudyObject.MACE_Ba_T*1e4),CM.MACE_Ba_T_RelErr,CM.StudyObject.ELossFlag,...
        maxE_el,ELossBinStep,RFBinStep);
else
    rf_path = [getenv('SamakPath'),sprintf('/inputs/CovMat/RF/LookupTables/PartVar/RFPartVar/')];
    rf_filename = sprintf('RFInfo_%u-Trials_%u-NIS_BField-%s_RhoXsection-%s_Eloss-%s_%s_%g-molPercm2_%.2gerr_xsection%gerr_%.2fT-Bmax_%.2fT-Bs_BT%gerr_%.0fG-Ba_%.2gerr_elossBinning-%.geV-%.2geVStep_RFbinStep-%.2f.mat',...
        CM.nTrials, CM.StudyObject.NIS,CM.SysEffect.RF_BF,CM.SysEffect.RF_RX,CM.SysEffect.RF_EL,...
        CM.StudyObject.TD,CM.StudyObject.WGTS_CD_MolPerCm2,...
        CM.WGTS_CD_MolPerCm2_RelErr,CM.ISXsection_RelErr,...
        CM.StudyObject.MACE_Bmax_T,CM.StudyObject.WGTS_B_T,CM.WGTS_B_T_RelErr,...
        floor(CM.StudyObject.MACE_Ba_T*1e4),CM.MACE_Ba_T_RelErr,...
        maxE_el,ELossBinStep,RFBinStep);
    if strcmp(CM.SysEffect.RF_EL,'ON')
        rf_filename = [rf_filename,sprintf('_%s',CM.StudyObject.ELossFlag)];
    end
   
end

rf_file = strcat(rf_path,rf_filename);
d = importdata(rf_file);

%% plot
f1 = figure('Renderer','opengl');
set(f1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
qu = 10;
[l,a] = boundedline(d.Te-d.qU(qu),d.RFmean(qu,:),10.*d.RFstd(qu,:));
l.Color = A.PlotColor; l.LineWidth=2; a.FaceColor = A.PlotColor; a.FaceAlpha = 0.2;
PrettyFigureFormat;
xlim([0 30]);
xlabel('surplus energy (eV)')
ylabel('transmission probability')
str_leg = '';savelabel = '';
if strcmp(CM.SysEffect.RF_RX,'ON')
str_cd = sprintf('\\Delta\\rhod\\sigma  = %.1f %% \n',CM.WGTS_CD_MolPerCm2_RelErr*100);
str_leg = [str_leg,str_cd];
savelabel = [savelabel,sprintf('%.2f-rhodsigmaErr',CM.WGTS_CD_MolPerCm2_RelErr*100)];
end
if strcmp(CM.SysEffect.RF_BF,'ON')
str_bf = sprintf('\\DeltaB_s     = %.1f %% \n\\DeltaB_a     = %.1f %% \n\\DeltaB_{max} = %.1f %%',...
               CM.WGTS_B_T_RelErr*100,CM.MACE_Ba_T_RelErr*100,CM.MACE_Bmax_T_RelErr*100);
str_leg = [str_leg,str_bf]; 
savelabel  = [savelabel,sprintf('_%.2f-BmaxErr',CM.WGTS_B_T_RelErr*100)];
end  

leg = legend([l,a],['response function ',A.ModelObj.TD],sprintf('1\\sigma band x 10'));
leg.Location= 'northwest';
legend boxoff
 a=annotation('textbox', [0.7 0.1 0.5 0.23], ...
                    'String', str_leg, ...
                    'EdgeColor', 'none', ...
                    'HorizontalAlignment', 'left');
                a.FontSize=20;a.FontWeight='bold';
                
% save plot
savepath = [getenv('SamakPath'),'knm1ana/knm1_systematics/plots/'];
savename = [savepath,sprintf('knm1_rf_%s_%s.png',RunList,savelabel)];

if ~exist(savepath,'dir')
    system(['mkdir ',savepath]);
end
print(gcf,savename,'-dpng','-r400');

%%
histogram(d.WGTS_CD_MolPerCm2_local)
PrettyFigureFormat;
%xlabel(sprintf('B_s (T)'));
xlabel(sprintf('\\rho d \\sigma'));
box off;
print(gcf,[savepath,'RXDist.png'],'-dpng','-r400');