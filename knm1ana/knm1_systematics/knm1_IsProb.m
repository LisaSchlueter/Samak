% Diagnosis plots for sys uncertainty of column density and inel. cross section -> covariance matrix
% since uncertainty on product: oly 1 uncertainty given -> has same effect (in formula always product)
%% settings
RunList = 'KNM1';
nTrials = 1000;
RecomputeFlag = 'OFF';
% Init Model Object and covariance matrix object
if ~exist('A','var')
    A = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','DataType','Twin','SysBudget',2);
end

SysEffects = struct('RF_RX','ON');
A.ComputeCM('SysEffects',SysEffects,'nTrials',nTrials,'BkgCM','OFF')
CM = A.FitCM_Obj;
%CM = CovarianceMatrix('StudyObject',A.ModelObj, 'nTrials',nTrials,'SysEffect',SysEffects,'RecomputeFlag',RecomputeFlag,...
 %   'ISXsection_RelErr',ISXsection_RelErr,'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr);
 CM.ComputeCM_RF;
 
 %% load inelastic scattering probabilities and plot uncertainties
 if strcmp(CM.SysEffect.RF_BF,'ON') && strcmp(CM.SysEffect.RF_RX,'ON')
     isp_dir =   [getenv('SamakPath'),sprintf('inputs/CovMat/RF/LookupTables/ISProb/')];
 else
     isp_dir =   [getenv('SamakPath'),sprintf('inputs/CovMat/RF/LookupTables/PartVar/ISProbPartVar/')];
 end
 
 isp_common = sprintf('ISProb-LookupTable_%u-Trials_%.0fNIS',CM.nTrials,CM.StudyObject.NIS+1);
 str_cd = sprintf('_%.5g-molPercm2',CM.StudyObject.WGTS_CD_MolPerCm2);
 str_is = sprintf('_IsX%.5gm2',CM.StudyObject.ISXsection);
 if strcmp(CM.SysEffect.RF_RX,'ON')
     str_cd = [str_cd,sprintf('-%.3gerr',CM.WGTS_CD_MolPerCm2_RelErr)];
     str_is = [str_is,sprintf('-%.3gerr',CM.ISXsection_RelErr)];
 end
 str_rx = [str_cd,str_is];
 
 str_bs   = sprintf('_%.3gBs',CM.StudyObject.WGTS_B_T);
 str_bmax = sprintf('_%.3gBmax',CM.StudyObject.MACE_Bmax_T);
 if strcmp(CM.SysEffect.RF_BF,'ON')
     str_bs    = [str_bs,sprintf('-%.3gerr',CM.WGTS_B_T_RelErr)];
     str_bmax  = [str_bmax,sprintf('-%.3gerr',CM.MACE_Bmax_T_RelErr)];
 end
 str_bf = [str_bs,str_bmax];
 
 isp_filename = [isp_common,str_rx,str_bf,'.mat'];
 
 f33 = figure('Renderer','opengl');
 set(f33, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);

d = importdata([isp_dir,isp_filename]);
x = 0:A.ModelObj.NIS;
b = bar(x,mean(d.Pis_m,2),'FaceColor',rgb('Silver'),'EdgeColor',rgb('DimGray'));
hold on;
e = errorbar(x,mean(d.Pis_m,2),std(d.Pis_m'),'s','LineWidth',2,'Color',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'));
hold off;
PrettyFigureFormat;
xlabel('scatterings');
ylabel('probability (%)')

CDsigmaRelErr = round(std(d.WGTS_CD_MolPerCm2_v.*d.ISXsection_v)/mean(d.WGTS_CD_MolPerCm2_v.*d.ISXsection_v)*100,1);
leg = legend([b,e],sprintf('\\rhod = %.3g mol/cm^2 \n\\sigma_{inel} = %.5g m^2',...
    A.ModelObj.WGTS_CD_MolPerCm2,A.ModelObj.ISXsection),sprintf('\\Delta \\rhod\\sigma = %.1f %% ',CDsigmaRelErr));
legend boxoff
set(gca,'YScale','log');

% save plot
savepath = [getenv('SamakPath'),'knm1ana/knm1_systematics/plots/'];
savename = [savepath,sprintf('knm1_isProb_%s_%.2f-rhodsigmaErr.png',RunList,CDsigmaRelErr)];

if ~exist(savepath,'dir')
    system(['mkdir ',savepath]);
end
print(gcf,savename,'-dpng','-r400');
%%
%fprintf('%.3f \n',mean(d.Pis_m,2)')
%fprintf('+/- %.3f \n',std(d.Pis_m'));


