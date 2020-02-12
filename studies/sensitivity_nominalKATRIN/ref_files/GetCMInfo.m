function GetCMInfo_RF(CMObj,SysBudget)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(CMObj.SysEffect.RF_BF,'ON') && strcmp(CMObj.SysEffect.RF_RX,'ON') && strcmp(CMObj.SysEffect.RF_EL,'ON')
    rf_path = sprintf('../../inputs/CovMat/RF/LookupTables/RF/');
    rf_filename = sprintf('RF-LookupTable_%u-Trials_%u-NIS_%s_%u-Runs_%g-molPercm2_%.2ferr_xsection%gerr_%.2f-Bmax_%.2f-Bs_BT%gerr.mat',...
        CMObj.nTrials, CMObj.StudyObject.NIS, CMObj.StudyObject.TD,CMObj.nRuns, CMObj.StudyObject.WGTS_CD_MolPerCm2,...
        CMObj.WGTS_CD_MolPerCm2_RelErr,CMObj.ISXsection_RelErr,...
        CMObj.StudyObject.MACE_Bmax_T,CMObj.StudyObject.WGTS_B_T,CMObj.WGTS_B_T_RelErr);
    rf_file = strcat(rf_path,rf_filename);
else
    rf_path = sprintf('../../inputs/CovMat/RF/LookupTables/PartVar/RFPartVar/');
    rf_filename = sprintf('RF-LookupTable_%u-Trials_%u-NIS_BField-%s_RhoXsection-%s_Eloss-%s_%s_%u-Runs_%g-molPercm2_%.2ferr_xsection%gerr_%.2f-Bmax_%.2f-Bs_BT%gerr.mat',...
        CMObj.nTrials, CMObj.StudyObject.NIS,CMObj.SysEffect.RF_BF,CMObj.SysEffect.RF_RX,CMObj.SysEffect.RF_EL,...
        CMObj.StudyObject.TD,CMObj.nRuns,CMObj.StudyObject.WGTS_CD_MolPerCm2,...
        CMObj.WGTS_CD_MolPerCm2_RelErr,CMObj.ISXsection_RelErr,...
        CMObj.StudyObject.MACE_Bmax_T,CMObj.StudyObject.WGTS_B_T,CMObj.WGTS_B_T_RelErr);
    rf_file = strcat(rf_path,rf_filename);
end

if exist(rf_file,'file')==2
    fprintf('--------------------------------------------------------------------------\n')
    fprintf('Loading Response Functions from File \n')
    RF_CM=importdata(rf_file); %Response Functions and Covariance Matrix
else 
    return
end

%% plots
f44 = figure(44);
set(f44, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 1]);

Systitle = '';
SysNames = fieldnames(CMObj.SysEffect);
for i=1:numel(SysNames)
   if strcmp(CMObj.SysEffect.(SysNames{i}),'ON')
Systitle = [Systitle,SysNames{i}];
   end
end
maintitle = sprintf('Variation of input Parameters - Systematic Effect %s',strrep(Systitle,'_',' '));
names = fieldnames(RF_CM);
for i=1:numel(names)-1
        subplot(3,2,i);
        histogram(RF_CM.(names{i}),'FaceColor',rgb('CadetBlue'));
        legend(sprintf('\\sigma_{rel}=%.2f %%',100*std(RF_CM.(names{i}))/mean(RF_CM.(names{i}))));
        PrettyFigureFormat;   
        format long;
        RF1_tmp = RF_CM.(names{i});
   if all(RF_CM.(names{i})==RF1_tmp(1)) %when all the same    
        xlabel([strrep(strrep(names{i},'_',' '),'local',''),' fixed']);
   else
        xlabel(strrep(strrep(names{i},'_',' '),'local',''))
    end
end
subplot(3,2,6);
RF_Samples = RF_CM.ResponseFunction(:,1:100,:);
RF_Samples(isnan(RF_Samples))=0;
RF_i = CMObj.StudyObject.RF;
for qu = 1:CMObj.StudyObject.nqU
plot(CMObj.StudyObject.Te-CMObj.StudyObject.qU(qu),squeeze(RF_Samples(qu,:,:))'-RF_i(:,qu));
xlabel('E_{kin} - qU eV');
 ylabel('RF_{Samples} - RF_{ini} (Transmission Probability %)');
PrettyFigureFormat;
legend boxoff;
end
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=18;a.FontWeight='bold';

save_plot = sprintf('CovMatInfo_ResponseFunction_%s_SysBudget%s',Systitle,SysBudget);
if ~exist('./CovMatInfo/plots/png/','dir')
    mkdir ./CovMatInfo/plots/png/
    mkdir ./CovMatInfo/plots/pdf/
    mkdir ./CovMatInfo/plots/fig/
end
print(f44,['./CovMatInfo/plots/png/',save_plot,'.png'],'-dpng');
publish_figurePDF(f44,['./CovMatInfo/plots/pdf/',save_plot,'.pdf']);
savefig(f44,['./CovMatInfo/plots/fig/',save_plot,'.fig']);

end

