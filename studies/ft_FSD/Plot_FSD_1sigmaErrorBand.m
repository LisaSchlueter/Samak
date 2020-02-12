% Configuration
nTrials = 1000;
Isotopologue = 'TT';
Mode = 'ground state';
PlotFontSize = 18;

% Load
switch Mode
    case 'ground state'
save_name = sprintf('./results/%s_FSDdiagnostic_GroundState_%.0fSamples.mat',Isotopologue,nTrials);
    case 'all'
save_name = sprintf('./results/%s_FSDdiagnostic_%.0fSamples.mat',Isotopologue,nTrials);
end

if exist(save_name,'file')
    load(save_name);
    close
else
    fprintf('file doesnt exist! \n');
    return
end


%% plot
 IndexFSDRelErr = 6;
 title_norm  = sprintf('normalization uncertainty 1%% \n');
 title_shapeGS = sprintf('bin-to-bin uncorrelated uncertainty (GS) %.1f%%', FSDShapeGS_RelErr(IndexFSDRelErr)*100);
 switch Mode
     case 'ground state'
         title_shapeES = '';
         xmax =  4.9;%A.ModelObj.TTexE(A.ModelObj.TTGSTh);
         save_name = sprintf('SampleFSD_GroundStates_%.2gGSerr_%.0fSamples',FSDShapeGS_RelErr(IndexFSDRelErr),nTrials);
     case 'all'
         title_shapeES = sprintf(' (ES) %.1f%%', FSDShapeES_RelErr(IndexFSDRelErr)*100);
         xmax = max(Energy);
         save_name = sprintf('%s_SampleFSD_AllStates_%.2gESerr_%.0fSamples',Isotopologue,FSDShapeES_RelErr(IndexFSDRelErr),nTrials);
 end
 mytitle = [title_norm,title_shapeGS,title_shapeES];


fig55 = figure('Name','FSD','Renderer','opengl');
set(fig55, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.7, 0.7]);
[l,p] = boundedline(Energy,FSD_P_norm_mean(:,IndexFSDRelErr)*100 ,FSD_P_norm_std(:,IndexFSDRelErr)*100);
l.LineWidth = 1.5; l.Color = rgb('CadetBlue');
p.FaceColor = rgb('CadetBlue'); p.FaceAlpha = 0.5;
set(gca,'XScale','log')
xlabel('excitation energy (eV)');
ylabel('probability (%)');
title(mytitle);
sample_leg = sprintf('1\\sigma error band');
theo_leg = sprintf('%s Theory (%s)',strrep(Isotopologue,'TT','T_2'),strrep(FSDCalc,'SAENZ','Saenz et. al'));
legend([l p],theo_leg, sample_leg,'Location','northeast'); legend boxoff;
PrettyFigureFormat;

xlim([0.5 xmax]);
ylim([min(FSD_P_theory*100) 8.2]);%max(FSD_P_norm_mean(:,IndexFSDRelErr)*100)+max(FSD_P_norm_std(:,IndexFSDRelErr))*100]);
set(gca,'FontSize',PlotFontSize);

save_dir = sprintf('./plots/FSD_diagnostic/');
print([save_dir,'png/',save_name,'.png'],'-dpng');
export_fig([save_dir,'fig/',save_name,'.fig']);
publish_figurePDF(fig55,[save_dir,'pdf/',save_name,'.pdf'])

%%
sqrt(varFSDvar_GS(IndexFSDRelErr))./meanFSDvar_GS(IndexFSDRelErr)*100;