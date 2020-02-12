% Configuration
nTrials = 10000;
Isotopologue = 'TT';
Mode = 'ground state';

%% Load
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
%% Title legend etc preparation
IndexFSDRelErr = 6;
PlotFontSize = 22;
title_norm  = sprintf('normalization uncertainty 1%% \n');
title_shapeGS = sprintf('bin-to-bin uncorrelated uncertainty (GS) %.1f%%',...
    FSDShapeGS_RelErr(IndexFSDRelErr)*100);
switch Mode
    case 'ground state'
        PlotVar =  FSDvar_GS(:,IndexFSDRelErr);
        maintitle=  sprintf('\\mu=%.3f eV^2, \\sigma=%.3feV^2',...
            meanFSDvar_GS(IndexFSDRelErr),sqrt(varFSDvar_GS(IndexFSDRelErr)));
        title_shapeES = '';
        sample_leg = sprintf('\\sigma_{FSD}=%.1f %%',sqrt(varFSDvar_GS(IndexFSDRelErr))./meanFSDvar_GS(IndexFSDRelErr)*100);
        save_name = sprintf('Variance_%s_%.1f_%.0fSamples',Mode,FSDShapeGS_RelErr(IndexFSDRelErr)*100,nTrials);
        adim =  [0.0 0.2 1 0];
    case 'all'
        PlotVar = FSDvar(:,IndexFSDRelErr);
        maintitle=  sprintf('\\mu=%.1f eV^2, \\sigma=%.1feV^2',...
            meanFSDvar(IndexFSDRelErr),sqrt(varFSDvar(IndexFSDRelErr)));
        title_shapeES = sprintf(' (ES) %.1f%%', FSDShapeES_RelErr(IndexFSDRelErr)*100);
        sample_leg = sprintf('\\sigma_{FSD}=%.1f %%',sqrt(varFSDvar(IndexFSDRelErr))./meanFSDvar(IndexFSDRelErr)*100);
        save_name = sprintf('Variance_%s_%.1f_%.0fSamples',Mode,FSDShapeES_RelErr(IndexFSDRelErr)*100,nTrials);
        adim = [0 0.2 0.97 0];
end
mytitle = [title_norm,title_shapeGS,title_shapeES];
%% Plot
close
fig67 = figure('Name','Variances','Renderer','opengl');

set(fig67, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.7, 0.7]);
    h1 = nhist(PlotVar,'color',rgb('CadetBlue'),'linewidth',0.01);
    l1 = legend(sample_leg,'Location','northwest');
    legend boxoff
    xlabel('variance \sigma^2 (eV²)'); ylabel('samples');
    title(mytitle);
    PrettyFigureFormat;
    set(gca,'FontSize',PlotFontSize);
    l1.FontSize = PlotFontSize+2;
 
a=annotation('textbox', adim, ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=PlotFontSize;a.FontWeight='bold';    
    
save_dir = sprintf('./plots/FSD_diagnostic/');

export_fig([save_dir,'png/',save_name,'.png']);
export_fig([save_dir,'fig/',save_name,'.fig']);
publish_figurePDF(fig67,[save_dir,'pdf/',save_name,'.pdf']);
%close;