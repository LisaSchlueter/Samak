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

%% Title etc
PlotFontSize = 22;
switch Mode
    case 'ground state'
        ShapeErr   = FSDShapeGS_RelErr*100;
        SigmaFSD = sqrt(varFSDvar_GS)./meanFSDvar_GS*100;
        title1 = sprintf('%s Ground State Sampling (%u Samples) \nNormalization Uncertainty %u %%',Isotopologue,nTrials, FSDNorm_RelErr*100);
        title2 = '';
        save_name = sprintf('%s_FSDVariance_GroundStates_%.0fSamples',Isotopologue,nTrials);
    case 'all'
        ShapeErr = FSDShapeES_RelErr*100;
        SigmaFSD = sqrt(varFSDvar)./meanFSDvar*100; 
        title1 = sprintf('%s Sampling (%u Samples) \nNormalization Uncertainty %u %%',Isotopologue,nTrials, FSDNorm_RelErr*100);
        title2 = sprintf(', Bin-to-Bin %.1f%% (GS)',FSDShapeGS_RelErr(1)*100);
        save_name = sprintf('%s_FSDVariance_%.0fSamples',Isotopologue,nTrials);
end

%% Variance in dependence of GS Error
fig681 = figure(681);
set(fig681, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.7, 0.7],'Renderer','opengl');
plot(ShapeErr,SigmaFSD,'-o','MarkerSize',8,'LineWidth',5,'Color',rgb('FireBrick'),'MarkerFaceColor',rgb('FireBrick'));
hold on;

PrettyFigureFormat;
ylim([min(SigmaFSD) max(SigmaFSD)])
xlim([min(ShapeErr) max(ShapeErr)]);
title([title1,title2]);
set(gca,'FontSize',PlotFontSize);
xlabel(sprintf('\\sigma_{Fluct} (%%)'),'FontSize',PlotFontSize+5);
ylabel('\sigma_{FSD} (%)','FontSize',PlotFontSize+5);
grid on;
switch Mode
    case 'ground state'
        plot(ShapeErr,ones(numel(ShapeErr),1),'--k','LineWidth',2);
    case 'all'
        plot(ShapeErr,2*ones(numel(ShapeErr),1),'--k','LineWidth',2);
        plot(ShapeErr,3*ones(numel(ShapeErr),1),'--k','LineWidth',2);
end

save_dir = sprintf('./plots/FSD_diagnostic/');

export_fig([save_dir,'png/',save_name,'.png']);
export_fig([save_dir,'fig/',save_name,'.fig']);
publish_figurePDF(fig681,[save_dir,'pdf/',save_name,'.pdf']);