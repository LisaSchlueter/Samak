% Load MultiRunAnalysis Object for KNM1
if ~exist('Real','var')
    options = {...
        'RunList','KNM1',...
        'exclDataStart',2,...
        'fixPar','1 5 6 7 8 9 10 11',...
        'chi2','chi2Stat',...
        'ELossFlag','KatrinD2',...
        'StackTolerance',1,...
        'NonPoissonScaleFactor',1,...
        'Debug','ON',...
        'AnaFlag','StackPixel',...
        'PixList',[]};
    Real=MultiRunAnalysis('DataType','Real',options{:});
end

if isempty('Real.SingleRunObj')
    Real.LoadSingleRunObj;
end

%% Loop over all KNM1 runs
% Compute IS probabilities
% Store Results in is Arrays
if isempty('is')
    for i=1:1:numel(Real.RunList)
        is(i,:) = Real.SingleRunObj{i}.ComputeISProb;
    end
end

%% Get Mean
mean_is = mean(is,1);

%% Get Stacked IS
stacked_is = Real.ModelObj.ComputeISProb;

%% Plot Distributions
figure('Units', 'pixels','Position', [0 0 1200 1200]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', 'KNM1 Run-wise Inelastic Scattering Probability Distributions', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=20;a.FontWeight='bold';

for i=1:8
    s(i) = subplot(4,2,i);
    h(i)=histogram(is(:,i),'FaceColor',rgb('IndianRed'),'EdgeAlpha',0.5);
    title(sprintf('%d-fold scattering',i-1));
    xlabel('%');
    ylabel('runs');
    
    hold on
    m(i)=line([mean_is(i) mean_is(i)],[0 max(h(i).Values)],'Color',rgb('Red'),'LineWidth',2);
    l(i)=line([stacked_is(i) stacked_is(i)],[0 max(h(i).Values)],'Color',rgb('SteelBlue'),'LineWidth',2);    
    hold off
    
    leg = legend([h(i) l(i) m(i)],sprintf('\\sigma / mean: %.2g %% \n min= %.3g %% \n max= %.3g %%',...
        100*std(is(:,i)/mean(is(:,i))),min(is(:,i)),max(is(:,i))),'stacked','mean','Location','northwest');
    leg.Color = 'none'; legend boxoff; leg.FontSize = 10;
    
    PrettyFigureFormat
end
export_fig(gcf,sprintf('plots/KNM1_RunWiseISprobabilities_Mean_Stacked'),'-m3');

%% Print Latex Table
t = PrintTable('KNM1 Inelastic Scattering Probabilities');
t.addRow('Number of Scattering','Mean Value','Stacked Value','Relative Difference (%)')
for i=1:8
t.addRow(sprintf('%d',i),...
    sprintf('%5.6f %',mean_is(i)),sprintf('%5.6f %',stacked_is(i)),...
    (mean_is(i)-stacked_is(i))/((mean_is(i)+stacked_is(i))/2)*100);
end
t.display;
%t.HasHeader = true;
t.Format = 'tex';
t.Caption = sprintf('KNM1 Inelastic Scattering Probabilities - %.0f runs',numel(Real.RunList)');
t.print;