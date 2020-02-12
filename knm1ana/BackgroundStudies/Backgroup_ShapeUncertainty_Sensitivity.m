%% settings
RunList               = 'KNM1';
exclDataStart         = 14;
NonPoissonScaleFactor = 1.0;
nTrials               = 1000;
RecomputeFlag         = 'OFF';
BkgShape              = 'ON';
SysEffects            = struct('RF_RX','OFF','RF_EL','OFF','RF_BF','OFF','BkgShape',BkgShape);
chi2                  = 'chi2CMShape';
            
%% Init Model Object and covariance matrix object
%if ~exist('Twin','var')
    Twin = MultiRunAnalysis('RunList',RunList,...
        'chi2',chi2,'DataType','Twin',...
        'exclDataStart',exclDataStart,...
        'fixPar',' 5 6 7 8 9 10 11',...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSDFlag','BlindingKNM1');
%end
qUmin = round(Twin.ModelObj.qU(Twin.exclDataStart));
qUmax = round(Twin.ModelObj.qU(end));
range = round(Twin.ModelObj.qU(Twin.exclDataStart)-Twin.ModelObj.Q_i);

Twin.ComputeCM(...
        'SysEffects',SysEffects,...
        'nTrials',nTrials,...
        'RecomputeFlag',RecomputeFlag,...
        'BkgCM',BkgShape);
    
%% Fit With ASIMOV + Stat-only
Twin.Fit;

%% plot of the correlation matrix - background shape cov matrix 
Twin.exclDataStart=14;
myMainTitle = sprintf('KATRIN - KNM1 Background Shape Uncertainty Correlation Matrix - %d runs - [%.0f - %0.f] eV',...
    numel(Twin.RunList),qUmin,qUmax);
maintitle   = myMainTitle;
savefile    = sprintf('Backgroup_ShapeUncertainty_Sensitivity%d_%.0feVbelowE0_2.png',...
    numel(Twin.RunList),abs(range));
fig1      = figure('Name','Background Shape Uncertainty Correlation Matrix','NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
a=annotation('textbox', [0 0.91 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=22;a.FontWeight='bold';
corplot((Twin.FitCM(Twin.exclDataStart:end,Twin.exclDataStart:end)));
set(gca,'XAxisLocation','top');
n = size(Twin.FitCM(Twin.exclDataStart:end,Twin.exclDataStart:end),1)+1;
set(gca,'XTick',.5+(1:2:(n-1)),'XTickLabel',round(Twin.RunData.qU(Twin.exclDataStart:2:Twin.exclDataStart+(n-1))-Twin.ModelObj.Q_i,1));
set(gca,'YTick',.5+(1:2:(n-1)),'YTickLabel',round(Twin.RunData.qU(Twin.exclDataStart:2:Twin.exclDataStart+(n-1))-Twin.ModelObj.Q_i,1));
set(gca,'xaxisLocation','bottom');
xlabel('eV');ylabel('eV')
xtickangle(45);
PrettyFigureFormat
savepath = [getenv('SamakPath'),'knm1ana/BackgroundStudies/plots/'];
savename = [savepath,savefile];
if ~exist(savepath,'dir')
    system(['mkdir ',savepath]);
end
print(gcf,savename,'-dpng','-r400');
