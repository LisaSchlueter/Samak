%% settings
RunList = 'KNM1';
NonPoissonScaleFactor=1.064;
nTrials = 1000;
RecomputeFlag = 'OFF';

% Init Model Object and covariance matrix object
if ~exist('A','var')
    A = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','DataType','Real','exclDataStart',14);
end
qUmin = round(A.ModelObj.qU(A.exclDataStart));
qUmax = round(A.ModelObj.qU(end));
range = round(A.ModelObj.qU(A.exclDataStart)-A.ModelObj.Q_i);
A.Fit;

SysEffects = struct('RF_RX','OFF','RF_EL','OFF','RF_BF','OFF');
CM = CovarianceMatrix('StudyObject',A.ModelObj, 'nTrials',nTrials,'SysEffect',SysEffects,...
    'RecomputeFlag',RecomputeFlag,'NonPoissonScaleFactor',NonPoissonScaleFactor);
CM.ComputeCM_Background;


%% plot of the variances
myMainTitle = sprintf('KATRIN - KNM1 Signal / Background Variances - %d runs - [%.0f - %0.f] eV',...
    numel(A.RunList),qUmin,qUmax);
maintitle   = myMainTitle;
savefile    = sprintf('BackgroundVariance_%d_%.0feVbelowE0_1.png',...
    numel(A.RunList),abs(range));
fig1      = figure('Name','KATRIN - KNM1 Signal / Background Variances','NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
data=stairs(A.RunData.qU(A.exclDataStart:end)-A.ModelObj.Q_i,...
    (A.RunData.TBDIS(A.exclDataStart:end))./(A.RunData.TBDIS(A.exclDataStart:end)),'LineWidth',4);
hold on
model=stairs(A.RunData.qU(A.exclDataStart:end)-A.ModelObj.Q_i,...
    A.ModelObj.TBDIS(A.exclDataStart:end)./(A.RunData.TBDIS(A.exclDataStart:end)),'LineWidth',4);
modelS=stairs(A.RunData.qU(A.exclDataStart:end)-A.ModelObj.Q_i,...
    (A.ModelObj.TBDIS(A.exclDataStart:end)-(A.ModelObj.BKG_RateSec.*A.RunData.qUfrac(A.exclDataStart:end).*A.RunData.TimeSec))./(A.RunData.TBDIS(A.exclDataStart:end)),'LineWidth',4);
modelB=stairs(A.RunData.qU(A.exclDataStart:end)-A.ModelObj.Q_i,...
    (A.ModelObj.BKG_RateSec.*A.RunData.qUfrac(A.exclDataStart:end).*A.RunData.TimeSec)./(A.RunData.TBDIS(A.exclDataStart:end)),'LineWidth',4);
shapeB=stairs(A.RunData.qU(A.exclDataStart:end)-A.ModelObj.Q_i,...
    (diag(CM.CovMat(A.exclDataStart:end,A.exclDataStart:end)))./(A.RunData.TBDIS(A.exclDataStart:end)),'LineWidth',4);
hold off
leg = legend([data model modelS modelB shapeB],...
    ' Data',...
    ' Model',...
    ' Tritium Model',...
    ' Background Rate Model',...
    ' Background Shape Model (correlated)'...
    );
leg.Location= 'southeast';
legend boxoff
set(gca,'YScale','log')
xlabel(sprintf('retarding potential - %.1f eV',A.ModelObj.Q_i));
ylabel('\sigma^2 / \sigma^2_{Data}')  
ylim([1e-3 2])
PrettyFigureFormat
set(gca,'FontSize',24);

% save plot
savepath = [getenv('SamakPath'),'knm1ana/BackgroundStudies/plots/'];
savename = [savepath,savefile];
if ~exist(savepath,'dir')
    system(['mkdir ',savepath]);
end
print(gcf,savename,'-dpng','-r400');


%% plot of the correlation matrix - background shape cov matrix 
A.exclDataStart=14;
myMainTitle = sprintf('KATRIN - KNM1 Background Shape Uncertainty Correlation Matrix - %d runs - [%.0f - %0.f] eV',...
    numel(A.RunList),qUmin,qUmax);
maintitle   = myMainTitle;
savefile    = sprintf('BackgroundVariance_%d_%.0feVbelowE0_2.png',...
    numel(A.RunList),abs(range));
fig1      = figure('Name','Background Shape Uncertainty Correlation Matrix','NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
a=annotation('textbox', [0 0.91 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=22;a.FontWeight='bold';
corplot((CM.CovMat(A.exclDataStart:end,A.exclDataStart:end)));
set(gca,'XAxisLocation','top');
n = size(CM.CovMat(A.exclDataStart:end,A.exclDataStart:end),1)+1;
set(gca,'XTick',.5+(1:2:(n-1)),'XTickLabel',round(A.RunData.qU(A.exclDataStart:2:A.exclDataStart+(n-1))-A.ModelObj.Q_i,1));
set(gca,'YTick',.5+(1:2:(n-1)),'YTickLabel',round(A.RunData.qU(A.exclDataStart:2:A.exclDataStart+(n-1))-A.ModelObj.Q_i,1));
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
