% Plot Covariance Matrix for TASR single systematic effect
% systematics setting
myEffect      = 'TASR';
RecomputeFlag = 'OFF';
nTrials       = 1000; 
SysBudget     = 31; % 31= knm2 preliminary input

% model setting
RunList   = 'KNM2_Prompt';
AnaFlag   = 'StackPixel';%'Ring'; % FPD segmentation

RunArg = {'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag',AnaFlag,...
    'chi2','chi2Stat',...      % switch on later by hand
    'RunList',RunList,...      
    'fixPar','E0 Norm Bkg',... % free parameter
    'DataType','Twin',...      % MC twins
    'TwinBias_Q',18573.56,...  % twin endpoint
    'SysBudget',SysBudget,...
    'RingMerge','Full'};        

CmArg = {'BkgCM','OFF',...
    'SysEffects',struct(myEffect,'ON'),...
    'RecomputeFlag',RecomputeFlag,...
    'nTrials',nTrials,...
    };
%% init model
M = MultiRunAnalysis(RunArg{:});

%% calculate covariance matrix
M.chi2 = 'chi2CMShape';
M.ComputeCM(CmArg{:});

%% display and save to plots
M.FitCM_Obj.PlotCM('qUWindowIndexMax',10,'saveplot','ON','Convergence','OFF');

%% Plot subrunwise activity fluctuations
savedir = [getenv('SamakPath'),'knm2ana/knm2_systematics/plots/'];
MakeDir(savedir);
[WGTS_TASR_RelErr,SubRunActivity,TASR_CorrMat] =  M.Get_DataDriven_RelErr_TASR;
nqU = M.ModelObj.nqU;
nruns = size(SubRunActivity,2);
x          = reshape(repmat(M.ModelObj.qU(:,1),[1,nruns]),[nqU*nruns,1]);
y          = reshape(SubRunActivity,[nqU*nruns,1]);
y = y-mean(y);
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
l = plot(x-18574,zeros(nqU*nruns,1),'-','LineWidth',2,'Color',rgb('Silver'));
hold on;
sg  = scatter(x-18574,y,'filled');
sg.MarkerFaceAlpha = 0.2;
sg.MarkerFaceColor = rgb('DodgerBlue');
xlabel('Retarding potential - 18574 (eV)');
ylabel('Rel. Tritium activity');
PrettyFigureFormat('FontSize',24);
xlim([-41,0]);
%ylim([0.995, 1.006]);
MeanErrorOfMean = mean(std(SubRunActivity,0,2)./sqrt(nruns));
leg = legend(sg,sprintf('%.0f scans , \\langle\\sigma / {\\surd %.0f}\\rangle = %.2e',...
   nruns,nruns,MeanErrorOfMean));
leg.EdgeColor = rgb('LightGray');
leg.LineWidth = 2;
leg.Location = 'northwest';

savename1 = [savedir,'TASR_RelActivityFluctSubrun.png'];
print(f1,savename1,'-dpng','-r450');
    
%% plot correlation plot
f2 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
[cp,~, c] = corplot(TASR_CorrMat);
PrettyFigureFormat('FontSize',24);
ax = gca;
ax.XTick      = [1.5,nqU+0.5];
ax.XTickLabel = {'-40','0'};
ax.YTick  = [];
set(gca,'XMinorTick','off');
ax.XAxisLocation = 'bottom';
xlabel('Retarding potential -18574 (eV)');
c.Label.String = 'Correlation';
c.Label.FontSize = ax.XLabel.FontSize;
pos = ax.XLabel.Position;
ax.XLabel.Position =  [pos(1),pos(2)-1,pos(3)];
savename2 = [savedir,'TASR_CorrMat.pdf'];
export_fig(f2,savename2);





