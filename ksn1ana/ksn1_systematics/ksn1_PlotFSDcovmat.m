% display FSD systematics for KSN1
% Lisa, May 2020
%% settings
SysBudget     = 22; % 22 = knm1 default
myEffect      = 'FSD';
RecomputeFlag = 'ON';
nTrials       =  501;
FSDFlag = 'SibilleFull';
%% collect model + covmat settings
RunArg = {'FSDFlag',FSDFlag,...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...      % switch on later by hand
    'RunList','KNM1',...
    'fixPar','E0 Norm Bkg',... % free parameter
    'DataType','Twin',...      % MC twins
    'TwinBias_Q',18573.7,...  % twin endpoint
    'SysBudget',SysBudget,...
    'RingMerge','Full',...
    'NonPoissonScaleFactor',1,...
    'AngularTFFlag','OFF'};

CmArg = {'BkgCM','OFF',...
    'RecomputeFlag',RecomputeFlag,...
    'nTrials',nTrials,...
    'SysEffects',struct(myEffect,'ON'),...
    };
%% init model
M = MultiRunAnalysis(RunArg{:});
M.chi2 = 'chi2CMShape';

%% label
plotdir = [getenv('SamakPath'),'ksn1ana/ksn1_systematics/plots/'];
MakeDir(plotdir);

%% calculate/load FSD cov mat
M.ComputeCM(CmArg{:});

%% plot 1 & 2 : covmat, correlation, convergence 
M.FitCM_Obj.PlotCM('qUWindowIndexMax',10,'saveplot','ON','Convergence',...
    'ON','Mode','Shape','PlotEffect',myEffect,'savedir',plotdir);
M.FitCM_Obj.PlotCorr('qUWindowIndexMax',10,'saveplot','ON',...
    'savedir',plotdir,'savename',myEffect);

%% plot 3: FSD with 1 sigma error band
d = importdata(M.FitCM_Obj.CovMatFile);
Iso = 'TT';
ProbSamples = squeeze(d.([Iso,'_P_norm']));
Energy  = d.obj.StudyObject.TTexE;
Prob    = mean(ProbSamples,2);
ProbErr = std(ProbSamples,0,2);

figureHandle = figure('Units','normalized','Position',[0.1,0.1,0.7,0.5]);
[line,area] = boundedline(Energy,Prob,ProbErr);
area.FaceAlpha = 1;
area.FaceColor = rgb('PowderBlue');
line.LineWidth = 2;
line.Color = rgb('DodgerBlue');
PrettyFigureFormat('FontSize',24);
xlabel('Excitation energy (eV)');
ylabel('Probability');
xlim([0 100]);
ylim([0 0.03]);
leg = legend([line,area],...
    sprintf('%s %s',Iso, FSDFlag),sprintf('1\\sigma error band'),...
    'EdgeColor',rgb('Silver'));

plotname = sprintf('%sFSDerrband_%s%s.png',plotdir,FSDFlag,Iso);
print(plotname,'-dpng','-r450');
