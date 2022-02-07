% KNM1 Alternative run lists
% save or load results

savedir = [getenv('SamakPath'),'knm1ana/knm1_AlternativeRunLists/results/'];
fcommon = [savedir,'knm1_AltRunList'];

% settings
range   = 40;         % 40eV range = 27 subruns
% Init Model Object and covariance matrix object
RunAnaArg = { 'chi2','chi2Stat',...
    'DataType','Real',...
    'fixPar','mNu E0 Norm Bkg',... free parameter
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1,...
    'minuitOpt','min ; minos',...
    'FSDFlag','Sibille0p5eV',...
    'ELossFlag','KatrinT2',...
    'SysBudget',22,...
    'AngularTFFlag','OFF'};

%% define run lists
RunLists = {'KNM1',...
    'KNM1upScan',...
    'KNM1downScan',...
    'KNM1-FirstHalfTime',...
    'KNM1-LastHalfTime',...
    'KNM1-MiddleHalfTime',...
    'KNM1-RhoD',...
    'KNM1-AntiRhoD',...
    'KNM1-TpurityLow',...
    'KNM1-TpurityHigh',...
    'KNM1_m183mvRW',...
    'KNM1_m149mvRW',...
    };

RunLists_Labels = {'All scans','Up scans','Down scans',...
                 'First third','Middle third','Last third',...
                 sprintf('\\sigma(\\rhod) < 1'), sprintf('\\sigma(\\rhod) > 1'),...
                 sprintf('High \\epsilon_T'), sprintf('Low \\epsilon_T'),...
                 sprintf('{\\itU}_{RW} = -183 mV'),...
                 sprintf('{\\itU}_{RW} = -149 mV'),...
                 };
%% fit run lists
nRuns_v = zeros(numel(RunLists),1);
mNuSq_v = zeros(numel(RunLists),1);
mNuSqErr_v = zeros(numel(RunLists),1);
chi2min_v = zeros(numel(RunLists),1);
p_v = zeros(numel(RunLists),1);
TimeSec_v = zeros(numel(RunLists),1);

for i=1:numel(RunLists)
    filename = sprintf('%s_%s.mat',fcommon,RunLists{i});
    if exist(filename,'file')
        fprintf('Load file %s\n',filename);
        d = importdata(filename);
        nRuns_v(i) = d.nRuns;
        mNuSq_v(i) = d.FitResult.par(1);
        mNuSqErr_v(i) =  d.FitResult.err(1);
        chi2min_v(i) =  d.FitResult.chi2min;
        p_v(i) = 1-chi2cdf(d.FitResult.chi2min,d.FitResult.dof);
        TimeSec_v(i) = d.TimeSec;
    else
        fprintf('Fit %s\n',RunLists{i});
        R = MultiRunAnalysis('RunList',RunLists{i},RunAnaArg{:});
        R.exclDataStart = R.GetexclDataStart(range);
        R.Fit;
        FitResult = R.FitResult;
        nRuns = R.nRuns;
        TimeSec = sum(R.RunData.TimeSec.*R.RunData.qUfrac(R.exclDataStart:end)); % measurement time in analysis interval
        save(filename,'FitResult','nRuns','TimeSec','RunAnaArg','range');
        fprintf('save file %s \n',filename);
    end
end

%% print results on screen
for i=1:numel(RunLists)
    fprintf('m2 = %.2f +- %.2f eV^2 , chi2min = %.1f, p = %.2f, %.0f runs, %.0f hours, %s \n',...
        mNuSq_v(i),mNuSqErr_v(i),chi2min_v(i),p_v(i),nRuns_v(i),TimeSec_v(i)./3600,RunLists{i});
end

%%
figureHandle = figure('Units','normalized','Position',[0.1,0.1,0.6,0.5]);
s1 = subplot(1,4,1:3);
plot(mNuSq_v(1).*ones(10,1),linspace(0,numel(RunLists)+1,10),':','LineWidth',2,'Color',rgb('Silver'));
hold on;
errorbar(mNuSq_v,1:numel(RunLists),mNuSqErr_v,'horizontal',...
    '.','MarkerSize',20,'LineStyle','none','Color',rgb('DodgerBlue'),'CapSize',0,'LineWidth',2);
PrettyFigureFormat('FontSize',24);
xlabel(sprintf('{\\itm}_\\nu^{ 2} (eV^{ 2})'));
ylabel('Scan selection');
set(gca,'YMinorTick','off');
yticks(1:numel(RunLists));
yticklabels(RunLists_Labels)
ax1 = gca;
xlim([-8.7 2.3]);
ylim([0.5 numel(RunLists)+0.5])

s2 = subplot(1,4,4);
area(linspace(0,0.05,10),numel(RunLists)+1.*ones(10,1),'FaceColor',rgb('Red'));%(0.05.*ones(10,1),linspace(0,numel(RunLists)+1,10),'-','MarkerSize',20,'LineWidth',2,'Color',rgb('Red'));
hold on;
plot(p_v,1:numel(RunLists),'.','MarkerSize',20,'LineWidth',2,'Color',rgb('DodgerBlue'));
PrettyFigureFormat('FontSize',24);
yticklabels('');
xlabel(sprintf('{\\it p}'));
ax2 = gca;
xlim([0 1])
ylim([0.5 numel(RunLists)+0.5])
linkaxes([s1,s2],'y');


%
ax1.Position(2) = 0.17;
ax2.Position(2) = ax1.Position(2);
ax2.Position(4) = ax1.Position(4);
ax1.Position(1) = 0.2;
ax2.Position(1) = 0.79;
 

 %% save plot
 pltname = [strrep(savedir,'results','plots'),'knm1_AltRunLists.pdf'];
 export_fig(gcf,pltname);