% KNM1 Alternative run lists
% save or load results
% plot results
chi2 = 'chi2Stat';
NP = 1.064;

savedir = [getenv('SamakPath'),'knm1ana/knm1_AlternativeRunLists/results/'];
fcommon = [savedir,sprintf('knm1_AltRunList_%s_NP%2g',chi2,NP)];

% settings
range   = 40;         % 40eV range = 27 subruns
% Init Model Object and covariance matrix object
RunAnaArg = { 'chi2',chi2,...
    'DataType','Real',...
    'fixPar','mNu E0 Norm Bkg',... free parameter
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',NP,...
    'minuitOpt','min ; minos',...
    'FSDFlag','Sibille0p5eV',...
    'ELossFlag','KatrinT2',...
    'SysBudget',22,...
    'AngularTFFlag','OFF'};

% define run lists

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

RunLists_Labels = {'Golden scans','Up scans','Down scans',...
                 'First third','Middle third','Last third',...
                 sprintf('\\sigma(\\rhod) < 1'), sprintf('\\sigma(\\rhod) > 1'),...
                 sprintf('Low \\epsilon_T'), sprintf('High \\epsilon_T'),...
                 sprintf('{\\itU}_{RW} = -183 mV'),...
                 sprintf('{\\itU}_{RW} = -149 mV'),...
                 };
%% fit or load run lists
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
        ftmp = importdata(filename);
        nRuns_v(i) = ftmp.nRuns;
        mNuSq_v(i) = ftmp.FitResult.par(1);
       % mNuSqErr_v(i) =  ftmp.FitResult.err(1);
        mNuSqErr_v(i) =  0.5.*(ftmp.FitResult.errPos(1)-ftmp.FitResult.errNeg(1));
        chi2min_v(i) =  ftmp.FitResult.chi2min;
        p_v(i) = 1-chi2cdf(ftmp.FitResult.chi2min,ftmp.FitResult.dof);
        TimeSec_v(i) = ftmp.TimeSec;
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

%% plot
RandList = 'ON';
LocalFontSize = 18;
figureHandle = figure('Units','normalized','Position',[0.1,0.1,0.6,0.5]);
s1 = subplot(1,4,1:3);

if strcmp(RandList,'ON')
    nRunList = 2000;
    chi2 = 'chi2Stat';
    NP = 1.064;
    savenameRand = [savedir,sprintf('RandomHalfRunList_Unblinded_%s_NP%2g_%.0ffits_merged.mat',chi2,NP,nRunList)];
    frand = importdata(savenameRand);
    mNuSq_rand =  cell2mat(cellfun(@(x) x.par(1),frand.FitResults,'UniformOutput',0));
    [l,a_rand] = boundedline(mean(mNuSq_rand).*ones(10,1),linspace(-1,numel(RunLists)+1,10),std(mNuSq_rand).*ones(10,1),'orientation','horiz');
    hold on
    l.LineStyle = 'none'; l.Color = rgb('DimGray'); l.LineWidth = 2;
    a_rand.FaceColor = rgb('Silver');
    a_rand.FaceAlpha = 0.5;   
end

%
pref = plot(mNuSq_v(1).*ones(10,1),linspace(-1,numel(RunLists)+2,10),':','LineWidth',2.5,'Color',rgb('LimeGreen'));
hold on;
e_fits = errorbar(mNuSq_v,1:numel(RunLists),mNuSqErr_v,'horizontal',...
    '.','MarkerSize',20,'LineStyle','none','Color',rgb('DodgerBlue'),'CapSize',0,'LineWidth',2);
e_bf = errorbar(mNuSq_v(1),1,mNuSqErr_v(1),'horizontal',...
    '.','MarkerSize',20,'LineStyle','none','Color',pref.Color,'CapSize',0,'LineWidth',2);
xlabel(sprintf('{\\itm}_\\nu^{ 2} (eV^{ 2})'));
ylabel('Scan selection');

PrettyFigureFormat('FontSize',LocalFontSize);
set(gca,'YMinorTick','off');
yticks(1:numel(RunLists));
yticklabels(RunLists_Labels)
ax1 = gca;
xlim([-8.8 2.3]);

if strcmp(RandList,'ON')
    leg = legend([a_rand],...
       sprintf('Random half scans: 1\\sigma band (%.0f samples)',nRunList));
    PrettyLegendFormat(leg,'alpha',0.8); leg.EdgeColor = 'none';
    leg.Location= 'northwest';
    leg.FontSize = LocalFontSize-2;
end

s2 = subplot(1,4,4);
area(linspace(0,0.05,10),numel(RunLists)+2.*ones(10,1),'FaceColor',rgb('Red'));%(0.05.*ones(10,1),linspace(0,numel(RunLists)+1,10),'-','MarkerSize',20,'LineWidth',2,'Color',rgb('Red'));
hold on;
plot(p_v,1:numel(RunLists),'.','MarkerSize',20,'LineWidth',2,'Color',rgb('DodgerBlue'));
plot(p_v(1),1,'.','MarkerSize',20,'LineWidth',2,'Color',pref.Color);

yticklabels('');
xlabel(sprintf('{\\it p}'));
ax2 = gca;
xlim([0 1])
ylim([0.5 numel(RunLists)+0.5])
PrettyFigureFormat('FontSize',LocalFontSize);
linkaxes([s1,s2],'y');

ax1.Position(2) = 0.17;
ax2.Position(2) = ax1.Position(2);
ax2.Position(4) = ax1.Position(4);
ax1.Position(1) = 0.2;
ax2.Position(1) = 0.79;
 
if strcmp(RandList,'ON')
    ylim([0 numel(RunLists)+2])
else
    ylim([0.5 numel(RunLists)+0.5])
end

 % save plot
 pltname = [strrep(savedir,'results','plots'),'knm1_AltRunLists.pdf'];
 export_fig(gcf,pltname);
 
 %% how many results lie within the 1sigma band from the random half?
mNuSq_low =  mean(mNuSq_rand)-std(mNuSq_rand);
mNuSq_up =  mean(mNuSq_rand)+std(mNuSq_rand);

LogIdx = mNuSq_v>=mNuSq_low & mNuSq_v<=mNuSq_up;
Frac = sum(LogIdx)./numel(LogIdx);
fprintf('%.0f out of %.0f (%.1f%%) fit results are within 1sigma\n',sum(LogIdx),numel(LogIdx),1e2*Frac);
 
%% calculate sigmas
Sigma = (mNuSq_v-mean(mNuSq_rand))./std(mNuSq_rand);


