% neutrino mass sensitivity as a function of bkg slope constraint
% method: gaussian randomization

% background slope (qU) systematics
% strategies: gaussian randomization

%% systematics setting
RecomputeFlag = 'OFF';
CovMatRecomputeFlag = 'ON';
MaxSlopeCpsPereV = (2:1:15).*1e-06;%5.2*1e-06;

savedir = [getenv('SamakPath'),'knm2ana/knm2_systematics/results/'];
MakeDir(savedir);

%% model setting
RunList   = 'KNM2_Prompt';
AnaFlag   = 'StackPixel'; % FPD segmentation

RunArg = {'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag',AnaFlag,...
    'chi2','chi2Stat',...      % switch on later by hand
    'RunList',RunList,...
    'fixPar','mNu E0 Norm Bkg',... % free parameter
    'DataType','Twin',...      % MC twins
    'TwinBias_Q',18573.56,...  % twin endpoint
    'SysBudget',34,...
    'RingMerge','Full',...
    'NonPoissonScaleFactor',1};
%% stat. only
savenameStat = sprintf('%smNuSqBkgSlopeSys_Stat_40eVrange.mat',savedir);
if exist(savenameStat,'file') && strcmp(RecomputeFlag,'OFF')
    dstat = importdata(savenameStat);
    FitResultsStat = dstat.FitResultsStat;
elseif ~exist(savenameStat,'file') %%
    M = MultiRunAnalysis(RunArg{:});
    M.exclDataStart = M.GetexclDataStart(40);
    M.Fit;
    FitResultsStat = M.FitResult;
    save(savenameStat,'FitResultsStat','RunArg');
end

mNuSysNeg  = zeros(numel(MaxSlopeCpsPereV),1);
mNuSysPos  = zeros(numel(MaxSlopeCpsPereV),1);
mNuSysMean = zeros(numel(MaxSlopeCpsPereV),1);
%% stat + syst
for i=1:numel(MaxSlopeCpsPereV)
    progressbar(i/numel(MaxSlopeCpsPereV))
    savename = sprintf('%sknm2_BKGsys_GaussOnly_%.1fmcpskeV.mat',savedir,MaxSlopeCpsPereV(i)*1e6);
    if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
        d = importdata(savename);
        mNuSysNeg(i) = sqrt(d.FitResultsCM.errNeg(1)^2-FitResultsStat.errNeg(1)^2);
        mNuSysPos(i) = sqrt(d.FitResultsCM.errPos(1)^2-FitResultsStat.errPos(1)^2);
        mNuSysMean(i) = 0.5*(mNuSysNeg(i)+mNuSysPos(i));
    else
        CmArg = {'BkgCM','ON',...
            'SysEffects',struct('BkgShape','ON'),...
            'nTrials',50000,...
            'MaxSlopeCpsPereV',MaxSlopeCpsPereV(i)};
        %% init model
        M = MultiRunAnalysis(RunArg{:});
        M.exclDataStart = M.GetexclDataStart(40);
        M.chi2 = 'chi2CMShape';
        %% stat. + syst: calculate covariance matrix: slope randn
        M.ComputeCM(CmArg{:},'BkgMode','Gauss');
        CMFracGauss      = M.FitCM_Obj.CovMatFrac;
        CMFracGaussShape = M.FitCM_Obj.CovMatFracShape;
        M.Fit;
        FitResultsCM = M.FitResult;
        
        save(savename,'FitResultsCM','CMFracGauss','CMFracGaussShape');
        
        mNuSysNeg(i) = sqrt(FitResultsCM.errNeg(1)^2-FitResultsStat.errNeg(1)^2);
        mNuSysPos(i) = sqrt(FitResultsCM.errPos(1)^2-FitResultsStat.errPos(1)^2);
        mNuSysMean(i) = 0.5*(mNuSysNeg(i)+mNuSysPos(i));
    end
end

 %% display results
 
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
x_plot = linspace(min(MaxSlopeCpsPereV),max(MaxSlopeCpsPereV),1e3);
plot(x_plot.*1e6,interp1(MaxSlopeCpsPereV,mNuSysMean,x_plot,'spline'),'-','LineWidth',2,'Color',rgb('DodgerBlue'));
hold on;
hFT = PlotLine(5.2,interp1(MaxSlopeCpsPereV,mNuSysMean,5.2.*1e-06,'spline'),'SlateGray','-');
hKNM12    = PlotLine(13.3,interp1(MaxSlopeCpsPereV,mNuSysMean,13.3.*1e-06,'spline'),'Crimson','-.');
hKNM12fix = PlotLine(9.8,interp1(MaxSlopeCpsPereV,mNuSysMean,9.8.*1e-06,'spline'),'LimeGreen','--');
hKNM2 = PlotLine(11.7,interp1(MaxSlopeCpsPereV,mNuSysMean,11.7.*1e-06,'spline'),'Orange',':');
leg = legend([hFT,hKNM12fix,hKNM2,hKNM12],'FT             =   5.2 mcps / keV',...
                                          'KNM1(+2) =   9.8 mcps / keV (baseline fixed)',...
                                          'KNM2       = 11.7 mcps / keV',...
                                          'KNM1(+2) = 13.3 mcps / keV',...
                                   'Location','southeast','EdgeColor',rgb('Silver'));
%plot(MaxSlopeCpsPereV.*1e6,mNuSysMean,'o','Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue'));
t = title('Background slope constraint','FontWeight','normal');
PrettyFigureFormat;
xlabel(sprintf('1\\sigma constraint (mcps / keV)'));
ylabel(sprintf('\\sigma({\\itm}_\\nu^2) (eV^2)'));
xlim([min(MaxSlopeCpsPereV) max(MaxSlopeCpsPereV)].*1e6)
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = sprintf('%sknm2_mNuSqBkgSlopeSys_40eV.png',plotdir);
print(f1,plotname,'-dpng','-r500');

%  mNuStatAsym = 0.5*(-FitResultsStat.errNeg(1)+FitResultsStat.errPos(1));
% fprintf('mnuSq sensitivity stat only        = %.3f eV^2 \n',mNuStatAsym);
% % fprintf('mnuSq sensitivity stat + syst only = %.3f eV^2  (gauss) \n',0.5*(-FitResultsCM.errNeg(1)+FitResultsCM.errPos(1)));
% fprintf(2,'mnuSq sensitivity syst only (asym) = %.3g eV^2  (gauss) \n',mNuSysMean);
% % %% plot cov mat
% %M.FitCM_Obj.PlotCM('qUWindowIndexMax',-20,'saveplot','ON');

function phandle = PlotLine(x,y,color,lstyle)
%phandle = plot(x.*ones(10,1),linspace(0,y,10),lstyle,'LineWidth',2,'Color',rgb(color));
phandle =plot(x,y,'o','Color',rgb(color),'MarkerFaceColor',rgb(color),'MarkerSize',8);
hold on;
plot(linspace(0,x,10),y.*ones(10,1),lstyle,'LineWidth',2,'Color',rgb(color));
end