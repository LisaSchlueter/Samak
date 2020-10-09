%% sensitivity on background slope from KNM2 data
%% settings
range = 40;
NonPoissonScaleFactor = 1.112;
MC_BkgSlope = [0,5,10,15,20].*1e-06;
SigmaBkgSlope = 14.3*1e-06;%4.74.*1e-06;%
SavePlot = 'ON';

%% init
mNuSqStat        = zeros(numel(MC_BkgSlope),1);
mNuSqBkgFree     = zeros(numel(MC_BkgSlope),1);
mNuSqBkgFreePull = zeros(numel(MC_BkgSlope),1);
mNuSqBkgFixCM    = zeros(numel(MC_BkgSlope),1);

mNuSqErrStat        = zeros(numel(MC_BkgSlope),1);
mNuSqErrBkgFree     = zeros(numel(MC_BkgSlope),1);
mNuSqErrBkgFreePull = zeros(numel(MC_BkgSlope),1);
mNuSqErrBkgFixCM    = zeros(numel(MC_BkgSlope),1);
%% load or calculate
savedir = [getenv('SamakPath'),'knm2ana/knm2_background/results/'];
MakeDir(savedir);

for i=1:numel(MC_BkgSlope)
    progressbar(i/numel(MC_BkgSlope));
    savename = sprintf('%sknm2_BkgSlopemNuSqBias_Bslope%.0fmcpsKeV_NPfactor%.2f_%.0feV_BkgSlopeSigma%.1f.mat',savedir,MC_BkgSlope(i)*1e6,NonPoissonScaleFactor,range,1e6*SigmaBkgSlope);
    if exist(savename,'file')
        d = importdata(savename);
        mNuSqStat(i)        = d.FitResults_mNuSqBkgFixStat.par(1);
        mNuSqBkgFree(i)     = d.FitResults_mNuSqBkgFree.par(1);
        mNuSqBkgFreePull(i) = d.FitResults_mNuSqBkgFreePull.par(1);
        mNuSqBkgFixCM(i)    = d.FitResults_BkgFixCM.par(1);
        
        mNuSqErrStat(i)        = 0.5*(-d.FitResults_mNuSqBkgFixStat.errNeg(1)+d.FitResults_mNuSqBkgFixStat.errPos(1));
        mNuSqErrBkgFree(i)     = 0.5*(-d.FitResults_mNuSqBkgFree.errNeg(1)+d.FitResults_mNuSqBkgFree.errPos(1));
        mNuSqErrBkgFreePull(i) = 0.5*(-d.FitResults_mNuSqBkgFreePull.errNeg(1)+d.FitResults_mNuSqBkgFreePull.errPos(1));
        mNuSqErrBkgFixCM(i)    = 0.5*(-d.FitResults_BkgFixCM.errNeg(1)+d.FitResults_BkgFixCM.errPos(1));
    else
        %% set up model
        RunAnaArg = {'RunList','KNM2_Prompt',...
            'fixPar','mNu E0 Bkg Norm',...
            'DataType','Twin',...
            'FSDFlag','BlindingKNM2',...
            'ELossFlag','KatrinT2A20',...
            'AnaFlag','StackPixel',...
            'chi2','chi2Stat',...
            'ROIFlag','Default',...
            'SynchrotronFlag','ON',...
            'AngularTFFlag','ON',...
            'ISCSFlag','Edep',...
            'TwinBias_Q',18573.56,...
            'NonPoissonScaleFactor',NonPoissonScaleFactor};
        T = MultiRunAnalysis(RunAnaArg{:});
        T.exclDataStart = T.GetexclDataStart(range);
        T.InitModelObj_Norm_BKG;
        T.ModelObj.ComputeTBDDS('BSlope_bias',MC_BkgSlope(i));
        T.ModelObj.ComputeTBDIS;
        TBDIS_i = T.ModelObj.TBDIS;
        
        T.RunData.TBDIS = TBDIS_i;
        
        %% reference: stat only + fixed bkg slope
        T.Fit;
        FitResults_mNuSqBkgFixStat = T.FitResult;
        %% fit with free nu-mass + free bkg slope + pull term
        T.pullFlag = 10;
        T.pulls = SigmaBkgSlope;
        T.fixPar = 'mNu E0 Bkg Norm BkgSlope';
        T.InitFitPar;
        T.Fit;
        FitResults_mNuSqBkgFreePull = T.FitResult;
        
        %% fit with free nu-mass + free bkg slope
        T.pullFlag = 99;
        T.fixPar = 'mNu E0 Bkg Norm BkgSlope';
        T.InitFitPar;
        T.Fit;
        FitResults_mNuSqBkgFree = T.FitResult;
        
        %% fit free nu-mass + fixed background slope + cov mat
        T.ModelObj.ComputeTBDDS('BSlope_bias',0);
        T.ModelObj.ComputeTBDIS;
        T.fixPar = 'mNu E0 Bkg Norm';
        T.InitFitPar;
        T.pullFlag = 99;
        
        % get cov. mat.
        T.chi2 = 'chi2CMShape';
        CmArg = {'BkgCM','ON',...
            'SysEffects',struct('BkgShape','ON'),...
            'nTrials',10000,...
            'MaxSlopeCpsPereV',SigmaBkgSlope};
        T.ComputeCM(CmArg{:},'BkgMode','Gauss');
        
        
        T.Fit;
        FitResults_BkgFixCM = T.FitResult;
        
        save(savename,'FitResults_BkgFixCM','FitResults_mNuSqBkgFree',...
            'FitResults_mNuSqBkgFixStat','FitResults_mNuSqBkgFreePull',...
            'range','RunAnaArg','range','NonPoissonScaleFactor')
        
        mNuSqStat(i)        = FitResults_mNuSqBkgFixStat.par(1);
        mNuSqBkgFree(i)     = FitResults_mNuSqBkgFree.par(1);
        mNuSqBkgFreePull(i) = FitResults_mNuSqBkgFreePull.par(1);
        mNuSqBkgFixCM(i)    = FitResults_BkgFixCM.par(1);
        
        mNuSqErrStat(i)        = 0.5*(-FitResults_mNuSqBkgFixStat.errNeg(1)+FitResults_mNuSqBkgFixStat.errPos(1));
        mNuSqErrBkgFree(i)     = 0.5*(-FitResults_mNuSqBkgFree.errNeg(1)+FitResults_mNuSqBkgFree.errPos(1));
        mNuSqErrBkgFreePull(i) = 0.5*(-FitResults_mNuSqBkgFreePull.errNeg(1)+FitResults_mNuSqBkgFreePull.errPos(1));
        mNuSqErrBkgFixCM(i)    = 0.5*(-FitResults_BkgFixCM.errNeg(1)+FitResults_BkgFixCM.errPos(1));
    end
end


%% plot results: neutrino mass bias
GetFigure;
CommonPlotArg = {'LineWidth',2,'MarkerSize',20};
pref = plot(linspace(-1,21,10),zeros(10,1),'-',CommonPlotArg{:},'Color',rgb('Silver'));
hold on;
pBfixStat = plot(1e6.*MC_BkgSlope,mNuSqStat,'.-.','Color',rgb('FireBrick'),CommonPlotArg{:});
hold on;
pBfree     = plot(1e6.*MC_BkgSlope,mNuSqBkgFree,'.--','Color',rgb('DodgerBlue'),CommonPlotArg{:});
pBfreePull = plot(1e6.*MC_BkgSlope,mNuSqBkgFreePull,'.-','Color',rgb('Orange'),CommonPlotArg{:});
pBfixCM    = plot(1e6.*MC_BkgSlope,mNuSqBkgFixCM,'.:','Color',rgb('ForestGreen'),CommonPlotArg{:});

xlabel(sprintf('MC data {\\itB} slope (mcps / keV)'));
ylabel(sprintf('\\Delta{\\itm}_\\nu^2 (eV^2)'));
PrettyFigureFormat('FontSize',22);
leg = legend([pBfixStat,pBfixCM,pBfree,pBfreePull],...
    sprintf('{\\itB} slope fixed to 0 (stat. only)'),...
    sprintf('{\\itB} slope fixed to 0 (stat. + cov. mat. \\sigma = %.1f mcps / keV %',SigmaBkgSlope*1e6),...
    sprintf('{\\itB} slope free  (stat. only)'),...
    sprintf('{\\itB} slope free  (stat. + pull term \\sigma = %.1f mcps / keV %',SigmaBkgSlope*1e6),...
    'EdgeColor',rgb('Silver'),'Location','northwest');
xlim([-1 21])
ylim([-0.01 0.2]);
if strcmp(SavePlot,'ON')
    plotdir = strrep(savedir,'results','plots');
    MakeDir(plotdir)
    plotname = sprintf('%sknm2_BkgSlopemNuSqBias.pdf',plotdir);
    export_fig(plotname);
    plotnamePNG = sprintf('%sknm2_BkgSlopemNuSqBias_BkgSlopeSigma%.1f.png',plotdir,BkgPullSigma*1e6);
    print(plotnamePNG,'-dpng','-r500')
    fprintf('save plot to %s \n',plotname)
end

%% plot syst. uncertainty on mNuSq
% comparison pull term with cov. mat.
GetFigure;
CommonPlotArg = {'LineWidth',2,'MarkerSize',20};
mNuSqSysCM   = sqrt(mNuSqErrBkgFixCM.^2-mNuSqErrStat.^2);
mNuSqSysPull = sqrt(mNuSqErrBkgFreePull.^2-mNuSqErrStat.^2);

pPull = plot(1e6.*MC_BkgSlope,mNuSqSysPull ,'.-','Color',rgb('Orange'),CommonPlotArg{:});
hold on;
pCM = plot(1e6.*MC_BkgSlope,mNuSqSysCM,'.:','Color',rgb('ForestGreen'),CommonPlotArg{:});


hold off
xlabel(sprintf('MC data {\\itB} slope (mcps / keV)'));
ylabel(sprintf('Syst. only. \\sigma({\\itm}_\\nu^2) (eV^2)'));
PrettyFigureFormat('FontSize',22);
leg = legend([pCM,pPull],...
    sprintf('{\\itB} slope fixed to 0 (stat. + cov. mat. \\sigma = %.1f mcps / keV %',SigmaBkgSlope*1e6),...
    sprintf('{\\itB} slope free  (stat. + pull term \\sigma = %.1f mcps / keV %',SigmaBkgSlope*1e6),...
    'EdgeColor',rgb('Silver'),'Location','northwest');
xlim([-1 21]);

if strcmp(SavePlot,'ON')
    plotdir = strrep(savedir,'results','plots');
    MakeDir(plotdir)
    plotname2 = sprintf('%sknm2_BkgSlopemNuSqBias_Sensitivity.pdf',plotdir);
    export_fig(plotname2);
    plotname2PNG = sprintf('%sknm2_BkgSlopemNuSqBias_Sensitivity_BkgSlopeSigma%.1f.png',plotdir,BkgPullSigma*1e6);
    print(plotname2PNG,'-dpng','-r500')
    fprintf('save plot to %s \n',plotname)
end