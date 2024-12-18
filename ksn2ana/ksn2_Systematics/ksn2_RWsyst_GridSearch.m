% ksn2 grid earch on MC with TBDIS_wgts+TBDIS_rw
%% settings that might change
Mode                  = 'Display';
chi2                  = 'chi2Stat';
DataType              = 'Twin';
nGridSteps            = 25;
range                 = 40;
NonPoissonScaleFactor = 1.112;
RandMC                = 2;
%% configure RunAnalysis object
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar','E0 Norm Bkg',...%free par
    'SysBudget',40,...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    'RadiativeFlag','ON',...
    'FSDFlag','KNM2_0p5eV',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2',chi2,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'FSD_Sigma',sqrt(0.0124+0.0025),...
    'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
    'TwinBias_Q',18573.7,...
    'PullFlag',99,...;%99 = no pull
    'BKG_PtSlope',3*1e-06,...
    'TwinBias_BKG_PtSlope',3*1e-06,...
    'DopplerEffectFlag','FSD'};
A = MultiRunAnalysis(RunAnaArg{:});
A.exclDataStart = A.GetexclDataStart(range);
%% calculate twin spectrum with TBDIS_wgts + TBDIS_rw
ScaleRW = 0.012/1.181;%0.58*1e-03;%5e-04;
E0ShifteV = +1.64; % eV
B_RW_Source = 1.23;%2.52;%
freePar = 'mNu E0 Norm Bkg';
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
savename = [savedir,sprintf('ksn2_RWsyst_Fit%s_%.0feV_RWE0shift%.3geV_ScaleRWRate%.3g.mat',strrep(freePar,' ',''),range,E0ShifteV,ScaleRW)];
if B_RW_Source~=2.52
    savename =  strrep(savename,'.mat',sprintf('_%.3gWGTS_B_T.mat',B_RW_Source));
end

if exist(savename,'file')
    d = importdata(savename);
    fprintf('load from file %s \n',savename)
    A.RunData.TBDIS = d.TBDIS_sum;
else
    fprintf(2,'file %s not found \n Calculate spectrum with ksn2_RWsyst_mNuShift \n',savename)
end
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'range',range,...
    'RandMC',RandMC,...
    'RandMC_TBDIS',d.TBDIS_sum};
S =  SterileAnalysis(SterileArg{:});
switch Mode
    case 'Compute'
        S.GridSearch;
    case 'Display'
        S.RandMC = '';
        S.RandMC_TBDIS = '';
        S.RunAnaObj.FSDFlag = 'KNM2_0p1eV';
        S.nGridSteps = 50;
        S.LoadGridFile;
        S.Interp1Grid('Maxm4Sq',38.2^2);
        preg = S.ContourPlot('HoldOn','OFF','Color',rgb('DodgerBlue'));
        
        S.RandMC = RandMC;
        S.RandMC_TBDIS = d.TBDIS_sum;
        S.RunAnaObj.FSDFlag = 'KNM2_0p5eV';
        S.nGridSteps = 25;
        S.LoadGridFile;
        S.Interp1Grid('Maxm4Sq',36^2);
        prw = S.ContourPlot('HoldOn','ON','Color',rgb('Orange'),'LineStyle','-.');
        
        leg = legend([preg,prw],sprintf('{\\itR}_{WGTS}'),sprintf('{\\itR}_{WGTS} + {\\itR}_{RW}'));
        PrettyLegendFormat(leg);
        leg.Title.String = 'Simulated data';
        leg.Title.FontWeight = 'normal';
        
        ylim([1 max(ylim)])
        xlim([5e-03 0.5])
        
        % save
        plotdir = strrep(savedir,'results','plots');
        plotname = [plotdir,sprintf('ksn2_RWsyst_GridSearch_E0shift%.3geV_Bsource%.3gT_ScaleRW%.3g.png',...
            E0ShifteV,B_RW_Source,ScaleRW)];
        print(gcf,plotname,'-dpng','-r300');
        fprintf('save plot to %s \n',plotname);
end
