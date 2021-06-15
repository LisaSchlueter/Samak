% investigate possible impact of new syst. effect: Tritium on rear wall
% tritium spectrum modifications
% 1) different response function (starting pos.)
% 2) Endpoint shifted
% 3) Signal normalization
% 4) No additional background
% 5) maybe different FSD

ScaleRW = 0.012/1.181;%0.58*1e-03;%5e-04;
E0ShifteV = +1.64; % eV
B_RW_Source = 1.23;%2.52;%

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
savename = [savedir,sprintf('ksn2_RWsyst_spectrum_E0shift%.3geV_Bsource%.3gT_ScaleRW%.3g.mat',E0ShifteV,B_RW_Source,ScaleRW)];

if exist(savename,'file')
    load(savename)
    fprintf('Load spectra from file %s \n',savename)
else

    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType','Twin',...
        'fixPar','E0 Norm Bkg',...%free par
        'SysBudget',40,...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'RadiativeFlag','ON',...
        'FSDFlag','KNM2_0p1eV',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'NonPoissonScaleFactor',1,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'TwinBias_Q',18573.7,...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',3*1e-06,...
        'TwinBias_BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
    A = MultiRunAnalysis(RunAnaArg{:});
       % init model and save init settings
    FitResults                 = A.InitModelObj_Norm_BKG;
    A.ModelObj.BKG_RateSec_i   = A.ModelObj.BKG_RateSec_i+FitResults.par(3);
    A.ModelObj.NormFactorTBDDS = A.ModelObj.NormFactorTBDDS.*(1+FitResults.par(4));
    NormFactorTBDDS_i          = A.ModelObj.NormFactorTBDDS;
    Bkg_i                      = A.ModelObj.BKG_RateSec_i;
    Q_i                        = A.ModelObj.Q_i;
    WGTS_CD_MolPerCm2_i        = A.ModelObj.WGTS_CD_MolPerCm2;
    
    % compute TBDIS for regular tritium spectrum (wgts)
    A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
    TBDIS_wgts             = A.ModelObj.TBDIS;
    
    % for refernence: no (additional) background
     A.ModelObj.BKG_RateSec_i = 0; 
     A.ModelObj.BKG_PtSlope_i = 0;
     A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
     TBDIS_wgts_noBkg =  A.ModelObj.TBDIS;
     
    %% init rear wall spectrum
     A.ModelObj.KTFFlag = 'RW_WGTSMACE';
     A.ModelObj.WGTS_B_T = B_RW_Source;
     A.ModelObj.ComputeNormFactorTBDDS;
     A.ModelObj.InitializeRF;
     
     NormFactorTBDDS_rw       = A.ModelObj.NormFactorTBDDS; % different norm factor, becauseof B-source
     A.ModelObj.Q_i           = Q_i + E0ShifteV;
     A.ModelObj.BKG_RateSec_i = 0; % no (additional) background
     A.ModelObj.BKG_PtSlope_i = 0;
     A.ModelObj.NormFactorTBDDS = NormFactorTBDDS_rw*ScaleRW; % scale to lower activity
     
     A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
     TBDIS_rw = A.ModelObj.TBDIS;
     
     % save
     MakeDir(savedir);
     qU = A.RunData.qU;
     TBDIS_Data = A.RunData.TBDIS;
     Time = A.RunData.qUfrac.*A.RunData.TimeSec;
     
      TBDIS_sum =   TBDIS_wgts+TBDIS_rw;
       
     save(savename,'qU','TBDIS_wgts_noBkg','TBDIS_sum','TBDIS_rw','TBDIS_wgts','E0ShifteV','Q_i','Time','NormFactorTBDDS_rw','NormFactorTBDDS_i');
end
 %% display
 PlotRwScale = 1e2;
 GetFigure;
 p_sum = plot(qU-18574,(TBDIS_wgts+PlotRwScale.*TBDIS_rw)./Time,'-','LineWidth',2.5,'Color',rgb('DodgerBlue'));
 hold on;
 p_wgts = plot(qU-18574,TBDIS_wgts./Time,':','LineWidth',2.5,'Color',rgb('Orange'));
 p_rw = plot(qU-18574,PlotRwScale.*TBDIS_rw./Time,'-.','LineWidth',2.5,'Color',rgb('FireBrick'));
 
 PrettyFigureFormat('FontSize',22); 
 xlabel('Retarding energy - 18574 (eV)');
ylabel('Rate');
xlim([-41 10]);
ylim([0.1,1e2])

leg = legend([p_sum,p_wgts,p_rw],sprintf('{\\itR}_{WGTS} + {\\itR}_{RW} \\times %.3g',PlotRwScale),...
    sprintf('{\\itR}_{WGTS}'),...
    sprintf('{\\itR}_{RW} \\times %.3g',PlotRwScale),'Location','northeast');
 PrettyLegendFormat(leg);
 leg.Title.String = 'Integral tritium spectrum';
 leg.Title.FontWeight = 'normal';
 set(gca,'YScale','log')
% save
plotdir = strrep(savedir,'results','plots');
plotname = [plotdir,sprintf('ksn2_RWsyst_spectrum_E0shift%.3geV_Bsource%.3gT_ScaleRW%.3g_PlotScaleRW%.3g.png',...
    E0ShifteV,B_RW_Source,ScaleRW,PlotRwScale)];
print(gcf,plotname,'-dpng','-r300');
fprintf('save plot to %s \n',plotname);

%% number of counts
fprintf('============================================================= \n')
fprintf('Sum Counts RW         = %.3g  (40 eV range - %.0f subruns) \n',sum(TBDIS_rw(11:end)),numel(qU(11:end)));
fprintf('Sum Counts WGTS - BKG = %.3g  (40 eV range - %.0f subruns) \n',sum(TBDIS_wgts_noBkg(11:end)),numel(qU(11:end)));
fprintf('Sum Counts WGTS + BKG = %.3g  (40 eV range - %.0f subruns) \n',sum(TBDIS_wgts(11:end)),numel(qU(11:end)));
fprintf('============================================================= \n')