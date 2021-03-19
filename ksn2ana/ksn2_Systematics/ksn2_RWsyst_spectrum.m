% investigate possible impact of new syst. effect: Tritium on rear wall
% tritium spectrum modifications
% 1) different response function (starting pos.)
% 2) Endpoint shifted
% 3) Signal normalization
% 4) No additional background
% 5) maybe different FSD

E0ShifteV = 1.5; % eV
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
savename = [savedir,sprintf('ksn2_RWsyst_spectrum_E0shift%.3geV.mat',E0ShifteV)];

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
     FitResults= A.InitModelObj_Norm_BKG;
     A.ModelObj.BKG_RateSec_i = A.ModelObj.BKG_RateSec_i+FitResults.par(3);
     A.ModelObj.NormFactorTBDDS = A.ModelObj.NormFactorTBDDS.*(1+FitResults.par(4));
     A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
     
     TBDIS_wgts = A.ModelObj.TBDIS;
     Q_i = A.ModelObj.Q_i;
     %% init rear wall spectrum
     A.ModelObj.KTFFlag = 'RW_WGTSMACE';
     A.ModelObj.InitializeRF;
     A.ModelObj.Q_i = Q_i+E0ShifteV;
     A.ModelObj.BKG_RateSec_i = 0; % no (additional) background
     A.ModelObj.BKG_PtSlope_i = 0;
     A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
   
     TBDIS_rw = A.ModelObj.TBDIS;
     
     % save
     MakeDir(savedir);
     qU = A.RunData.qU;
     TBDIS_Data = A.RunData.TBDIS;
     Time = A.RunData.qUfrac.*A.RunData.TimeSec;
     save(savename,'qU','TBDIS_Data','TBDIS_rw','TBDIS_wgts','E0ShifteV','Q_i','Time');
end
 %% display
Scale = 1e3;

 GetFigure;
 p_sum = plot(qU-18574,(TBDIS_rw.*Scale/1e3+TBDIS_wgts)./Time,'-','LineWidth',2.5,'Color',rgb('DodgerBlue'));
 hold on;
 p_wgts = plot(qU-18574,TBDIS_wgts./Time,':','LineWidth',2.5,'Color',rgb('Orange'));
 p_rw = plot(qU-18574,(TBDIS_rw.*Scale/1e3)./Time,'-.','LineWidth',2.5,'Color',rgb('FireBrick'));
 
 PrettyFigureFormat('FontSize',22);
 xlabel('Retarding energy - 18574 (eV)');
ylabel('Rate');
xlim([-41 137]);
if Scale==1
    ylim([1e-04,1e2])
    ScaleStr = '';
else
ylim([0.1,1e2])
ScaleStr = sprintf(' \\times %.3g',Scale);
end
leg = legend([p_sum,p_wgts,p_rw],sprintf('RF_{WGTS} + RF_{RW}%s',ScaleStr),'RF_{WGTS}',...
    sprintf('RF_{RW}%s',ScaleStr),'Location','northeast');
 PrettyLegendFormat(leg);
 set(gca,'YScale','log')
%% save
plotdir = strrep(savedir,'results','plots');
plotname = [plotdir,sprintf('ksn2_RWsyst_spectrum_Scale%.3g_E0Shift%.3geV.png',Scale,E0ShifteV)];
print(gcf,plotname,'-dpng','-r300');
fprintf('save plot to %s \n',plotname);