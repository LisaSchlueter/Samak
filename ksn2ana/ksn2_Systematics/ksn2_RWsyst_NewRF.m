% investigate possible impact of new syst. effect: Tritium on rear wall
% response function needs modification: no integration over starting position z
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
WGTS_B_T = 1.23;
savename = [savedir,sprintf('ksn2_RWsyst_NewRW_Bs%.3gT.mat',WGTS_B_T)];

if exist(savename,'file')
    load(savename)
    fprintf('Load RF from file %s \n',savename)
else
    %% configure RunAnalysis object
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
    
    RF_WGTS = A.ModelObj.RF;
    %%
    A.ModelObj.KTFFlag = 'RW_WGTSMACE';
    A.ModelObj.WGTS_B_T = WGTS_B_T;
    A.ModelObj.InitializeRF;
    RF_RW = A.ModelObj.RF;
    
    E = A.ModelObj.Te;
    qU = A.ModelObj.qU;
    
    WGTS_CD_MolPerCm2 = A.ModelObj.WGTS_CD_MolPerCm2;
    MACE_Bmax_T = A.ModelObj.MACE_Bmax_T;
    MACE_Ba_T = A.ModelObj.MACE_Ba_T;
    WGTS_B_T = A.ModelObj.WGTS_B_T;
    MakeDir(savedir);
    save(savename,'E','qU','RF_RW','RF_WGTS','WGTS_CD_MolPerCm2','MACE_Bmax_T','MACE_Ba_T','WGTS_B_T');
end
%% display
nqU = 20;
GetFigure;
p_wgts = plot(E-qU(nqU),RF_WGTS(:,nqU),'-','LineWidth',2);
hold on;
p_rw = plot(E-qU(nqU),RF_RW(:,nqU),'-.','LineWidth',2);
PrettyFigureFormat('FontSize',22);
xlabel('Surplus energy (eV)');
ylabel('Transmission');
xlim([-10 100]);
leg = legend([p_wgts,p_rw],'RF_{WGTS}','RF_{RW}','Location','northwest');
PrettyLegendFormat(leg);
% save
plotdir = strrep(savedir,'results','plots');
plotname = [plotdir,'ksn2_RWsyst_NewRF.png'];
print(gcf,plotname,'-dpng','-r300');
fprintf('save plot to %s \n',plotname);