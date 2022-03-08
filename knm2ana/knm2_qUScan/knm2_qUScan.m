% Uniform fit on KNM2 data/twins
% March 2022, Lisa

chi2 = 'chi2CMShape';
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
SigmaSq =  0.0124+0.0025;
BKG_PtSlope = 3*1e-06;
ELossFlag = 'KatrinT2A20';
FSDFlag = 'KNM2_0p1eV';
SysBudget = 40;

RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType','Real',...
    'fixPar','mNu E0 Bkg Norm',...
    'RadiativeFlag','ON',...
    'DopplerEffectFlag','FSD',...
    'minuitOpt','min ; minos',...
    'FSDFlag',FSDFlag,...
    'ELossFlag',ELossFlag,...
    'SysBudget',SysBudget,...
    'AnaFlag','StackPixel',...
    'chi2',chi2,...
    'TwinBias_Q',18573.7,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'FSD_Sigma',sqrt(SigmaSq),...
    'TwinBias_FSDSigma',sqrt(SigmaSq),...
    'BKG_PtSlope',BKG_PtSlope,...
    'TwinBias_BKG_PtSlope',BKG_PtSlope};
D = MultiRunAnalysis(RunAnaArg{:});
range = 40;               % fit range in eV below endpoint        
D.exclDataStart = D.GetexclDataStart(range); % find correct data, where to cut spectrum
%% extra label for qU-scan with "wrong" energy loss function
saveStr = ['_',ELossFlag];
%% 
% plotting
D.chi2 = 'chi2CMShape';
D.NonPoissonScaleFactor = 1.112;
D.ComputeCM;
[parqUCM, errqUCM, chi2qUCM, dofqUCM,eCM] = D.qUScan('qURange',[90 22],'RecomputeFlag','OFF',...
        'saveplot','ON','ErrorBarScaling',1,'saveStr',saveStr,'HoldOn','ON1','CorrMean','OFF');
    
D.chi2 = 'chi2Stat';
D.NonPoissonScaleFactor = 1;
D.ComputeCM_StatPNP;
[parqU, errqU, chi2qU, dofqU,e1] = D.qUScan('qURange',[90 22],'RecomputeFlag','OFF',...
        'saveplot','ON','ErrorBarScaling',1,'saveStr',saveStr,'HoldOn','ON');
   
leg = legend([e1,eCM],'Stat. only','Stat. and syst.',...
    'EdgeColor',rgb('Silver'),'Location','southwest','FontSize',get(gca,'FontSize'));
plotname = [getenv('SamakPath'),...
    sprintf('knm2ana/knm2_qUScan/plots/knm2_qUScanStatSyst_mNuSq_%s_%s.png',...
    FSDFlag,ELossFlag)];
print(plotname,'-dpng','-r350');