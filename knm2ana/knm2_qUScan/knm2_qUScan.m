% Uniform fit on KNM2 data/twins
% January 2020, Lisa
ELossFlag = 'KatrinT2A20';%'Aseev'; 
chi2 = 'chi2Stat';
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
 SigmaSq =  0.0124+0.0025;
 
RunAnaArg = {'RunList','KNM2_Prompt',...  % define run number -> see GetRunList
    'fixPar','mNu E0 Bkg Norm',...         % free Parameter !!
    'DataType','Real',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag',ELossFlag,...         % energy loss function     ( different parametrizations available)
    'minuitOpt','min ; minos',...
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2',chi2,...              % statistics only
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'TwinBias_Q',18573.7,...
    'DopplerEffectFlag','FSD',...
    'FSD_Sigma',sqrt(SigmaSq),...
    'TwinBias_FSDSigma',sqrt(SigmaSq),...
    'SysBudget',38};
    
    %% build object of MultiRunAnalysis class
D = MultiRunAnalysis(RunAnaArg{:});
   
%% modify some parameters in your analysis
range = 40;               % fit range in eV below endpoint        
D.exclDataStart = D.GetexclDataStart(range); % find correct data, where to cut spectrum
%% extra label for qU-scan with "wrong" energy loss function
%if ~strcmp(ELossFlag,'KatrinT2')
    saveStr = ['_',ELossFlag];
%else
%    saveStr = '';
%end
%% 
% plotting
D.chi2 = 'chi2CMShape';
D.NonPoissonScaleFactor = 1.112;
[parqUCM, errqUCM, chi2qUCM, dofqUCM,eCM] = D.qUScan('qURange',[90 22],'RecomputeFlag','OFF',...
        'saveplot','ON','ErrorBarScaling',1,'saveStr',saveStr,'HoldOn','ON1','CorrMean','OFF');
    
D.chi2 = 'chi2Stat';
D.NonPoissonScaleFactor = 1;
[parqU, errqU, chi2qU, dofqU,e1] = D.qUScan('qURange',[90 22],'RecomputeFlag','OFF',...
        'saveplot','ON','ErrorBarScaling',1,'saveStr',saveStr,'HoldOn','ON');
   
leg = legend([e1,eCM],'Stat. only','Stat. and syst.',...
    'EdgeColor',rgb('Silver'),'Location','southwest','FontSize',get(gca,'FontSize'));
plotname = [getenv('SamakPath'),'knm2ana/knm2_qUScan/plots/knm2_qUScanStatSyst_mNuSq.png'];
print(plotname,'-dpng','-r350');