% old
R = MultiRunAnalysis('RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
    'chi2','chi2Stat',...                 % uncertainties: statistical or stat + systematic uncertainties
    'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
    'fixPar','E0 Norm Bkg',...            % free Parameter!!
    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
    'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
    'FSDFlag','Sibille0p5eV',...               % final state distribution
    'ELossFlag','KatrinT2',...            % energy loss function
    'SysBudget',22,...
    'DopplerEffectFlag','FSD_Knm1',...
    'TwinBias_Q',18573.73); 
%%
savedir = [getenv('SamakPath'),'tritium-data/mat/TwinKNM1/'];
for i=1:numel(R.RunList)
savename = sprintf('%sTwin%.0f_E018573.73eV.mat',savedir,R.RunList(i));
StartTimeStamp = R.SingleRunData.StartTimeStamp;
qU = d.qU(2:end,:);
qUfrac = d.qUfrac(2:end,:);
WGTS_MolFrac_TT_SubRun = WGTS_MolFrac_TT_SubRun(2:end);
WGTS_MolFrac_DT_SubRun = WGTS_MolFrac_DT_SubRun(2:end);
WGTS_MolFrac_HT_SubRun = WGTS_MolFrac_HT_SubRun(2:end);
end
%SR = RunAnalysis('RunNr',51410,'DataType','Twin','NonPoissonScaleFactor',1.064);