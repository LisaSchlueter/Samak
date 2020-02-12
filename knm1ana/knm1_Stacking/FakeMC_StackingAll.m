% investigate effect of different qU, qUfrac, Time, rhod, isotopologues on stacking with fake MC
% test stacking qU correction

%% settings
CleanUp = 'OFF';
exclDataStart = 14;
InitFile = @ref_FakeRun_KNM1_StackingStudyAll;

%% Get qU uncertainty from real data
Real = MultiRunAnalysis('RunList','KNM1','DataType','Real');
qU_RelErr = std(Real.SingleRunData.qU,0,2)./mean(Real.SingleRunData.qU,2);
qUfrac_RelErr = std(Real.SingleRunData.qUfrac,0,2)./mean(Real.SingleRunData.qUfrac,2);
Time_RelErr   = std(Real.SingleRunData.TimeSec)./mean(Real.SingleRunData.TimeSec);
WGTS_CD_MolPerCm2_RelErr = std(Real.SingleRunData.WGTS_CD_MolPerCm2)./mean((Real.SingleRunData.WGTS_CD_MolPerCm2));
WGTS_MolFrac_TT_RelErr = std(Real.SingleRunData.WGTS_MolFrac_TT)./mean((Real.SingleRunData.WGTS_MolFrac_TT));
WGTS_MolFrac_HT_RelErr = std(Real.SingleRunData.WGTS_MolFrac_HT)./mean((Real.SingleRunData.WGTS_MolFrac_HT));
WGTS_MolFrac_DT_RelErr = std(Real.SingleRunData.WGTS_MolFrac_DT)./mean((Real.SingleRunData.WGTS_MolFrac_DT));
%% compute fake MC runs
nRuns = Real.nRuns;
savedir = [getenv('SamakPath'),'tritium-data/mat/',extractAfter(func2str(InitFile),'ref_'),'/'];
savename = [extractAfter(func2str(InitFile),'ref_'),'_',num2str(nRuns),'.mat'];

if ~exist([savedir,savename],'file') || strcmp(CleanUp,'ON')
    
    MC = McRunGenerator('McFlag','Fake','InitFile',InitFile,'nRuns',nRuns,'qU_RelErr',mean(qU_RelErr),'qUfrac_RelErr',mean(qUfrac_RelErr),...
        'TimeSec_RelErr',Time_RelErr,'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr,'WGTS_MolFrac_TT_RelErr',WGTS_MolFrac_TT_RelErr,...
        'WGTS_MolFrac_DT_RelErr',WGTS_MolFrac_DT_RelErr,'WGTS_MolFrac_HT_RelErr',WGTS_MolFrac_HT_RelErr,'VarDist','Uniform');
    
    MC.CleanUpFakeRuns;
    MC.ComputeFakeRun;
end
%%
M = MultiRunAnalysis('RunList',1:nRuns,'DataType','Fake','FakeStudyName',extractAfter(func2str(InitFile),'ref_'),'fixPar','5 6 7 8 9 10 11');
M.exclDataStart = 1;%exclDataStart;
M.Fit; M.PlotFit;
%M.RhoDScan;
% %%
%  Mcorr = MultiRunAnalysis('RunList',1:nRuns,'DataType','Fake','FakeStudyName',extractAfter(func2str(InitFile),'ref_'),'fixPar','5 6 7 8 9 10 11',...
%                          'StackqUCorrFlag','ON','StackTolerance',10);
%  Mcorr.Fit;
