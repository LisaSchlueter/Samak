% save my contour to .dat file for other fitters to read

% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';
nGridSteps = 30;
range = 40;
freePar = 'mNu E0 Norm Bkg';
NH = 'OFF';
%% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar',freePar,...%free par
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
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range,...
    'LoadGridArg',{'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'}};
S = SterileAnalysis(SterileArg{:});
%%
S.InterpMode = 'spline';
if ~contains(freePar,'mNu') && strcmp(NH,'OFF') 
    GridArg_i = S.LoadGridArg;
    S.LoadGridArg = {GridArg_i{:},'Extsin2T4','ON'};
    S.LoadGridFile(S.LoadGridArg{:});
    S.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',33^2);
    S.FindBestFit;
    S.FindBestFit('Mode','Imp');
    chi2_ref = S.chi2_ref;
    mNuSq_bf = S.mNuSq_bf;
    sin2T4_bf = S.sin2T4_bf;
    mNu4Sq_bf = S.mNu4Sq_bf;
    S.LoadGridArg = GridArg_i;
end

S.LoadGridFile(S.LoadGridArg{:});
S.Interp1Grid('Maxm4Sq',40^2);

if ~contains(freePar,'mNu') && strcmp(NH,'OFF') 
    S.chi2_ref  = chi2_ref;
    S.mNuSq_bf  = mNuSq_bf ;
    S.sin2T4_bf = sin2T4_bf;
    S.mNu4Sq_bf = mNu4Sq_bf;
else
    S.FindBestFit;
    S.FindBestFit('Mode','Imp');
    
    S.LoadGridFile(S.LoadGridArg{:});
    S.Interp1Grid('Maxm4Sq',40^2);
    
end

%
S.ContourPlot('NullHypothesis',NH,'ReCalcBF','OFF'); 

%% export contour

savedir = sprintf('%sksn2ana/ksn2_FitterCrossCheck/results/',getenv('SamakPath'));
MakeDir(savedir);

savename =  sprintf('%sKSN2_contour_Samak%s_%s_40eV_%s',...
    savedir,DataType,strrep(freePar,' ',''),strrep(chi2,'chi2',''));
if strcmp(NH,'ON')
  savename =  sprintf('%sKSN2_contour_Samak%s_%s_40eV_%s_NH',...
    savedir,DataType,strrep(freePar,' ',''),strrep(chi2,'chi2',''));  
end
Write2Txt('filename',savename,...
    'Format','dat','variable',[S.sin2T4_contour;S.mNu4Sq_contour],'nCol',2,'variableName','sinT4Sq m4Sq');

%% Write best fit
if strcmp(DataType,'Real')  && strcmp(NH,'OFF')
    Write2Txt('filename',[savename,'_bf'],...
    'Format','dat','variable',[S.chi2_bf; S.sin2T4_bf; S.mNu4Sq_bf;  S.mNuSq_bf ; ],'nCol',4,'variableName','chi2min sinT4Sq m4Sq mNuSq');
fprintf('Best fit: sin2T4 = %.2f m4Sq = %.2feV^2 mNuSq = %.2feV^2 chi2min = %.2f \n',S.sin2T4_bf,S.mNu4Sq_bf,S.mNuSq_bf,S.chi2_bf);
end
%% test plot
d = importdata([savename,'.dat']);
plot(d.data(:,1),d.data(:,2));
if strcmp(DataType,'Real')  && strcmp(NH,'OFF')
    dbf = importdata([savename,'_bf.dat']);
    hold on;
    plot(dbf.data(2),dbf.data(3),'x','LineWidth',2);
end
set(gca,'XScale','log');
set(gca,'YScale','log');
xlim([1e-03,0.5]);
ylim([1 40^2]);
xlabel('sint4^2'); ylabel('m4^2');
