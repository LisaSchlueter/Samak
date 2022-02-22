% save my contour to .dat file for other fitters to read

% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';
nGridSteps = 50;
range = 40;
freePar = 'E0 Norm Bkg';

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
    'LoadGridArg',{'mNu4SqTestGrid',5,'Extsin2T4','OFF','ExtmNu4Sq','ON'}};
S = SterileAnalysis(SterileArg{:});
%%
S.InterpMode = 'spline';

S.LoadGridArg = {S.LoadGridArg{:}};
S.LoadGridFile(S.LoadGridArg{:});

m4Sq   = S.mNu4Sq;
sin2t4 = S.sin2T4;
chi2   = S.chi2';

close all
GetFigure;
contour(sin2t4,m4Sq,chi2-min(min(chi2)),[5.99,5.99],'LineWidth',2);
set(gca,'XScale','log');
set(gca,'YScale','log');

S.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',38.2^2);
S.FindBestFit;
chi2_ref = S.chi2_ref;

m4Sq_interp   = S.mNu4Sq;
sin2t4_interp = S.sin2T4;
chi2_interp   = S.chi2;

hold on;

contour(sin2t4_interp,m4Sq_interp,chi2_interp-chi2_ref,[5.99,5.99],'LineWidth',2);


%% export chi2 map

savedir = sprintf('%sksn2ana/ksn2_Export/results/',getenv('SamakPath'));
MakeDir(savedir);

savename =  sprintf('%sKSN2_chi2map.h5', savedir);

system(sprintf('rm %s',savename));

h5create(savename,'/chiSqmin',[1,1]);
h5write(savename, '/chiSqmin', chi2_ref);

h5create(savename,'/chiSq',[50,50]);
h5write(savename, '/chiSq', chi2);

h5create(savename,'/m4Sq',[50,50]);
h5write(savename, '/m4Sq', m4Sq);

h5create(savename,'/sint4Sq',[50,50]);
h5write(savename, '/sint4Sq', sin2t4);

% h5create(savename,'/chiSq_interp',[1e3,1e3]);
% h5write(savename, '/chiSq_interp', chi2_interp);
% 
% h5create(savename,'/m4Sq_interp',[1e3,1e3]);
% h5write(savename, '/m4Sq_interp', m4Sq_interp);
% 
% h5create(savename,'/sint4Sq_interp',[1e3,1e3]);
% h5write(savename, '/sint4Sq_interp', sin2t4_interp);

h5disp(savename)

%%
savenameData =  sprintf('%sKSN2_Data.h5', savedir);
system(sprintf('rm %s',savenameData));

qU    = A.RunData.qU(A.exclDataStart:end);
Counts   = A.RunData.TBDIS(A.exclDataStart:end);
Time_sec = A.RunData.qUfrac(A.exclDataStart:end).*A.RunData.TimeSec;
CovMat = A.FitCMShape(A.exclDataStart:end,A.exclDataStart:end);

h5create(savenameData,'/RetardingEnergy_eV',[28,1]);
h5write(savenameData, '/RetardingEnergy_eV', qU);

h5create(savenameData,'/Counts',[28,1]);
h5write(savenameData, '/Counts', Counts);

h5create(savenameData,'/Time_sec',[28,1]);
h5write(savenameData, '/Time_sec', Time_sec);

h5create(savenameData,'/CovarianceMatrix',[28,28]);
h5write(savenameData, '/CovarianceMatrix', CovMat);

h5disp(savenameData);
%%
Write2Txt('filename',savename,...
    'Format','dat',...
    'variable',[S.sin2T4_contour;S.mNu4Sq_contour],'nCol',2,'variableName','sinT4Sq m4Sq');

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
