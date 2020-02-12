FirstRun = 51442;
LastRun = 1e9;
StepSize = 25;
RunExcList = [51691,51641,51861,51896,51897,51899,51918];
Version = 'RunSummary-Durable2a-fpd00';
fixE0 = 'OFF';

savedir = [getenv('SamakPath'),'knm1ana/knm1Twins/results/'];
if strcmp(fixE0,'ON')
    savefile = sprintf('FitResults_TwinStackingRunList_%.0f-%.0g_%.0fRunSteps_fixE0.mat',FirstRun,LastRun,StepSize);
else
    savefile = sprintf('FitResults_TwinStackingRunList_%.0f-%.0g_%.0fRunSteps.mat',FirstRun,LastRun,StepSize);
end

if exist([savedir,savefile],'file')
load([savedir,savefile]);
else
% Read All KNM1 HD5 File
tmp = dir([getenv('SamakPath'), '/tritium-data/hdf5/',GetDataSet(FirstRun), '/*.h5']);
h5list = arrayfun(@(x) x.name,tmp,'UniformOutput',0);
h5list =  extractBefore(h5list,'.h5');
h5list = str2double(extractAfter(h5list,Version));
h5list(isnan(h5list)) = []; % %delete everything that has a different version

% Truncate:
h5list(h5list<FirstRun)=[];
h5list(h5list>LastRun)=[];
h5list(ismember(h5list,RunExcList))=[];
HDF5readallruns('h5runlist',h5list,'reConvert','OFF','DataSet',GetDataSet(FirstRun)); %looks for unconverted runs and converts if needed

RunList=sort(h5list);
nRuns = numel(RunList);

nIteration = floor(nRuns/StepSize);

mNuSq = zeros(nIteration+1,1);
E0    = zeros(nIteration+1,1);
chi2min = zeros(nIteration+1,1);
qUfracStd = zeros(nIteration+1,40);
if strcmp(fixE0,'ON')
    fixPar = '2 5 6 7 8 9 10 11';
else
    fixPar = '5 6 7 8 9 10 11';
end

for i=1:nIteration
    progressbar(i/nIteration)
    M = MultiRunAnalysis('RunList',RunList(((i-1)*StepSize+1):i*StepSize)','DataType','Twin',...
        'Twin_SameCDFlag','ON','Twin_SameIsotopFlag','ON','Twin_SameqUFlag','ON','Twin_SameqUfracFlag','OFF');
    M.exclDataStart = 14;
    M.fixPar = fixPar;
    M.Fit;
    mNuSq(i) = M.FitResult.par(1);
    E0(i)    = M.FitResult.par(2);
    chi2min(i) = M.FitResult.chi2min;
    qUfracStd(i,:) = std(M.SingleRunData.qUfrac');
end

M = MultiRunAnalysis('RunList',RunList((nIteration*StepSize+1):end)','DataType','Twin',...
     'Twin_SameCDFlag','ON','Twin_SameIsotopFlag','ON','Twin_SameqUFlag','ON','Twin_SameqUfracFlag','OFF');
M.exclDataStart = 14;
M.fixPar = fixPar;
M.Fit;
mNuSq(end) = M.FitResult.par(1);
E0(end)    = M.FitResult.par(2);
chi2min(end) = M.FitResult.chi2min;
qUfracStd(end,:) = std(M.SingleRunData.qUfrac');

save([savedir,savefile],'mNuSq','chi2min','qUfracStd','M')
end
%%
x = qUfracStd(:,4);
%x = StepSize:StepSize:((nIteration+1)*StepSize);
plot(x,mNuSq,'x');
PrettyFigureFormat;
ylabel('nu mass shift (eV^2)')

%p = plot(x,qUfracStd')

%%
plot(E0,mNuSq,'x')
PrettyFigureFormat;
xlabel('endpoint shift (eV)');
ylabel('m_\nu shift (eV^)');
