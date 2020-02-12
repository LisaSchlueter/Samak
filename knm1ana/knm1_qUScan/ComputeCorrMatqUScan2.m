function ComputeCorrMatqUScan2(obj,varargin)
RunList = 'KNM1';
DataType = 'Twin';
chi2 = 'chi2CMShape';
SysBudet = 22;
fixPar = '5 6 7 8 9 10 11';
firstPoint = 2;
lastPoint = 20;
nSamples = 100;
FSDFlag ='SibilleFull'; 
T = MultiRunAnalysis('RunList',RunList,...
    'chi2',chi2,'DataType',DataType,...
    'fixPar',fixPar,...
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1.064,...
    'minuitOpt','min ; migrad',...
    'FSDFlag',FSDFlag,...
    'ELossFlag','KatrinT2',...
    'SysBudget',22);
T.ComputeCM;
TBDIS = mvnrnd(T.RunData.TBDIS',T.FitCMShape,nSamples)';
%%
progressbar('Fit twins with stat. fluctuations')
T.fixPar = '5 6 7 8 9 10 11'; % free nu mass

parqU = zeros(11,numel(firstPoint:lastPoint),nSamples);
errqU = zeros(11,numel(firstPoint:lastPoint),nSamples);
chi2qU = zeros(numel(firstPoint:lastPoint),nSamples);
dofqU = zeros(numel(firstPoint:lastPoint),nSamples);

for i=1:nSamples
    progressbar(i/nSamples)
    
    T.RunData.TBDIS = TBDIS(:,i);
    T.RunData.TBDISE = sqrt(TBDIS(:,i));
   
    [parqU(:,:,i), errqU(:,:,i), chi2qU(:,i), dofqU(:,i)] = ...
    T.qUScan('firstPoint',firstPoint,'lastPoint',lastPoint,'saveplot','OFF','RecomputeFlag','OFF','CorrMean','OFF');
end

savedir = [getenv('SamakPath'),'knm1ana/knm1_qUScan/results/'];
system(['mkdir -p ',savedir]);
savefile = [savedir,sprintf('CorrMat_qUScan_%s_%.0f-%.0f_%s_%.0f_%s.mat',RunList,firstPoint,lastPoint,chi2,nSamples)];

save(savefile,'parqU','errqU','chi2qU','dofqU');

CorrMatE0 = corr(squeeze(parqU(2,:,:))');
CorrMatmNuSq = corr(squeeze(parqU(1,:,:))');
CorrMatB = corr(squeeze(parqU(3,:,:))');
CorrMatN = corr(squeeze(parqU(4,:,:))');

save(savefile,'CorrMatmNuSq','CorrMatE0','CorrMatB','CorrMatN','-append');
end


