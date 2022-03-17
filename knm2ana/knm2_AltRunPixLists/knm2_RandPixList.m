% Uniform fit on KNM2 data
% Random half of all golden runs
% March 2022, Lisa
nRandPixList = 2000;
freePar = 'mNu E0 Bkg Norm';
DataType = 'Real';
RunList = 'KNM2_Prompt';
range = 40;               % fit range in eV below endpoint
FSDFlag = 'KNM2_0p1eV';
BKG_PtSlope = 3*1e-06;
savedir = [getenv('SamakPath'),'knm2ana/knm2_AltRunPixLists/results/'];
savename = sprintf('%sknm2_PixListRandHalf_%s_%s_%.0feV_%.0ffits_%s_BkgPT%.2gmucpsPers.mat',...
            savedir,DataType,strrep(freePar,' ',''),range,nRandPixList,FSDFlag,BKG_PtSlope*1e6);
       
% draw random pixels
if exist(savename,'file')
    load(savename)
else
  SigmaSq =  0.0124+0.0025;
    PixList = cell(nFits,1);
    FitResult = cell(nFits,1);
 
    RunAnaArg = {'RunList',RunList,...  % define run number -> see GetRunList
        'fixPar',freePar,...         % free Parameter !!
        'DataType',DataType,...              % Real, Twin or Fake
        'FSDFlag',FSDFlag,...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2','chi2Stat',...              % statistics only
        'NonPoissonScaleFactor',1.112,...
        'TwinBias_Q',18573.7,...
        'DopplerEffectFlag','FSD',...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'BKG_PtSlope',BKG_PtSlope};
    
    M = MultiRunAnalysis(RunAnaArg{:});
    M.exclDataStart = M.GetexclDataStart(range);
    
    PixList_gold = M.PixList; % regular golden pixel list
    nPix = numel(M.PixList);
    nPixHalf     = ceil(nPix/2);           % take only half of all runs
    
    D = copy(repmat(M,nRandPixList,1));
    D = reshape(D,numel(D),1);
    
    FitResults = cell(nRandPixList,1);
    PixLists   = cell(nRandPixList,1);
    
    parfor i=1:nRandPixList
        RandIndex    = randperm(nPix);         % randomly permute runlist indices
        PixLists{i}   = PixList_gold(RandIndex(1:nPixHalf));
        
        D(i).PixList =  PixLists{i};
        D(i).StackRuns('CutOnSC','OFF','SCsigma',3,'CutOnFitSingleRuns','OFF');
        
        D(i).SimulateStackRuns;
        D(i).exclDataStart = D(i).GetexclDataStart(40);
        D(i).Fit; 
        
        FitResults{i} = D(i).FitResult;
    end
    
    save(savename,'FitResults','PixLists','M','RunAnaArg');
    fprintf('save to %s \n',savename);
    
%     Q_i     = D.ModelObj.Q_i;
%     RunList = D.RunList;
%     PixList  = cell2mat(PixList);
%     E0       = cellfun(@(x) x.par(2),FitResult);
%     E0Err    = cellfun(@(x) x.err(2),FitResult);
%     mNuSq    = cellfun(@(x) x.par(1),FitResult);
%     mNuSqErr = cellfun(@(x) x.err(1),FitResult);
%     
%     MakeDir(savedir);
%     save(savename,'PixList','FitResult','RunAnaArg','Q_i','RunList',...
%         'E0','E0Err','mNuSq','mNuSqErr');
end

%%
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
meanE0 =wmean(E0,1./E0Err.^2);
h1 =histogram(E0-meanE0,'FaceColor',rgb('SkyBlue'),'FaceAlpha',1);
xlabel(sprintf('{\\itE}_0^{fit} - \\langle{\\itE}_0^{fit}\\rangle (eV)'));
ylabel('Occurrence');
PrettyFigureFormat('FontSize',22)
t = title(sprintf('%.0f stacked runs - uniform fit to %.0f random pixels',numel(RunList),size(PixList,2)),...
    'FontWeight','normal');
leg = legend(sprintf('\\langle{\\itE}_0^{fit}\\rangle = %.2f eV , \\sigma =  %.0f meV',meanE0+Q_i,1e3*std(E0)));
leg.EdgeColor = rgb('Silver');


