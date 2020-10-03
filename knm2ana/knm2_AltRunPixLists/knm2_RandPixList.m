% Uniform fit on KNM2 data
% Random half of all golden runs
% March 2020, Lisa
nFits = 500;
freePar = 'mNu E0 Bkg Norm';
DataType = 'Real';
RunList = 'KNM2_Prompt';
range = 40;               % fit range in eV below endpoint
savedir = [getenv('SamakPath'),'knm2ana/knm2_AltRunPixLists/results/'];
savename = sprintf('%sknm2_PixListRandHalf_%s_%s_%.0feV_%.0ffits.mat',...
            savedir,DataType,strrep(freePar,' ',''),range,nFits);
% draw random pixels
PixList_def = GetPixList('Knm2');
nPix = numel(PixList_def);

if exist(savename,'file')
    load(savename)
else
  SigmaSq =  0.0124+0.0025;
    PixList = cell(nFits,1);
    FitResult = cell(nFits,1);
   % E0 = knm2FS_GetE0Twins('SanityPlot','OFF');
    RunAnaArg = {'RunList',RunList,...  % define run number -> see GetRunList
        'fixPar',freePar,...         % free Parameter !!
        'DataType',DataType,...              % Real, Twin or Fake
        'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2','chi2Stat',...              % statistics only
        'NonPoissonScaleFactor',1.112,...
        'TwinBias_Q',18573.7,...
        'DopplerEffectFlag','FSD',...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq)};
    
    for i=1:nFits
        RandIndex    = randperm(nPix);         % randomly permute runlist indices
        nPixHalf     = ceil(nPix/2);           % take only half of all runs
        PixList{i}   = PixList_def(RandIndex(1:nPixHalf));
        
        %% build object of MultiRunAnalysis class
        D = MultiRunAnalysis(RunAnaArg{:},'PixList',PixList{i});
        
%         if strcmp(DataType,'Twin')
%             Sigma = std(E0);
%             FSDArg = {'SanityPlot','OFF','Sigma',Sigma};
%             D.ModelObj.LoadFSD(FSDArg{:});
%             D.ModelObj.ComputeTBDDS; D.ModelObj.ComputeTBDIS;
%         end   
        D.exclDataStart = D.GetexclDataStart(range); % find correct data, where to cut spectrum     
        %% Fit -> fit results are in property: A.FitResult
        D.InitModelObj_Norm_BKG('RecomputeFlag','ON');
        D.Fit;
        
        % save
        FitResult{i} = D.FitResult;
    end
    Q_i     = D.ModelObj.Q_i;
    RunList = D.RunList;
    PixList  = cell2mat(PixList);
    E0       = cellfun(@(x) x.par(2),FitResult);
    E0Err    = cellfun(@(x) x.err(2),FitResult);
    mNuSq    = cellfun(@(x) x.par(1),FitResult);
    mNuSqErr = cellfun(@(x) x.err(1),FitResult);
    
    MakeDir(savedir);
    save(savename,'PixList','FitResult','RunAnaArg','Q_i','RunList',...
        'E0','E0Err','mNuSq','mNuSqErr');
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


