% Fit fraction of T- ion to KNM2 data
% Fit three RW setting periods seperately
DataSet = 'KNM2';
DataType   = 'Real';
ELossFlag  = 'KatrinT2';
AnaFlag    = 'StackPixel'; % uniform FPD
chi2       = 'chi2Stat';
range = 90;
RecomputeFlag = 'OFF';

switch DataSet
    case 'KNM1'
        RunList = 'KNM1';
        FSDFlag    = 'Sibille0p5eV';
        fixPar = 'mNu E0 Bkg Norm FracTm'; % free parameter
        NonPoissonScaleFactor = 1.064;
        if range==40
            exclDataStart = 14; % 2==90eV, 14==40eV
        elseif range==90
            exclDataStart = 2; % 2==90eV, 14==40eV
        end
    case 'KNM2'
        RunList  = 'KNM2_Prompt';
        FSDFlag  = 'BlindingKNM2';
        fixPar = 'E0 Bkg Norm FracTm'; % free parameter
        if range==40
            exclDataStart = 11;
        elseif range==90
            exclDataStart = 1;
        end
end

RunAnaArg = {'DataType',DataType,...
    'FSDFlag',FSDFlag,...
    'ELossFlag',ELossFlag,...
    'AnaFlag',AnaFlag,...
    'chi2',chi2,...
    'exclDataStart',exclDataStart,...
    'minuitOpt','min;minos',...
    'fixPar',fixPar};

%%labeling
savedir = [getenv('SamakPath'),'knm2ana/knm2_Plasma/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('knm2_FitTminusIon_%s_%s_%s_%.0feV.mat',DataSet,chi2,strrep(fixPar,' ',''),range)];

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
else
    % set up model
    M = MultiRunAnalysis(RunAnaArg{:},'RunList',RunList);%
    if strcmp(DataSet,'KNM1')
        M.NonPoissonScaleFactor = NonPoissonScaleFactor;
    end
    % fit
    M.Fit;
    FitResult =  M.FitResult;
    save(savename,'RunAnaArg','FitResult');
    
    %% plot
    M.PlotFit('saveplot','pdf',...
        'ErrorBarScaling',1,...
        'YLimRes',[-3,3],...
        'Colors','RGB',...
        'FitResultsFlag','OFF',...
        'qUDisp','Rel');
end


