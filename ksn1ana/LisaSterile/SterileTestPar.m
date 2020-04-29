% test parallel grid search for sterile analysis
% Lisa, April 2020
savedir = [getenv('SamakPath'),'ksn1ana/LisaSterile/results/'];
MakeDir(savedir);
range = 95;
nGridSteps = 50;
chi2 = 'chi2CMShape';
if strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor=1.064;
else
    NonPoissonScaleFactor=1;
end
DataType = 'Twin';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'ON';
if strcmp(SmartGrid,'ON')
    AddSin2T4 = 0.1;
   extraStr = sprintf('_SmartGrid%.0e',AddSin2T4);
else
    extraStr = '';
end
savefile = sprintf('%sKSN1_GridSearch_%s_%s_%s_%.0feVrange_%s_%.0fnGrid%s.mat',...
    savedir,RunList,DataType,strrep(freePar,' ',''),range,chi2,nGridSteps,extraStr);
if exist(savefile,'file')
    load(savefile)
else
    if range<=40
        FSDFlag = 'Sibille0p5eV';
    elseif range>40
        FSDFlag = 'SibilleFull';
    end
        
    RunAnaArg = {'RunList',RunList,...
        'fixPar',freePar,...
        'DataType',DataType,...
        'FSDFlag',FSDFlag,...
        'ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel',...
        'chi2',chi2,...
        'ROIFlag','Default',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF',...
        'ISCSFlag','Edep',...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'TwinBias_Q',18573.73};
    T = MultiRunAnalysis(RunAnaArg{:});
    T.exclDataStart = T.GetexclDataStart(range);
    
    if strcmp(DataType,'Real')
        T.fixPar = [freePar,'mnu4Sq , sin2T4'];
        T.InitFitPar;
        T.pullFlag = 9;
    end
    
    T.Fit;
    FitResults_ref = T.FitResult;
    chi2_ref       = T.FitResult.chi2min;
    
    if strcmp(DataType,'Real')
        T.fixPar = freePar;
        T.InitFitPar;
        T.pullFlag = 99;
    end
    %% define grid
    switch SmartGrid
        case 'OFF'
            sin2T4      = logspace(-3,log10(0.5),nGridSteps); %linspace(0.001,0.5,nGridSteps)
            mnu4Sq      = logspace(-1,log10((range+5)^2),nGridSteps)';
            mnu4Sq      = repmat(mnu4Sq,1,nGridSteps);
            sin2T4      = repmat(sin2T4,nGridSteps,1);
        case 'ON'
           [mnu4Sq,sin2T4] = GetSmartKsn1Grid('range',range,'ConfLevel',95,...
                'nGridSteps',nGridSteps,'AddSin2T4',AddSin2T4,...
                'SanityPlot','ON');
    end
    %% make copy of models for parallel computing
    D = copy(repmat(T,nGridSteps.*nGridSteps,1));

    %% scan over msq4-sin2t4
    D              = reshape(D,numel(D),1);
    chi2Grid       = zeros(nGridSteps*nGridSteps,1);
    FitResultsGrid = cell(nGridSteps*nGridSteps,1);
    mnu4Sq_Grid    = reshape(mnu4Sq',nGridSteps*nGridSteps,1);
    sin2T4_Grid    = reshape(sin2T4',nGridSteps*nGridSteps,1);
    
    parfor i= 1:(nGridSteps*nGridSteps)
        D(i).SimulateStackRuns;
        D(i).ModelObj.SetFitBiasSterile(mnu4Sq_Grid(i),sin2T4_Grid(i));
        D(i).Fit
        chi2Grid(i) = D(i).FitResult.chi2min;
        FitResultsGrid{i} = D(i).FitResult;
    end
    
    chi2       = reshape(chi2Grid,nGridSteps,nGridSteps);
    FitResults = reshape(FitResultsGrid,nGridSteps,nGridSteps);
    %% save
    save(savefile,'chi2_ref','FitResults_ref','RunAnaArg',...
                 'chi2','mnu4Sq','sin2T4','FitResults');
    fprintf('save file to %s \n',savefile)
    
    try
        %% get contour at X sigma
        [mnu4Sq_contour95, sin2T4_contour95] = GetSterileContour(mnu4Sq,sin2T4,chi2,chi2_ref,95);
        save(savefile,'mnu4Sq_contour95','sin2T4_contour95','-append')
         [mnu4Sq_contour90, sin2T4_contour90] = GetSterileContour(mnu4Sq,sin2T4,chi2,chi2_ref,90);
        save(savefile,'mnu4Sq_contour90','sin2T4_contour90','-append')
    catch
    end
end

% %% get contour at X sigma
% CL = 95;
% [mnu4Sq_contour, sin2T4_contour] = GetSterileContour(mnu4Sq,sin2T4,chi2,chi2_ref,CL);
% %%
% GetFigure;
% plot(sin2T4_contour,mnu4Sq_contour,'.-','LineWidth',2,'MarkerSize',20);
% set(gca,'YScale','log');
% set(gca,'XScale','log');
% xlim([1e-02 1])
% PrettyFigureFormat;
% xlabel('|U_{e4}|^2');
% ylabel(sprintf('{\\itm}_4 (eV^2)'));
% 
% plotdir = strrep(savedir,'results','plots');
% MakeDir(plotdir);
% plotname = strrep(strrep(savefile,'results','plots'),'.mat','.png');
% print(plotname,'-dpng','-r450');