
% KNM1 Multi ring fit Results
% Golden Run List
% Golden Pixel List
%% settings
RunList               = 'KNM1';
exclDataStart         = 14; % 40eV range = 27 subruns
chi2                  = 'chi2CMShape';
fixPar = 'mNu E0 Bkg Norm qU mTSq';
pullFlag = [4,8];
RecomputeFlag = 'OFF';
    RunArg = {'RunList',RunList,...
        'chi2','chi2Stat','DataType','Real',...
        'exclDataStart',exclDataStart,...
        'fixPar',fixPar,...% free Parameter !
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',1.064,...
        'minuitOpt','min ; minos',...
        'FSDFlag','Sibille0p5eV',...
        'ELossFlag','KatrinT2',...
        'SysBudget',22,...
        'AnaFlag','Ring',...
        'RingMerge','Full',...
        'chi2',chi2,...
        'pullFlag',pullFlag};
    % Init Model Object and covariance matrix object
    
    Real = MultiRunAnalysis(RunArg{:});
    
    savedir = [getenv('SamakPath'),'knm1ana/knm1_MultiRingFit/results/'];
    MakeDir(savedir);
    savename = [savedir,sprintf('knm1_MultiRingFit_%s_%s_%s_pull%.0f_%.0frange_RingMerge%s.mat',...
        RunList,chi2,strrep(fixPar,' ',''),pullFlag,exclDataStart,Real.RingMerge)];
    if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
        d = importdata(savename,'FitResult');
        FitResult = d.FitResult;
        Real.FitResult = FitResult;
        Real.FitResult.par(2) = FitResult.par(2)-Real.ModelObj.Q_i;
        if isfield(d,'ScanResults')
            ScanResults = d.ScanResults;
        end
    else
        Real.Fit;
        FitResult = Real.FitResult;
        FitResult.par(2) = FitResult.par(2)+Real.ModelObj.Q_i;
        RingList = Real.RingList;
        RingPixList = Real.RingPixList;
        nPixRing = cell2mat(cellfun(@(x) numel(x),RingPixList,'UniformOutput',0));
        save(savename,'FitResult','RingPixList','RingList','nPixRing','RunArg');
        
        % ScanResults = Real.GetAsymFitError('Mode','Smart','Parameter','mNu','ParScanMax',1.5);
        % save(savename,'ScanResults','-append');
    end
    %%
    FitResult.par(1)
    FitResult.chi2min
%     %%
%     Real.PlotFit('Ring',2,...
%     'saveplot','pdf',...
%     'ErrorBarScaling',50,...
%     'YLimRes',[-2.2,2.9],...
%     'Colors','RGB',...
%     'FitResultsFlag','OFF',...
%     'qUDisp','Rel');
    %% plots
% % plot ringwise fit results: normal x,y plot
% try
%     Real.PlotFitMultiRing('PlotPar','mTSq','linFit','OFF','savePlot','ON');
%     t = 1;
% catch
%     fprintf('Cannot plot this, because mTSq was fixed \n')
%     t=0;
% end
% Real.PlotFitMultiRing('PlotPar','qU','linFit','OFF','savePlot','ON');
% Real.PlotFitMultiRing('PlotPar','Norm','linFit','OFF','savePlot','ON');
% Real.PlotFitMultiRing('PlotPar','Bkg','linFit','OFF','savePlot','ON');
% 
% % plot ringwise fitresult: FPDviewer
% plotdir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/MultiRingFit/',Real.DataSet)];
% MakeDir(plotdir);
% 
% if t==1
%     [plotHandlemTSq, ~]  = PlotRingWiseFitPar(Real,'PlotPar','mTSq');
%     plotname1= [plotdir,sprintf('FPDViewer_MultiRing%s_%s_%s_fixPar%s_%s.pdf',...
%         Real.RingMerge,Real.RunData.RunName,Real.chi2,strrep(strrep(Real.fixPar,'fix ',''),' ;',''),'mTSq')];
%     export_fig(plotHandlemTSq,plotname1,'-painters');
% end
% 
% [plotHandleqU, ~]  = PlotRingWiseFitPar(Real,'PlotPar','qU');
% plotname2= [plotdir,sprintf('FPDViewer_MultiRing%s_%s_%s_fixPar%s_%s.pdf',...
%     Real.RingMerge,Real.RunData.RunName,Real.chi2,strrep(strrep(Real.fixPar,'fix ',''),' ;',''),'qU')];
% export_fig(plotHandleqU,plotname2,'-painters');
% 
% 
% 
% 
% 
