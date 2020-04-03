% KNM2 twin multi ring fit
% March 2020, Lisa

RunList = 'KNM2_Prompt';
E0 = knm2FS_GetE0Twins('SanityPlot','OFF');
range = 40;
chi2 = 'chi2Stat';
pullFlag = 4;
freePar = 'mNu E0 Norm Bkg';
DataType = 'Real';
RingMerge = 'Full';

savedir = [getenv('SamakPath'),'knm2ana/knm2_MultiRingFit/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('knm2_MultiRingFit_%s_%s_%s_%s_pull%.0f_%.0feVrange_RingMerge%s.mat',...
    DataType, RunList,chi2,strrep(freePar,' ',''),pullFlag,range,RingMerge)];

if exist(savename,'file')
    load(savename)
    MR.FitResult = FitResults;
else
    % read data and set up model
    RunArg = {'RunList',RunList,...
        'chi2',chi2,...
        'DataType',DataType,...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag','Sibille0p5eV',...
        'ELossFlag','KatrinT2',...
        'SysBudget',34,...
        'AnaFlag','Ring',...
        'RingMerge',RingMerge,...
        'chi2',chi2,...
        'pullFlag',pullFlag,...
        'TwinBias_Q',E0,...
        'ROIFlag','14keV',...
        'MosCorrFlag','OFF',...
        'NonPoissonScaleFactor',1};
    
    MR = MultiRunAnalysis(RunArg{:});
    MR.exclDataStart = MR.GetexclDataStart(range);
    if ~strcmp(chi2,'chi2Stat')
        MR.NonPoissonScaleFactor = 1.112;
        MR.SetNPfactor; % convert to right dimension (if multiring)
        MR.chi2 = chi2;
        MR.ComputeCM;
    end
    MR.ModelObj.RFBinStep = 0.02;
    MR.ModelObj.InitializeRF;
    MR.Fit('SaveFit','ON','CATS','ON');
    FitResults = MR.FitResult;
    save(savename,'FitResults','RunArg','MR');
    
    %% fit with broadening of RF + broadening/shift of FSD
    TimeSec = zeros(3,1);
    TimeSec(1) = sum(MR.SingleRunData.TimeSec(1:171));
    TimeSec(2) = sum(MR.SingleRunData.TimeSec(172:268));
    TimeSec(3) = sum(MR.SingleRunData.TimeSec(269:361));
    MultiWeights = TimeSec./sum(TimeSec);
    MultiPos = [E0(1),E0(end-120),E0(end)]';
    MultiPosRel = repmat(MultiPos-wmean(MultiPos,MultiWeights),1,MR.nRings);
    Sigma = repmat(std(E0),3,MR.nRings);
    FSDArg = {'MultiPos',MultiPosRel,'MultiWeight',MultiWeights,...
        'SanityPlot','ON','Sigma',Sigma};
    MR.ModelObj.LoadFSD(FSDArg{:});
    MR.ModelObj.RFBinStep = 0.002; % finer binning to resolve effect
    MR.ModelObj.InitializeRF;
    MR.Fit('SaveFit','OFF');
    FitResults_imp = MR.FitResult;
    save(savename,'FitResults_imp','FSDArg','E0','-append');
end
%% result
fprintf('--------------------------------------\n')
fprintf('mNuSq = %.4f (%.3f +%.3f) eV^2  (ref) \n',FitResults.par(1),FitResults.errNeg(1),FitResults.errPos(1))
fprintf('mNuSq = %.4f (%.3f +%.3f) eV^2  (imp) \n',FitResults_imp.par(1),FitResults_imp.errNeg(1),FitResults_imp.errPos(1))
fprintf('Average mNuSq sensitivity = %.3f eV^2 \n',0.5*(-FitResults_imp.errNeg(1)+FitResults_imp.errPos(1)))
fprintf('--------------------------------------\n')
fprintf('E0    = %.0e (+-%.2f) eV (ref) \n',FitResults.par(2)+MR.ModelObj.Q_i-mean(E0),FitResults.err(2))
fprintf('E0    = %.0e (+-%.2f) eV (imp) \n',FitResults_imp.par(2)+MR.ModelObj.Q_i-mean(E0),FitResults_imp.err(2))
fprintf('--------------------------------------\n')

% MR.PlotFitMultiRing('PlotPar','qU','linFit','ON','savePlot','ON','Blind','OFF');
% plotdir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/MultiRingFit/',MR.DataSet)];
% freePar = ConvertFixPar('freePar',MR.fixPar,'nPar',MR.nPar,'nPixels',numel(MR.RunData.MACE_Ba_T),'Mode','Reverse');
% plotname = [plotdir,sprintf('FPDViewer_MultiRing%s_%s_%s_freePar%s_%s.pdf',...
%     MR.RingMerge,MR.RunData.RunName,MR.chi2,freePar,'qU')];
% PlotRingWiseFitPar(MR,'SaveAs',plotname); %FPD viewer
% %
% %%
%MR.FitCM_Obj.PlotCM('qUWindowIndexMax',10,'saveplot','ON','PlotEffect','total',...
%        'Convergence','OFF','CovMatInput',MR.FitCMFracShape,...
%       'savedir',[getenv('SamakPath'),'knm2ana/knm2_systematics/plots/'],'savename','FS_MultiRing');
%   
%   
  