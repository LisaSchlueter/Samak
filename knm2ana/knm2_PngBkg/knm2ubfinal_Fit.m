% unblinded fit with penning trap background slope
range     = 40;
freePar   = 'mNu E0 Bkg Norm';
chi2      = 'chi2CMShape';
DataType  = 'Real';
AnaFlag   = 'StackPixel';
RingMerge = 'Full';
DopplerEffectFlag = 'FSD';
BKG_PtSlope = 3*1e-06;
TwinBias_BKG_PtSlope = 3*1e-06;
FSDFlag   = 'KNM2';

PullFlag = 99;%99;%6;%[7,24]; %24 = 3.0 mucps/s

if strcmp(AnaFlag,'Ring')
    if strcmp(RingMerge,'Full')
        AnaStr = AnaFlag;
        SysBudget = 41;
    else
        AnaStr = sprintf('Ring%s',RingMerge);
        SysBudget = 40;
    end
else
    SysBudget = 40;
    AnaStr = AnaFlag;
end

savedir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savename = sprintf('%sknm2ubfinal_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s.mat',...
    savedir,BKG_PtSlope*1e6,DataType,range,strrep(freePar,' ',''),chi2,AnaStr,FSDFlag);

if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
elseif strcmp(chi2,'chi2Stat+')
    NonPoissonScaleFactor = 1.112;
    chi2 = 'chi2Stat';
end

if ~strcmp(chi2,'chi2Stat')
    savename = strrep(savename,'.mat',sprintf('_SysBudget%.0f.mat',SysBudget));
end

if strcmp(DataType,'Twin')
    savename = strrep(savename,'.mat',sprintf('_TwinBpng-%.1fmucpsPers.mat',1e6*TwinBias_BKG_PtSlope));
end

if any(PullFlag~=99)
    if numel(PullFlag)==2
        savename = strrep(savename,'.mat',sprintf('_pull%.0f_%.0f.mat',PullFlag(1),PullFlag(2)));
    else
        savename = strrep(savename,'.mat',sprintf('_pull%.0f.mat',PullFlag));
    end
end

if exist(savename,'file')
    load(savename,'FitResult','RunAnaArg','A');
else
    SigmaSq =  0.0124+0.0025;
    
    
    
    if strcmp(RingMerge,'None') && strcmp(chi2,'chi2CMShape') && strcmp(AnaFlag,'Ring')
        chi2tmp = 'chi2Stat';
    else
        chi2tmp  = chi2;
    end
    
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag',FSDFlag,...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',SysBudget,...
        'AnaFlag',AnaFlag,...
        'chi2',chi2tmp,...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'RingMerge',RingMerge,...
        'PullFlag',PullFlag,...;%99 = no pull
        'BKG_PtSlope',BKG_PtSlope,...
        'TwinBias_BKG_PtSlope',TwinBias_BKG_PtSlope,...
        'DopplerEffectFlag',DopplerEffectFlag};
    A = MultiRunAnalysis(RunAnaArg{:});
    %%
    A.exclDataStart = A.GetexclDataStart(range);
    
    if strcmp(DataType,'Twin')
        A.ModelObj.RFBinStep = 0.01;
        A.ModelObj.InitializeRF;
    end
    
    if  contains(freePar,'BkgPTSlope') && contains(freePar,'BkgSlope')  && strcmp(chi2,'chi2CMShape')
        A.ComputeCM('BkgPTCM','OFF','BkgCM','OFF');
    elseif  contains(freePar,'BkgSlope') && strcmp(chi2,'chi2CMShape')
        A.ComputeCM('BkgPTCM','ON','BkgCM','OFF');
    elseif contains(freePar,'BkgPTSlope') && strcmp(chi2,'chi2CMShape')
        A.ComputeCM('BkgPTCM','OFF','BkgCM','ON');
    end
    
    if strcmp(RingMerge,'None') && strcmp(chi2,'chi2CMShape') && strcmp(AnaFlag,'Ring')
        A.chi2 = chi2;
        A.ComputeCM('SysEffects',  struct(...
            'RF_EL','OFF',...   % Response Function(RF) EnergyLoss
            'RF_BF','OFF',...   % RF B-Fields
            'RF_RX','OFF',...   % Column Density, inel cross ection
            'FSD','ON',...
            'TASR','ON',...
            'TCoff_RAD','OFF',...
            'TCoff_OTHER','ON',...
            'DOPoff','OFF',...
            'Stack','OFF',...
            'FPDeff','ON',...
            'LongPlasma','ON'),...
            'BkgPTCM','ON',...
            'BkgCM','ON');
    else
        
    end
    A.Fit;
    FitResult = A.FitResult;
    MakeDir(savedir);
    save(savename,'FitResult','RunAnaArg','A','SigmaSq')
end
%%

%A.PlotFit;
fprintf('============================================\n');

fprintf(2,'m_nu^2 = %.3f + %.2f %.2f eV^2, ',FitResult.par(1),FitResult.errPos(1),FitResult.errNeg(1))
fprintf('mean err = %.2f eV^2 \n',(FitResult.errPos(1)-FitResult.errNeg(1))/2)
if strcmp(A.AnaFlag,'Ring')
    RingIdx = 2*A.nRings+9:3*A.nRings+8;
    E0eff    = FitResult.par(2)+A.ModelObj.Q_i+FitResult.par(RingIdx);
    E0effErr = sqrt(FitResult.err(2).^2+FitResult.err(RingIdx).^2);
    for r=1:A.nRings
        fprintf('E_0 Ring %.0f = %.3f + %.3f eV  \n',r,E0eff(r),E0effErr(r))
    end
    MeanE0eff = wmean(E0eff,1./E0effErr.^2);
    MeanE0effErr = sqrt(sum(E0effErr.^2))./sqrt(A.nRings);
    fprintf(2,'E_0 Mean    = %.3f + %.3f eV  \n',MeanE0eff,MeanE0effErr)
else
    fprintf('E_0 = %.3f + %.3f eV  \n',FitResult.par(2)+A.ModelObj.Q_i,FitResult.err(2))
end

if strcmp(A.AnaFlag,'Ring')
    for r=1:A.nRings
        fprintf('B Ring %.0f = %.1f + %.1f mcps  (%.0f pixel)\n',...
            r,(FitResult.par(2+r)+A.ModelObj.BKG_RateSec_i(r)).*1e3,1e3.*FitResult.err(2+r),numel(A.RingPixList{r}))
    end
    fprintf(2,'Sum B: %.1f +- %.1f mcps (%.0f pixel) \n',sum((FitResult.par(2+(1:A.nRings))+A.ModelObj.BKG_RateSec_i).*1e3),...
        1e3.*sqrt(sum(FitResult.err(2+(1:1:A.nRings)).^2)),numel(A.PixList))   
else
    fprintf('B  = %.1f + %.1f mcps  \n',...
        (FitResult.par(3)+A.ModelObj.BKG_RateSec_i).*1e3,1e3.*FitResult.err(3));
end
if strcmp(A.AnaFlag,'Ring')
    for r=1:A.nRings
        fprintf('N Ring %.0f = %.3f + %.3f  \n',r,(FitResult.par(2+A.nRings+r)+1),FitResult.err(2+A.nRings+r));
    end   
    fprintf(2,'Mean N: %.3f +- %.3f \n',mean((FitResult.par(6+(1:1:A.nRings)))+1),...
        sqrt(sum(FitResult.err(6+(1:1:A.nRings)).^2))./sqrt(A.nRings))
else
    fprintf('N  = %.3f + %.3f mcps  \n',...
        (FitResult.par(4)+1),FitResult.err(4));
end
fprintf('chi2 = %.3f (%.0f dof), p = %.3f  \n',FitResult.chi2min,FitResult.dof,1-chi2cdf(FitResult.chi2min,FitResult.dof));
fprintf('============================================\n ');

%%

Plot = 'OFF';
A.ErrorBarScaling = 50;
if strcmp(Plot,'ON')
    if strcmp(AnaFlag,'StackPixel')
        A.PlotFit('LabelFlag','FinalKNM1',...
            'saveplot','pdf',...
            'ErrorBarScaling',50,...
            'YLimRes',[-2.2,2.9],...
            'Colors','RGB',...
            'DisplayStyle','PRL',...
            'FitResultsFlag','OFF',...
            'qUDisp','Abs',...
            'TickDir','Out');
    else
        if A.nRings<=4
            A.PlotResidualsMultiRing('saveplot','pdf','YLimRes',[-3 3],'DisplayMTD','OFF','Colors','BW','XLims',[-42 138]);
        else
            A.PlotResidualsMultiRing('saveplot','pdf','YLimRes',[-3.5 3.5],'DisplayMTD','OFF','Colors','BW','XLims',[-42 138]);
            
        end
        A.PlotSpectrumMultiRing('SavePlot','ON','DisPlayStyle','Rel');
        A.PlotFitMultiRing('PlotPar','Bkg','savePlot','ON','linFitFlag','ON','RefLine','ON');
        A.PlotFitMultiRing('PlotPar','Norm','savePlot','ON','linFitFlag','ON','RefLine','ON');
        A.PlotFitMultiRing('PlotPar','E0eff','savePlot','ON','linFitFlag','ON','RefLine','ON');
    end
end
%  A.FitCM_Obj.PlotCM('qUWindowIndexMax',10,'qUWindowIndexMin',40,'saveplot',...
%         'ON','Convergence','OFF','CovMatInput',A.FitCMFracShape,'PlotEffect','total',...
%%         %         'savename','KNM2_UB1_MultiRingFull');
% % A.FitCM_Obj.PlotCorr('qUWindowIndexMax',10,'qUWindowIndexMin',90,'saveplot',...
% % 'ON','CovMatInput',A.FitCMFracShape,...
% % 'savename',sprintf('KNM2_Final_%s',AnaFlag));

%% some statistics: runs test on residuals
Residuals_abs  = A.RunData.TBDIS(A.exclDataStart:end,:)-A.ModelObj.TBDIS(A.exclDataStart:end,:);
DiagErr = sqrt(reshape(diag(A.FitCMShape),A.ModelObj.nqU,A.nRings));
Residuals_norm = Residuals_abs./DiagErr(A.exclDataStart:end,:);

MaxDeviation =  max(Residuals_norm); % maximal sigma for each rings

% runstest for each ring
pRuns = zeros(A.nRings,1);
for i=1:A.nRings
[~,pRuns(i)] = runstest(Residuals_norm(:,i));
fprintf('runs test ring %.0f, p = %.2f\n',i,pRuns(i));
end






