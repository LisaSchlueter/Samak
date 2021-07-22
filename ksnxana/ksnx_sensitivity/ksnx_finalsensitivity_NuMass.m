% KSNX KATRIN final sensitivity
% Lisa, May 2020
% 750 days
% 130 mcps background
% TDR like systematics
% neutrino mass sensitivity


%% settings for runanalysis
mNu4Sq_true = 0;
FixmNuSq_all = (mNu4Sq_true-0.5):0.1:(mNu4Sq_true+0.5);
nGridSteps = 30;
DataType   = 'Fake';
chi2       = 'chi2Stat';
RecomputeFlag = 'OFF';
savedir = [getenv('SamakPath'),'ksnxana/ksnx_sensitivity/results/'];

savefileCombi = sprintf('%sksnx_%s_%s_NuMassSensitivityGridSearch_CombiSTEREO-%s_mNuTrue%.2feV2.mat',...
    savedir,DataType,chi2,'OFF',mNu4Sq_true);

if exist(savefileCombi,'file') && strcmp(RecomputeFlag,'OFF')
    load(savefileCombi);
else
    if mNu4Sq_true==0
        FakeInitFile = @ref_KSNX_KATRIN_Final2;
    elseif mNu4Sq_true==0.1
        FakeInitFile = @ref_KSNX_KATRIN_Final_mNu0p1;
    elseif mNu4Sq_true==0.2
        FakeInitFile = @ref_KSNX_KATRIN_Final_mNu0p2;
    elseif   mNu4Sq_true==0.5
        FakeInitFile = @ref_KSNX_KATRIN_Final_mNu0p5;
    end
    freePar    = 'E0 Norm Bkg';
    
    % init
    chi2min   =  zeros(numel(FixmNuSq_all),1);
    sin2T4_bf =  zeros(numel(FixmNuSq_all),1);
    mNu4Sq_bf =  zeros(numel(FixmNuSq_all),1);
    
    for i=1:numel(FixmNuSq_all)
        
        RunAnaArg = {'RunNr',1,...
            'fixPar',freePar,...
            'DataType','Fake',...
            'FSDFlag','KNM2_0p5eV',...
            'ELossFlag','KatrinT2A20',...
            'AnaFlag','StackPixel',...
            'chi2','chi2Stat',...
            'SynchrotronFlag','ON',...
            'AngularTFFlag','ON',...
            'ISCSFlag','Edep',...
            'TwinBias_Q',18573.73,...
            'SysBudget',66,...
            'pullFlag',99,...
            'NonPoissonScaleFactor',1,...
            'FakeInitFile',FakeInitFile,...
            'PixList',1:136,...
            'DopplerEffectFlag','FSD'};
        
        T = RunAnalysis(RunAnaArg{:});
        T.chi2 = chi2;
        %% settings sterile class
        SterileArg = {'RunAnaObj',T,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
            'nGridSteps',nGridSteps,...
            'SmartGrid','OFF',...
            'RecomputeFlag','OFF',...
            'SysEffect','all',...
            'RandMC','OFF',...
            'range',40,...
            'LoadGridArg',{'ExtmNu4Sq','ON','mNu4SqTestGrid',2}};
        
        S = SterileAnalysis(SterileArg{:});
        %    S.GridSearch(S.LoadGridArg{:},'FixmNuSq',FixmNuSq_all(i));
        T.ModelObj.mnuSq_i = FixmNuSq_all(i);
        S.LoadGridFile(S.LoadGridArg{:},'FixmNuSq',FixmNuSq_all(i));
        S.InterpMode = 'Mix';
        S.Interp1Grid('MaxM4Sq',38^2)
        S.GridPlot('BestFit','ON','Contour','ON','SavePlot','OFF',...
            'ExtraStr',sprintf('_mNuSq%.2geV2',FixmNuSq_all(i))); 
      %  close all
        
        sin2T4_bf(i) = S.sin2T4_bf;
        mNu4Sq_bf(i) = S.mNu4Sq_bf;
        
        if S.chi2_bf<S.chi2_Null
            chi2min(i) = S.chi2_bf;
        else
            chi2min(i) = S.chi2_Null;
        end
    end
    
    %%
    T.fixPar = ConvertFixPar('freePar','mNu E0 Norm Bkg','nPar',T.nPar);
    T.ModelObj.SetFitBias(0);
    T.ModelObj.SetFitBiasSterile(0,0);
    T.ModelObj.mnuSq_i = 0;
    T.Fit;
    FitResult_NoSterile_mNuFree  = T.FitResult;
    %%
    T.fixPar = ConvertFixPar('freePar','E0 Norm Bkg','nPar',T.nPar);
    T.ModelObj.SetFitBias(0);
    T.ModelObj.SetFitBiasSterile(0,0);
    T.ModelObj.mnuSq_i = 0;
    T.Fit;
    FitResult_NoSterile_mNu0  = T.FitResult;
    %%
    save(savefileCombi,'chi2min','mNu4Sq_bf','sin2T4_bf','FixmNuSq_all',...
        'FitResult_NoSterile_mNuFree','FitResult_NoSterile_mNu0');
    
   
end

%% interpolate
x = linspace(min(FixmNuSq_all),max(FixmNuSq_all),1e3);
y = interp1(FixmNuSq_all,chi2min,x,'spline');
% best fit
chi2bf      = min(y);
Idx_bf      = find(y==chi2bf);
mNuSq_bf    = x(Idx_bf);
% uncertainties
mNuSqDown = interp1(y(x<mNuSq_bf),x(x<mNuSq_bf),chi2bf+1,'spline');
mNuSqUp   = interp1(y(x>mNuSq_bf),x(x>mNuSq_bf),chi2bf+1,'spline');
mNuSqErrDown = mNuSq_bf-mNuSqDown;
mNuSqErrUp = mNuSqUp-mNuSq_bf;
mNuSqErr   = 0.5.*(mNuSqErrDown+mNuSqErrUp);
%%
dof = FitResult_NoSterile_mNuFree.dof-2;

GetFigure;
plot(x,y,'LineWidth',2);
hold on;
pNone = plot(NaN,NaN,'w.');
PrettyFigureFormat('FontSize',22);
leg = legend(sprintf('3\\nu+1 model: {\\itm}_\\nu^2 = %.2f ^{-%.2g}_{+%.3f} eV^2',mNuSq_bf,mNuSqErrDown,mNuSqErrUp),...
             sprintf('3\\nu   model: {\\itm}_\\nu^2 = %.2f ^{- %.3f}_{+%.3f} eV^2',FitResult_NoSterile_mNuFree.par(1),-FitResult_NoSterile_mNuFree.errNeg(1),FitResult_NoSterile_mNuFree.errPos(1)),...
            'Location','northwest');
PrettyLegendFormat(leg);
xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
ylabel(sprintf('\\chi^2 (%.0f dof)',dof));
ylim([-0.5 5])
xlim([min(x) max(x)])
