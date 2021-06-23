range = 40;
DataType = 'Real';
freePar = 'mNu E0 Norm Bkg';
FinalSensitivity = 'ON';
pullFlag = 15;
switch pullFlag
    case 99
        pullStr = '';
    case 15
        pullStr = sprintf('_mNuSqPull1eV2');
end
savedir = [getenv('Samak3.0'),'ksn1ana/ksn1_NuBetaBeta/results/'];
if strcmp(FinalSensitivity,'OFF')
    savefile = sprintf('%sContour_%s_%s_%.0feV%s.mat',savedir,DataType,strrep(freePar,' ',''),range,pullStr);
else
    savefile = sprintf('%sContour_FinalSensitivity.mat',savedir);
end
%%
if exist(savefile,'file')
    load(savefile)
else
    RunAnaArg = {'RunList','KNM1',...
        'fixPar',freePar,...
        'DataType',DataType,...
        'FSDFlag','SibilleFull',...
        'ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'ROIFlag','Default',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF',...
        'ISCSFlag','Edep',...
        'TwinBias_Q',18573.73,...
        'SysBudget',24,...
        'pullFlag',pullFlag,...
        'NonPoissonScaleFactor',1};
    
    T = MultiRunAnalysis(RunAnaArg{:});
    T.chi2 = 'chi2CMShape';
    %% settings sterile class
    SterileArg = {'RunAnaObj',T,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'nGridSteps',50,...
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',...
        'range',range};
    
    S = SterileAnalysis(SterileArg{:});
    CL = 0.95;
    S.DeltaChi2 = GetDeltaChi2(CL,2);
    %%
    if strcmp(FinalSensitivity,'OFF')
        if contains(freePar,'mNu') && pullFlag==99
            S.InterpMode = 'lin';
        end
    elseif strcmp(FinalSensitivity,'ON')
        S.RunAnaObj.DataType = 'Fake';
        S.RunAnaObj.DataSet = 'Knm2';
        FakeInitFile = @ref_KSNX_KATRIN_Final;
        S.RunAnaObj.FakeInitFile = FakeInitFile;
        S.RunAnaObj.fixPar = 'E0 Norm Bkg';
        S.RunAnaObj.InitFitPar;
        S.RunAnaObj.pullFlag = 99;
        S.RunAnaObj.AngularTFFlag = 'ON';
        S.RunAnaObj.ELossFlag = 'KatrinT2A20';
        S.RunAnaObj.SysBudget = 66;
    end
    
    S.LoadGridFile('CheckSmallerN','ON');
    S.Interp1Grid('RecomputeFlag','ON');
    
    [M,pHandle]= contour(S.sin2T4,S.mNu4Sq,S.chi2-S.chi2_ref,...
        [S.DeltaChi2 S.DeltaChi2],'Color',rgb('Red'));
    
    mNu4Sq    = M(2,2:end);
    sinT4Sq = M(1,2:end);
    
    plot(sinT4Sq,mNu4Sq,'LineWidth',2);
    
    S.FindBestFit;
    mNu4Sq_bf = S.mNu4Sq_bf;
    sinT4Sq_bf = S.sin2T4_bf;
    
    MakeDir(savedir);
    save(savefile,'mNu4Sq','sinT4Sq','mNu4Sq_bf','sinT4Sq_bf');
end
plot(sinT4Sq,mNu4Sq,'LineWidth',2,'Color',rgb('Red'));
set(gca,'YScale','log');
set(gca,'XScale','log');
