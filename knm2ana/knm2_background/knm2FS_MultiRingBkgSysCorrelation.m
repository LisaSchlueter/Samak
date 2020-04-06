% KNM2 twin multi ring fit
% March 2020, Lisa

RunList = 'KNM2_Prompt';
E0 = knm2FS_GetE0Twins('SanityPlot','OFF');
range = 40;
pullFlag = 4;
freePar = 'mNu E0 Norm Bkg';
DataType = 'Twin';
RingMerge = 'Full';

MaxSlopeCpsPereV = 5.2*-06;
savedir = [getenv('SamakPath'),'knm2ana/knm2_systematics/results/'];
MakeDir(savedir);

CorrCoeff = (0:0.2:1);

for i=1:numel(CorrCoeff)
    savename = [savedir,sprintf('knm2_MultiRingFit_BkgSys_Constrain%.3gCpsPerEv_%s_%s_%s_pull%.0f_%.0feVrange_RingMerge%s_CorrCoeff%.2f.mat',...
        MaxSlopeCpsPereV,DataType, RunList,strrep(freePar,' ',''),pullFlag,range,RingMerge,CorrCoeff(i))];
    if exist(savename,'file')
        load(savename)
    else
        % read data and set up model
        RunArg = {'RunList',RunList,...
            'chi2','chi2Stat',...
            'DataType',DataType,...
            'fixPar',freePar,...
            'RadiativeFlag','ON',...
            'minuitOpt','min ; minos',...
            'FSDFlag','BlindingKNM2',...
            'ELossFlag','KatrinT2',...
            'SysBudget',33,...
            'AnaFlag','Ring',...
            'RingMerge',RingMerge,...
            'chi2','chi2Stat',...
            'pullFlag',pullFlag,...
            'TwinBias_Q',E0,...
            'ROIFlag','14keV',...
            'MosCorrFlag','OFF',...
            'NonPoissonScaleFactor',1};
        
        MR = MultiRunAnalysis(RunArg{:});
        MR.exclDataStart = MR.GetexclDataStart(range);
        MR.ModelObj.RFBinStep = 0.02;
        MR.ModelObj.InitializeRF;
        Sigma = repmat(std(E0),3,MR.nRings);
        FSDArg = {'SanityPlot','OFF','Sigma',Sigma};
        MR.ModelObj.LoadFSD(FSDArg{:});
        
        % stat only
        savenameStat = strrep(savename,'.mat','_StatOnly');
        if exist(savenameStat,'file')
            load(savenameStat)
        else  
            MR.Fit;
            FitResultStat = MR.FitResult;
            save(savenameStat,'FitResultStat','RunArg','MR','FSDArg','E0');
        end
        % compute CM
        MR.chi2 = 'chi2CMShape';
        MR.NonPoissonScaleFactor = 1.112;
        MR.SetNPfactor; % convert to right dimension (if multiring)
        BkgRingCorrCoeff = CorrCoeff(i);
        MR.ComputeCM('SysEffects',struct('FSD','OFF'),'BkgCM','ON',...
            'MaxSlopeCpsPereV',MaxSlopeCpsPereV,'BkgRingCorrCoeff',BkgRingCorrCoeff);
        BkgCovMatFrac      = MR.FitCM_Obj.CovMatFrac;
        BkgCovMatFracShape = MR.FitCM_Obj.CovMatFracShape;
        BkgCovMat          = MR.FitCM_Obj.CovMat;
        
        MR.NonPoissonScaleFactor = 1;
        MR.SetNPfactor; % convert to right dimension (if multiring)
        
        MR.Fit;
        FitResultBkgCM = MR.FitResult;
        save(savename,'FitResultBkgCM','FitResultStat','RunArg','MR','FSDArg','E0','BkgRingCorrCoeff',...
            'BkgCovMatFracShape','BkgCovMatFrac','BkgCovMat');
    end
    
    mNuStat = 0.5*(-FitResultStat.errNeg(1)+FitResultStat.errPos(1));
    mNuCM   = 0.5*(-FitResultBkgCM.errNeg(1)+FitResultBkgCM.errPos(1));
    mNuSys =  sqrt(mNuCM^2-mNuStat^2);
    fprintf('mnuSq sensitivity stat only        = %.3f eV^2 \n',mNuStat);
    fprintf('mnuSq sensitivity stat + syst only = %.3f eV^2 \n',mNuCM);
    fprintf('mnuSq sensitivity syst only = %.3f eV^2 \n',mNuSys);
    
end
%%



