% KNM2 twin multi ring fit
% March 2020, Lisa

RunList   = 'KNM2_Prompt';
E0        = 18573.7;%knm2FS_GetE0Twins('SanityPlot','OFF');
range     = 40;
freePar   = 'mNu E0 Norm Bkg';
DataType  = 'Twin';
RingMerge = 'Full';
CorrCoeff = 0; %1
MaxSlopeCpsPereV = [4.6, 5.4, 5.5, 4.2].*1e-06; % [2.2, 3.0, 3.1, 1.9].*1e-06:
Mode      = 'Gauss'; 
RecomputeCMFlag = 'ON';
RecomputeFlag = 'ON';

savedir = [getenv('SamakPath'),'knm2ana/knm2_systematics/results/'];
MakeDir(savedir);

savename = [savedir,sprintf('knm2_MultiRingFitBkgSys_Constraint%.3gmcpsPerKeV_nc%.0f_%s_%s_%s_%.0feVrange_RingMerge%s_%.0fCorrCoeff_%s.mat',...
    norm(MaxSlopeCpsPereV).*1e6,numel(MaxSlopeCpsPereV),DataType,RunList,strrep(freePar,' ',''),range,RingMerge,CorrCoeff,Mode)];
if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
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
        'AnaFlag','Ring',...
        'RingMerge',RingMerge,...
        'chi2','chi2Stat',...
        'TwinBias_Q',E0,...
        'ROIFlag','Default',...
        'MosCorrFlag','OFF',...
        'NonPoissonScaleFactor',1};
    
    MR = MultiRunAnalysis(RunArg{:});
    MR.exclDataStart = MR.GetexclDataStart(range);
    MR.ModelObj.RFBinStep = 0.02;
    MR.ModelObj.InitializeRF;
    if numel(E0)>1
        Sigma = repmat(std(E0),3,MR.nRings);
        FSDArg = {'SanityPlot','OFF','Sigma',Sigma};
        MR.ModelObj.LoadFSD(FSDArg{:});
    end
    MR.Fit;
    FitResultStat = MR.FitResult;
    
    if numel(MaxSlopeCpsPereV)>1
        ScalingOpt = 99; % no scaling
    elseif CorrCoeff==1
        ScalingOpt = 1;  % scale with statistics
    elseif CorrCoeff==0
        ScalingOpt = 2;  % scale with stat. uncertainties
    end
    
    % compute CM
    MR.chi2 = 'chi2CMShape';
    CMArg = {'SysEffects',struct('FSD','OFF'),'BkgCM','ON',...
        'MaxSlopeCpsPereV',MaxSlopeCpsPereV,'BkgRingCorrCoeff',CorrCoeff,...
        'BkgScalingOpt',ScalingOpt,'BkgMode',Mode};
    MR.ComputeCM(CMArg{:},'RecomputeFlag',RecomputeCMFlag);
    
    MR.Fit;
    FitResultCM = MR.FitResult;
    
    CovMatFrac      = MR.FitCM_Obj.CovMatFrac;
    CovMatFracShape = MR.FitCM_Obj.CovMatFracShape;
    CovMat          = MR.FitCM_Obj.CovMat;
    CovMatFile      = MR.FitCM_Obj.CovMatFile;
    d = importdata(CovMatFile);
    Slopes      = d.Slopes;
    SlopesExcl  = d.SlopesExcl;
    save(savename,'FitResultCM','FitResultStat','RunArg',...
                   'MaxSlopeCpsPereV','ScalingOpt',...
                   'CovMatFrac','CovMatFracShape','CovMat',...
                   'Slopes','SlopesExcl','CovMatFile');
end

%%
mNuStat = 0.5*(-FitResultStat.errNeg(1)+FitResultStat.errPos(1));
mNuCM   = 0.5*(-FitResultCM.errNeg(1)+FitResultCM.errPos(1));
mNuSys =  sqrt(mNuCM^2-mNuStat^2);
fprintf('mnuSq sensitivity stat only        = %.3f eV^2 \n',mNuStat);
fprintf('mnuSq sensitivity stat + syst only = %.3f eV^2 \n',mNuCM);
fprintf('mnuSq sensitivity syst only = %.3f eV^2 \n',mNuSys);





