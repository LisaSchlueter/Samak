
InitFile = @ref_FakeRun_KNM2_CD84_308x2hours;
range = 40;

% save label
savedir = [getenv('SamakPath'),'knm2ana/knm2_systematics/results/'];
MakeDir(savedir);

PixList = [1:28,29:29+27,65:65+27,101:101+27];
RunAnaArg = {'RunNr',1,...
    'DataType','Fake',...
    'FakeInitFile',InitFile,...
    'fixPar','mNu E0 Bkg Norm',....%free Par
    'NonPoissonScaleFactor',1,...
    'AnaFlag','Ring',...
    'RingMerge','Full',...
    'PixList',PixList,...
    'ELossFlag','KatrinT2',...
    'ISCSFlag','Edep'};

MR = RunAnalysis(RunAnaArg{:});
MR.exclDataStart= MR.GetexclDataStart(range);

BkgRingCorrCoeff = [1,0];
ScalingOpt       = [1,2];
Mode = 'SlopeFit';
MaxSlopeCpsPereV = 5.2*1e-06;

% init
mNuSqErr        = zeros(numel(BkgRingCorrCoeff)+1,1);
CovMatFracShape = cell(numel(BkgRingCorrCoeff),1);
CovMatFrac      = cell(numel(BkgRingCorrCoeff),1);
CovMat          = cell(numel(BkgRingCorrCoeff),1);
FitResultBkgCM = cell(numel(BkgRingCorrCoeff),1);

% stat. only
savenameStat = sprintf('%sknm2_MRStatOnly_FakeRun%s.mat',savedir,func2str(InitFile));
if exist(savenameStat,'file')
    d = importdata(savenameStat);
    mNuSqErr(1) = d.mNuSqErr(1);
else
    MR.Fit;
    FitResultStat = MR.FitResult;
    mNuSqErr(1) = 0.5*(abs(FitResultStat.errNeg)+FitResultStat.errPos);
    save(savenameStat,'mNuSqErr','FitResultStat','RunAnaArg');
end

% stat + syst.
MR.chi2 = 'chi2CMShape';

for i=1:numel(BkgRingCorrCoeff)
    savename = sprintf('%sMultiRingBkgSys_%s_%.0fCorr_%.2fMaxSlope_%.0fScaling_%s.mat',...
        savedir,func2str(InitFile),BkgRingCorrCoeff(i),MaxSlopeCpsPereV*1e6,ScalingOpt(i),Mode);
    if exist(savename,'file')
        d = importdata(savename);
        FitResultBkgCM{i}  = d.FitResultCM;
        CovMatFracShape{i} = d.BkgCovMatFracShape;
        CovMatFrac{i}      = d.BkgCovMatFrac;
        CovMat{i}          = d.BkgCovMat;
        mNuSqErr(i+1)      = 0.5*(abs(d.FitResultCM.errNeg(1))+d.FitResultCM.errPos(1));
    else
        MR.ComputeCM('SysEffect',struct('FSD','OFF'),'MaxSlopeCpsPereV',MaxSlopeCpsPereV,...
            'BkgRingCorrCoeff',BkgRingCorrCoeff(i),...
            'BkgScalingOpt',ScalingOpt(i),'BkgMode',Mode,...
            'PlotSaveCM','ON');
        MR.Fit;
        FitResultCM        =  MR.FitResult;
        BkgCovMatFrac      = MR.FitCM_Obj.CovMatFrac;
        BkgCovMatFracShape = MR.FitCM_Obj.CovMatFracShape;
        BkgCovMat          = MR.FitCM_Obj.CovMat;
        BkgCovMatFile      = MR.FitCM_Obj.CovMatFile;
        
        save(savename,'mNuSqErr','FitResultCM','BkgCovMatFracShape','BkgCovMatFrac','BkgCovMat','BkgCovMatFile');
        
        FitResultBkgCM{i} = MR.FitResult;
        CovMatFracShape{i} = BkgCovMatFracShape;
        CovMatFrac{i}      = BkgCovMatFrac;
        CovMat{i}          = BkgCovMat;
        mNuSqErr(i+1)      = 0.5*(abs(MR.FitResult.errNeg(1))+MR.FitResult.errPos(1));
    end
end


