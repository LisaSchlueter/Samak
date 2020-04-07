% KNM2 twin multi ring fit
% March 2020, Lisa

RunList = 'KNM2_Prompt';
E0 = knm2FS_GetE0Twins('SanityPlot','OFF');
range = 40;
pullFlag = 4;
freePar = 'mNu E0 Norm Bkg';
DataType = 'Twin';
RingMerge = 'Full';
ScalingOpt = 2;
MaxSlopeCpsPereV = 5.2*1e-06;
savedir = [getenv('SamakPath'),'knm2ana/knm2_systematics/results/'];
MakeDir(savedir);

CorrCoeff       = 0;%0:0.2:1;%0.9;%(0:0.2:1);
mNuSqErr        = zeros(numel(CorrCoeff)+1,1);
CovMatFracShape = cell(numel(CorrCoeff),1);
CovMatFrac      = cell(numel(CorrCoeff),1);

for i=1:numel(CorrCoeff)
    if ScalingOpt==2
        extraStr = '_Scaling2';
    else
        extraStr = '';
    end
    
    savename = [savedir,sprintf('knm2_MultiRingFit_BkgSys_Constrain%.3gCpsPerEv_%s_%s_%s_pull%.0f_%.0feVrange_RingMerge%s_CorrCoeff%.2f%s.mat',...
        MaxSlopeCpsPereV,DataType, RunList,strrep(freePar,' ',''),pullFlag,range,RingMerge,CorrCoeff(i),extraStr)];
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
        
    % stat only
    savenameStat = [savedir,'knm2_MRStatOnly.mat'];
    if  mNuSqErr(1)~=0
    elseif exist(savenameStat,'file')
        d = importdata(savenameStat);
        mNuSqErr(1) =0.5*(-d.FitResultStat.errNeg(1)+d.FitResultStat.errPos(1)); % d.FitResultStat.err(1);%
        fprintf('load from file %s \n',savenameStat)
    else
        % read data and set up model
        MR = MultiRunAnalysis(RunArg{:});
        MR.exclDataStart = MR.GetexclDataStart(range);
        MR.ModelObj.RFBinStep = 0.02;
        MR.ModelObj.InitializeRF;
        Sigma = repmat(std(E0),3,MR.nRings);
        FSDArg = {'SanityPlot','OFF','Sigma',Sigma};
        MR.ModelObj.LoadFSD(FSDArg{:});
        MR.Fit;
        FitResultStat = MR.FitResult;
        save(savenameStat,'FitResultStat','RunArg','MR','FSDArg','E0');
    end
    
    if exist(savename,'file') 
        d = importdata(savename);
        mNuSqErr(i+1)   = 0.5*(-d.FitResultBkgCM.errNeg(1)+d.FitResultBkgCM.errPos(1)); % d.FitResultBkgCM.err(1);
        CovMatFracShape{i} = d.BkgCovMatFracShape;
        CovMatFrac{i}      = d.BkgCovMatFrac;
        fprintf('load from file %s \n',savename)
    else
        if ~exist('MR','var')
            MR = MultiRunAnalysis(RunArg{:});
            MR.exclDataStart = MR.GetexclDataStart(range);
            MR.ModelObj.RFBinStep = 0.02;
            MR.ModelObj.InitializeRF;
            Sigma = repmat(std(E0),3,MR.nRings);
            FSDArg = {'SanityPlot','OFF','Sigma',Sigma};
            MR.ModelObj.LoadFSD(FSDArg{:});
        end
    
        
        % stat only
        savenameStat = [savedir,'knm2_MRStatOnly.mat'];
        if exist(savenameStat,'file')
            d = importdata(savenameStat);
            mNuSqErr(1) = 0.5*(-d.FitResultStat.errNeg(1)+d.FitResultStat.errPos(1));
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
            'MaxSlopeCpsPereV',MaxSlopeCpsPereV,'BkgRingCorrCoeff',BkgRingCorrCoeff,...
            'BkgScalingOpt',ScalingOpt);
        BkgCovMatFrac      = MR.FitCM_Obj.CovMatFrac;
        BkgCovMatFracShape = MR.FitCM_Obj.CovMatFracShape;
        BkgCovMat          = MR.FitCM_Obj.CovMat;
        
        MR.NonPoissonScaleFactor = 1;
        MR.SetNPfactor; % convert to right dimension (if multiring)
        MR.ComputeCM('SysEffects',struct('FSD','OFF'),'BkgCM','ON',...
            'MaxSlopeCpsPereV',MaxSlopeCpsPereV,'BkgRingCorrCoeff',BkgRingCorrCoeff,...
            'BkgScalingOpt',ScalingOpt);
        
        MR.Fit;
        FitResultBkgCM = MR.FitResult;
        save(savename,'FitResultBkgCM','RunArg','MR','FSDArg','E0','BkgRingCorrCoeff',...
            'BkgCovMatFracShape','BkgCovMatFrac','BkgCovMat');
    end
    
    
%     mNuStat = 0.5*(-FitResultStat.errNeg(1)+FitResultStat.errPos(1));
%     mNuCM   = 0.5*(-FitResultBkgCM.errNeg(1)+FitResultBkgCM.errPos(1));
%     
%     mNuSys =  sqrt(mNuCM^2-mNuStat^2);
%     fprintf('mnuSq sensitivity stat only        = %.3f eV^2 \n',mNuStat);
%     fprintf('mnuSq sensitivity stat + syst only = %.3f eV^2 \n',mNuCM);
%     fprintf('mnuSq sensitivity syst only = %.3f eV^2 \n',mNuSys);
    
end
%%
%plot(CorrCoeff,sqrt(mNuSqErr(2:end).^2-mNuSqErr(1).^2));
sqrt(mNuSqErr.^2-mNuSqErr(1)^2)

f1 = figure('Units','normalized','Position',[0.1,0.1,0.6,0.5]);
for i=1:numel(CorrCoeff)
    subplot(ceil(numel(CorrCoeff)/3),3,i);
    imagesc(CovMatFrac{i});
    pbaspect([1 1 1])
    colorbar
    title(sprintf('\\rho = %.1f',CorrCoeff(i)),'FontWeight','normal','FontSize',get(gca,'FontSize'));
    ax = gca;
    mypos = ax.Position;
    
    PrettyFigureFormat;
    set(gca,'XMinorTick','off');
    set(gca,'YMinorTick','off');
    set(gca,'LineWidth',0.5);
    ax.Position = [mypos(1) mypos(2)-0.03 mypos(3:4)+0.03];
    xticks([]);
    yticks([]);
end
%%
f2 = figure('Units','normalized','Position',[0.1,0.1,0.6,0.5]);
for i=1:numel(CorrCoeff)
    subplot(ceil(numel(CorrCoeff)/3),3,i);
    imagesc(corrcoef(CovMatFrac{i}));
    pbaspect([1 1 1])
    colorbar
    title(sprintf('\\rho = %.1f',CorrCoeff(i)),'FontWeight','normal','FontSize',get(gca,'FontSize'));
    ax = gca;
    mypos = ax.Position;
    
    PrettyFigureFormat;
    set(gca,'XMinorTick','off');
    set(gca,'YMinorTick','off');
    set(gca,'LineWidth',0.5);
    ax.Position = [mypos(1) mypos(2)-0.03 mypos(3:4)+0.03];
    colormap(parula)
    xticks([]);
    yticks([]);
end