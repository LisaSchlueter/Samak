%
% KNM1 Final fit Results
% Uniform Fit
% Golden Run List
% Golden Pixel List
% Last Updated: 23/03/2021

%% settings
RecomputeFlag         = 'OFF';
Plots                 = 'OFF';
DataType              = 'Real';
range                 = 40;
chi2                  = 'chi2CMShape';
freePar               = 'mNu E0 Norm Bkg';
FSDFlag               = 'KNM2_0p1eV';%'Sibille0p5eV';%
ELossFlag             = 'KatrinT2A20';%A20';
AngTF                 = 'ON';
SysBudget             = 200;%200;
BKG_PtSlope           = -2.2*1e-06;
BkgCmMode             = '';%'Gauss';

if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.064;
end

% label
savedir = [getenv('SamakPath'),'knm12Combi/knm1ReAnalysis/results/'];
savefile = sprintf('%sknm1_UniformFit_ReAna_%s_%s_NP%.4g_%s_%.0feV_%s_%s_AngTF%s.mat',...
    savedir,DataType,chi2,NonPoissonScaleFactor,strrep(freePar,' ',''),range,FSDFlag,ELossFlag,AngTF);

if strcmp(chi2,'chi2CMShape')
    savefile = strrep(savefile,chi2,sprintf('%s_SysBudget%.0f',chi2,SysBudget));
    if strcmp(BkgCmMode,'Gauss')
        savefile = strrep(savefile,'.mat','_BkgCmGauss.mat');
    end
end
if BKG_PtSlope~=0
    savefile = strrep(savefile,'.mat',sprintf('_BkgPtSlope%.3gmuCpsS.mat',1e6.*BKG_PtSlope));
end


if exist(savefile,'file') && strcmp(RecomputeFlag,'OFF')
    load(savefile);
    fprintf('load from file %s \n',savefile)
else
    
    if contains(FSDFlag,'Sibille')
        DopplerEffectFlag = 'OFF'; % already in FSD
    else
        DopplerEffectFlag = 'FSD'; % need to be added in FSD
    end
    
    %% Init Model Object and covariance matrix object
    Real = MultiRunAnalysis('RunList','KNM1',...
        'chi2',chi2,...
        'DataType',DataType,...
        'fixPar',freePar,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'SysBudget',SysBudget,...
        'minuitOpt','min ; minos',...
        'FSDFlag',FSDFlag,...
        'ELossFlag',ELossFlag,...
        'AngularTFFlag',AngTF,...
        'SynchrotronFlag','ON',...
        'RadiativeFlag','ON',...
        'DopplerEffectFlag',DopplerEffectFlag,...
        'BKG_PtSlope',BKG_PtSlope);
    
    Real.exclDataStart =    Real.GetexclDataStart(range);
    
    if strcmp(chi2,'chi2Stat')
        Real.InitModelObj_Norm_BKG('RecomputeFlag','ON');
    else
        Real.InitModelObj_Norm_BKG('RecomputeFlag','ON');
        if strcmp(BkgCmMode,'Gauss')
            Real.ComputeCM('BkgMode','Gauss');
        else
            Real.ComputeCM;
        end
    end
    
    Real.Fit;
    
    FitResult = Real.FitResult;
    MakeDir(savedir);
    save(savefile,'FitResult','Real');
end
%%
fprintf('m^2 = %.3f (%.3f +%.3f) eV^2 \n',FitResult.par(1),FitResult.errNeg(1),FitResult.errPos(1))
fprintf('E0 = %.2f  +- %.2f eV \n',FitResult.par(2)+Real.ModelObj.Q_i,FitResult.err(2))
fprintf('B  = %.1f (%.1f +%.1f) mcps \n',(FitResult.par(3)+Real.ModelObj.BKG_RateSec_i)*1e3,1e3*FitResult.errNeg(3),1e3*FitResult.errPos(3))
fprintf('chi^2 = %.1f (%.0f dof) \n',FitResult.chi2min,FitResult.dof)
if strcmp(Plots,'ON')
    Real.PlotFit('LabelFlag','data','saveplot','pdf','ErrorBarScaling',1,'YLimRes',[-2.2,2.9],'Colors','RGB','DisplayStyle','PRL');
end

