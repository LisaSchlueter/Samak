% sanity plots for stacking uncertainty
% different Mode to choose from: 
% 1. Effect of Stacking on neutrino mass using twin runs
% 2. show variation of qU und qUfrac in covariance matrix

Mode = 2;

%% 1. Effect on neutrino mass
if Mode==1
    if ~exist('M','var')
        M = MultiRunAnalysis('RunList','KNM1','chi2','chi2Stat','fixPar','5 6 7 8 9 10 11',...
            'DataType','Twin','FSDFlag','Sibille0p5eV',...
            'minuitOpt','min;minos');
    end
    
    %stat only
    exclDataStart = [2,14,17];
    mNuSq = zeros(numel(exclDataStart),1);  mNuSqErr = zeros(numel(exclDataStart),1);
    E0    = zeros(numel(exclDataStart),1);  E0Err = zeros(numel(exclDataStart),1);
    M.chi2 = 'chi2Stat';
    for i=1:numel(exclDataStart)
        M.exclDataStart=exclDataStart(i);
        M.Fit;
        mNuSq(i)    = M.FitResult.par(1);   mNuSqErr(i) = M.FitResult.err(1);
        E0(i)       = M.FitResult.par(2);   E0Err(i) = M.FitResult.err(2);
    end
    
    % stat + stacking CM
    M.chi2 = 'chi2CMShape';
    M.ComputeCM('SysEffects',struct('Stack','ON'),'BkgCM','OFF')
    mNuSqCM = zeros(numel(exclDataStart),1);  mNuSqErrCM = zeros(numel(exclDataStart),1);
    E0CM    = zeros(numel(exclDataStart),1);  E0ErrCM = zeros(numel(exclDataStart),1);
    for i=1:numel(exclDataStart)
        M.exclDataStart=exclDataStart(i);
        M.Fit;
        mNuSqCM(i)    = M.FitResult.par(1);   mNuSqErrCM(i) = M.FitResult.err(1);
        E0CM(i)       = M.FitResult.par(2);   E0ErrCM(i) = M.FitResult.err(2);
    end
end
%%
if Mode==2
    if ~exist('M','var')
        M = MultiRunAnalysis('RunList','KNM1','chi2','chi2Stat','fixPar','5 6 7 8 9 10 11',...
            'DataType','Twin','FSDFlag','Sibille0p5eV',...
            'minuitOpt','min;migrad');
    end
    
% %     M.Get_DataDriven_RelErr_qU('Debug','ON')
 %     M.ComputeCM('SysEffects',struct('Stack','ON'),'BkgCM','OFF');
    M.FitCM_Obj.RecomputeFlag ='ON';
    M.FitCM_Obj.nTrials = 1000;
      M.FitCM_Obj.ComputeCM_Stacking('SanityPlot','ON');
%      M.FitCM_Obj.PlotCM('savePlot','ON');

c = corplot(M.FitCM_Obj.CovMatFracShape(2:end,2:end));
PrettyFigureFormat
pbaspect([1 1 1])
set(gca,'xtick',[5,35])
set(gca,'ytick',[])
qUmin = sprintf('qU_{min} = E_0-%.0fV',abs(M.FitCM_Obj.StudyObject.qU(2)-M.FitCM_Obj.StudyObject.Q_i));
qUmax = sprintf('qU_{max} = E_0-%.0fV',abs(M.FitCM_Obj.StudyObject.qU(end)-M.FitCM_Obj.StudyObject.Q_i));

set(gca,'xticklabel',{qUmin,qUmax}); set(gca,'yticklabel',[]);
set(gca,'FontSize',16)
print(gcf,'./plots/StackCM_CorrPlot.png','-dpng','-r450');
end

if Mode==3
    M.FitCM_Obj.ComputeCM_Stacking;
    M.FitCM_Obj.PlotCM;
   
end
