NScan=30;
fitter='minuit';
savename = './RelicNuBkg/Misc/KRN2ExhaustiveSystematics_matlabCorrected.mat';
savename2= './RelicNuBkg/Misc/KRN2narrowScan.mat';
if exist(savename,'file')
    load(savename);

A = MultiRunAnalysis('RunList','KNM2_Prompt',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                            'chi2','chi2CMShape',...                 % uncertainties: statistical or stat + systematic uncertainties
                            'DataType','Real',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                            'fixPar','mNu E0 Norm Bkg eta',...                   % free Parameter!!
                            'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                            'NonPoissonScaleFactor',1.112,...     % background uncertainty are enhanced
                            'fitter',fitter,...
                            'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                            'FSDFlag','KNM2',...          % final state distribution
                            'ELossFlag','KatrinT2A20',...            % energy loss function
                            'SysBudget',40,...                    % defines syst. uncertainties -> in GetSysErr.m;
                            'DopplerEffectFlag','FSD',...
                            'Twin_SameCDFlag','OFF',...
                            'Twin_SameIsotopFlag','OFF',...
                            'SynchrotronFlag','ON',...
                            'AngularTFFlag','ON',...
                            'TwinBias_Q',18573.7,...
                            'TwinBias_mnuSq',0,...
                            'FSD_Sigma',sqrt(0.0124+0.0025),...
                            'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
                            'BKG_PtSlope',3*1e-06,...
                            'TwinBias_BKG_PtSlope',3*1e-06);
B = MultiRunAnalysis('RunList','KNM2_Prompt',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                            'chi2','chi2CMShape',...                 % uncertainties: statistical or stat + systematic uncertainties
                            'DataType','Real',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                            'fixPar','mNu E0 Norm Bkg',...                   % free Parameter!!
                            'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                            'NonPoissonScaleFactor',1.112,...     % background uncertainty are enhanced
                            'fitter',fitter,...
                            'minuitOpt','min ; imp',...         % technical fitting options (minuit)
                            'FSDFlag','KNM2',...          % final state distribution
                            'ELossFlag','KatrinT2A20',...            % energy loss function
                            'SysBudget',40,...                    % defines syst. uncertainties -> in GetSysErr.m;
                            'DopplerEffectFlag','FSD',...
                            'Twin_SameCDFlag','OFF',...
                            'Twin_SameIsotopFlag','OFF',...
                            'SynchrotronFlag','ON',...
                            'AngularTFFlag','ON',...
                            'TwinBias_Q',18573.7,...
                            'TwinBias_mnuSq',0,...
                            'FSD_Sigma',sqrt(0.0124+0.0025),...
                            'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
                            'BKG_PtSlope',3*1e-06,...
                            'TwinBias_BKG_PtSlope',3*1e-06);
A.exclDataStart = A.GetexclDataStart(40);
B.exclDataStart = B.GetexclDataStart(40);
%A.exclDataStop  = A.ModelObj.nqU-1;
%B.exclDataStop  = B.ModelObj.nqU-1;
else                      
SysEffectsList = categorical({'Total','Stat','FSD','RF','TASR','Stack','FPD','TC','Bkg','Plasma','PT','NP'});
SysEffectsList = reordercats(SysEffectsList,{'Total','Stat','FSD','RF','TASR','Stack','FPD','TC','Bkg','Plasma','PT','NP'});
SysEffects     = categories(SysEffectsList);
Chi2Profiles   = zeros(NScan,numel(SysEffectsList));
ErrorBand      = linspace(-9e10,9e10,NScan);

for i=1:numel(SysEffectsList)
    if ~(strcmp(SysEffects{i},'Total') || strcmp(SysEffects{i},'Bkg') || strcmp(SysEffects{i},'NP') || strcmp(SysEffects{i},'PT') || strcmp(SysEffects{i},'Stat'))
        A.ComputeCM('SysEffects',struct(SysEffects{i},'ON'),'BkgCM','OFF','BkgPtCM','OFF');
        B.ComputeCM('SysEffects',struct(SysEffects{i},'ON'),'BkgCM','OFF','BkgPtCM','OFF');
    elseif strcmp(SysEffects{i},'Bkg')
        A.ComputeCM('SysEffects',struct(),'BkgCM','ON','BkgPtCM','OFF');
        B.ComputeCM('SysEffects',struct(),'BkgCM','ON','BkgPtCM','OFF');
    elseif strcmp(SysEffects{i},'NP')
        A.NonPoissonScaleFactor=1.112;
        B.NonPoissonScaleFactor=1.112;
        A.ComputeCM('SysEffects',struct(),'BkgCM','OFF','BkgPtCM','OFF');
        B.ComputeCM('SysEffects',struct(),'BkgCM','OFF','BkgPtCM','OFF');
    elseif strcmp(SysEffects{i},'PT')
        A.ComputeCM('SysEffects',struct(),'BkgCM','OFF','BkgPtCM','ON');
        B.ComputeCM('SysEffects',struct(),'BkgCM','OFF','BkgPtCM','ON');
    elseif strcmp(SysEffects{i},'Stat')
        A.NonPoissonScaleFactor=1;
        B.NonPoissonScaleFactor=1;
        A.ComputeCM('SysEffects',struct(),'BkgCM','OFF','BkgPtCM','OFF');
        B.ComputeCM('SysEffects',struct(),'BkgCM','OFF','BkgPtCM','OFF');
    end
    A.Fit;
    ResultArray(i)=A.FitResult;
    for j=1:NScan
        B.ModelObj.eta_i=ResultArray(1).par(18)*1e10+ErrorBand(j);
        B.ModelObj.ComputeNormFactorTBDDS;
        B.ModelObj.ComputeTBDDS;
        B.ModelObj.ComputeTBDIS;
        B.Fit;
        Chi2Profiles(j,i)=B.FitResult.chi2min;
    end
end
save(savename,'ResultArray','SysEffectsList','Chi2Profiles','ErrorBand');
end
EtaVec = ResultArray(1).par(18)*1e10+ErrorBand;
Error  = zeros(1,numel(SysEffectsList));
Syst   = zeros(1,numel(SysEffectsList)-2);
for i=1:numel(SysEffectsList)
   etaBest    = EtaVec(Chi2Profiles(:,i)==min(Chi2Profiles(:,i)));
   HighChi2   = EtaVec(Chi2Profiles(:,i)>min(Chi2Profiles(:,i))+1);
   LowChi2    = find(Chi2Profiles(:,i)<min(Chi2Profiles(:,i))+1);
   etaLowVec  = HighChi2(HighChi2<etaBest);
   etaHighVec = HighChi2(HighChi2>etaBest);
   etaLow     = ((EtaVec(LowChi2(1))-etaLowVec(end))./(Chi2Profiles(LowChi2(1),i)-Chi2Profiles(LowChi2(1)-1,i))).*...
       (1+min(Chi2Profiles(:,i)))-Chi2Profiles(LowChi2(1),i).*((EtaVec(LowChi2(1))-etaLowVec(end))./...
       (Chi2Profiles(LowChi2(1),i)-Chi2Profiles(LowChi2(1)-1,i)))+EtaVec(LowChi2(1));
   etaHigh     = ((etaHighVec(1)-EtaVec(LowChi2(end)))./(Chi2Profiles(LowChi2(end)+1,i)-Chi2Profiles(LowChi2(end),i))).*...
       (1+min(Chi2Profiles(:,i)))-Chi2Profiles(LowChi2(end)+1,i).*((etaHighVec(1)-EtaVec(LowChi2(end)))./...
       (Chi2Profiles(LowChi2(end)+1,i)-Chi2Profiles(LowChi2(end),i)))+etaHighVec(1);
   Error(i)    = (etaHigh-etaLow)/2;
end
for i=1:numel(SysEffectsList)-2
    Syst(i)=sqrt(Error(i+2).^2-Error(2).^2);
end
bar(SysEffectsList,[Error(1) Error(2) Syst]);

SaveDir='./RelicNuBkg/Plots/Animation/';
ErrorBar_narrow = linspace(-9.4e10,-5.4e10,NScan);%linspace(ResultArray(1).par(18)*1e10,-5.4e10,NScan);
%B.ModelObj.eta_i=0;
B.Fit;
times = B.ModelObj.qUfrac*B.ModelObj.TimeSec;
qU    = B.ModelObj.qU; qU    = qU-B.TwinBias_Q;
qULimit = -40;
% Spectrum no relics 
IS  = B.ModelObj.TBDIS; 
YIs = IS./times;
DIS = B.RunData.TBDIS;
DIS = DIS./times;
% YIs=YIs(qU>qULimit);
% DIS=DIS(qU>qULimit);
% times=times(qU>qULimit);
% qU=qU(qU>qULimit);
for i=1:NScan
    B.ModelObj.eta_i=ErrorBar_narrow(i);
    B.ModelObj.ComputeNormFactorTBDDS;
    B.ModelObj.ComputeTBDDS;
    B.ModelObj.ComputeTBDIS;
    B.Fit;
    ResultArray_narrow(i)=B.FitResult;
    % Spectrum relics
    YI = B.ModelObj.TBDIS; 
    %YI=YI(qU>qULimit);
    YI = YI./times;
    
    % Error bar
    err  = (diag(sqrt(B.FitCMShape)));
    %err=err(qU>qULimit);
    err  = err./times;
    err  = err./YI;


    RSP  = (YI./YIs);
    RSPd = RSP;
    LocalFontSize = 20;
    prlG = [81 126 102]/255;
    prlB = [50 148 216]/255;
    FitStyleArg = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',4,'MarkerEdgeColor',rgb('Black')};
    fig_i=figure('Renderer','Painters');
    set(fig_i, 'Units', 'normalized', 'Position', [0.001, 0.001,0.45, 0.45]);
    % Plot
    hr1 = plot(qU,RSP,'color',rgb('Salmon'),'LineWidth',3,'LineStyle','-');
    hold on;
    hr2 = plot(qU,ones(1,numel(qU)),'color',prlB,'LineWidth',3,'LineStyle',':');
    hr3 = errorbar(qU,(DIS./YIs),err,FitStyleArg{:},'CapSize',0);
    yl2=ylabel('Ratio');
    katrinsim   = sprintf('\\eta=0');
    sterilemod  = sprintf('Best fit: \\eta=%.2g',B.ModelObj.eta);
    hl=legend([hr2 hr1],{katrinsim,sterilemod},'Location','best','box','off');
    hl.NumColumns=2;
    hl.FontSize = LocalFontSize-2;

    xlim([-40 7]);
    ylim([0.983 1.012]);

    PRLFormat;
    set(gca,'FontSize',LocalFontSize);
    set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
    set(get(gca,'YLabel'),'FontSize',LocalFontSize+4); 
    hl.Position(2) = 0.333;
    ax2 = gca;
    SaveName = sprintf('LocalFit_eta%.2g_%i.pdf',ErrorBar_narrow(i),i);
    export_fig(fig_i, [SaveDir,SaveName]);
end
save(savename2,'ResultArray_narrow');