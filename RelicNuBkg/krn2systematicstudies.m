NScan=20;
A = MultiRunAnalysis('RunList','KNM2_Prompt',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                            'chi2','chi2CMShape',...                 % uncertainties: statistical or stat + systematic uncertainties
                            'DataType','Real',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                            'fixPar','mNu E0 Norm Bkg eta',...                   % free Parameter!!
                            'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                            'NonPoissonScaleFactor',1.112,...     % background uncertainty are enhanced
                            'fitter','minuit',...
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
                            'fitter','minuit',...
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
A.exclDataStart = A.GetexclDataStart(40);
B.exclDataStart = B.GetexclDataStart(40);
                        
SysEffectsList = categorical({'Total','Stat','FSD','RF','TASR','Stack','FPD','TC','Bkg','Plasma','PT','NP'});
SysEffectsList = reordercats(SysEffectsList,{'Total','Stat','FSD','RF','TASR','Stack','FPD','TC','Bkg','Plasma','PT','NP'});
SysEffects     = categories(SysEffectsList);
Chi2Profiles   = zeros(NScan,numel(SysEffectsList));
ErrorBand      = linspace(-6e10,6e10,NScan);

for i=1:numel(SysEffectsList)
    if ~(strcmp(SysEffects{i},'Total') || strcmp(SysEffects{i},'Bkg') || strcmp(SysEffects{i},'NP') || strcmp(SysEffects{i},'PT'))
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
%     for j=1:NScan
%         B.ModelObj.eta_i=A.FitResult.par(18)*1e10+ErrorBand(j);
%         B.ModelObj.ComputeNormFactorTBDDS;
%         B.ModelObj.ComputeTBDDS;
%         B.ModelObj.ComputeTBDIS;
%         B.Fit;
%         Chi2Profiles(j,i)=B.FitResult.chi2min;
%     end
end
save('./RelicNuBkg/Misc/KRN2ExhaustiveSystematics.mat','ResultArray','SysEffectsList','Chi2Profiles','ErrorBand');