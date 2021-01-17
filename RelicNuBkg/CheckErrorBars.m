matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Chi2Scans/')];
savename=[matFilePath,sprintf('RelicChi2Scan_Twin_BiasmnuSq%g_Syst%s_range%g_%s_[0 %g]_%s.mat',1,'ON',40,'KNM1',5e11,'mNu E0 Norm Bkg')];
load(savename);

F = MultiRunAnalysis('RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
     'chi2','chi2CMShape',...                 % uncertainties: statistical or stat + systematic uncertainties
     'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
     'fixPar','mNu E0 Norm Bkg',...                   % free Parameter!!
     'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
     'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
     'fitter','minuit',...
     'minuitOpt','min ; minos',...         % technical fitting options (minuit)
     'FSDFlag','SibilleFull',...          % final state distribution
     'ELossFlag','KatrinT2',...            % energy loss function
     'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
     'DopplerEffectFlag','FSD',...
     'Twin_SameCDFlag','OFF',...
     'Twin_SameIsotopFlag','OFF',...
     'SynchrotronFlag','ON',...
     'AngularTFFlag','OFF',...
     'TwinBias_Q',18573.73,...
     'TwinBias_mnuSq',1);

F.exclDataStart = F.GetexclDataStart(40);

steps = 11;
 
eta=(0:(Netabins-1))*((etafactor*10^(etarange))/(Netabins-1));
ProblemPoints = [eta(10)];%[eta(2) eta(5) eta(8) eta(10)];
ProblemMasses = [mnuSq(10)];%[mnuSq(2) mnuSq(5) mnuSq(8) mnuSq(10)];
ScanRange = [linspace(-0.8,0.8,steps); linspace(-mnuSq_err(5),mnuSq_err(5),steps); linspace(-mnuSq_err(8),mnuSq_err(8),steps); linspace(-mnuSq_err(10),mnuSq_err(10),steps)];
Chi2mnu = ones(numel(ProblemMasses),numel(ScanRange(1,:)));

for i=1:numel(ProblemPoints)
    F.ModelObj.eta_i = ProblemPoints(i);
    for j=1:numel(ScanRange(1,:))
        F.ModelObj.mnuSq_i = ProblemMasses(i)+ScanRange(i,j);
        F.ModelObj.ComputeNormFactorTBDDS;
        F.ModelObj.ComputeTBDDS;
        F.ModelObj.ComputeTBDIS;
        F.Fit;
        Chi2mnu(i,j) = F.FitResult.chi2min;
    end
end

hold on;
for i=1:numel(ScanRange(:,1))
    plot(ScanRange(i,:),Chi2mnu(i,:));
end
hold off;