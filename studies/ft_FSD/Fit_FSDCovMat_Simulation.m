% script to do fits to Simulation with FSD covmat for Stacked Pixel
%for First Tritium
% Plan: 1. simulate spectrum, where sys effect is slightly different
%       2. add stat flutuations -> cholesky
%       % fit with statistics: does effect influence endpoint?
%       % fit with systematics included: is fit improved? 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General settings
myRunList = {'StackCD100all','StackCD100up','StackCD100down'}; 
ringCutFlag = 'ex2';
% Systematics settings for FIT
nTrials = 1000;
RecomputeFlag = 'OFF';
StackCM = 'OFF';
mySysEffects = struct('FSD','ON','RF_EL','OFF','RF_RX','OFF','RF_BF','OFF','TASR','ON'); % sys effects used in fit with CM
%Fit settings
fixPar = '1 5 6';
chi2Flag = 'chi2CM';
% Simulation settings
StatFluct = 'ON'; % add stat fluctuation 
SysShift = 'ON';  % add systematic shift
nuMass = 0;         %Neutrino mass of simulation model
nSim = 10;         %Number of simulated spectra
mySysEffectsSim =  struct('FSD','ON','RF_EL','ON','RF_RX','ON','RF_BF','ON','TASR','ON'); % sys effects in MC data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize RunAnalysis Model
A = MultiRunAnalysis('RunList',myRunList{1},'AnaFlag','StackPixel',...
    'fixPar',fixPar,'pullFlag',3);

% Init other variables
TBDISsim       = zeros(A.ModelObj.nqU,nSim);
FitResults_par = zeros(6,nSim);
FitResults_err = zeros(6,nSim);
Fitchi2min     = zeros(nSim,1);

% Initialize Simulation 
S = ref_RunAnalysis(myRunList{1},ringCutFlag,'','ISCS','Theory',...
                    'recomputeRF','OFF','mnuSq_i',nuMass);              
S.ComputeTBDDS;
S.ComputeTBDIS;
TBDIS_i = S.TBDIS;

% Add Fluctuation: Stat and/or Sys
if strcmp(SysShift,'ON')
    A.ComputeCM('SysEffects',mySysEffectsSim,'nTrials',nTrials,...
        'RecomputeFlag',RecomputeFlag,'InitNormFit','OFF',...
        'StackCM',StackCM);
    CMFrac = A.FitCM_Obj.CovMatFrac;
    TBDIS_i = mvnrnd(S.TBDIS,S.TBDIS.*CMFrac.*S.TBDIS',1);
end
if strcmp(StatFluct,'ON')
    TBDISsim = TBDIS_i' + sqrt(TBDIS_i)'.*randn(S.nqU,nSim);
else
    TBDISsim = repmat(S.TBDIS,1,nSim);
end

% Get Covariance Matrix for Fit (if needed)
A.chi2 = chi2Flag;
if ~strcmp(A.chi2,'chi2Stat')
    % check if mySysEffects and mySysEffectsSim  are the same
    names = fieldnames(mySysEffects);
    tmp =0;
    for i=1:numel(names) 
        if  strcmp(mySysEffects.(names{i}),mySysEffectsSim.(names{i}))
            tmp = tmp +1;
        end
    end
    % if not the same
    if tmp~=numel(names) || strcmp('SysShift','OFF')
        A.ComputeCM('SysEffects',mySysEffects,'nTrials',nTrials,...
            'RecomputeFlag',RecomputeFlag,'InitNormFit','OFF',...
            'StackCM',StackCM);
    end
end

% Do the Fit
A.RunData.qU    = S.qU;
for i=1:nSim
A.RunData.TBDIS = TBDISsim(:,i);
A.Fit;
FitResults_par(:,i) = A.FitResult.par;
FitResults_err(:,i) = A.FitResult.err;
Fitchi2min(i)       = A.FitResult.chi2min;
end
%% Display Results
f100 = figure(100);
set(f100, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.8, 1.2]);
maintitle= 'main title '; 
        a=annotation('textbox', [0 0.91 1 0.1], ...
            'String', maintitle, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
subplot(2,3,1);
histogram(FitResults_par(1,:),'FaceColor',rgb('CadetBlue'));
xlabel('\nu-mass² (eV²)');
ylabel('');
legend(sprintf('mean = %.2f eV² \n std = %.2f eV⁴',mean(FitResults_par(1,:)),std(FitResults_par(1,:))));
PrettyFigureFormat;
set(gca,'FontSize',16);

subplot(2,3,2);
histogram(FitResults_par(2,:),'FaceColor',rgb('CadetBlue'));
xlabel(['E_0 - ',sprintf('%.1f (eV)',A.ModelObj.Q_i)]);
ylabel('');
legend(sprintf('mean = %.2e eV \n std = %.2e eV²',mean(FitResults_par(2,:)),std(FitResults_par(2,:))));
PrettyFigureFormat;
set(gca,'FontSize',16);

subplot(2,3,4);
histogram(FitResults_err(1,:),'FaceColor',rgb('CadetBlue'));
xlabel('\nu-mass² uncertainty (eV²)');
ylabel('');
legend(sprintf('mean = %.2f eV² \n std = %.2f eV⁴',mean(FitResults_err(1,:)),std(FitResults_err(1,:))));
PrettyFigureFormat;
set(gca,'FontSize',16);

subplot(2,3,5);
histogram(FitResults_err(2,:),'FaceColor',rgb('CadetBlue'));
xlabel('E_0 uncertainty');
ylabel('');
legend(sprintf('mean = %.2e eV \n std = %.2e eV²',mean(FitResults_err(2,:)),std(FitResults_err(2,:))));
PrettyFigureFormat;
set(gca,'FontSize',16);

subplot(2,3,3)%chi2 distribution
h1 = histogram(Fitchi2min,'FaceColor',rgb('CadetBlue'),'Normalization','probability');
hold on;
x = min(h1.BinLimits):max(h1.BinLimits);
chi2Dist = chi2pdf(x,A.FitResult.dof);
p1 = plot(x,chi2Dist*h1.BinWidth,'LineWidth',3,'color',rgb('DarkRed'));
xlabel(['\chi²',sprintf('/ %u dof',A.FitResult.dof)]);
ylabel('');
legend([h1,p1],sprintf('mean = %.2f \n std = %.2f',mean(Fitchi2min),std(Fitchi2min)),'chi2 pdf');
PrettyFigureFormat;
set(gca,'FontSize',16);

subplot(2,3,6)
pvalues = chi2cdf(Fitchi2min,A.FitResult.dof);
histogram(pvalues,'FaceColor',rgb('CadetBlue'));
xlabel('p-value');
ylabel('');
legend(sprintf('mean = %.2f \n std = %.2f',mean(pvalues),std(pvalues)));
PrettyFigureFormat;
set(gca,'FontSize',16);


