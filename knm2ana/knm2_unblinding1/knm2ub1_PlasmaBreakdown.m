% Calculte and/or plot systematics breakdown for KNM2 unblinding 1


file  = [getenv('SamakPath'),'inputs/CovMat/LongPlasma/CM/LongPlasma_KNM2Prompt_5000Trials_0meVEloss_99000meVElossErr_0meVSigma_127meVSigmaErr_0CorrCoeff_Troitsk+.mat'];
d = importdata(file);


f1 = figure('Units','normalized','Position',[0.1,0.1,0.4,0.45]);
h1 = histogram(d.MACE_Var_v,'Normalization','probability','FaceColor',rgb('DodgerBlue'));
PrettyFigureFormat;
xlabel(sprintf('\\sigma^2 (eV^2)'));
ylabel('Frequency');
set(gca,'LineWidth',1);
legend(sprintf('std = %.4f eV^2',std(d.MACE_Var_v)));
legend box off

%%
f2 = figure('Units','normalized','Position',[0.1,0.1,0.4,0.45]);
h1 = histogram(d.is_EOffset_v,'Normalization','probability','FaceColor',rgb('DodgerBlue'));
PrettyFigureFormat;
xlabel(sprintf('\\Delta_{10} (eV)'));
ylabel('Frequency');
set(gca,'LineWidth',1);
legend(sprintf('std = %.4f eV',std(d.is_EOffset_v)));
legend box off

%%
Var_i          = 0.0149;%0.0149; % in eV^2
VarErr         = 0.0161; % in eV^2
PlasmaPar1     = randn(5e4,1).*VarErr+Var_i; % broadening
DeltaMax       = sqrt(abs(PlasmaPar1))./1.3;
PlasmaPar2     = rand(5e4,1).*2.*DeltaMax-DeltaMax;

std(PlasmaPar2)