%
% Compute TBDIS Covariance Matrix
% May 2018 Run Like Simulation
%
% Thierry Lasserre
% 21/03/2018
%

%% Init
spec = ref_tfcm();
nTrials = 200;

% Display
spec.DisplayTDBInfo;
spec.DisplayWGTSMACEInfo;

% Compute Nominal Spectra
spec.ComputeTBDDS
spec.ComputeTBDIS

% Inelastic Scattering Covariance Matrix
% spec.ComputeISProbMatrix('nTrials',nTrials,'Display','OFF','NIS',11,...
%     'MACE_Bmax_T_RelErr',1e-2,...
%     'WGTS_B_T_RelErr',1e-2,...
%     'WGTS_CD_MolPerCm2_RelErr',10e-2,...
%     'ISXsection_RelErr',1e-2);

%% TF Covariance Matrix
spec.ComputeRFCovMatrix('nTrials',nTrials,'Debug','OFF','NIS',11,...
    'Emin',-5,'Emax',200,...
    'MACE_Bmax_T_RelErr',1e-2,...
    'WGTS_B_T_RelErr',1e-2,...
    'MACE_Ba_T_RelErr',1e-2);

%% TBDIS WGTS+MACE Covariance Matrix
spec.ComputeTBDISCovarianceMatrix('nTrials',nTrials,...
    'StatFluct','OFF','Debug','OFF');

%% Display
cm=importdata('WGTSMACE_CovMat_200Trials_Simulation_Sec_86400.mat');
figure(2)
syst=errorbar(spec.qU-spec.Q,cm.TBDIS_av,sqrt(diag(cm.WGTSMACE_CovMat)));
hold on
stat=errorbar(spec.qU-spec.Q,cm.TBDIS_av,sqrt(cm.TBDIS_av),'Color','Red');
legend([stat,syst],'stat','syst');
xlabel('T_e - E_0 eV');
ylabel('TBD Integral Spectrum');
titlestr = sprintf('KATRIN - %g m/cm^2 - 1 day - 10 percent-ish uncertainties',spec.WGTS_CD_MolPerCm2);
title(titlestr);
hold off
%set(gca,'yscale','log')
PrettyFigureFormat

figure(3)
errorbar(spec.qU-spec.Q,cm.TBDIS_av./cm.TBDIS_av,sqrt(diag(cm.WGTSMACE_CovMat))./cm.TBDIS_av);
hold on
errorbar(spec.qU-spec.Q,cm.TBDIS_av./cm.TBDIS_av,sqrt(cm.TBDIS_av)./cm.TBDIS_av,'Color','Red');
hold off
xlabel('T_e - E_0 eV');
ylabel('Ratio');
PrettyFigureFormat

