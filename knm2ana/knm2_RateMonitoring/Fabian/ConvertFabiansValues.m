
%% uniform
savedir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/'];
savenameBU = sprintf('%sSamakKNM2_NPBroadeningsInRW123uniform_mV.mat',savedir);
savenameSU = sprintf('%sSamakKNM2_ShiftDriftInRW123uniform_mVperDay.mat',savedir);

dBU = importdata(savenameBU);
dSU = importdata(savenameSU);

NPfactor    = dBU.NP';
NPfactorErr = dBU.NP_E'.*1e-03; 
Sigma       = dBU.Broadening'.*1e-03;                   % in eV
SigmaErr    = dBU.Broadening_E'.*1e-03; 
Shift       = dSU.CrateEquivalentmVAverage'.*1e-03;     % in eV
ShiftErr    = dSU.CrateEquivalentmVAverageE'.*1e-03;    % in eV
Drift       = dSU.SlopeEquivalent_mV';                  % in mV / day
DriftErr    = dSU.SlopeErrorEquivalent_mV';             % in mV / day

savenameU = sprintf('%sknm2_RManalysis_Uniform.mat',savedir);
save(savenameU,'NPfactor','NPfactorErr','Sigma','SigmaErr','Shift','ShiftErr','Drift','DriftErr');

%% ring-wise
clear all
savedir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/'];
savenameBR = sprintf('%sSamakKNM2_NPBroadeningsInRW123PSR1234_mV.mat',savedir);
savenameSR = sprintf('%sSamakKNM2_ShiftDriftInRW123PSR1234_mVperDay.mat',savedir);

dBR = importdata(savenameBR);
dSR = importdata(savenameSR);

NPfactor    = dBR.NP;
NPfactorErr = dBR.NP_E;
Sigma       = dBR.Broadening.*1e-03;                   % in eV
SigmaErr    = dBR.Broadening_E.*1e-03;
Shift       = dSR.CrateEquivalentmVAverage'.*1e-03;     % in eV
ShiftErr    = dSR.CrateEquivalentmVAverageE'.*1e-03;    % in eV
Drift       = dSR.SlopeEquivalent_mV';                  % in mV / day
DriftErr    = dSR.SlopeErrorEquivalent_mV';             % in mV / day

savenameR = sprintf('%sknm2_RManalysis_Ring.mat',savedir);
save(savenameR,'NPfactor','NPfactorErr','Sigma','SigmaErr','Shift','ShiftErr','Drift','DriftErr');



