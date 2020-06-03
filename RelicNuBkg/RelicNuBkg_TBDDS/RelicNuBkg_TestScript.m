% test script
% can be deleted later

% ini TBD object with Design Report (TDR) properties:
A = ref_RelicNuBkg_TDR('DopplerEffectFlag','OFF');

% Compute diff spec.
A.ComputeTBDDS;

A.PlotTBDDS;
