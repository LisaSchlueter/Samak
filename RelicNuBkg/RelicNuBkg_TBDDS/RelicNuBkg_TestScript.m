% test script
% can be deleted later

% init TBD object with Design Report (TDR) properties:
A = ref_RelicNuBkg_TDR;

% Compute diff spec.
A.ComputeTBDDS;
A.PlotTBDDS;