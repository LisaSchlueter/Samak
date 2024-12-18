%--------------------------------------------------------------------------%
% Simple reference script to se if all computers give the same fit results
% Single Run
% With Covariance Matrix
% Exclude 2 outer Rings
% Fixed neutrino mass (and FSD probabilities)
%% ------------------------------------------------------------------------%
A = MultiRunAnalysis('RunList',40668,'AnaFlag','StackPixel',...
    'chi2','chi2Stat','ringCutFlag','ex2','fixPar','1 5 6');
A.InitializeCM('DataDriven','ON','WGTS_CD_MolPerCm2_RelErr',0.05);
A.ComputeCM;
A.Fit;
A.PlotFit;


%% ------------------------------------------------------------------------%
% Results Lisa (29.06.2018)
%  m^2       = 0 +/- 0 eV^2
%  m         = 0 +/- 0 eV
% - - - - - - - - - - - - - - - - - - - - - - - 
%  (E0)eff  = 18573.4 +/- 0.875794 eV
% - - - - - - - - - - - - - - - - - - - - - - - 
%  B   = 331.019 +/- 10.6102 mcps
%  B total   = 331.019 +/- 10.6102 mcps
% - - - - - - - - - - - - - - - - - - - - - - - 
%  N average = 0.950984 +/- 0.0108764
% - - - - - - - - - - - - - - - - - - - - - - - 
%  Chi2/dof  = 6.66203/23
%  p-value  = 0.999639
%--------------------------------------------------------------------------%
% Results Thierry  (29.06.2018)
%   m^2       = 0 � 0 eV^2
%   m         = 0 � 0 eV
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   (E0)eff  = 18573.4 � 0.872629 eV
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   B   = 331.006 � 10.6089 mcps
%   B total   = 331.006 � 10.6089 mcps
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   N average = 0.952265 � 0.0108399
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   Chi2/dof  = 6.62428/23
%   p-value  = 0.999656
%--------------------------------------------------------------------------%
% Results Pablo (29.06.2018)
%   m^2       = 0 � 0 eV^2
%   m         = 0 � 0 eV
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   (E0)eff  = 18573.4 � 0.876284 eV
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   B   = 330.984 � 10.6132 mcps
%   B total   = 330.984 � 10.6132 mcps
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   N average = 0.9509 � 0.010925
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   Chi2/dof  = 6.5949/23
%   p-value  = 0.999669

