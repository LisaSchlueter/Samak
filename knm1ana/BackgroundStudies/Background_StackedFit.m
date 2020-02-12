%
% Investigate KNM1 Non-Poissonian Component of Background
% T. Lasserre
% Last Modified: 19/6/2019
%

%%
% Read Real Data
% Fit Staxcked Spectrum - Focusing on Background Region
% Neutrino Mass is set to zero
%
if ~exist('Real','var')
    RunList = 'KNM1';
    Real = MultiRunAnalysis('RunList',RunList,...
        'DataType','Real',...
        'minuitOpt','min;minos','fitter','minuit',...
        'exclDataStart',30,...
        'fixPar','1 5 6 7 8 9 10 11');
    qUmin = round(Real.ModelObj.qU(Real.exclDataStart));
    qUmax = round(Real.ModelObj.qU(end));
    range = round(Real.ModelObj.qU(Real.exclDataStart)-Real.ModelObj.Q_i);
end

for InitRange=14:1:35
    Real.exclDataStart=InitRange;
    Real.Fit
end

