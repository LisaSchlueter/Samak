qURange  = [95,20];
RunList = 'KNM1';
DataType = 'Real';
chi2 = 'chi2CMShape';
SysBudet = 22;
fixPar = 'mNu E0 Norm Bkg';% free par
FSDFlag       = 'SibilleFull';
ELossFlag     = 'KatrinT2A20';
AngularTFFlag = 'ON';

M = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat',...
    'DataType',DataType,...
    'fixPar',fixPar,...
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1.064,...
    'minuitOpt','min ; migrad',...
    'FSDFlag',FSDFlag,...
    'ELossFlag',ELossFlag,...
    'SysBudget',22,...
    'AngularTFFlag',AngularTFFlag);

M.chi2 = chi2;
%% systematics
if strcmp(chi2,'chi2CMShape')
    if strcmp(fixPar,'5 6 7 8 11')
        M.ComputeCM('FSDNorm_RelErr',0);
    else
        M.ComputeCM;
    end
end
%%
if strcmp(AngularTFFlag,'ON')
    AddsaveStr = sprintf('_AngTF_%s',ELossFlag);
else
    AddsaveStr = sprintf('_%s',ELossFlag);
end

[parqU, errqU, chi2qU, dofqU] = ...
    M.qUScan('qURange',qURange,...
    'saveplot','ON',...
    'RecomputeFlag','OFF',...
    'CorrMean','OFF',...
    'HoldOn','OFF',...
    'RelFlag','OFF',...
    'saveStr',AddsaveStr,...
    'RefLine',40);

