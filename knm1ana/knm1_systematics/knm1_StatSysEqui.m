% fit neutrino mass at different ranges (-> qU Scan)
% once stat. only
% once with systematics
% search for equilibrium of statistical and systematic uncertainties
% calculate, plot in different script
qURange  = [95,20];
DataType = 'Real';
chi2 = 'chi2CMShape';
NP = 1.064;

% create mini short cut file for plotting
savedir = [getenv('SamakPath'),'knm1ana/knm1_qUScan/results/'];
savename = sprintf('%sknm1_qUScan_Mini_%s_%s_NP%.4g_%.0feV_to_%.0feV.mat',...
    savedir,DataType,chi2,NP,qURange(1),qURange(2));

if exist(savename,'file')
    load(savename)
else
    SysBudet = 22;
    fixPar = 'mNu E0 Norm Bkg';% free par
    FSDFlag       = 'SibilleFull';
    ELossFlag     = 'KatrinT2';%A20';
    AngularTFFlag = 'OFF';
    
    M = MultiRunAnalysis('RunList','KNM1',...
        'chi2',chi2,...
        'DataType',DataType,...
        'fixPar',fixPar,...
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',NP,...
        'minuitOpt','min ; minos',...
        'FSDFlag',FSDFlag,...
        'ELossFlag',ELossFlag,...
        'SysBudget',22,...
        'AngularTFFlag',AngularTFFlag);
    
    %%
    if strcmp(AngularTFFlag,'ON')
        AddsaveStr = sprintf('_AngTF_%s',ELossFlag);
    else
        AddsaveStr = sprintf('_%s',ELossFlag);
    end
    
    %%
    [parqU, errqU, chi2qU, dofqU,~,~,mNuSqErr] = ...
        M.qUScan('qURange',qURange,...
        'saveplot','OFF',...
        'RecomputeFlag','ON',...
        'CorrMean','OFF',...
        'HoldOn','OFF',...
        'RelFlag','OFF',...
        'saveStr',AddsaveStr,...
        'RefLine',40,...
        'PlotResult','OFF',...
        'randMC','OFF');
    %
        save(savename,'parqU','errqU','chi2qU','dofqU','M','mNuSqErr');
end


