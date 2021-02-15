%% Settings
TwinmnuSq = -1;
mNuBins   = 9;
mNuUpper  = 4;   %eV

A = RelicNuAnalysis('Params','KNM1');

matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/')];
savename = [matFilePath,sprintf('SensVsMnuSq_%s_[0_%g]_mnuSqfix.mat',A.Params,mNuUpper)];

if exist(savename,'file')
    load(savename,'ScanPoints','Sensitivities')
else
    Sensitivities = linspace(0,mNuUpper,mNuBins);
    ScanPoints = linspace(4,mNuUpper,mNuBins);

    for i=1:numel(ScanPoints)
        A.Chi2Twin('Recompute','ON',...
            'TwinBias_mnuSq',TwinmnuSq,...
            'fitPar','E0 Norm Bkg',...
            'etarange',11,...
            'etafactor',6.001,...
            'NetaBins',2,...
            'DeltaChi2',2.71,...
            'Plot','OFF',...
            'mnuSqfix',ScanPoints(i));
        Sensitivities(i) = A.etaSensitivity;
    end

    save(savename,'ScanPoints','Sensitivities');
end

fig1=figure(1);
plot(ScanPoints,Sensitivities,'LineWidth',2);
xlabel('m_{\nu}^{2} (eV^{2})','FontSize',12);
ylabel('\eta','FontSize',12);
PrettyFigureFormat;