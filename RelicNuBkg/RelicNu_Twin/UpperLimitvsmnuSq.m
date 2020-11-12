%% Settings
mNuBins  = 10;
mNuUpper = 2;   %eV

A = RelicNuDebug('Params','KNM1');

matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/')];
savename = [matFilePath,sprintf('SensVsMnuSq_%s_[0_%g].mat',A.Params,mNuUpper)];

if ~exist(savename,'file')
    load(savename,'ScanPoints','Sensitivities')
else
    Sensitivities = 0:mNuBins;
    ScanPoints = [(0:mNuBins)*mNuUpper./mNuBins 5];

    for i=1:numel(ScanPoints)
        A.Chi2Twin('TwinBias_mnuSq',ScanPoints(i),...
            'etarange',11,...
            'etafactor',5.001,...
            'NetaBins',2,...
            'Plot','OFF');
        Sensitivities(i) = A.etaSensitivity;
    end

    save(savename,'ScanPoints','Sensitivities');
end

fig1=figure(1);
plot(ScanPoints,Sensitivities,'LineWidth',2);
xlabel('m_{\nu}^{2} (eV^{2})','FontSize',12);
ylabel('\eta','FontSize',12);
PrettyFigureFormat;