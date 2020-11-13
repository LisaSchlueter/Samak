%% Settings
mNuBins  = 10;
mNuUpper = 2;   %eV

A = RelicNuDebug('Params','TDR');

matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/')];
savename = [matFilePath,sprintf('SensVsMnuSq_TDR_[0_%g].mat',mNuUpper)];

if exist(savename,'file')
    load(savename,'ScanPoints','Sensitivities')
else
    Sensitivities = 0:(mNuBins+1);
    ScanPoints = [(0:mNuBins)*mNuUpper./mNuBins 5];

    for i=1:numel(ScanPoints)
        A.Chi2Fake('RunNr',10,...
            'etafactor',5.001,...
            'etarange',9,...
            'Init_Opt',{'mNuSq_i',ScanPoints(i)},...
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