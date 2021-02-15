%% Settings
mNuBins  = 5;
mNuUpper = 4;   %eV²

A = RelicNuAnalysis('Params','KNM2');

matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/')];
savename = [matFilePath,sprintf('SensVsmnuSQ_KNM2_[1_%g].mat',mNuUpper)];

if exist(savename,'file')
    load(savename,'ScanPoints','Sensitivities')
else
    Sensitivities = zeros(1,mNuBins);
    ScanPoints = linspace(0,mNuUpper,mNuBins);

    for i=1:numel(ScanPoints)
        A.Chi2Fake('Recompute','ON',...
            'RunNr',10,...
            'range',40,...
            'etafactor',3,...
            'etarange',11,...
            'fitPar','mNu E0 Norm Bkg',...
            'Syst','ON',...
            'Init_Opt',{'mnuSq_i',ScanPoints(i)},...
            'NetaBins',2,...
            'Plot','OFF');
        Sensitivities(i) = A.etaSensitivity;
    end

    save(savename,'ScanPoints','Sensitivities');
end

fig1=figure(1);
plot(ScanPoints./2018850,Sensitivities,'LineWidth',2);
xlabel('Measurement Time (multiple of KNM1)','FontSize',12);
ylabel('\eta','FontSize',12);
PrettyFigureFormat;