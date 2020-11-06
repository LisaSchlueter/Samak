%% Settings
mNuBins  = 10;
mNuUpper = 2;   %eV

A = RelicNuDebug('Params','TDR');

matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/')];
savename = [matFilePath,sprintf('SensVsMnuSq_%s_[0_%g].mat',A.Params,mNuUpper)];

if exist(savename,'file')
    load(savename,'mNuBins','mNuUpper','Sensitivities')
else
    Sensitivities = 0:mNuBins;

    for i=1:mNuBins+1
        A.Chi2Fake('RunNr',10,...
            'Init_Opt',{'BKG_RateAllFPDSec',0.4,'mNuSq_i',mNuUpper*(i-1)/mNuBins,'TimeSec',3*365.242*24*3600},...
            'NetaBins',2,...
            'Plot','OFF');
        Sensitivities(i) = A.etaSensitivity;
    end

    save(savename,'mNuBins','mNuUpper','Sensitivities');
end

fig1=figure(1);
plot((0:mNuBins)*mNuUpper./mNuBins,Sensitivities,'LineWidth',2);
xlabel('m_{\nu}^{2} (eV^{2})','FontSize',12);
ylabel('\eta','FontSize',12);
PrettyFigureFormat;