%% Settings
mNuBins  = 10;
mNuUpper = 2;   %eV

A = RelicNuDebug('Params','KNM1');
Sensitivities = 1:mNuBins;

for i=0:mNuBins
    A.Chi2Twin('TwinBias_mnuSq',mNuUpper*i/mNuBins);
    Sensitivities(i) = A.etaSensitivity;
end

matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/')];
savename = sprintf('SensVsMnuSq_%s_[0 %g].mat',A.Params,mNuUpper);

fig1=figure(1);
plot((0:mNuBins)*mNuUpper,Sensitivities,'LineWidth',2);
xlabel('\eta','FontSize',12);
ylabel('\chi^2','FontSize',12);
PrettyFigureFormat;