% test response function broadening
MACE_Sigma = 1e-3.*[1,5,10,15,20:10:100,125,150,200];

savedir = [getenv('SamakPath'),'knm2ana/knm2_ResponseFuncton/results/'];
MakeDir(savedir);
savename = sprintf('%sknm2_broadenRF_min%.0fmeV_max%.0fmeV.mat',savedir,min(MACE_Sigma)*1e3,max(MACE_Sigma)*1e3);

if exist(savename,'file')
    load(savename)
else
RunAnaArg = {'RunList','KNM2_Prompt',...
    'fixPar','mNu E0 Norm Bkg',...
    'SysBudget',31,... % 31= knm2 preliminary input
    'FSDFlag','BlindingKNM2',...
    'fitter','minuit',...
    'NonPoissonScaleFactor',1,...
    'DataType','Twin',...
    'TwinBias_Q',18573.70};

A = MultiRunAnalysis(RunAnaArg{:});
A.exclDataStart = A.GetexclDataStart(40);
TBDIS_ref       = A.ModelObj.TBDIS;
A.RunData.TBDIS = TBDIS_ref;
A.Fit;
mNuref = A.FitResult.par(1);
%%

mNuSq = zeros(numel(MACE_Sigma),1);

for i=1:numel(MACE_Sigma)
A.ModelObj.MACE_Sigma = MACE_Sigma(i);
A.ModelObj.InitializeRF;
A.Fit;
mNuSq(i) = A.FitResult.par(1);
end

save(savename,'mNuSq','mNuSqref','RunAnaArg','MACE_Sigma');

end


%% plot

f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
plot(MACE_Sigma*1e3,mNuSq,'LineWidth',2);
hold on;
plot(MACE_Sigma*1e3,2.*MACE_Sigma.^2,'--','LineWidth',2);
PrettyFigureFormat('FontSize',24);
xlabel(sprintf('\\sigma (meV)'));
ylabel(sprintf('\\Delta{\\it m}_\\nu^2 (eV^2)'));
leg = legend('Response function broadening',sprintf('\\Delta{\\it m}_\\nu^2 = 2\\sigma^2'),...
    'Location','northwest');
leg.EdgeColor = rgb('Silver');



