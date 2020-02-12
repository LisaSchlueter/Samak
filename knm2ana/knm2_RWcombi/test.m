% Calculate neutrino mass squared shift as a function of FSD broadening width
RunList = 'KNM2_Prompt';
FSDFlag = 'BlindingKNM2';
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
chi2 = 'chi2Stat';
DataType = 'Twin';
range = 40;
Twin_SameqUFlag = 'ON';

CommonArg = {'RunList',RunList,...
    'FSDFlag',FSDFlag,...
    'ELossFlag',ELossFlag,...
    'AnaFlag',AnaFlag,...
    'chi2',chi2,...
    'Twin_SameqUFlag',Twin_SameqUFlag,....
    'DataType','Twin',...
    'TwinBias_Q',18573.7,...
    'fixPar','mNu E0 Bkg Norm'};

savedir = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];
savename = [savedir,sprintf('knm2_RWcombi_FSDSigma2mNuSq.mat')];

if exist(savename,'file')
    load(savename);
else
    
    Sigma = [0.01,0.025,0.05,0.075,0.1,0.15,0.18,0.2];
    FitResult = cell(numel(Sigma),1);
    
    M = MultiRunAnalysis(CommonArg{:});
    M.exclDataStart = M.GetexclDataStart(range);
    
    for i=1:numel(Sigma)
        progressbar(i/numel(Sigma));
        M.ModelObj.LoadFSD('Sigma',Sigma(i));
        M.ModelObj.ComputeTBDDS;
        M.ModelObj.ComputeTBDIS;
        
        M.Fit;
        FitResult{i} = M.FitResult;
    end
    
    save(savename,'FitResult','Sigma','CommonArg');
end

%% plot
mNuSq = cellfun(@(x) x.par(1),FitResult);
plot(Sigma.*1e3,mNuSq.*1e3);
hold on;
plot(Sigma.*1e3,(2.*Sigma.^2).*1e3)
PrettyFigureFormat;
xlabel(sprintf('\\sigma_G (meV)'));
ylabel(sprintf('{\\itm}_\\nu^2 (meV^2)'));