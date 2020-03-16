% Impact of energy loss function convolution range on mnu-bias
% KNM2 twins
% March 2020, Lisa

savedir = [getenv('SamakPath'),'knm2ana/knm2_ElossRangeBinning/results/'];
savename = sprintf('%sknm2_ElossRangeBinningTest.mat',savedir);
RecomputeFlag = 'OFF';

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
else
    RunAnaArg = {'RunList','KNM2_Prompt',...  % define run number -> see GetRunList
        'fixPar','mNu E0 Bkg Norm',...         % free Parameter !!
        'DataType','Twin',...              % Real, Twin or Fake
        'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2','chi2Stat',...              % statistics only
        'NonPoissonScaleFactor',1,...
        'MosCorrFlag','OFF',...
        'TwinBias_Q',18573.7,...
        'ROIFlag','14keV',...
        'DopplerEffectFlag','FSD'};
    
    %% build object of MultiRunAnalysis class
    D = MultiRunAnalysis(RunAnaArg{:});
    range = 40;               % fit range in eV below endpoint
    D.exclDataStart = D.GetexclDataStart(range); % find correct data, where to cut spectrum
    
    %% Fit -> fit results are in property: A.FitResult
    ELossRanges = [500,1000,4000,9288];
    FitResults = cell(numel(ELossRanges),1);
    
    for i=1:numel(ELossRanges)
        D.ModelObj.recomputeRF = 'ON';
        D.ModelObj.RF = CalcRF(D.ModelObj,ELossRanges(i));
        D.Fit;
        FitResults{i} = D.FitResult;
    end
    %%
    mNuSq = cell2mat(cellfun(@(x) x.par(1),FitResults,'UniformOutput',false));
    MakeDir(savedir);
    save(savename,'FitResults','mNuSq','ELossRanges');
end

%% plot
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
plot(ELossRanges,mNuSq-mean(mNuSq),'-o','LineWidth',2,'MarkerFaceColor',rgb('DodgerBlue'));
xlabel('Energy-loss range (eV)');
ylabel(sprintf('{\\itm}_\\nu^2 - \\langle{\\itm}_\\nu^2\\rangle(eV)'));
PrettyFigureFormat('FontSize',24);
MakeDir(strrep(savedir,'results','plots'));
plotname = strrep(strrep(savename,'results','plots'),'.mat','.pdf');
export_fig(f1,plotname);
%%
function out =  CalcRF(obj,ELossRange)
tfloc = @obj.ComputeRF; %function handle of response function
parTe = obj.Te;
parqU = obj.qU;
parRF = zeros(obj.nTe,obj.nqU);

for ii = 1:obj.nqU
    parRF(:,ii) = tfloc(parTe,parqU(ii),'ELossRange',ELossRange);
end
out = parRF;
end



