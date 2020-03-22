range = 40;               % fit range in eV below endpoint
RunListName = 'KNM2_Prompt';
is_E0Offsets = (-0.3:0.05:0.3);

savedir = [getenv('SamakPath'),'knm2ana/knm2_ElossShift/results/'];
savename = sprintf('%sknm2_ElossShift_%s_nOffsets%.0f_min%.0fmeV_max%.0fmeV_%.0feVrange.mat',...
    savedir,RunListName,numel(is_E0Offsets),min(is_E0Offsets)*1e3,max(is_E0Offsets)*1e3,range);
if exist(savename,'file')
    load(savename)
else

E0 = knm2FS_GetE0Twins('SanityPlot','OFF');
RunAnaArg = {'RunList',RunListName,...  % define run number -> see GetRunList
    'fixPar','mNu E0 Bkg Norm',...         % free Parameter !!
    'DataType','Twin',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'MosCorrFlag','OFF',...
    'TwinBias_Q',E0,...% 18573.73 = default settings, 18574= const ISX, 18575= same Te for all runs
    'ROIFlag','14keV',...    % 18573.73 == new RF binning without interp
    'DopplerEffectFlag','FSD'};%,...

D = MultiRunAnalysis(RunAnaArg{:});
D.exclDataStart = D.GetexclDataStart(range); % find correct data, where to cut spectrum
D.ModelObj.RFBinStep = 0.02;
D.ModelObj.InitializeRF;
RunList = D.RunList;
%% broaden FSDs
Sigma = std(E0);
FSDArg = {'SanityPlot','ON','Sigma',Sigma};
D.ModelObj.LoadFSD(FSDArg{:});
D.ModelObj.ComputeTBDDS; D.ModelObj.ComputeTBDIS;

D.Fit; % test fit
%%
FitResults = cell(numel(is_E0Offsets),1); 

for i=1:numel(is_E0Offsets)
    progressbar(i/numel(is_E0Offsets));
    D.ModelObj.is_EOffset = is_E0Offsets(i);
    D.ModelObj.InitializeRF;
    D.Fit;
    FitResults{i} = D.FitResult;
end

MakeDir(savedir);
save(savename,'FitResults','is_E0Offsets','FSDArg','RunAnaArg','E0','RunList','range');
end

%% plot
mNuSq = cell2mat(cellfun(@(x) x.par(1),FitResults,'UniformOutput',0));

fel = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
plot(is_E0Offsets,zeros(numel(mNuSq),1),'-k','LineWidth',2);
hold on;
plot(is_E0Offsets,mNuSq,'LineWidth',3,'Color',rgb('DodgerBlue'));
PrettyFigureFormat('FontSize',24);
xlabel('Energy-loss shift (eV)');
ylabel(sprintf('\\Delta{\\itm}_\\nu^2 (eV^{ 2})'));
title(sprintf('%.0f eV range',range),'FontWeight','normal','FontSize',get(gca,'FontSize'));

%save
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = strrep(strrep(savename,'results','plots'),'.mat','.pdf');
export_fig(fel,plotname);
fprintf('save plot to %s\n',plotname)