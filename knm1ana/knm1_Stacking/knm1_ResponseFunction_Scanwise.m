% scan-wise response functio for KNM1 (twins)
savedir = [getenv('SamakPath'),'knm1ana/knm1_Stacking/results/'];
savename = [savedir,sprintf('ResponseFunction_Scanwise_KNM1.mat')];

if exist(savename,'file')
    load(savename)
else 
    RunList = 'KNM1';
    range   = 40;         % 40eV range = 27 subruns
    Real = MultiRunAnalysis('RunList',RunList,...
        'chi2','chi2Stat',...
        'DataType','Twin',...
        'fixPar','mNu E0 Norm Bkg',... free parameter
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',1.064,...
        'minuitOpt','min ; minos',...
        'FSDFlag','Sibille',...
        'ELossFlag','KatrinT2',...
        'SysBudget',22,...
        'AngularTFFlag','OFF',...
        'DopplerEffectFlag','OFF',...
        'SynchrotronFlag','ON',...
        'TwinBias_Q',18573.70);
    Real.exclDataStart = Real.GetexclDataStart(range);
    Real.LoadSingleRunObj
    Te_Stack = Real.ModelObj.Te;
    qU27_Stack = Real.ModelObj.qU(27);
    RF_Stack = Real.ModelObj.RF(:,27);
    
    % response function for one illustrative retarding energy
    RF_single_raw = cellfun(@(x) x.RF(:,27)',Real.SingleRunObj,'UniformOutput',false);
    Te_single = cellfun(@(x) x.Te,Real.SingleRunObj,'UniformOutput',false);
    qU27_single = cell2mat(cellfun(@(x) x.qU(27),Real.SingleRunObj,'UniformOutput',false));
    %% interpolate at the same TE, do NOT correct for different qU(27)!!!
    RF_single = zeros(Real.nRuns,numel(Te_Stack));
    for r = 1:Real.nRuns
        RF_single(r,:) = interp1(Te_single{r},RF_single_raw{r},Te_Stack);
    end
    
    save(savename,'Te_Stack','qU27_Stack','RF_Stack','Te_single','RF_single_raw','qU27_single','RF_single');

end


%%
GetFigure;
[l,a] = boundedline(Te_Stack-qU27_Stack,mean(RF_single)',std(RF_single)');
hold on;
l.Color = rgb('Orange'); l.LineStyle = '-'; l.LineWidth = 2;
a.FaceAlpha = 1; a.FaceColor = rgb('SkyBlue');
pstack = plot(Te_Stack-qU27_Stack,RF_Stack,'r:','LineWidth',2);
xlim([-0.15,3.3]);
ylim([0 0.82]);

xlabel(sprintf('{\\itE}_{kin} - \\langle{\\itqU}\\rangle (eV)'));
ylabel('Transmission probability');
PrettyFigureFormat('FontSize',24);
leg2 = legend([a,l,pstack],sprintf('Scan-wise: 1\\sigma band'),....
    sprintf('\\langleScan-wise \\rangle'),...
    sprintf('Stacked scan'),'Location','northwest'); 
leg2.Title.String = 'Response function'; leg2.Title.FontWeight = 'normal'; 
PrettyLegendFormat(leg2);
leg2.FontSize = get(gca,'FontSize');

plotdir = strrep(savedir,'results','plots');
plotname = [plotdir,sprintf('ResponseFunction_Scanwise_KNM1.pdf')];
export_fig(plotname);

