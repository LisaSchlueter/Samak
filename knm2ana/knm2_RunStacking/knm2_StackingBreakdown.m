% --------------------------------------------------------------------------------
% Goal:
% estimate systematic  bias of run stacking for KNM2
% - investigate influence of column density, qU, isotopologues
%
% Method:
% - calculate twins with different settings
% - estimate neutrino mass bias
%
% Lisa Schlüter
% March 2020
% --------------------------------------------------------------------------------
%%  common settings
SavePlot = 'ON';
RecomputeFlag = 'ON';
range = 40; % in eV below E0
RunList = 'KNM2_Prompt';
fitter = 'minuit';
CommonArg = {'RunList',RunList,...
    'fixPar','mNu E0 Bkg Norm',...
    'DataType','Twin',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'TwinBias_Q',18573.70,...
    'NonPoissonScaleFactor',1,...
    'fitter',fitter};

switch fitter
    case 'matlab'
        fitStr = '_matlab';
    case 'minuit'
        fitStr = '';
end
TwinOpt = {'Default','Twin_SameqUFlag','Twin_SameqUfracFlag','Twin_SameTime',...
           'Twin_SameMTD','Twin_SameCDFlag','Twin_SameIsotopFlag','Sameall'};
mNuSq    = zeros(numel(TwinOpt),1);
mNuSqErr = zeros(numel(TwinOpt),1);%
%TwinBias_Time
%%
for i=1:numel(TwinOpt)
    %% load if possible
    savedir = [getenv('SamakPath'),'knm2ana/knm2_RunStacking/results/'];
    savename = [savedir,sprintf('knm2_RunStacking_%s_%.0feV_TwinOpt-%s%s.mat',RunList,range,TwinOpt{i},fitStr)];
    
    if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
        d = importdata(savename);
        mNuSq(i)    = d.FitResult.par(1);
        mNuSqErr(i) = d.FitResult.err(1);
        if isfield(d,'Runs')
            Runs = d.Runs;
        end
    else
        if strcmp(TwinOpt{i},'Default')
            A = MultiRunAnalysis(CommonArg{:}); % default
        elseif strcmp(TwinOpt{i},'DefaultImp')
            A = MultiRunAnalysis(CommonArg{:}); % default
            RFsigma = std(A.SingleRunData.qU,0,2);
            A.ModelObj.MACE_Sigma = mean(RFsigma);
            A.ModelObj.InitializeRF('RebinMode','Integral');
        elseif strcmp(TwinOpt{i},'Twin_SameMTD')
            A = MultiRunAnalysis(CommonArg{:},...% all flags on
                'Twin_SameqUFlag','ON',...
                'Twin_SameqUfracFlag','ON');
        elseif strcmp(TwinOpt{i},'Sameall')
            A = MultiRunAnalysis(CommonArg{:},...% all flags on
                'Twin_SameqUFlag','ON',...
                'Twin_SameCDFlag','ON',...
                'Twin_SameIsotopFlag','ON',...
                'Twin_SameqUfracFlag','ON');
        elseif strcmp(TwinOpt{i},'Twin_SameTime')
              A = MultiRunAnalysis(CommonArg{:},'TwinBias_Time',7416.0277); 
              % all slow control default, but all runs have the same total measurement time
        else
            A = MultiRunAnalysis(CommonArg{:},TwinOpt{i},'ON');
        end
        A.exclDataStart = A.GetexclDataStart(range);              % get correct range
        A.Fit;                                                    % fit
        FitResult = A.FitResult;                                  % save
        
        MakeDir(savedir);
        Runs = A.RunList;
        save(savename,'FitResult','CommonArg','range','Runs');
        
        mNuSq(i)   = FitResult.par(1);
        mNuSqErr(i) = FitResult.err(1);
    end
end

%% plot

f1 = figure('Units','normalized','Position',[0.1,0.1,0.8,0.5]);
x = 1:numel(TwinOpt);
pref = plot(linspace(-1,10,numel(x)),zeros(numel(x),1),'-k','LineWidth',2);
hold on;

switch range
    case 40
        PlotArg = {'--d','Color',rgb('Orange'),'MarkerFaceColor',rgb('Orange')};
        legStr  = '40 eV range';
    case 90
        PlotArg = {':o','Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue')};
        legStr = '90 eV range';
end

p1 = errorbar(x,mNuSq,zeros(numel(x),1),PlotArg{:},...
    'MarkerSize',9,'LineWidth',pref.LineWidth,'CapSize',0);
leg = legend(p1,legStr);
xticks(x);
leg.EdgeColor = rgb('Silver');
%xticklabels(cellfun(@(x) strrep(strrep(strrep(x,'Twin_',''),'Flag',''),'Same','Same '),TwinOpt,'UniformOutput',0));

xticklabels({'default',sprintf('same {\\itqU}'),sprintf('same {\\itqU}frac'),'same Time','same MTD',...
    sprintf('same {\\it\\rhod\\sigma}'),sprintf('same TT,HT,DT'),'all same'});

xlabel('Twin species')
ylabel(sprintf('\\Delta{\\itm}_\\nu^2 (eV^2)'));
t = title(sprintf('%.0f KNM2 twin runs - common E_0 = %.2f eV',numel(Runs),18573.70));
t.FontWeight = 'normal';
PrettyFigureFormat('FontSize',22);
set(gca,'XMinorTick','off');
xlim([x(1)-0.2 x(end)+0.2])
%ylim([-14 6]*1e-3);

if strcmp(SavePlot,'ON')
    grid on
savedirplot = strrep(savedir,'results','plots');
MakeDir(savedirplot);
savename = sprintf('%sknm2_StackingBreakdown_%.0feV%s.pdf',savedirplot,range,fitStr);
export_fig(f1,savename);
print(f1,strrep(savename,'.pdf','.png'),'-dpng','-r300')
fprintf('save plot to %s \n',savename);
end





