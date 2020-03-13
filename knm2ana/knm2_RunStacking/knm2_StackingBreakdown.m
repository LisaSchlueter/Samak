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
range = 40; % in eV below E0
RunList = 'KNM2_Prompt';
CommonArg = {'RunList',RunList,...
    'fixPar','mNu E0 Bkg Norm',...
    'DataType','Twin',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'TwinBias_Q',18573.70,...
    'NonPoissonScaleFactor',1};

TwinOpt = {'Default','Twin_SameqUFlag','Twin_SameCDFlag','Twin_SameIsotopFlag','Sameall'};
mNuSq = zeros(numel(TwinOpt),1);

for i=1:numel(TwinOpt)
    %% load if possible
    savedir = [getenv('SamakPath'),'knm2ana/knm2_RunStacking/results/'];
    savename = [savedir,sprintf('knm2_RunStacking_%s_%.0feV_TwinOpt-%s.mat',RunList,range,TwinOpt{i})];
    
    if exist(savename,'file')
        d = importdata(savename);
        mNuSq(i) = d.FitResult.par(1);
    else
        if strcmp(TwinOpt{i},'Default')
            A = MultiRunAnalysis(CommonArg{:}); % default
        elseif strcmp(TwinOpt{i},'Sameall')
              A = MultiRunAnalysis(CommonArg{:},...% all flags on
                  'Twin_SameqUFlag','ON',...
                  'Twin_SameCDFlag','ON',...
                  'Twin_SameIsotopFlag','ON',...
                  'Twin_SameqUfracFlag','ON'); 
        else
            A = MultiRunAnalysis(CommonArg{:},TwinOpt{i},'ON');
        end
        A.exclDataStart = A.GetexclDataStart(range);              % get correct range
        A.Fit;                                                    % fit
        FitResult = A.FitResult;                                  % save
        
        MakeDir(savedir);
        save(savename,'FitResult','CommonArg','range');
        
        mNuSq(i) = FitResult.par(1);
    end
end

%% plot
figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
x = 1:numel(TwinOpt);
pref = plot(x,zeros(numel(x),1),'-k','LineWidth',2);
hold on;
p1 = errorbar(x,mNuSq,zeros(numel(x),1),':o',...
    'LineWidth',pref.LineWidth,'Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue'),'MarkerSize',8);
xticks(x);
xticklabels(cellfun(@(x) strrep(strrep(strrep(x,'Twin_',''),'Flag',''),'Same','Same '),TwinOpt,'UniformOutput',0));
xlabel('Twin species')
ylabel(sprintf('\\Delta{\\itm}_\\nu^2 (eV^2)'));
PrettyFigureFormat('FontSize',22);





