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

TwinOpt = {'Default','Twin_SameqUFlag','Twin_SameCDFlag','Twin_SameIsotopFlag'};

for i=1:numel(TwinOpt)
    %% load if possible
    savedir = [getenv('SamakPath'),'knm2ana/knm2_RunStacking/results/'];
    savename = [savedir,sprintf('knm2_RunStacking_%s_%.0feV_TwinOpt-%s.mat',RunList,range,TwinOpt{i})];
    
    if exist(savename,'file')
       d = importdata(savename);
    else
        if strcmp(TwinOpt{i},'Default')
            A = MultiRunAnalysis(CommonArg{:}); % default
        else
            A = MultiRunAnalysis(CommonArg{:},TwinOpt{i},'ON');
        end
        A.exclDataStart = A.GetexclDataStart(range);              % get correct range
        A.Fit;                                                    % fit
        FitResult = A.FitResult;                                  % save
        
        MakeDir(savedir);
        save(savename,'FitResult','CommonArg','range');
    end
end

