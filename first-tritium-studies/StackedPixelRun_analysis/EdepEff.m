function EdepEff(varargin)
% Fit Run List Run-by-Run and study E0 distribution
% With and Without Efficiency Correction
%  - Stacked-Pixels
%  - Stacked-Runs
%
% Input:
%  - RunList         0 (for all) or list as [run1 ... runX]
%  - chi2:           'chi2Stat', 'chi2CM', 'chi2CMFrac','chi2CMShape'...
%  - fitter:         'minuit', 'matlab'
%  - exclDataStart:  7=400eV, 9=200eV, ...
%  - displayHist:   'ON','OFF'
%  - displayFit:    'ON','OFF'
%  - fixPar:        '1 5 6' (1=mass, 2=E0, 3=B, 4=B, 5=Pgs, 6=1-Pgs)
%  - ringCutFlag:   'ex2','ex2b','OFF'
%  - AnaFlag:       'StackPixel'
% Output:
% -
%
% T. Lasserre, Pablo Morales, Lisa Schlueter
% Last Update: June 24 2018
%

close all
unix('rm -rf plots/tmp');
mkdir  plots/tmp

% Parser
p = inputParser;
p.addParameter('RunList',0);
p.addParameter('chi2','chi2Stat',@(x)ismember(x,{'chi2Stat', 'chi2CM', 'chi2CMFrac','chi2CMShape', 'chi2P','chi2Nfix'}));
p.addParameter('fitter','minuit',@(x)ismember(x,{'minuit','matlab'}));
p.addParameter('exclDataStart',7,@(x)isfloat(x));
p.addParameter('displayHist','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('displayFit','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('fixPar','1 5 6',@(x)ischar(x));
p.addParameter('ringCutFlag','ex2',@(x)ismember(x,{'ex2','ex2b','OFF'}));
p.addParameter('AnaFlag','StackPixel');
p.parse(varargin{:});

RunList       = p.Results.RunList;
chi2          = p.Results.chi2;
fitter        = p.Results.fitter;
exclDataStart = p.Results.exclDataStart;
displayHist   = p.Results.displayHist;
displayFit    = p.Results.displayFit;
fixPar        = p.Results.fixPar;
ringCutFlag   = p.Results.ringCutFlag;
AnaFlag       = p.Results.AnaFlag;

index =  [1602  1402  1204  1002   804   602   402   304   202   177   154   129   104    92    82    72    62    52    45    32    22    15     2    -6   -18   -28];

%%
if RunList == 0
    [ n, ~ , RunList ] = GetRunList( '../../tritium-data/hdf5/','2f-fpd00*.h5',1,'string');
    Remove1=RunList==40773;RunList(Remove1)=[]; % Sterile Neutrino Run - High Count Rate
    Remove2=RunList==40806;RunList(Remove2)=[]; % Sterile Neutrino Run - High Count Rate
    Remove3=RunList==40995;RunList(Remove3)=[]; % Sterile Neutrino Run - High Count Rate
    Remove4=RunList==40761;RunList(Remove4)=[]; % Stability Run - No Scan
    Remove5=RunList==40762;RunList(Remove5)=[]; % Stability Run - No Scan
    Remove6=RunList==40807;RunList(Remove6)=[]; % Anomaleous DT uncertainty
    Remove7=RunList==40769;RunList(Remove7)=[]; % Sterile Neutrinos? 31 qU
    Remove7=RunList==40770;RunList(Remove7)=[]; % Sterile Neutrinos? 31 qU
    Remove8=RunList==40805;RunList(Remove8)=[]; % Sterile Neutrinos? 18 qU
    RunList=RunList(RunList<41003);
    RunList=RunList(RunList>40537);
    RunList1=RunList(RunList<40693);
    RunList2=RunList(RunList>40975);
    RunList=[RunList1 RunList2]; % 100% CD
    RunList=RunList(RunList>40538 & RunList<40693);
    disp(RunList);
end

%% Fit All RunList
for i=1:numel(RunList)
    fprintf(2,'processing run %g ...\n',RunList(i))
    close all
    FT_edepOFF(i) = MultiRunAnalysis('RunList',RunList(i),'chi2',chi2,'exclDataStart',exclDataStart,...
        'fixPar','1 5 6','ringCutFlag','ex2','AnaFlag','StackPixel','DataEffCorr','OFF');
    close all
    FT_edepON(i)  = MultiRunAnalysis('RunList',RunList(i),'chi2',chi2,'exclDataStart',exclDataStart,...
        'fixPar','1 5 6','ringCutFlag','ex2','AnaFlag','StackPixel','DataEffCorr','ROI+PileUp');
    
    if ~strcmp(chi2,'chi2Stat')
        FT_edepOFF(i).InitializeCM('DataDriven','ON','WGTS_CD_MolPerCm2_RelErr',0.05);
        FT_edepOFF(i).ComputeCM;
        FT_edepON(i).InitializeCM('DataDriven','ON','WGTS_CD_MolPerCm2_RelErr',0.05);
        FT_edepON(i).ComputeCM;
    end
    FT_edepOFF(i).Fit;
    FT_edepON(i).Fit;
    switch displayFit
        case 'ON'
            FT_edepOFF(i).PlotFit('saveplot','OFF');
            FT_edepON(i).PlotFit('saveplot','OFF');
    end
end

%% Gather Results
DeltaE0 = [];
DeltaChi2=[];
PvalueON=[];
PvalueOFF=[];
r='';
for i=1:numel(RunList)
    % Endpoint
    DeltaE0     = [DeltaE0    FT_edepON(i).FitResult.par(2)-FT_edepOFF(i).FitResult.par(2)];
    % p-value
    DeltaChi2   = [DeltaChi2  FT_edepON(i).FitResult.chi2min,-FT_edepOFF(i).FitResult.chi2min];
    PvalueON    = [PvalueON  chi2pvalue(FT_edepON(i).FitResult.chi2min,FT_edepON(i).FitResult.dof)];
    PvalueOFF   = [PvalueOFF chi2pvalue(FT_edepOFF(i).FitResult.chi2min,FT_edepOFF(i).FitResult.dof)];
    
    % Labels
    r    = [r num2str(RunList(i))];
    r    = string(r);
end


switch displayHist
    case 'ON'
        myMainTitle=[sprintf('Impact of Energy Dependent Efficiency Correction - ROI+Pile-Up')];
        
        %% Fit Parameters
        maintitle=myMainTitle;
        savefile=sprintf('plots/AllRunsROIPU.pdf');
        x   = linspace(1,numel(RunList),numel(RunList));
        fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 2000]);
        a=annotation('textbox', [0 0.91 1 0.1], ...
            'String', maintitle, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';
        subplot(2,1,1)
        hoff=errorbar(x,PvalueOFF,PvalueOFF.*0,'-s','LineStyle','none','MarkerSize',10,...
            'MarkerEdgeColor','red','MarkerFaceColor','red');
        hold on
        hon=errorbar(x,PvalueON,PvalueON.*0,'-s','LineStyle','none','MarkerSize',10,...
            'MarkerEdgeColor','blue','MarkerFaceColor','blue');
        line([0.5,numel(RunList)+0.5],[0.05 0.05]);
        hold off
        legend([hon hoff],'With efficiency correction','Without efficiency correction')
        xticks(x)
        xticklabels(r)
        ylabel('p-value')
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        set(gca, 'YScale', 'log')
        sp1 = sprintf('<\\Delta\\chi^2>=%.2f \\pm %.2f (std)',mean(DeltaChi2),std(DeltaChi2));title(sp1)
        grid on
        set(gca,'FontSize',16);
        PrettyFigureFormat
        subplot(2,1,2)
        bar(x,DeltaE0,'facecolor',rgb('CadetBlue'));
        xticks(x)
        xticklabels(r)
        ylabel('\Delta(E_0) (eV) ON-OFF')
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        sp1 = sprintf('<E0>=%.2f eV \\pm %.2f eV (std)',mean(DeltaE0),std(DeltaE0));title(sp1)
        grid on
        set(gca,'FontSize',16);
        PrettyFigureFormat
        %        publish_figurePDF(gcf,savefile);
        
        
end

end

