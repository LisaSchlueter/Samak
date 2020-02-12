function FitResults = FitPlot_FitParameterEvolution(varargin)
% Run-by-Run Fit (w, w/o systematics) to study fit parameter evolution
% Input:
%  -  MultiRunAnalysis Object
%  - Parameter      'E0','N','B', which evolution shall be studied
%  - displayHist:   'ON','OFF' 
%  - displayFit:    'ON','OFF' 

% Output:
% - Fit Results
% - Fit and Plot distribution of E0-18573.7 eV
%   with chisquare distribution and E0 histogram
% 
% T. Lasserre, Pablo Morales, Lisa Schlueter
% Last Update: September 10, 2018
%
if exist('./plots/tmp/','dir')~=7 %if folder doesn't exist
mkdir  plots/tmp
end

% Parser
p = inputParser;
p.addParameter('RunObj','',@(x)isa(x,'MultiRunAnalysis'));
p.addParameter('displayHist','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('displayFit','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('Parameter','E0',@(x)ismember(x,{'E0','N','B'}));
p.addParameter('FitResultsFlag','ON',@(x)ismember(x,{'ON','OFF'})); % for public plots
p.addParameter('mergePDF','ON',@(x)ismember(x,{'ON','OFF'})); %merge all plots into 1 pdf and remove single plots
p.addParameter('saveplot','ON',@(x)ismember(x,{'ON','OFF','png','pdf','eps'}));
p.addParameter('SysFlag','all',@(x)ismember(x,{'all','corr'})); %all systematics or only correlated systematics (-> all but TASR)

p.parse(varargin{:});

RunObj         = p.Results.RunObj;
displayHist    = p.Results.displayHist;
displayFit     = p.Results.displayFit;
Parameter      = p.Results.Parameter;
FitResultsFlag = p.Results.FitResultsFlag;
mergePDF       = p.Results.mergePDF;
saveplot       = p.Results.saveplot;
SysFlag        = p.Results.SysFlag;

RunList = RunObj.SingleRunData.Runs(RunObj.SingleRunData.Select_all);
%% Fit All RunList
FT = cell(numel(RunList),1);
parfor i=1:numel(RunList)
    FT{i} = RunAnalysis('RunNr',RunList(i),'chi2',RunObj.chi2,'exclDataStart',RunObj.exclDataStart,...
        'fixPar',RunObj.fixPar,'ringCutFlag',RunObj.ringCutFlag,'AnaFlag',RunObj.AnaFlag,'DataEffCorr',RunObj.DataEffCorr);
    if strcmp(FT{i}.chi2,'chi2CM') && strcmp(SysFlag,'corr') %take all sys - TASR
        [CM_TASR, ~, ~, ~] = FT{i}.ComputeMultiPixelCM('CMFrac',FT{i}.FitCM_Obj.MultiCovMatFrac.CM_TASR);
        %[CM_RF, ~, ~, ~]   = FT{i}.ComputeMultiPixelCM('CMFrac',FT{i}.FitCM_Obj.MultiCovMatFrac.CM_RF);
        FT{i}.FitCM =   FT{i}.FitCM - CM_TASR; %- CM_RF;
    end
    FT{i}.Fit;
    switch displayFit
        case 'ON'
            FT{i}.PlotFit('saveplot','ON');
    end
end
%% Gather Results
    % Time
    T     = cellfun(@(x) x.ModelObj.TimeSec ,FT);
    TC    = cumsum(T); % cumulated time
    % Endpoint
    E0     = cellfun(@(x) x.FitResult.par(2) ,FT);
    E0Err  = cellfun(@(x) x.FitResult.err(2) ,FT);
    % Background
    B     = cellfun(@(x) x.ModelObj.BKG_RateAllFPDSec+x.FitResult.par(3),FT);
    BErr  = cellfun(@(x) x.FitResult.err(3),FT);
    % Normalizationn  
    N     = cellfun(@(x) x.FitResult.par(4),FT);
    NErr  = cellfun(@(x) x.FitResult.err(4),FT);
    % p-value
    pValue = cellfun(@(x) chi2pvalue(x.FitResult.chi2min,x.FitResult.dof),FT);
    Chi2   =  cellfun(@(x) x.FitResult.chi2min,FT);
    dof    = cellfun(@(x) x.FitResult.dof,FT);
    % RhoD
    rhoD     = cellfun(@(x) x.ModelObj.WGTS_CD_MolPerCm2,FT);
    rhoDErr  = cellfun(@(x) x.ModelObj.WGTS_CD_MolPerCm2*0,FT);
    % DT
    dt      = cellfun(@(x) x.ModelObj.WGTS_MolFrac_DT,FT);
    dtErr   = cellfun(@(x) x.ModelObj.WGTS_MolFracRelErr_DT,FT);
    
    % Labels
    r    = string(RunList);
    
    FitResults = struct('E0',E0,'E0Err',E0Err,'N',N,'NErr',NErr,'B',B,'BErr',BErr,...
        'pValue',pValue,'chi2min',Chi2,'dof',dof);
    x   = linspace(1,numel(RunList),numel(RunList));

%% Plot

switch displayHist
    case 'ON'
        index =  [1602  1402  1204  1002   804   602   402   304   202   177   154   129   104    92    82    72    62    52    45    32    22    15     2    -6   -18   -28];
        LabelSize = 16;
        if numel(RunList)>50 && numel(RunList)<80
            LabelSize = 10;
        elseif numel(RunList)>80
            LabelSize = 8;
        end
        maintitle=[sprintf('KATRIN First Tritium Samak Fit: Fit Parameter Evolution (%s) %s , %g',...
            RunObj.ringCutFlag,RunObj.chi2,index(RunObj.exclDataStart)),'eV below E0'];
        %% ALL Fit Parameters: overview
        savefile=sprintf('plots/tmp/KATRIN_FT_AllPixels_Samak_%s_%geVbelowE0-1_EvolutionOverview',RunObj.chi2,index(RunObj.exclDataStart));
        x   = linspace(1,numel(RunList),numel(RunList));
        fig1 = figure('Name','VFT Fits','Renderer','opengl','NumberTitle','off','rend','painters','pos',[10 10 1400 2000]);
        a=annotation('textbox', [0 0.91 1 0.1], ...
            'String', maintitle, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        a.FontSize=18;a.FontWeight='bold';
        
        s1 = subplot(3,1,1); % Endpoint
        bar(x,E0,'facecolor',rgb('CadetBlue'));
        hold on
        errorb(x,E0,E0Err,'ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1)
        hold off
        xticks(x)
        xticklabels(r)
        ylabel(sprintf('E_0-%.1f (eV)',FT{1}.ModelObj.Q_i))
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        switch FT{1}.chi2
            case 'chi2Stat'
                sp1 = sprintf('<E_0>=%.2f eV \\pm %.2f eV (std over runs) \\pm %.2f eV (standard error of mean)',...
                    FT{1}.ModelObj.Q_i+mean(E0),std(E0),mean(E0Err)/sqrt(numel(E0Err)));
            case {'chi2CM','chi2CMShape'}
                sp1 = sprintf('<E_0>=%.2f eV \\pm %.2f eV (std over runs)',...
                    FT{1}.ModelObj.Q_i+mean(E0),std(E0));
        end
        title(sp1);
        grid on
        set(gca,'FontSize',LabelSize)
        PrettyFigureFormat
        set(get(gca,'XAxis'),'FontSize',LabelSize);
        set(get(gca,'XLabel'),'FontSize',12);
        
        s2 = subplot(3,1,2); % Background
        bar(x,B*1e3,'facecolor',rgb('CadetBlue'));
        hold on
        errorb(x,B*1e3,BErr*1e3,'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1)
        hold off
        xticks(x)
        xticklabels(r)
        ylabel('Background (mcps)')
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        sp2 = sprintf('<B>=%.0f mcps \\pm %.0f mcps (std over runs) \\pm %.0f mcps (standard error of mean)',...
            mean(B*1e3),std(B*1e3),mean(BErr*1e3)/sqrt(numel(BErr)));title(sp2)
        grid on
        set(gca,'FontSize',LabelSize)
        PrettyFigureFormat
        set(get(gca,'XAxis'),'FontSize',LabelSize);
        set(get(gca,'XLabel'),'FontSize',12);
        
        s3 = subplot(3,1,3); % Normalization
        bar(x,N+1,'facecolor',rgb('CadetBlue'));
        hold on
        errorb(x,N+1,NErr,'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1)
        hold off
        xticks(x)
        xticklabels(r)
        ylabel('Normalization')
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        sp2 = sprintf('<N>=%.2f \\pm %.1e (std over runs) \\pm %.1e (standard error of mean)',...
            mean(N+1),std(N+1),mean(NErr)/sqrt(numel(NErr)));title(sp2)
        grid on
        set(gca,'FontSize',LabelSize);
        PrettyFigureFormat
        set(get(gca,'XAxis'),'FontSize',LabelSize);
        set(get(gca,'XLabel'),'FontSize',12);
        
        linkaxes([s1,s2,s3],'x');
        
        if strcmp(saveplot,'ON') || strcmp(saveplot,'pdf')  
        publish_figurePDF(gcf,[savefile,'.pdf']);
        else
            export_fig([savefile,sprintf('.%s',saveplot)]);
        end
        
        %% Parameter of Interest: Evolution with respect to mean and Chi2/p-values distribution
       if strcmp(FitResultsFlag,'ON')
        maintitle=[sprintf('KATRIN First Tritium Samak %s Fit: %s Evolution (%s) %s , %g',...
            strrep(RunObj.chi2,'chi2',''),Parameter,RunObj.ringCutFlag,RunObj.chi2,index(RunObj.exclDataStart)),'eV below E0'];
       else
            maintitle=[sprintf('KATRIN First Tritium %s Fit: %s Evolution (%g',...
            strrep(RunObj.chi2,'chi2',''),Parameter,index(RunObj.exclDataStart)),'eV below E0)'];
       end
        savefile=sprintf('plots/tmp//KATRIN_FT_AllPixels_Samak_%s_%geVbelowE0_%s_Evolution',RunObj.chi2,index(RunObj.exclDataStart),Parameter);
        x   = linspace(1,numel(RunList),numel(RunList));
        fig1 = figure('Name','VFT Fits','Renderer','opengl','pos',[10 10 1400 2000]);
        a=annotation('textbox', [0 0.88 1 0.1], ...
            'String', maintitle, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        a.FontSize=18;a.FontWeight='bold';
        
        subplot(2,2,[1 2])
        bar(x,FitResults.(Parameter)-mean(FitResults.(Parameter)),'facecolor',rgb('CadetBlue'));
        hold on
        errorb(x,FitResults.(Parameter)-mean(FitResults.(Parameter)),FitResults.([Parameter,'Err']),'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1)
        hold off
        xticks(x)
        xticklabels(r)
        ylabel(sprintf('%s-<%s> (eV)',Parameter,Parameter)); 
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        if strcmp(FitResultsFlag,'ON')
            switch FT{1}.chi2
                case 'chi2Stat'
        sp1 = sprintf('<%s>=%.2f eV \\pm %.2f eV (std over runs) \\pm %.2f eV (standard error of mean)',...
            Parameter,FT{1}.ModelObj.Q_i+mean(FitResults.(Parameter)),std(FitResults.(Parameter)),mean(FitResults.([Parameter,'Err']))/sqrt(numel(FitResults.(Parameter))));
                case 'chi2CM'
                   sp1 = sprintf('<%s>=%.2f eV \\pm %.2f eV (std over runs)',...
            Parameter,FT{1}.ModelObj.Q_i+mean(FitResults.(Parameter)),std(FitResults.(Parameter)));  
            end
        else
            sp1 = '';
        end
        title(sp1);
        grid on
        set(gca,'FontSize',LabelSize)
        PrettyFigureFormat
        set(get(gca,'XAxis'),'FontSize',LabelSize);
        set(get(gca,'XLabel'),'FontSize',20);
        set(get(gca,'YAxis'),'FontSize',20);
        
        subplot(2,2,3)
        h=histfit(Chi2,8,'gamma');
        h(1).FaceColor = rgb('CadetBlue');
        ylabel('runs')
        xlabel('\chi^2'); 
        sp3 = sprintf('<p-value>=%.2f \\pm %.2f (std)',mean(pValue),std(pValue));title(sp3)
        grid on   
        PrettyFigureFormat
        set(gca,'FontSize',20);
        
        subplot(2,2,4)
        h=histfit((FitResults.(Parameter)-mean(FitResults.(Parameter)))./FitResults.([Parameter,'Err']),5,'normal');
        ParStd = std(FitResults.(Parameter)-mean(FitResults.(Parameter)));
        EoLeg = sprintf('\\sigma = %.2f eV',ParStd);
        h(1).FaceColor = rgb('CadetBlue');
        ylabel('runs') 
        xlabel([sprintf('%s-<%s> / ',Parameter,Parameter) ,'\sigma']); 
        grid on
        PrettyFigureFormat
        set(gca,'FontSize',20);
        title(EoLeg);
        
        if strcmp(saveplot,'ON') || strcmp(saveplot,'pdf')
            publish_figurePDF(gcf,[savefile,'.pdf']);
        else
            export_fig([savefile,sprintf('.%s',saveplot)]);
        end
        
        %% Build pdf file with all fits 
        switch mergePDF
            case 'ON'
                PATH = getenv('PATH');
                setenv('PATH', [PATH ':/usr/local/bin/']);
                cd plots/tmp
                command = 'gs -sDEVICE=pdfwrite -sOutputFile="run.pdf" -dNOPAUSE -dEPSCrop -c "<</Orientation 0>> setpagedevice" -f KATRIN_FT_AllPixels_Samak_*.pdf -c quit';
                unix(command);
                unix('rm KATRIN_FT_AllPixels_Samak_*.pdf');
                mycommand1 = sprintf('mv run.pdf ../FT_FitParameterEvolution%s-%s-%ubelowE0-%s-samak.pdf',Parameter,RunObj.chi2,index(RunObj.exclDataStart),datestr(now,'dd-mmm-yyyy'));
                unix(mycommand1);
                mycommand2 = sprintf('open ../FT_FitParameterEvolution%s-%s-%ueVbelowE0-%s-samak.pdf',Parameter,RunObj.chi2,index(RunObj.exclDataStart),datestr(now,'dd-mmm-yyyy'));
                unix(mycommand2);
                cd ../..
                unix('rm -rf plots/tmp');
        end
end
end