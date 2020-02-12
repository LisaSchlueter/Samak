Mode = 'Real';
RunList = 'KNM1';
chi2 = 'chi2Stat';
RecomputeFlag = 'OFF';

if strcmp(RunList,'KNM1')
    WGTS_CD_MolPerCm2_V = 0.75:0.02:1.25;
    exclDataStartAll = 2:1:11;
    ELossFlag = 'KatrinT2';
    FSDFlag = 'SibilleFull';
    NonPoissonScaleFactor = 1.064;
    fixPar = '5 6 7 8 9 10 11';
elseif strcmp(RunList,'FTpaper')
    WGTS_CD_MolPerCm2_V = 0.85:0.01:1.15;
    exclDataStartAll = 1:1:11;
    ELossFlag = 'Abdurashitov';
    FSDFlag = 'SAENZ';
    NonPoissonScaleFactor = 1.0;
    fixPar = '1 5 6 7 8 9 10 11';
end

nIt = numel(exclDataStartAll);

switch chi2
    case {'chi2CM','chi2CMShape'}
        SysEffects            = struct('TASR','ON','FSD','ON','RF_RX','OFF','RF_EL','ON','RF_BF','ON','BkgShape','ON','TCoff_RAD','ON','TCoff_OTHER','ON','FPDeff','ON','Stack','ON');
        BkgCM                 = 'ON';
end

ModelArg = {'RunList',RunList,'chi2',chi2,...
    'minuitOpt','min  ; migrad ',...
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'FSDFlag',FSDFlag,...
    'ELossFlag',ELossFlag,...
    'SysBudget',22,...
    'fixPar',fixPar};

        switch Mode
            case 'Real'
                A = MultiRunAnalysis('DataType','Real',ModelArg{:});
            case 'Twin'
                A = MultiRunAnalysis('DataType','Twin',ModelArg{:});
            case 'Twin_CDsame'
                A = MultiRunAnalysis('DataType','Twin','Twin_SameCDFlag','ON',ModelArg{:});
        end
        
        %save/load
        switch RunList
            case 'KNM1'
                savedir = [getenv('SamakPath'),'knm1ana/knm1_ColumnDensity/results/'];
            case 'FTpaper'
                savedir = [getenv('SamakPath'),'first-tritium-studies/ft_ColumnDensity/results/'];
        end
        MakeDir(savedir);
        
        savefile = [savedir,sprintf('RDScan_%s_%s_%s_Fit%.0f-%.0feV_RD%.1f-%.2f-%.1f.mat',RunList,Mode,chi2,18575-A.ModelObj.qU(exclDataStartAll(1)),...
            18575-A.ModelObj.qU(exclDataStartAll(end)),min(WGTS_CD_MolPerCm2_V),WGTS_CD_MolPerCm2_V(2)-WGTS_CD_MolPerCm2_V(1),max(WGTS_CD_MolPerCm2_V))];
        
        if exist(savefile,'file') && strcmp(RecomputeFlag,'OFF')
            load(savefile,'CD_bestfit','CDErr_bestfit_up','CDErr_bestfit_low','chi2RD','parRD');
        else
            
            % init
            CD_bestfit        = zeros(nIt,1); %Fit Value for column density
            CDErr_bestfit_up  = zeros(nIt,1);
            CDErr_bestfit_low = zeros(nIt,1);
            chi2RD            = zeros(nIt,numel(WGTS_CD_MolPerCm2_V));
            parRD             = zeros(nIt,11,numel(WGTS_CD_MolPerCm2_V));
            errRD             = zeros(nIt,11,numel(WGTS_CD_MolPerCm2_V));
            
            progressbar('column density scans ..')
            for i=1:nIt
                progressbar(i/nIt)
                A.exclDataStart = exclDataStartAll(i);
                [parRD(i,:,:), errRD(i,:,:), chi2RD(i,:), ~, CD_bestfit(i)] = ...
                    A.RhoDScan('WGTS_CD_MolPerCm2_V',WGTS_CD_MolPerCm2_V,'saveplot','OFF');
                close all;
                CDErr_bestfit_low(i) = A.rhoDlowerUnc;
                CDErr_bestfit_up(i) = A.rhoDupperUnc;
            end
            save(savefile,'CD_bestfit','CDErr_bestfit_up','CDErr_bestfit_low','chi2RD','parRD','errRD');
        end
        
        %%
        chi2Both = 'ON';
        if strcmp(chi2Both,'ON')
            d = importdata(strrep(savefile,'chi2Stat','chi2CMShape'));
        end
        close all
        fig12345 = figure(12345);
        A.GetPlotColor;
        set(fig12345, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.7]);
        %round(A.ModelObj.qU(exclDataStartAll)-18575)
        x = A.ModelObj.qU(exclDataStartAll)-18575;
        [l1,a1] = boundedline(linspace(min(x)-2,max(x)+4,10),100.*ones(10,1),5.*ones(10,1));%,...
        l1.LineStyle = '--'; l1.Color = rgb('DimGray'); l1.LineWidth=3;
        a1.FaceColor = rgb('DimGray'); a1.FaceAlpha=0.3;
        hold on;
        switch A.DataType
            case 'Real'
                PlotColor2 = rgb('DarkCyan');
            case 'Twin'
               PlotColor2 = rgb('IndianRed');
        end
        e1 = errorbar(x,100.*CD_bestfit./A.ModelObj.WGTS_CD_MolPerCm2,...
            100.*(CD_bestfit./A.ModelObj.WGTS_CD_MolPerCm2-CDErr_bestfit_low./A.ModelObj.WGTS_CD_MolPerCm2),...
            100.*(CD_bestfit./A.ModelObj.WGTS_CD_MolPerCm2-CDErr_bestfit_up./A.ModelObj.WGTS_CD_MolPerCm2),...
            'o','MarkerSize',10,'MarkerFaceColor',A.PlotColor,...
            'Color',PlotColor2,'LineWidth',3);
        hold on;
        xmax = max(x);
        e1.CapSize = 0;
        if strcmp(chi2Both,'ON')
            e2 = errorbar(x+1,100.*d.CD_bestfit./A.ModelObj.WGTS_CD_MolPerCm2,...
                100.*(d.CD_bestfit./A.ModelObj.WGTS_CD_MolPerCm2-d.CDErr_bestfit_low./A.ModelObj.WGTS_CD_MolPerCm2),...
                100.*(d.CD_bestfit./A.ModelObj.WGTS_CD_MolPerCm2-d.CDErr_bestfit_up./A.ModelObj.WGTS_CD_MolPerCm2),...
                'o','MarkerSize',10,'MarkerFaceColor',rgb('DarkGoldenRod'),...
                'Color',rgb('GoldenRod'),'LineWidth',3);
            xmax = max(x)+1;
            e2.CapSize=0;
        end
        hold off;
        PrettyFigureFormat;
        xlabel('Lower fit boundary below E_0 (eV)');
        ylabel('Column density (%)');
        switch chi2
            case 'chi2Stat'
                leg = legend([l1,a1,e1],sprintf('Expected value = %.3g mol/cm^2',A.ModelObj.WGTS_CD_MolPerCm2),...
                    '5% Uncertainty band','Column density fit (stat. only)','Location','southwest');
            case {'chi2CM','chi2CMShape'}
                leg = legend([l1,a1,e1],sprintf('Expected value = %.3g mol/cm^2',A.ModelObj.WGTS_CD_MolPerCm2),...
                    '5% Uncertainty band','Column density fit (stat+sys)','Location','southwest');
        end
        
        if strcmp(chi2Both,'ON')
            leg = legend([l1,a1,e1,e2],sprintf('expected value = %.3g mol/cm^2',A.ModelObj.WGTS_CD_MolPerCm2),...
                '5% uncertainty band','column density fit (stat)','column density fit (stat+sys)',...
                'Location','southwest');
        end
        legend boxoff
        set(gca,'FontSize',20);
        %xlim([x(2) xmax])
        xlim([min(x)-2 xmax+2])
        %grid on;
        
        % save
        if ~exist(strrep(savedir,'results','plots'),'dir')
            system(['mkdir ',strrep(savedir,'results','plots')]);
        end
        
        savename_plot = strrep(strrep(savefile,'results','plots'),'.mat','.png');
        if strcmp(chi2Both,'ON')
            savename_plot = strrep(strrep(savename_plot,'.png','_Both.png'),[chi2,'_'],'');
        end
        
        print(fig12345,savename_plot,'-dpng','-r450');
       % publish_figurePDF(fig12345,strrep(savename_plot,'.png','.pdf'));

