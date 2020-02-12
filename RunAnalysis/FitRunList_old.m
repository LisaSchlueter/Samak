function FitRunList(varargin)
close all
unix('rm -rf plots/tmp');
mkdir  plots/tmp

% Parser
p = inputParser;
p.addParameter('RunList',0);
p.addParameter('chi2','chi2Stat',@(x)ismember(x,{'chi2Stat', 'chi2CM', 'chi2CMFrac','chi2CMShape', 'chi2P','chi2Nfix'}));
p.addParameter('fitter','minuit',@(x)ismember(x,{'minuit','matlab'}));
p.addParameter('exclDataStart',1,@(x)isfloat(x));
p.addParameter('displayHist','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('displayFit','OFF',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});

RunList       = p.Results.RunList;
chi2          = p.Results.chi2;
fitter        = p.Results.fitter;
exclDataStart = p.Results.exclDataStart;
displayHist   = p.Results.displayHist;
displayFit    = p.Results.displayFit;


index =  [1602  1402  1204  1002   804   602   402   304   202   177   154   129   104    92    82    72    62    52    45    32    22    15     2    -6   -18   -28];

%% 
if RunList == 0
[ n, datafiles , RunList ] = GetRunList( '../../tritium-data/hdf5/','*.h5',1,'string');
end

%% Fit All RunList
for i=1:numel(RunList)
    close all
    FT(i) = RunAnalysis('RunNr',RunList(i),'fitter',fitter,'chi2',chi2,'exclDataStart',exclDataStart);
    FT(i).Fit; 
    switch displayFit
        case 'ON'
            FT(i).PlotFit('saveplot','ON');
    end
%      FT(i).RhoDScan('saveplot','ON');
%      close all
%      FT(i).qUScan('saveplot','ON','exclnPoints',17);
%      close all
end

%% Gather Results
T = []; E0 = [];E0Err=[];B=[];BErr=[];N=[];NErr=[];Chi2=[];
r='';rhoD = []; rhoDErr = []; dt=[]; dtErr=[];
for i=1:numel(RunList)
    % Endpoint
    T     = [T    FT(i).ModelObj.TimeSec];
    TC    = cumsum(T);
    % Endpoint
    E0    = [E0    FT(i).FitResult.par(2)];
    E0Err = [E0Err FT(i).FitResult.err(2)];
    % Background
    B    = [B    FT(i).ModelObj.BKG_RateAllFPDSec+FT(i).FitResult.par(3)];
    BErr = [BErr FT(i).FitResult.err(3)];
    % Normalizationn
    N    = [N    FT(i).FitResult.par(4)];
    NErr = [NErr FT(i).FitResult.err(4)];
    % Chi2
    Chi2 = [Chi2 chi2pvalue(FT(i).FitResult.chi2min,FT(i).FitResult.dof)];
    % RhoD
    rhoD    = [rhoD FT(i).ModelObj.WGTS_CD_MolPerCm2];
    rhoDErr = [rhoDErr FT(i).ModelObj.WGTS_CD_MolPerCm2*0];
    % DT
    dt      = [dt FT(i).ModelObj.WGTS_MolFrac_DT];
    dtErr   = [dtErr FT(i).ModelObj.WGTS_MolFracRelErr_DT];
    % Labels
    r    = [r num2str(RunList(i))];
    r    = string(r);
end
x   = linspace(1,numel(RunList),numel(RunList));

%% Plot
switch displayHist
    case 'ON'
        
        %% Fit Parameters
        maintitle=sprintf('KATRIN  First Tritium - All pixels but 2 ext Rings - Samak Fit: %s , %g eV below E_0',chi2,index(exclDataStart));
        savefile=sprintf('plots/tmp//KATRIN_FT_AllPixels_Samak_%s_%geVbelowE0-1.pdf',chi2,index(exclDataStart));
        x   = linspace(1,numel(RunList),numel(RunList));
        fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
        a=annotation('textbox', [0 0.91 1 0.1], ...
            'String', maintitle, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';
        subplot(2,1,1)
        bar(x,E0,'facecolor',rgb('CadetBlue'));
        hold on
        errorb(x,E0,E0Err,'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1)
        hold off
        xticks(x)
        xticklabels(r)
        ylabel('E_0-18575 (eV)')
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        sp1 = sprintf('<E0>=%.2f eV \\pm %.2f eV (std)',mean(E0),std(E0-18575));title(sp1)
        grid on
        set(gca,'FontSize',16);
        PrettyFigureFormat
        subplot(2,1,2)
        bar(x,B,'facecolor',rgb('CadetBlue'));
        hold on
        errorb(x,B,BErr,'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1)
        hold off
        xticks(x)
        xticklabels(r)
        ylabel('Background (cps)')
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        sp2 = sprintf('<B>=%.2f cps \\pm %.2f cps (std)',mean(B),std(B));title(sp2)
        grid on
        set(gca,'FontSize',16);
        PrettyFigureFormat
        PrettyFigureFormat
        publish_figurePDF(gcf,savefile);
        
        maintitle=sprintf('KATRIN  First Tritium - All pixels but 2 ext Rings - Samak Fit: %s , %g eV below E_0',chi2,index(exclDataStart));
        savefile=sprintf('plots/tmp//KATRIN_FT_AllPixels_Samak_%s_%geVbelowE0-2.pdf',chi2,index(exclDataStart));
        x   = linspace(1,numel(RunList),numel(RunList));
        fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
        a=annotation('textbox', [0 0.91 1 0.1], ...
            'String', maintitle, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';
        subplot(2,1,1)
        bar(x,N,'facecolor',rgb('CadetBlue'));
        hold on
        errorb(x,N,NErr,'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1)
        hold off
        xticks(x)
        xticklabels(r)
        ylabel('1+Normalization')
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        sp3 = sprintf('<1+N>=%.2f \\pm %.2f  (std)',mean(N),std(N));title(sp3)
        grid on
        set(gca,'FontSize',16);
        PrettyFigureFormat
        subplot(2,1,2)
        bar(x,Chi2,'facecolor',rgb('CadetBlue'))
        hold on
        line([min(x)-0.5,max(x)+0.5],[0.05 0.05],'Color','Red','LineWidth',1,'LineStyle','--')
        hold off
        xticks(x)
        xticklabels(r)
        ylabel('p-value')
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        sp3 = sprintf('<p-value>=%.2f \\pm %.2f (std)',mean(Chi2),std(Chi2));title(sp3)
        set(gca,'yscale','log');
        grid on
        set(gca,'FontSize',16);
        PrettyFigureFormat
        publish_figurePDF(gcf,savefile);
        
        %% Slow Control Parameters
        maintitle=sprintf('KATRIN  First Tritium - Slow Control Parameters');
        savefile=sprintf('plots/tmp//KATRIN_FT_AllPixels_Samak_%s_%geVbelowE0_3.pdf',chi2,index(exclDataStart));
        x   = linspace(1,numel(RunList),numel(RunList));
        fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
        a=annotation('textbox', [0 0.91 1 0.1], ...
            'String', maintitle, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';
        subplot(2,1,1)
        errorbar(x,rhoD./mean(rhoD),rhoDErr./mean(rhoD),'ks','MarkerSize',12,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        xticks(x)
        xticklabels(r)
        ylabel('\rho d / <\rho d>')
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        sp1 = sprintf('<\rho d>=%g mol/cm^2 \\pm %g eV (std)',mean(rhoD),std(rhoD));title(sp1)
        grid on
        set(gca,'FontSize',16);
        PrettyFigureFormat
        subplot(2,1,2)
        errorbar(x,dt./mean(dt),dtErr,'ks','MarkerSize',12,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        hold off
        xticks(x)
        xticklabels(r)
        ylabel('[DT]/<[DT]>')
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        sp2 = sprintf('<DT>=%g  \\pm %.g (std)',mean(dt),std(dt));title(sp2)
        grid on
        set(gca,'FontSize',16);
        PrettyFigureFormat
        publish_figurePDF(gcf,savefile);
        
         %% Data Taking Time Parameters
        maintitle=sprintf('KATRIN  First Tritium - Time Control Parameters');
        savefile=sprintf('plots/tmp//KATRIN_FT_AllPixels_Samak_%s_%geVbelowE0_4.pdf',chi2,index(exclDataStart));
        x   = linspace(1,numel(RunList),numel(RunList));
        fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
        a=annotation('textbox', [0 0.91 1 0.1], ...
            'String', maintitle, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';
        subplot(2,1,1)
        bar(x,T,'facecolor',rgb('CadetBlue'));
        xticks(x)
        xticklabels(r)
        ylabel('Time (sec)')
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        sp1 = sprintf('<Run Time>=%g seconds \\pm %g eV (std)',mean(T),std(T));title(sp1)
        grid on
        set(gca,'FontSize',16);
        PrettyFigureFormat
        subplot(2,1,2)
        bar(x,TC./86400,'facecolor',rgb('CadetBlue'));
        hold off
        xticks(x)
        xticklabels(r)
        ylabel('Cumulative Time (day)')
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        sp2 = sprintf('Total Time %g days',TC(end)./86400);title(sp2)
        grid on
        set(gca,'FontSize',16);
        PrettyFigureFormat
        publish_figurePDF(gcf,savefile);
        
        %% Fit Parameters Verus rhoD
        maintitle=sprintf('KATRIN  First Tritium - All pixels but 2 ext Rings - Samak Fit: %s , %g eV below E_0',chi2,index(exclDataStart));
        savefile=sprintf('plots/tmp//KATRIN_FT_AllPixels_Samak_%s_%geVbelowE0_5.pdf',chi2,index(exclDataStart));
        x   = linspace(1,numel(RunList),numel(RunList));
        fig1 = figure('Name','VFT Fits Verus Column Density','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
        a=annotation('textbox', [0 0.91 1 0.1], ...
            'String', maintitle, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';
        subplot(2,2,1)
        errorbar(rhoD,E0,E0Err,'ks','MarkerSize',12,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        ylabel('E_0-18575 (eV)')
        xlabel('\rho d (mol/cm^2)');xtickangle(45);
        PrettyFigureFormat
        subplot(2,2,2)
        errorbar(rhoD,B,BErr,'ks','MarkerSize',12,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        ylabel('B (cps)')
        xlabel('\rho d (mol/cm^2)');xtickangle(45);
        PrettyFigureFormat
        subplot(2,2,3)
        errorbar(rhoD,N,NErr,'ks','MarkerSize',12,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        ylabel('N')
        xlabel('\rho d (mol/cm^2)');xtickangle(45);
        PrettyFigureFormat
        subplot(2,2,4)
        errorbar(rhoD,Chi2,Chi2.*0,'ks','MarkerSize',12,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        ylabel('\chi^2')
        xlabel('\rho d (mol/cm^2)');xtickangle(45);
        set(gca,'yscale','log');
        PrettyFigureFormat
        publish_figurePDF(gcf,savefile);
        
        %% Fit Parameters Verus DT
        maintitle=sprintf('KATRIN  First Tritium - All pixels but 2 ext Rings - Samak Fit: %s , %g eV below E_0',chi2,index(exclDataStart));
        savefile=sprintf('plots/tmp//KATRIN_FT_AllPixels_Samak_%s_%geVbelowE0_6.pdf',chi2,index(exclDataStart));
        x   = linspace(1,numel(RunList),numel(RunList));
        fig1 = figure('Name','VFT Fits Verus [DT]','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
        a=annotation('textbox', [0 0.91 1 0.1], ...
            'String', maintitle, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';
        subplot(2,2,1)
        errorbar(dt,E0,E0Err,'ks','MarkerSize',12,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        ylabel('E_0-18575 (eV)')
        xlabel('[DT]');xtickangle(45);
        PrettyFigureFormat
        subplot(2,2,2)
        errorbar(dt,B,BErr,'ks','MarkerSize',12,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        ylabel('B (cps)')
        xlabel('[DT]');xtickangle(45);
        PrettyFigureFormat
        subplot(2,2,3)
        errorbar(dt,N,NErr,'ks','MarkerSize',12,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        ylabel('N')
        xlabel('[DT]');xtickangle(45);
        PrettyFigureFormat
        subplot(2,2,4)
        errorbar(dt,Chi2,Chi2.*0,'ks','MarkerSize',12,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        ylabel('\chi^2')
        xlabel('[DT]');xtickangle(45);
        set(gca,'yscale','log');
        PrettyFigureFormat
        publish_figurePDF(gcf,savefile);
   
%% Build pdf file with all fits
PATH = getenv('PATH');
setenv('PATH', [PATH ':/usr/local/bin/']);
cd plots/tmp
command = 'gs -sDEVICE=pdfwrite -sOutputFile="run.pdf" -dNOPAUSE -dEPSCrop -c "<</Orientation 0>> setpagedevice" -f KATRIN_FT_AllPixels_Samak_*.pdf -c quit';
unix(command);
unix('rm KATRIN_FT_AllPixels_Samak_*.pdf');
mycommand1 = sprintf('mv run.pdf ../FTruns-%s-exclDataStart%g-%s-samak.pdf',chi2,exclDataStart,datestr(now,'dd-mmm-yyyy'));
unix(mycommand1);
mycommand2 = sprintf('open ../FTruns-%s-exclDataStart%g-%s-samak.pdf',chi2,exclDataStart,datestr(now,'dd-mmm-yyyy'));
unix(mycommand2);
cd ../..
unix('rm -rf plots/tmp');

end
end
