%
% Samak KATRIN Simulation and Analysis
% 
% Display Slow Control Parameters extracted from RS
% Input:
% - RunList:      0 for all, or a vector like [run1 run2 ...]
% - displayHist:  display histograms ON/OFF
% - PlotExt    :  export plot in png or pdf
% Output:
% - Excel Table
% 
% T. Lasserre, Pablo Morales, Lisa Schlueter
% Last Update: June 24 2018
%

function TexTable = DisplaySlowControlRS(varargin)
close all
unix('rm -rf plots/tmp');
mkdir  plots/tmp

% Parser
p = inputParser;
p.addParameter('RunList',0);
p.addParameter('displayHist','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('PlotExt','png',@(x)ismember(x,{'png','pdf'}));

p.parse(varargin{:});

RunList       = p.Results.RunList;
displayHist   = p.Results.displayHist;
PlotExt       = p.Results.PlotExt;

%% 
if RunList == 0
%[ n, ~ , RunList ] = GetRunList( '../tritium-data/hdf5/','2f-fpd*.h5',1,char('2f-fpd0040649.h5','40769.h5','40770.h5','40771.h5','40772.h5'));
[ n, ~ , RunList ] = GetRunList( '../../tritium-data/mat/','*ex2.mat',1,'string');
RunList=RunList(RunList>40538);             % First Good Run FT
RunList=RunList(RunList<40693);            % 40538-40693 List of 100%CD 3h runs
Remove1=RunList==40773;RunList(Remove1)=[]; % Sterile Neutrino Run - High Count Rate
Remove2=RunList==40806;RunList(Remove2)=[]; % Sterile Neutrino Run - High Count Rate
Remove3=RunList==40995;RunList(Remove3)=[]; % Sterile Neutrino Run - High Count Rate
Remove4=RunList==40761;RunList(Remove4)=[]; % Stability Run - No Scan
Remove5=RunList==40762;RunList(Remove5)=[]; % Stability Run - No Scan
Remove6=RunList==40807;RunList(Remove6)=[]; % Anomaleous DT uncertainty
Remove7=RunList==40769;RunList(Remove7)=[]; % Sterile Neutrinos? 31 qU
Remove7=RunList==40770;RunList(Remove7)=[]; % Sterile Neutrinos? 31 qU
Remove8=RunList==40805;RunList(Remove8)=[]; % Sterile Neutrinos? 18 qU
Remove9=RunList==40649;RunList(Remove9)=[]; % Firmware Test Run
RunList=RunList(RunList<41033);             % above are fading out scans (41035) or deep scans (41038,41039,41045:41056)  
RunList=[40531  40540  40541  40542  40543  40603  40604  40611  40613  40667  40668  40669  40670  40671  40672  40673  40674  40675  40676  40677  40678  40679  40680  40681  40682  40683  40684  40685  40686  40687  40688  40689  40690  40691  40692  40693  40976  40977  40979  40980  40982  40983  40985  40986  40988  40989  40991  40992  40994  40995  40997  41002  41005  41007  41008  41010  41011  41013  41014  41016  41017  41019  41020  41022  41023  41025  41026  41028  41029  41031  41032];
RunList=[40667  40668  40669  40670  40671  40672  40673  40674  40675  40676  40677  40678  40679  40680  40681  40682  40683  40684  40685  40686  40687  40688  40689  40690  40691  40692  40693];
RunList=[...% 100% column density - up scans (== from -1600eV to -1800eV)
                40541, 40543, 40604, 40611, 40613, 40667, 40669, 40671, 40673, ...
                40675, 40677, 40679, 40681, 40683, 40685, 40687, 40689, 40691, 40693, ...
                40977, 40980, 40983, 40986, 40989, 40992, 40995, 41002, 41005, 41008, ...
                41011, 41014, 41019, 41022,  41026,  41029, 41032 ...
                40531, 40540, 40542, 40603, 40668, 40670, 40672, 40674, ...
                40676, 40678, 40680, 40682, 40684, 40686, 40688, 40690, 40692, 40976, ...
                40979, 40982, 40985, 40988, 40991, 40994, 40997, 41007, 41010, 41013, ...
                41016, 41017, 41020, 41023, 41025,  41028, 41031];
end
disp(RunList);

LabelSize = 12;
if numel(RunList)>50 && numel(RunList)<80
    LabelSize = 10;
elseif   numel(RunList)>80
    LabelSize = 8;
end
%% Fit All RunList
for i=1:numel(RunList)
    fprintf(2,'processing run %g ...\n',RunList(i))
    close all
    FT(i) = RunAnalysis('RunNr',RunList(i));
end

%% Gather Results
T = []; r='';rhoD = []; rhoDErr = []; tbdis=[];
dt=[]; dtErr=[];tt=[]; ttErr=[]; ht=[]; htErr=[];

for i=1:numel(RunList)
    % Time 
    T     = [T    FT(i).RunData.TimeSec];
    TC    = cumsum(T);
    % TBDIS
    tbdis = [tbdis sum(FT(i).RunData.TBDIS)];
    % RhoD
    rhoD    = [rhoD FT(i).ModelObj.WGTS_CD_MolPerCm2];
    rhoDErr = [rhoDErr std(FT(i).RunData.WGTS_CD_MolPerCm2_SubRun)];
    % DT
    dt      = [dt FT(i).ModelObj.WGTS_MolFrac_DT];
    dtErr   = [dtErr mean(FT(i).RunData.WGTS_MolFrac_DT_SubRun_error_mean)];
    % TT
    tt      = [tt FT(i).ModelObj.WGTS_MolFrac_TT];
    ttErr   = [ttErr mean(FT(i).RunData.WGTS_MolFrac_TT_SubRun_error_mean)];
    % HT
    ht      = [ht FT(i).ModelObj.WGTS_MolFrac_HT];
    htErr   = [htErr mean(FT(i).RunData.WGTS_MolFrac_HT_SubRun_error_mean)];
    % Labels
    r    = [r num2str(RunList(i))];
    r    = string(r);
   
end
    % Latex Table
    unix('touch FT_SC_Table.xls; rm -rf FT_SC_Table.xls');
    TexTable = table(r',T',rhoD',rhoDErr'./rhoD'*100,dt',mean(dtErr,1)'*100,tt',mean(ttErr,1)'*100,ht',mean(htErr,1)'*100,tbdis',tbdis'./T','VariableNames',...
        {'run','seconds ','rhoD','rhoDerrorRel','DT','DTerrorRel','TT','TTerror','HT','HTerrorRel','Counts','Rate'});
    writetable(TexTable,'FT_SC_Table.xls');
    x   = linspace(1,numel(RunList),numel(RunList));

%% Time / Number of electrons
%maintitle=sprintf('KATRIN  First Tritium - Time / Events Control Parameters');
maintitle=sprintf('KATRIN  First Tritium - Time / Events Control Parameters');
        savefile=sprintf('plots/tmp/KATRIN_FT_SC_Time_Events');
        x   = linspace(1,numel(RunList),numel(RunList));
        fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
        a=annotation('textbox', [0 0.91 1 0.1], ...
            'String', maintitle, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';
        s1 = subplot(3,1,1);
        bar(x,tbdis,'facecolor',rgb('CadetBlue'));
        xticks(x)
        xticklabels(r)
        ylabel('Electrons (sum)')
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        sp1 = sprintf('<Electrons>=%g seconds \\pm %g  (std)',mean(tbdis),std(tbdis));title(sp1)
        grid on
        set(gca,'FontSize',LabelSize)
        PrettyFigureFormat
        set(get(gca,'XAxis'),'FontSize',LabelSize);
        set(get(gca,'XLabel'),'FontSize',12);
        
        s2=  subplot(3,1,2);
        bar(x,TC./86400,'facecolor',rgb('CadetBlue'));
        hold off
        xticks(x)
        xticklabels(r)
        ylabel('Cumulative Time (day)')
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        sp2 = sprintf('Total Time %g days',TC(end)./86400);title(sp2)
        grid on
        set(gca,'FontSize',LabelSize)
        PrettyFigureFormat
        set(get(gca,'XAxis'),'FontSize',LabelSize);
        set(get(gca,'XLabel'),'FontSize',12);
        
        s3=  subplot(3,1,3);
        bar(x,T./3600,'facecolor',rgb('CadetBlue'));
        hold off
        xticks(x)
        xticklabels(r)
        ylabel('Time (hour)')
        xlabel('run');xtickangle(45);
        xlim([0.5 numel(RunList)+0.5])
        grid on
        set(gca,'FontSize',LabelSize)
        PrettyFigureFormat
        set(get(gca,'XAxis'),'FontSize',LabelSize);
        set(get(gca,'XLabel'),'FontSize',12);
        linkaxes([s1,s2, s3],'x');
       
        switch PlotExt
            case 'pdf'
        publish_figurePDF(gcf,[savefile '.pdf']);
            case 'png'
        saveas(gcf,[savefile '.png']);
        end
        
%% Column Density
maintitle=sprintf('KATRIN  First Tritium - Slow Control Parameters');
savefile=sprintf('plots/tmp/KATRIN_FT_SC_rhoD1');
x   = linspace(1,numel(RunList),numel(RunList));
fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 900]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

s1 = subplot(2,1,1);
errorbar(x,rhoD,rhoDErr,'ks','MarkerSize',12,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
xticks(x)
xticklabels(r)
ylabel('\rho d')
xlabel('run');xtickangle(45);
xlim([0.5 numel(RunList)+0.5])
sp1 = sprintf('<\\rho d>=%g mol/cm^2 \\pm %g eV (std)',mean(rhoD),std(rhoD));title(sp1)
grid on
PrettyFigureFormat
set(get(gca,'XAxis'),'FontSize',LabelSize);
set(get(gca,'XLabel'),'FontSize',12);

s2 = subplot(2,1,2);
errorbar(x,rhoDErr./rhoD*100,rhoDErr*0,'ks','MarkerSize',12,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
hold off
xticks(x)
xticklabels(r)
ylabel('\rho d rel. uncertainty (%)')
xlabel('run');xtickangle(45);
xlim([0.5 numel(RunList)+0.5])
sp2 = sprintf('<\\rho d uncertainty>=%g  \\pm %.g (std)\n',mean(rhoDErr),std(rhoDErr));
sp22 =  sprintf('<\\rho d rel uncertainty>=%.3f %% \\pm %.3f %% (std)',mean(rhoDErr)/mean(rhoD)*100,std(rhoDErr)/mean(rhoD)*100);
title([sp2,sp22]);
grid on
PrettyFigureFormat
set(get(gca,'XAxis'),'FontSize',LabelSize);
set(get(gca,'XLabel'),'FontSize',12);
 linkaxes([s1,s2],'x');
 
switch PlotExt
            case 'pdf'
        publish_figurePDF(gcf,[savefile '.pdf']);
            case 'png'
        saveas(gcf,[savefile '.png']);
end
        
% Column Density Distributions
maintitle=sprintf('KATRIN  First Tritium - Slow Control Parameters');
savefile=sprintf('plots/tmp/KATRIN_FT_SC_rhoD2');
x   = linspace(1,numel(RunList),numel(RunList));
fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

subplot(2,1,1);
histogram(rhoD,'Facecolor',rgb('CadetBlue'),'FaceAlpha',.5,'EdgeColor','none');
ylabel('runs')
xlabel('\rho d');xtickangle(45);
sp1 = sprintf('<\\rho d>=%g mol/cm^2 \\pm %g eV (std)',mean(rhoD),std(rhoD));title(sp1)
grid on
set(gca,'FontSize',16);
PrettyFigureFormat

subplot(2,1,2);
histogram(rhoDErr./rhoD,'Facecolor',rgb('CadetBlue'),'FaceAlpha',.5,'EdgeColor','none');
hold off
xlabel('\rho d rel. uncertainty (%)');xtickangle(45);
sp2 = sprintf('<\\rho d rel. uncertainty>=%g  \\pm %.g (std)',mean(rhoDErr./rhoD),std(rhoDErr./rhoD));title(sp2)
grid on
set(gca,'FontSize',16);
PrettyFigureFormat

switch PlotExt
            case 'pdf'
        publish_figurePDF(gcf,[savefile '.pdf']);
            case 'png'
        saveas(gcf,[savefile '.png']);
end
        
% DT
maintitle=sprintf('KATRIN  First Tritium - Slow Control Parameters');
savefile=sprintf('plots/tmp/KATRIN_FT_SC_DT1');
x   = linspace(1,numel(RunList),numel(RunList));
fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

s1 = subplot(2,1,1);
errorbar(x,dt,mean(dtErr,1),'ks','MarkerSize',12,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
xticks(x)
xticklabels(r)
ylabel('[DT]')
xlabel('run');xtickangle(45);
xlim([0.5 numel(RunList)+0.5])
sp1 = sprintf('<[DT]>=%g \\pm %g  (std)',mean(dt),std(dt));title(sp1)
grid on
set(gca,'FontSize',16);
PrettyFigureFormat
set(get(gca,'XAxis'),'FontSize',LabelSize);
set(get(gca,'XLabel'),'FontSize',12);

s2 = subplot(2,1,2);
errorbar(x,mean(dtErr,1)./mean(dt)*100,mean(dtErr,1).*0,'ks','MarkerSize',12,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
hold off
xticks(x)
xticklabels(r)
ylabel('[DT] rel. uncertainty (%)>')
xlabel('run');xtickangle(45);
xlim([0.5 numel(RunList)+0.5])
sp2 = sprintf('<DT>=%g  \\pm %g (std)',mean(dt),std(dt));title(sp2)
grid on
set(gca,'FontSize',16);
PrettyFigureFormat
set(get(gca,'XAxis'),'FontSize',LabelSize);
set(get(gca,'XLabel'),'FontSize',12);
linkaxes([s1,s2],'x');
switch PlotExt
            case 'pdf'
        publish_figurePDF(gcf,[savefile '.pdf']);
            case 'png'
        saveas(gcf,[savefile '.png']);
end
        
% DT Distributions
maintitle=sprintf('KATRIN  First Tritium - Slow Control Parameters');
savefile=sprintf('plots/tmp/KATRIN_FT_SC_DT2');
x   = linspace(1,numel(RunList),numel(RunList));
fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

subplot(2,1,1);
histogram(dt,'Facecolor',rgb('CadetBlue'),'FaceAlpha',.5,'EdgeColor','none');
ylabel('runs')
xlabel('[DT]');xtickangle(45);
sp1 = sprintf('<[DT]>=%g \\pm %g ',mean(dt),std(dt));title(sp1)
grid on
set(gca,'FontSize',16);
PrettyFigureFormat

subplot(2,1,2);
histogram(mean(dtErr,1)./dt,'Facecolor',rgb('CadetBlue'),'FaceAlpha',.5,'EdgeColor','none');
hold off
xlabel('[DT] rel. uncertainty (%)');xtickangle(45);
%sp2 = sprintf('<[DT] rel. uncertainty>=%g  \\pm %.g (std)',mean(dtErr./dt),std(dtErr./dt));title(sp2)
grid on
set(gca,'FontSize',16);
PrettyFigureFormat
switch PlotExt
            case 'pdf'
        publish_figurePDF(gcf,[savefile '.pdf']);
            case 'png'
        saveas(gcf,[savefile '.png']);
        end
% HT
maintitle=sprintf('KATRIN  First Tritium - Slow Control Parameters');
savefile=sprintf('plots/tmp/KATRIN_FT_SC_HT');
x   = linspace(1,numel(RunList),numel(RunList));
fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

s1 = subplot(2,1,1);
errorbar(x,ht,mean(htErr,1),'ks','MarkerSize',12,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
xticks(x)
xticklabels(r)
ylabel('[HT]')
xlabel('run');xtickangle(45);
xlim([0.5 numel(RunList)+0.5])
sp1 = sprintf('<[HT]>=%g \\pm %g  (std)',mean(ht),std(ht));title(sp1)
grid on
set(gca,'FontSize',16);
PrettyFigureFormat
set(get(gca,'XAxis'),'FontSize',LabelSize);
set(get(gca,'XLabel'),'FontSize',12);

s2 = subplot(2,1,2);
errorbar(x,mean(htErr,1)./mean(ht)*100,mean(htErr,1).*0,'ks','MarkerSize',12,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
hold off
xticks(x)
xticklabels(r)
ylabel('[HT] rel. uncertainty (%)>')
xlabel('run');xtickangle(45);
xlim([0.5 numel(RunList)+0.5])
%sp2 = sprintf('<HT>=%g  \\pm %g (std)',mean(ht),std(ht));title(sp2)
grid on
set(gca,'FontSize',16);
PrettyFigureFormat
set(get(gca,'XAxis'),'FontSize',LabelSize);
set(get(gca,'XLabel'),'FontSize',12);
linkaxes([s1,s2],'x');
switch PlotExt
            case 'pdf'
        publish_figurePDF(gcf,[savefile '.pdf']);
            case 'png'
        saveas(gcf,[savefile '.png']);
        end
% TT
maintitle=sprintf('KATRIN  First Tritium - Slow Control Parameters');
savefile=sprintf('plots/tmp/KATRIN_FT_SC_TT');
x   = linspace(1,numel(RunList),numel(RunList));
fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

s1 = subplot(2,1,1);
errorbar(x,tt,mean(ttErr,1),'ks','MarkerSize',12,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
xticks(x)
xticklabels(r)
ylabel('[TT]')
xlabel('run');xtickangle(45);
xlim([0.5 numel(RunList)+0.5])
sp1 = sprintf('<[TT]>=%g \\pm %g  (std)',mean(tt),std(tt));title(sp1)
grid on
set(gca,'FontSize',16);
PrettyFigureFormat
set(get(gca,'XAxis'),'FontSize',LabelSize);
set(get(gca,'XLabel'),'FontSize',12);

s2 = subplot(2,1,2);
errorbar(x,mean(ttErr,1)./mean(tt)*100,mean(ttErr,1).*0,'ks','MarkerSize',12,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
hold off
xticks(x)
xticklabels(r)
ylabel('[TT] rel. uncertainty (%)>')
xlabel('run');xtickangle(45);
xlim([0.5 numel(RunList)+0.5])
%sp2 = sprintf('<TT>=%g  \\pm %g (std)',mean(tt),std(tt));title(sp2)
grid on
set(gca,'FontSize',16);
PrettyFigureFormat
set(get(gca,'XAxis'),'FontSize',LabelSize);
set(get(gca,'XLabel'),'FontSize',12);
linkaxes([s1,s2],'x');

switch PlotExt
            case 'pdf'
        publish_figurePDF(gcf,[savefile '.pdf']);
            case 'png'
        saveas(gcf,[savefile '.png']);
        
end
        
% Combine Plots
%% Build pdf file with all fits
switch PlotExt
            case 'pdf'
      PATH = getenv('PATH');
setenv('PATH', [PATH ':/usr/local/bin/']);
cd plots/tmp
command = 'gs -sDEVICE=pdfwrite -sOutputFile="run.pdf" -dNOPAUSE -dEPSCrop -c "<</Orientation 0>> setpagedevice" -f KATRIN_FT_*.pdf -c quit';
unix(command);
unix('rm KATRIN_FT__*.pdf');
mycommand1 = sprintf('mv run.pdf ../KATRIN_FT_%s-samak.pdf',datestr(now,'dd-mmm-yyyy'));
unix(mycommand1);
mycommand2 = sprintf('open ../KATRIN_FT_%s-samak.pdf',datestr(now,'dd-mmm-yyyy'));
unix(mycommand2);
cd ../..
%unix('rm -rf plots/tmp');
end

