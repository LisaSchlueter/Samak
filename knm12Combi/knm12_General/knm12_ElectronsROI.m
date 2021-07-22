% plot electrons as a function of live time in knm1 and knm2
savedir = sprintf('%sknm12Combi/knm12_General/results/',getenv('SamakPath'));
savename = sprintf('%sknm12_Data.mat',savedir);

if ~exist(savename,'file')
    M2 = MultiRunAnalysis('RunList','KNM2_Prompt','NonPoissonScaleFactor',1);
    M1 = MultiRunAnalysis('RunList','KNM1','NonPoissonScaleFactor',1);
    
    StackData1 = M1.RunData;
    Data1      = M1.SingleRunData;
    
    StackData2 = M2.RunData;
    Data2      = M2.SingleRunData;
    
    MakeDir(savedir);
    save(savename,'Data1','Data2','StackData1','StackData2');
else
    load(savename);
end
%% plot
Idx1 = 13; Idx2 = 11;
Time   = [Data1.StartTimeStamp, Data2.StartTimeStamp];
Counts1 = cumsum(sum(Data1.TBDIS(Idx1:end,:)));
Counts2 = cumsum(sum(Data2.TBDIS(Idx2:end,:)))+Counts1(end);

GetFigure;
p1 = plot(Data1.StartTimeStamp,Counts1,'LineWidth',3.5,'Color',rgb('Orange'));
hold on;
p2 = plot(Data2.StartTimeStamp-135,Counts2,'LineWidth',3.5,'Color',rgb('ForestGreen'));
PRLFormat;
ylabel('Cumulative electrons in ROI');
xlabel('');
xticks([mean(Data1.StartTimeStamp),mean(Data2.StartTimeStamp)-135]);
xticklabels({'April - May 2019','September - November 2019'});
set(gca,'FontSize',24);
%%
t1.delete;
t2.delete;
get(gca,'FontSize')
t1 = text(mean(Data1.StartTimeStamp)-10,6e6,...%mean(Counts1)*1.3,...
    sprintf('1^{st} campaign'),...%{\\itm}_\\nu
    'FontSize',ax.XLabel.FontSize,'FontName',get(gca,'FontName'),...
    'Color',p1.Color,'Rotation',0);
t2 = text(mean(Data2.StartTimeStamp)-145,6e6,...%mean(Counts2),...
    sprintf('2^{nd} campaign'),...
    'FontSize',ax.XLabel.FontSize,'FontName',get(gca,'FontName'),...
    'Color',p2.Color,'Rotation',0);
get(gca,'FontSize');

ylim([0 6.5e6])

pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);
pltname = sprintf('%sknm12_electronsROI.pdf',pltdir);
export_fig(pltname);
