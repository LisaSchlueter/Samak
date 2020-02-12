if ~exist('SPAll','var')
    BuildTable(1);
    load('BuildTableWorkspace')
end




E0_ring_s = SPix.TableShortStat(:,2) + E0_i - E_R_s;
E0_ring_m = SPix.TableMedStat(:,2) + E0_i - E_R_m;

 
FPDViewer([E0_ring_s;nan;nan])
title('First Tritium - Ring-wise Short Stat. E0 - <E_0> (eV)')
% publish_figurePDF(gcf,'plots/FPDRingShortStat.pdf')
export_fig('plots/FPDRingShortStat.png','-png')

FPDViewer([E0_ring_m;nan;nan])
title('First Tritium - Ring-wise Med. Stat. E0 - <E_0> (eV)')
% publish_figurePDF(gcf,'plots/FPDRingMedStat.pdf')
export_fig('plots/FPDRingMedStat.png','-png')


xRing = 1:11;
ylim = [-1 1];

figring = figure(8);
set(figring,'units','normalized','position',[0.001 0.05 0.9*(9/16) 0.7]);
% PLOT 1
subplot(2,2,1)
E0Ring_s = SPix.TableShortStat(:,2) + E0_i - E_R_s;
E0Ring_s_e = SPix.TableShortStat(:,8);
hold on
errorb(xRing,E0Ring_s,E0Ring_s_e);
h1 = plot(xRing,E0Ring_s,'sk','MarkerFaceColor','black');
plot(xRing,zeros(1,length(xRing)),'Color','red','LineStyle','-','LineWidth',3);
hold off
suptitle('First Tritium - fitted E_0 - <E_0> per Ring Stacked Runs (eV)');
axis([min(xRing)-1,max(xRing)+1,ylim]);
xticks(xRing)
%xticklabels([])
ylabel('Short (E_0 - 200 eV)');
PrettyFigureFormat;

% PLOT 2
subplot(2,2,2)
E0Ring_s = SPix.TableMedStat(:,2) + E0_i - E_R_m;
E0Ring_s_e = SPix.TableMedStat(:,8);
hold on
errorb(xRing,E0Ring_s,E0Ring_s_e);
h1 = plot(xRing,E0Ring_s,'sk','MarkerFaceColor','black');
plot(xRing,zeros(1,length(xRing)),'Color','red','LineStyle','-','LineWidth',3);
hold off
axis([min(xRing)-1,max(xRing)+1,ylim]);
xticks(xRing)
%xticklabels([])
ylabel('Med. (E_0 - 200 eV)');
PrettyFigureFormat;

% PLOT 3
subplot(2,2,3)

chi2_s = SPix.TableShortStat(:,13); dof_s = SPix.TableShortStat(1,14);
hold on
[~,Nbins,Wbins] = nhist(chi2_s);
binWidth = Wbins(2) - Wbins(1);
x_N = linspace(min(chi2_s),max(chi2_s),...
    length(Nbins)*50);
plot(x_N,length(chi2_s)*binWidth*chi2pdf(x_N,dof_s),'LineWidth',4);
hold off
xlabel('Short (E_0 - 200 eV)');

% PLOT 4
subplot(2,2,4)
chi2_s = SPix.TableMedStat(:,13); dof_s = SPix.TableMedStat(1,14);
hold on
[~,Nbins,Wbins] = nhist(chi2_s);
binWidth = Wbins(2) - Wbins(1);
x_N = linspace(min(chi2_s),max(chi2_s),...
    length(Nbins)*50);
plot(x_N,length(chi2_s)*binWidth*chi2pdf(x_N,dof_s),'LineWidth',4);
hold off
xlabel('Med. (E_0 - 200 eV)');

export_fig('plots/RingStat.png','-png');

