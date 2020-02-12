if ~exist('SPAll','var')
    BuildTable(1);
    load('BuildTableWorkspace')
end




E0_ring_s = R.TableShortStat(:,2) + E0_i - E_R_s;
E0_ring_m = R.TableMedStat(:,2) + E0_i - E_R_m;

 
FPDViewer([E0_ring_s;nan;nan])
title('First Tritium - Ring-wise Short Stat. E0 - <E_0> (eV)')
publish_figurePDF(gcf,'plots/FPDRingShortStat.pdf')

FPDViewer([E0_ring_m;nan;nan])
title('First Tritium - Ring-wise Med. Stat. E0 - <E_0> (eV)')
publish_figurePDF(gcf,'plots/FPDRingMedStat.pdf')


xRing = 1:11;
ylim = [-1 1];

figring = figure(8);
set(figring,'units','normalized','position',[0.001 0.05 0.9*(9/16) 0.7]);
% PLOT 1
subplot(2,2,1)
E0Ring_s = R.TableShortStat(:,2) + E0_i - E_R_s;
E0Ring_s_e = R.TableShortStat(:,8);
hold on
errorb(xRing,E0Ring_s,E0Ring_s_e);
h1 = plot(xRing,E0Ring_s,'sk','MarkerFaceColor','black');
plot(xRing,zeros(1,length(xRing)),'Color','red','LineStyle','-','LineWidth',3);
hold off
axis([min(xRing)-1,max(xRing)+1,ylim]);
xticks(xRing)
%xticklabels([])
ylabel('Short (E_0 - 200 eV)');
PrettyFigureFormat;

% PLOT 2
subplot(2,2,2)
E0Ring_s = R.TableMedStat(:,2) + E0_i - E_R_m;
E0Ring_s_e = R.TableMedStat(:,8);
hold on
errorb(xRing,E0Ring_s,E0Ring_s_e);
h1 = plot(xRing,E0Ring_s,'sk','MarkerFaceColor','black');
plot(xRing,zeros(1,length(xRing)),'Color','red','LineStyle','-','LineWidth',3);
hold off
axis([min(xRing)-1,max(xRing)+1,ylim]);
xticks(xRing)
%xticklabels([])
ylabel('Med. (E_0 - 400 eV)');
PrettyFigureFormat;

% PLOT 3
subplot(2,2,3)

chi2_s = R.TableShortStat(:,13); dof_s = R.TableShortStat(1,14);
hold on
[~,Nbins,Wbins] = nhist(chi2_s);
binWidth = Wbins(2) - Wbins(1);
x_N = linspace(min(chi2_s),max(chi2_s),...
    length(Nbins)*50);
plot(x_N,length(chi2_s)*binWidth*chi2pdf(x_N,dof_s),'LineWidth',4);
hold off
xlabel('\chi^2 short (E_0 - 200 eV)');

% PLOT 4
subplot(2,2,4)
chi2_s = R.TableMedStat(:,13); dof_s = R.TableMedStat(1,14);
hold on
[~,Nbins,Wbins] = nhist(chi2_s);
binWidth = Wbins(2) - Wbins(1);
x_N = linspace(min(chi2_s),max(chi2_s),...
    length(Nbins)*50);
plot(x_N,length(chi2_s)*binWidth*chi2pdf(x_N,dof_s),'LineWidth',4);
hold off
xlabel('\chi^2 med. (E_0 - 400 eV)');

suptitle('First Tritium - effective fitted E_0 - <E_0> per Ring Stacked Runs (eV)');


export_fig('plots/RingStat.pdf','-pdf');

