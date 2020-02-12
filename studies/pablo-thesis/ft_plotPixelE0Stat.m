if ~exist('SPAll','var')
    BuildTable(1);
    load('BuildTableWorkspace')
end

PixelList = [1:112,114:124];
FillerPixels = nan(148,1);


E0_pixel_s = SPix.TableShortStat(:,2) + E0_i - E_SPix_s;
E0_pixel_m = SPix.TableMedStat(:,2) + E0_i - E_SPix_m;

FillerPixels(PixelList)=E0_pixel_s;
E0_pixel_s = FillerPixels;
FPDViewer([E0_pixel_s])
title('First Tritium - pixel-wise Short Stat. E0 - <E_0> (eV)')
publish_figurePDF(gcf,'plots/FPDpixelShortStat.pdf')
%export_fig('plots/FPDpixelShortStat.pdf','-pdf')

FillerPixels(PixelList)=E0_pixel_m;
E0_pixel_m = FillerPixels;
FPDViewer([E0_pixel_m])
title('First Tritium - pixel-wise Med. Stat. E0 - <E_0> (eV)')
publish_figurePDF(gcf,'plots/FPDpixelMedStat.pdf')
%export_fig('plots/FPDpixelMedStat.pdf','-pdf')


xpixel = 1:length(PixelList);
ylim = [-3 3];

figpixel = figure(8);
set(figpixel,'units','normalized','position',[0.001 0.05 0.9*(9/16) 0.7]);
% PLOT 1
subplot(2,2,1)
E0pixel_s = SPix.TableShortStat(:,2) + E0_i - E_SPix_s;
E0pixel_s_e = SPix.TableShortStat(:,8);
hold on
errorb(xpixel,E0pixel_s,E0pixel_s_e);
h1 = plot(xpixel,E0pixel_s,'sk','MarkerFaceColor','black');
plot(xpixel,zeros(1,length(xpixel)),'Color','red','LineStyle','-','LineWidth',3);
hold off
axis([min(xpixel)-1,max(xpixel)+1,ylim]);
xticks(xpixel(1:20:end))
xticklabels(string(PixelList(1:20:end)-1))
ylabel('Short (E_0 - 200 eV)');
PrettyFigureFormat;

% PLOT 2
subplot(2,2,2)
E0pixel_s = SPix.TableMedStat(:,2) + E0_i - E_SPix_m;
E0pixel_s_e = SPix.TableMedStat(:,8);
hold on
errorb(xpixel,E0pixel_s,E0pixel_s_e);
h1 = plot(xpixel,E0pixel_s,'sk','MarkerFaceColor','black');
plot(xpixel,zeros(1,length(xpixel)),'Color','red','LineStyle','-','LineWidth',3);
hold off
axis([min(xpixel)-1,max(xpixel)+1,ylim]);
xticks(xpixel(1:20:end))
xticklabels(string(PixelList(1:20:end)-1))
ylabel('Med. (E_0 - 400 eV)');
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
xlabel('\chi^2 short (E_0 - 200 eV)');

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
xlabel('\chi^2 med. (E_0 - 400 eV)');
suptitle('First Tritium - effective fitted E_0 - <E_0> per pixel Stacked Runs (eV)');

export_fig('plots/pixelStat.pdf','-pdf');

