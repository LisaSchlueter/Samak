%clear;

% RFS = load('C:\Users\navar\Documents\samak2.0\inputs\ResponseFunction\samakRF_5e+17cm2_NIS11_Bm6T_Bs3.6T_Ba3G_multipix.mat');

RFS.RFip(RFS.E<0,:,:) = 0;
RF = squeeze(RFS.RFip(:,5,:));
Te = RFS.E;
Ba = RFS.MACE_Ba_TallPixels_save;

Bamean = mean(Ba);
[minBa,minID] = min(abs(Ba-Bamean));
%plot(abs(Ba-Bamean))

%RFmean = mean(RF,2);

RFmean = RF(:,minID);
RFRatio = RF./RFmean;
plot(Te,RF);
%xticks(1:40);
%xlim([0 40])
%ylim([-inf inf])
title('Ratio RF from all pixels to mean RF')
xlabel('qU - E_0 (eV)')
ylabel('Ratio')

PrettyFigureFormat;