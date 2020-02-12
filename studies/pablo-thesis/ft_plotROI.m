load('Corrections/ROI_Corrections')

rangePlot = 238:402;

plot(qU_correction(rangePlot),pixel_correction(rangePlot,1),'LineWidth',3)

title('ROI Correction Pixel 1');
ylabel('Correction Factor');
xlabel('retarding potential (kV)')

axis([min(qU_correction(rangePlot)) max(qU_correction(rangePlot)) -inf inf])

PrettyFigureFormat;

export_fig('ROICorrection.pdf','-pdf')

