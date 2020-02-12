% Plot for of FPD ROI versus FPD coverage
% FPD Coverage versus ROI
% Sanshiro data, 31/10/2018
FPD_File = load('FPDROICoverage.txt');
FPD_ROIlow =FPD_File(:,1); % lowerlimit ROI cut
FPD_Cov = @(x) interp1(FPD_File(:,1),FPD_File(:,2),x);

f1 = figure('Renderer','opengl');
set(f1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
plot(FPD_ROIlow,FPD_Cov(FPD_ROIlow),'color',rgb('IndianRed'),'LineWidth',4);
xlim([14,32]);
PrettyFigureFormat;
xlabel('ROI lower limit (keV)');
ylabel('FPD coverage')
grid on;
save_name = './plots/FPDcoverage-ROIcut';
print(f1,[save_name,'.png'],'-dpng');
publish_figurePDF(f1,[save_name,'.pdf']);