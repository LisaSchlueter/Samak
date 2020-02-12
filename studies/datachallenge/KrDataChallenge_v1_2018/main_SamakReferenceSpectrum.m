addpath(genpath('../../../Samak2.0'));
orange=1/255*[255,140,0];

Obj_KrDC= InitKATRIN_KrDataChallenge_benchmark();
Obj_KrDC.ComputeKrDS();
Obj_KrDC.ComputeKrIS();

RefSpec = [Obj_KrDC.qU , Obj_KrDC.KrIS , Obj_KrDC.KrISE ];
save('../../krypton-data/DataChallenge18V1/SamakReferenceSpec.txt','RefSpec','-ascii');

%figure('Name', 'Reference Spectrum');
figure(11);
s1= subplot(2,1,1)
gdiff = plot(Obj_KrDC.Te*1e-3, Obj_KrDC.KrDS, '-', 'Color', orange,'LineWidth', 2);
title('Samak Reference Spectrum');
xlim([30.4670 30.475]);
ylim([0 40]);
xlabel('Te (keV)','FontSize', 14);
ylabel('$\frac{\textrm{\textbf{d}}\dot{\textrm{\textbf{N}}}}{\textrm{\textbf{d T}}_\textrm{\textbf{e}}} \left(\frac{\textbf{cps}}{\textbf{0.2 eV}}\right)$', 'Interpreter','latex', 'FontSize', 14);
grid on;
diffleg = sprintf('Benchmark Parameters:\n A = 60 cps \n E = 30.4725 keV \n W = 1.19 eV \n O = 10 cps \n T = 100 K');
legend(diffleg,'Location', 'Northwest');
s2= subplot(2,1,2);
gint = errorbar(Obj_KrDC.qU*1e-03, Obj_KrDC.KrIS,Obj_KrDC.KrISE, 's','MarkerSize',2, 'Color', orange);
hold on;
offline = line([Obj_KrDC.qU(1)*1e-03, Obj_KrDC.qU(end)*1e-03], [Obj_KrDC.KrIS(end),Obj_KrDC.KrIS(end)],'LineStyle','--', 'LineWidth', 1.5);
hold off;
xlim([30.4670 30.475]);
xlabel('qU (keV)','FontSize', 14);
ylabel('$\dot{\textrm{\textbf{N}}}$ \textbf{(cps)}','Interpreter','latex','FontSize', 14);
intleg = sprintf('TD = flat (0.2eV steps)\nTime = 30s per qU');
intleg = legend(intleg, 'Offset');
grid on;
%export_fig ./plots/KrDataChallenge_SamakReference_new.png