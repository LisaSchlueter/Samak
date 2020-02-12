% rings
RK=[1:12]-0.2;
RF=[1:12];
RS=[1:12]+0.2;

% KaFit
E090_K    = [ 1.857313706316000   1.857323087214000   1.857337878555000   1.857351764326000   1.857342287781000   1.857343457951000   1.857341694033000   1.857354053137000 1.857352495437000   1.857367158391000   1.857351794755000   1.857369927661000]*1e4;
E090Err_K = [0.111945181119900   0.064799501938510   0.064994857668000   0.064258041808740   0.065237308753340   0.065312416900320   0.065585120165080   0.065551230096850 0.072916519504030   0.076863630190290   0.082308535823370   0.134556957215400];
E040_K    = [18573.14222955 18573.24813136 18573.46415022 18573.55991591 18573.47811784 18573.46110452 18573.38650465 18573.68527718 18573.5345524 18573.6535621 18573.64037174 18573.69633328];
E040Err_K = [0.18733     0.10885     0.10973     0.10837     0.10972     0.11007     0.11071     0.11102     0.12345     0.13099     0.14042     0.23201];
E090m_K     = wmean(E090_K,1./E090Err_K.^2);
E040m_K     = wmean(E040_K,1./E040Err_K.^2);
%E090_K    = E090_K-E0m_K;

% Fitrium
E090_F    = [18573.2566      18573.3235      18573.4539      18573.5502      18573.4744      18573.4483      18573.4517      18573.5623      18573.5057      18573.6627      18573.4988      18573.6596];
E090Err_F = [0.10556    0.064067    0.063852    0.063517    0.064331    0.064517      0.0647    0.064827    0.071945    0.075828     0.08088     0.13281];
E040_F    = [18573.285      18573.4021      18573.5922      18573.6573      18573.5724      18573.5235      18573.4714      18573.7382      18573.5801      18573.6909      18573.6429      18573.7567];
E040Err_F = [0.18423     0.10746     0.10812     0.10698     0.10868     0.10851      0.1085     0.10905     0.12095     0.12799      0.1373     0.22535];
E090m_F     = wmean(E090_F,1./E090Err_F.^2);
E040m_F     = wmean(E040_F,1./E040Err_F.^2);
%E090_F    = E090_F-E0m_F;

% Samak
E090_S    = [18573.181        18573.27       18573.399       18573.553       18573.444       18573.479       18573.468       18573.571       18573.561       18573.704       18573.547       18573.677];
E090Err_S = [0.134       0.077       0.077       0.076       0.077       0.077       0.077       0.077       0.085       0.089       0.096       0.154];
E040_S    = [18573.245       18573.317       18573.516       18573.631       18573.558       18573.557       18573.476       18573.797       18573.617       18573.697       18573.758        18573.59];
E040Err_S = [0.237       0.137       0.138       0.136       0.138       0.137       0.137       0.138       0.153       0.161       0.172       0.281];
E090m_S     = wmean(E090_S,1./E090Err_S.^2);
E040m_S     = wmean(E040_S,1./E040Err_S.^2);
%E090_S    = E090_S-E0m_S;

% Reference Scale
E090ref     = mean([E090m_K E090m_F E090m_S]);
E040ref     = mean([E040m_K E040m_F E040m_S]);

%%%
figure(1)
set(gcf, 'Position',  [100, 100, 1600, 1000])
MyMarkerSize=18;
hk=errorbar(RK,E090_K-E090ref,E090Err_K,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2);
hold on
hf=errorbar(RF,E090_F-E090ref,E090Err_F,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',2);
hs=errorbar(RS,E090_S-E090ref,E090Err_S,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('Amethyst'),'LineWidth',2) ;
hold off
ylabel('effective endpoint (eV) - Constant');
xlabel('ring');
leg = legend([hk,hf,hs],'KaFit','Fitrium','Samak','Location','northwest','FontSize',24);
leg.Color = 'none'; legend boxoff;
grid on
PrettyFigureFormat;
set(gca,'FontSize',20);
title('Ring-wise Fit - 90eV Range');

figure(2)
set(gcf, 'Position',  [100, 100, 1600, 1000])
hk=errorbar(RK,E040_K-E040ref,E040Err_K,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2) ;
hold on
hf=errorbar(RF,E040_F-E040ref,E040Err_F,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',2);
hs=errorbar(RS,E040_S-E040ref,E040Err_S,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('Amethyst'),'LineWidth',2);
hold off
ylabel('effective endpoint (eV)- Constant');
xlabel('ring');
leg = legend([hk,hf,hs],'KaFit','Fitrium','Samak','Location','northwest','FontSize',24);
leg.Color = 'none'; legend boxoff;
PrettyFigureFormat;
set(gca,'FontSize',20);
title('Ring-wise Fit - 40eV Range');
grid on

%%
figure(3)
set(gcf, 'Position',  [100, 100, 1600, 1000])
hk=plot(R,E090_K-E090_F,'LineWidth',5) ;
hold on
hf=plot(R,E090_K-E090_S,'LineWidth',5) ;
hs=plot(R,E090_F-E090_S,'LineWidth',5) ;
hold off
ylabel('effective endpoint difference (eV)');
xlabel('ring');
leg = legend([hk,hf,hs],'KaFit-Fitrium','KaFit-Samak','Samak-Fitrium','Location','northwest','FontSize',24);
leg.Color = 'none'; legend boxoff;
PrettyFigureFormat;
set(gca,'FontSize',20);
title('Ring-wise Fit - 90eV Range');
grid on

%%
figure(4)
set(gcf, 'Position',  [100, 100, 1600, 1000])
hk=plot(R,E040_K-E040_F,'LineWidth',5) ;
hold on
hf=plot(R,E040_K-E040_S,'LineWidth',5) ;
hs=plot(R,E040_F-E040_S,'LineWidth',5) ;
hold off
ylabel('effective endpoint difference (eV)');
xlabel('ring');
leg = legend([hk,hf,hs],'KaFit-Fitrium','KaFit-Samak','Samak-Fitrium','Location','northwest','FontSize',24);
leg.Color = 'none'; legend boxoff;
PrettyFigureFormat;
set(gca,'FontSize',20);
title('Ring-wise Fit - 90eV Range');
grid on
