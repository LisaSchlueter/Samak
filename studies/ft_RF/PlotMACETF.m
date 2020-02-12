
A = ref_TBD_NominalKATRIN('nTeBinningFactor',100,'WGTS_B_T',3.6,'MACE_Bmax_T',6,'MACE_Ba_T',7e-4);
qUIndex = 9;
TheoRes = 7e-4/6*A.qU(qUIndex)*1.2;

TF1 = A.ComputeMaceTF(A.Te,A.qU(qUIndex),'pixel',1,'SynchrotronLoss',1e-99,'gammaFacApprox',1.);
TF2 = A.ComputeMaceTF(A.Te,A.qU(qUIndex),'pixel',1,'SynchrotronLoss',1e-99,'gammaFacApprox',1.018);
TF3 = A.ComputeMaceTF(A.Te,A.qU(qUIndex),'pixel',1,'SynchrotronLoss',34e-3,'gammaFacApprox',1.);
TF4 = A.ComputeMaceTF(A.Te,A.qU(qUIndex),'pixel',1,'SynchrotronLoss',34e-3,'gammaFacApprox',1.018);

figure(1)
subplot(3,1,[1 2])
h1=plot(A.Te-A.qU(qUIndex),TF1,'LineWidth',2,'Color',rgb('FireBrick'));
hold on
h2=plot(A.Te-A.qU(qUIndex),TF2,'LineWidth',2,'Color',rgb('Orange'));
h3=plot(A.Te-A.qU(qUIndex),TF3,'LineWidth',2,'Color',rgb('CadetBlue'));
h4=plot(A.Te-A.qU(qUIndex),TF4,'LineWidth',2,'Color',rgb('DarkBlue'));
hold off
legend([h1 h2 h3 h4],'1) Non Relativistic - No synchrotron','2) Relativistic - No synchrotron','3) Non Relativistic - Synchrotron','4) Relativistic - Synchrotron');
xlim([-0.1 TheoRes])
PrettyFigureFormat;
%xlabel('electron surplus energy E - qU (eV)');
ylabel('transmission probability');
set(gca,'FontSize',16);
ylim([-0.05, 1.2]);
subplot(3,1,3)
h1=plot(A.Te-A.qU(qUIndex),(TF1-TF1),'LineWidth',2,'Color',rgb('FireBrick'));
hold on
h2=plot(A.Te-A.qU(qUIndex),(TF2-TF1),'LineWidth',2,'Color',rgb('Orange'));
h3=plot(A.Te-A.qU(qUIndex),(TF3-TF1),'LineWidth',2,'Color',rgb('CadetBlue'));
h4=plot(A.Te-A.qU(qUIndex),(TF4-TF1),'LineWidth',2,'Color',rgb('DarkBlue'));
hold off
legend([h1 h2 h3 h4],'1-1','2-1','3-1','4-1');
xlim([-0.1 TheoRes])
PrettyFigureFormat;
xlabel('electron surplus energy E - qU (eV)');
ylabel('Transmission difference');
set(gca,'FontSize',16);