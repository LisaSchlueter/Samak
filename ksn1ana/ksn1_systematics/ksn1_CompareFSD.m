% plot fitrium + samak fsd

%% load samak fsd
fsamak = [getenv('SamakPath'),'inputs/FSD/FSD_KNM1_T2_Doppler0p5eV_FullRange.txt'];
fS = importdata(fsamak);
EnergyS = fS(:,1);
ProbS   = fS(:,2);

%% load fitrium fsd
fFitrium = [getenv('SamakPath'),'ksn1ana/ksn1_systematics/results/FSD_Fitrium_Doppler.txt'];
fF = importdata(fFitrium);
EnergyF = fF(:,1);
ProbF   = fF(:,2);
%% plot
GetFigure;
pS = plot(EnergyS,ProbS,'-','LineWidth',2);
hold on;
pF = plot(EnergyF,ProbF,'-.','LineWidth',2);
PrettyFigureFormat('FontSize',22);

xlabel('Excitation energy (eV)');
ylabel('Probability');
xlim([0 100])
legend([pF,pS],'Fitrium FSD','Samak FSD','EdgeColor',rgb('Silver'));

%% find start rebinning
EnergyS(find(diff(EnergyS)>0.2,1,'first'))



