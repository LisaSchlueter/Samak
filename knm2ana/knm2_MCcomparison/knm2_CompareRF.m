% KNM2 response function comparison

%% Load RF Samak
savenameS = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/Knm2_SamakRF.mat'];
dSamak = importdata(savenameS);
qU = dSamak.qU;
TeS = dSamak.Te-qU; % energy
RFS = dSamak.RF; % transmission probability

%% Load KaFit
savenameK = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/Knm2_KaFitRF.dat'];
dKaFit = importdata(savenameK);
TeK = dKaFit(:,1); % energy
RFK = dKaFit(:,2); % transmission probability

%% Load Fitrium
savenameF = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/Knm2_FitriumRF.dat'];
RFF = ones(numel(TeS),1);

%% plot ratio 
f2 = figure('Units','normalized','Position',[0,0.1,0.5,0.5]);
pK = plot(TeS,RFS./RFK,'-','LineWidth',2.5,'Color',rgb('DodgerBlue'));
hold on;
pF = plot(TeS,RFS./RFS,'-.','LineWidth',2.5,'Color',rgb('GoldenRod'));
xlabel(sprintf('Energy - %.0f (eV)',qU));
ylabel('Rel. probability diff.');
PrettyFigureFormat('FontSize',22);
leg = legend([pK,pF],'Samak / KaFit','Samak / Fitrium');
leg.EdgeColor = rgb('Silver');
xlim([min(TeS),max(TeS)]);
%% plot difference
f2 = figure('Units','normalized','Position',[0.5,0.1,0.5,0.5]);
pK = plot(TeS,RFS-RFK,'-','LineWidth',2.5,'Color',rgb('DodgerBlue'));
hold on;
pF = plot(TeS,RFS-RFS,'-.','LineWidth',2.5,'Color',rgb('GoldenRod'));
xlabel(sprintf('Energy - %.0f (eV)',qU));
ylabel('Probability diff.');
PrettyFigureFormat('FontSize',22);
leg = legend([pK,pF],'Samak - KaFit','Samak - Fitrium');
leg.EdgeColor = rgb('Silver');
xlim([min(TeS),max(TeS)]);