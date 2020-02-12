% Energy Loss Function Samples for Covariance Matrix Computation
% Plot and Save Plot
ELossFlag = 'Abdurashitov'; %
PlotFontSize = 22;
if strcmp(ELossFlag,'Abdurashitov')
    nScatmax = 8;
elseif strcmp(ELossFlag,'Aseev')
    nScatmax = 10;
end

mypath = '../../inputs/CovMat/RF/LookupTables/ELoss/';
try
d = importdata([mypath,sprintf('ELossFunctions_%s_%.0fNIS_1000Trials.mat',ELossFlag,nScatmax)]);
catch
    fprintf(2,'File doesnt exist \n');
    return
end
%%
R = MultiRunAnalysis('RunList','KNM1_m149mvRW','chi2','chi2Stat','ElossFlag','Abdurashitov');
maxE = (R.ModelObj.Q_i-R.ModelObj.qUmin)*5; minE=-maxE; NbinE = (maxE-minE)/0.1;
E = linspace(minE,maxE,NbinE);
nSamples = 200;
l = cell(10,1);
c = colormap(jet(10));
if strcmp(ELossFlag,'Abdurashitov')
    nScatmax = 8;
elseif strcmp(ELossFlag,'Aseev')
    nScatmax = 10;
end

close;
f33 = figure('Renderer','opengl');
set(f33, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
pnone = plot(E,NaN*zeros(numel(E),1),'Color',rgb('White')); hold on;
pnone2 = plot(E,NaN*zeros(numel(E),1),'Color',rgb('White')); hold on;
pnone3 = plot(E,NaN*zeros(numel(E),1),'Color',rgb('White')); hold on;
for nScat=1:nScatmax   
ELossFun1_all = zeros(nSamples,numel(E));
for i =1:nSamples
ELossFun1 = d{nScat,i};
ELossFun1_all(i,:) = ELossFun1(E);
end
meanELoss =mean(ELossFun1_all); 
%maxDiff = max(abs(meanELoss-ELossFun1_all));
maxDiff = std(ELossFun1_all);
if nScat==1
  [l{nScat}, a] = boundedline(E,meanELoss./5,maxDiff./5);
elseif nScat==2
      [l{nScat}, a] = boundedline(E,meanELoss./2,maxDiff./2);
else
    [l{nScat}, a] = boundedline(E,meanELoss,maxDiff);
end
l{nScat}.Color= c(nScat,:); a.FaceColor = c(nScat,:);a.FaceAlpha = 0.5;
xlim([0 125]);
PrettyFigureFormat;
xlabel(['energy loss ',char(949),' (eV)']);
%title('Energy Loss Functions - Samples for Covariance Matrix')
ylabel('probability density (eV^{-1})')
ylim([0 0.055]);
end
%%
set(gca,'FontSize',PlotFontSize);
switch ELossFlag
    case 'Aseev'       
ELossInfo = sprintf('Energy Loss Uncertainties from \nEur. Phys. J. D 10, 39–52 (Aseev et. al.)');
   AnnotationPos = [0.48 0.35 1 0.1]; %[0.53 0.14 1 0.1]; for 500eV range
    case 'Abdurashitov'
ELossInfo = sprintf('Energy Loss Uncertainties from \nPhys.Part.Nucl.Lett. 14 (2017) \nno.6, 892-899 (Abdurashitov et al)');
  AnnotationPos = [0.58 0.55 1 0.1];
end
leg = legend([l{:}],string(1:10)); legend boxoff
leg.NumColumns = 2;
leg.Orientation = 'horizontal'; leg.Location = 'northeast';
leg.Title.String = 'number of scatterings';
% a=annotation('textbox', AnnotationPos, ...
%     'String', ELossInfo, ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'left');
% a.FontSize=16;a.FontWeight='bold';

%x1 = [0.3,0.165];
x1 = [0.35,0.21];
y1 = [0.83, 0.83];
y2 = [0.6, 0.6];
%x2 = [0.3,0.2];
x2 = [0.35,0.29];
a2= annotation(f33,'textarrow',x1,y1,'String',sprintf('scaled with 1/5'));
a3= annotation(f33,'textarrow',x2,y2,'String',sprintf('scaled with 1/2'));
a2.FontSize = 16;a2.FontWeight='bold';
a3.FontSize = 16;a3.FontWeight='bold';
a2.LineWidth=1.5; a3.LineWidth=1.5;
%% save 
if ~exist('../ft_CovarianceMatrices/plots/EnergyLossFunctions/','dir')
    mkdir ../ft_CovarianceMatrices/plots/EnergyLossFunctions/
end
save_name = sprintf('../ft_CovarianceMatrices/plots/EnergyLossFunctions/ELossFunction_%s_CovMatSamples',ELossFlag);
print(f33,[save_name,'.png'],'-dpng','-r450');
%publish_figurePDF(f33,[save_name,'.pdf']);

%% plot short range too
f34 = figure('Renderer','opengl');
set(f34, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
l = cell(2,1);
c = colormap(jet(2));
for nScat=1:2
ELossFun1_all = zeros(nSamples,numel(E));
for i =1:nSamples
ELossFun1 = d{nScat,i};
ELossFun1_all(i,:) = ELossFun1(E);
end
meanELoss =mean(ELossFun1_all); 
%maxDiff = max(abs(meanELoss-ELossFun1_all));
maxDiff = std(ELossFun1_all);
if nScat==1
  [l{nScat}, a] = boundedline(E,(meanELoss*100)./5,(maxDiff./5)*100);
elseif nScat==2
      [l{nScat}, a] = boundedline(E,(meanELoss*100)./2,(maxDiff./2)*100);
end
l{nScat}.Color= c(nScat,:); a.FaceColor = c(nScat,:);a.FaceAlpha = 0.5;
ylim([0 5]);
PrettyFigureFormat;
xlabel(['energy loss ',char(949),' (eV)']);
ylabel('probability density (eV^{-1})');
end
set(gca,'FontSize',PlotFontSize);
xlim([8 35]);
leg = legend([l{:}],string(1:2)); legend boxoff
leg.NumColumns = 2;
leg.Orientation = 'horizontal'; leg.Location = 'northeast';
leg.Title.String = 'number of scatterings';
x1 = [0.33,0.29];
y1 = [0.5, 0.5];
y2 = [0.33, 0.33];
switch ELossFlag
    case 'Abdurashitov'
        x2 = [0.64,0.60];
    case 'Aseev'
       x2 = [0.61,0.57];
end

a2= annotation(f34,'textarrow',x1,y1,'String',sprintf('scaled with 1/5'));
a3= annotation(f34,'textarrow',x2,y2,'String',sprintf('scaled with 1/2'));
a2.FontSize = 20;a2.FontWeight='bold';
a3.FontSize = 20;a3.FontWeight='bold';

save_name = sprintf('../ft_CovarianceMatrices/plots/EnergyLossFunctions/ELossFunction_%s_CovMatSamples_60eV',ELossFlag);
print(f34,[save_name,'.png'],'-dpng','-r450');
%publish_figurePDF(f34,[save_name,'.pdf']);
%% Display correlation Matrix (Abdurashitov et. al)
load('../../inputs/ELossFunction/ELossFunction_CorrMat.mat');

f35 = figure('Renderer','opengl');
set(f35, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.5, 0.7]);
corplot(ELossCorrMat_Average);
xticklabels({'N_1','N_2','W_1','W_2','P_1','P_2'});
yticklabels({'N_1','N_2','W_1','W_2','P_1','P_2'});
PrettyFigureFormat;
c = colorbar; colormap
set(gca,'FontSize',22);
c.Label.String = 'correlation coefficient';
c.Label.FontSize = 22;
save_name = sprintf('../ft_CovarianceMatrices/plots/EnergyLossFunctions/ELossFunction_CorrMatAverage');
print(f35,[save_name,'.png'],'-dpng','-r450');
%publish_figurePDF(f35,[save_name,'.pdf']);