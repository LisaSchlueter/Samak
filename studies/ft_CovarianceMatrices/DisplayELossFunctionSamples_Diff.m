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
d_Aseev = importdata([mypath,sprintf('ELossFunctions_%s_%.0fNIS_1000Trials.mat','Aseev',10)]);
d_Abdu = importdata([mypath,sprintf('ELossFunctions_%s_%.0fNIS_1000Trials.mat','Abdurashitov',8)]);
catch
    fprintf(2,'File doesnt exist \n');
    return
end
%%
R = MultiRunAnalysis('RunList','KNM1_m149mvRW','chi2','chi2Stat','ElossFlag','Abdurashitov');
maxE = (R.ModelObj.Q_i-R.ModelObj.qUmin)*5; minE=-maxE; NbinE = (maxE-minE)/0.1;
E = linspace(minE,maxE,NbinE);
nSamples = 200;
l1 = cell(2,1);
l2 = cell(2,1);
%c1 = colormap(jet(2));
%c2 = [rgb('DarkRed'); rgb('IndianRed')];
nScatmax = 2;
close;
f36 = figure('Renderer','opengl');
set(f36, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
for nScat=1:2
ELossFun1_all_Aseev = zeros(nSamples,numel(E));
ELossFun1_all_Abdu = zeros(nSamples,numel(E));
for i =1:nSamples
ELossFun1_Aseev = d_Aseev{nScat,i};
ELossFun1_all_Aseev(i,:) = ELossFun1_Aseev(E);
ELossFun1_Abdu = d_Abdu{nScat,i};
ELossFun1_all_Abdu(i,:) = ELossFun1_Abdu(E);
end
meanELoss_Aseev =mean(ELossFun1_all_Aseev); 
maxDiff_Aseev = std(ELossFun1_all_Aseev);
meanELoss_Abdu =mean(ELossFun1_all_Abdu); 
maxDiff_Abdu = std(ELossFun1_all_Abdu);
if nScat==1
  [l1{nScat}, a1] = boundedline(E,(meanELoss_Aseev*100)./5,(maxDiff_Aseev*100)./5);
  [l2{nScat}, a2] = boundedline(E,(meanELoss_Abdu*100)./5,(maxDiff_Abdu*100)./5);
elseif nScat==2
      [l1{nScat}, a1] = boundedline(E,(meanELoss_Aseev*100)./2,(maxDiff_Aseev./2)*100);
     [l2{nScat}, a2] = boundedline(E,(meanELoss_Abdu*100)./2,(maxDiff_Abdu./2)*100);
end
l1{nScat}.Color= rgb('CadetBlue'); a1.FaceColor = rgb('CadetBlue');a1.FaceAlpha = 0.5;
l1{nScat}.LineWidth = 2;
l2{nScat}.Color= rgb('IndianRed'); a2.FaceColor = rgb('IndianRed');a2.FaceAlpha = 0.5;
l2{nScat}.LineWidth = 2;
end
%%
xlim([8 35]);
ylim([0 5.3]);
PrettyFigureFormat;
xlabel(['energy loss ',char(949),' (eV)']);
ylabel('probability density (eV^{-1})')
set(gca,'FontSize',PlotFontSize);
leg = legend([a1, a2],'Aseev et al.','Abdurashitov et. al');
legend boxoff
%%
x1 = [0.33,0.29];
y1 = [0.5, 0.5];
y2 = [0.33, 0.33];
switch ELossFlag
    case 'Abdurashitov'
        x2 = [0.64,0.60];
    case 'Aseev'
       x2 = [0.61,0.57];
end

a2= annotation(f36,'textarrow',x1,y1,'String',sprintf('scaled with 1/5'));
a3= annotation(f36,'textarrow',x2,y2,'String',sprintf('scaled with 1/2'));
a2.FontSize = 20;a2.FontWeight='bold';
a3.FontSize = 20;a3.FontWeight='bold';
%% save 
if ~exist('../ft_CovarianceMatrices/plots/EnergyLossFunctions/','dir')
    mkdir ../ft_CovarianceMatrices/plots/EnergyLossFunctions/
end
save_name = sprintf('../ft_CovarianceMatrices/plots/EnergyLossFunctions/ELossFunction_Diff_CovMatSamples');
print(f36,[save_name,'.png'],'-dpng','-r450');
%publish_figurePDF(f36,[save_name,'.pdf']);



