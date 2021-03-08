 startup;
% Recompute = 'OFF';
% B=RelicNuAnalysis('Params','KNM1');
 matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Misc/')];
% B.SystBreakdown('TwinBias_mnuSq',0,'Plot','OFF','Recompute',Recompute);
% B.SystBreakdown('TwinBias_mnuSq',1,'Plot','OFF','Recompute',Recompute);
% B.SystBreakdown('TwinBias_mnuSq',5,'Plot','OFF','Recompute',Recompute);
% % B.SystBias('TwinBias_mnuSq',0,'Recompute',Recompute,'Plot','ON','DeltaChi2',1);
% % B.SystBias('TwinBias_mnuSq',1,'Recompute',Recompute,'Plot','ON','DeltaChi2',1);
% % B.SystBias('TwinBias_mnuSq',5,'Recompute',Recompute,'Plot','ON','DeltaChi2',1);
load([matFilePath,'SensitivityBreakdown_KNM1_mnuSq0_range40_RealData.mat']);
YD=Y;
load([matFilePath,'SensitivityBreakdown_KNM1_mnuSq0_range40.mat']);
Y0=Y;
load([matFilePath,'SensitivityBreakdown_KNM1_mnuSq1_range40.mat']);
Y1=Y;
load([matFilePath,'SensitivityBreakdown_KNM1_mnuSq5_range40.mat']);
Y5=Y;
load([matFilePath,'SystEffectOnSensitivity_KNM1_mnuSq0.mat']);
Bias0=Y;Bias0(3)=Bias(4);Bias0(4)=Bias(1);Bias0(5)=Bias(2);Bias0(6)=Bias(3);Bias0(7)=Bias(6);Bias0(8)=Bias(5);Bias0(9)=Bias(7);
load([matFilePath,'SystEffectOnSensitivity_KNM1_mnuSq1.mat']);
Bias1=Y;Bias1(3)=Bias(4);Bias1(4)=Bias(1);Bias1(5)=Bias(2);Bias1(6)=Bias(3);Bias1(7)=Bias(6);Bias1(8)=Bias(5);Bias1(9)=Bias(7);
load([matFilePath,'SystEffectOnSensitivity_KNM1_mnuSq5.mat']);
Bias5=Y;Bias5(3)=Bias(4);Bias5(4)=Bias(1);Bias5(5)=Bias(2);Bias5(6)=Bias(3);Bias5(7)=Bias(6);Bias5(8)=Bias(5);Bias5(9)=Bias(7);

% fig1=figure(1);
% bar(X,[Y0;Y1;Y5]);
% ax=gca;
% ax.YScale='log';
% ylabel('\eta');
% legend('m_{\nu}^{2}=0 eV^{2}','m_{\nu}^{2}=1 eV^{2}','m_{\nu}^{2}=5 eV^{2}');
% legend boxoff;
% PrettyFigureFormat;
% 
% fig10=figure(10);
% Y0rel = Y0./Y0(1);
% Y1rel = Y1./Y1(1);
% Y5rel = Y5./Y5(1);
% bar(X,[Y0rel;Y1rel;Y5rel]);
% ax=gca;
% ax.YScale='log';
% legend('m_{\nu}^{2}=0 eV^{2}','m_{\nu}^{2}=1 eV^{2}','m_{\nu}^{2}=5 eV^{2}');
% legend boxoff;
% ylabel('Fraction of total');
% PrettyFigureFormat;

fig100=figure(100);
X = categorical({'Total','Statistical','Final-state distribution','Response function','Scan fluctuations','Stacking','Detector efficiency','Theoretical corrections','Bkg slope','Bkg rate'});
X = reordercats(X,{'Response function','Theoretical corrections','Scan fluctuations','Detector efficiency','Final-state distribution','Stacking','Bkg slope','Bkg rate','Statistical','Total'});
bsingle = cell(numel(X),1);

PlotColor = {rgb('White'),rgb('Navy'),rgb('GoldenRod'),rgb('PowderBlue'),...
                rgb('CadetBlue'),rgb('DarkOrange'),rgb('FireBrick'),rgb('DarkSlateGray'),...
                rgb('YellowGreen'),rgb('Magenta'),...
                rgb('SeaGreen'),rgb('DodgerBlue')};

hold on;
for i=1:numel(X)
    bsingle{i}  = barh(X(i),YD(i));
    btmp= bsingle{i};
    btmp.LineStyle = 'none';
    btmp.FaceColor = PlotColor{i+1}; btmp.LineStyle ='none';
end
%B=barh(X,Y1,'green');
ax=gca;
ax.XScale='log';
xlabel('\eta');

% SysUncertainties = bsingle{1}.YData;
% a=annotation('textbox', [0.14 0.1 1 0.12], ...
%     'String', SysUncertainties, ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'left');
% a.FontSize=16;a.FontWeight='normal';

effs = ones(1,numel(bsingle));
vals = ones(1,numel(bsingle));
labels = string(zeros(10,1));
for i=1:numel(bsingle)
    effs(i) = bsingle{i}.YEndPoints;
    vals(i) = bsingle{i}.XEndPoints;
    labels(i)=sprintf('%.2g',bsingle{i}.YData);
end
text(effs,vals,labels,'HorizontalAlignment','left',...
   'VerticalAlignment','middle','FontSize',16);
xlim([1e9,1e12]);
PrettyFigureFormat;

% fig2=figure(2);
% bar(X,[Y0;Bias0]);
% ylabel('\eta');
% legend('m_{\nu}^{2}=0 ev^{2}, CM','m_{\nu}^{2}=0 ev^{2}, envelope');
% legend boxoff;
% PrettyFigureFormat;
% 
% fig3=figure(3);
% bar(X,[Y1;Bias1]);
% ylabel('\eta');
% legend('m_{\nu}^{2}=1 ev^{2}, CM','m_{\nu}^{2}=1 ev^{2}, envelope');
% legend boxoff;
% PrettyFigureFormat;
% 
% fig4=figure(4);
% bar(X,[Y5;Bias5]);
% ylabel('\eta');
% legend('m_{\nu}^{2}=5 ev^{2}, CM','m_{\nu}^{2}=5 ev^{2}, envelope');
% legend boxoff;
% PrettyFigureFormat;