%
% KNM1 Final fit Results
% Uniform Fit
% Golden Run List
% Golden Pixel List
% 
% Last Updated: 19/09/2019
%
%
%% settings
RunList               = 'KNM1';
exclDataStart         = 14; % 27 subruns
RecomputeFlag         = 'OFF';
BkgCM                 = 'ON';

%% Init Model Object and covariance matrix object
Real = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Real',...
    'exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11',...
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1.064,...
    'minuitOpt','min ; minos',...
    'FSDFlag','Sibille0p5eV',... 
    'ELossFlag','KatrinT2',...
    'SysBudget',22);
%% CM Shape
Real.chi2='chi2CMShape'; Real.ComputeCM('BkgCM',BkgCM);
Real.Fit('CATS','OFF');
Real.PlotFit('saveplot','pdf','ErrorBarScaling',1,'YLimRes',[-2.2,2.9],'Colors','RGB','DisplayStyle','PRL');
Rsysshape = Real.FitResult;

%% Scan For the fit chi2 profile on m2 marginalizing over all parameters
Real.fixPar = 'fix 1 ;fix 5 ;fix 6 ;fix 7 ;fix 8 ;fix 9 ;fix 10 ;fix 11 ;';
mass_counter=0;
for mass=-2:0.1:+2
    mass_counter         = mass_counter+1;
    Real.i_mnu           = mass;
    Real.Fit('CATS','OFF')
    myscan(mass_counter,:) = [Real.FitResult.par(1) Real.FitResult.chi2min];
end

%%
figure(1)
subplot(1,2,1)
p=plot(myscan(:,1),smooth(myscan(:,2)),'Color',rgb('DarkBlue'),'LineWidth',4);
hold on
e=myscan(:,1);
ep=e(e>=-1);
em=e(e<=-1);
p1=plot(ep,21.4228+(ep+1).^2./0.89^2,'Color',rgb('DarkGreen'),'LineWidth',2,'LineStyle','--');
p2=plot(em,21.4228+(em+1).^2./1.06^2,'Color',rgb('IndianRed'),'LineWidth',2,'LineStyle','--');
hold off
xlabel(sprintf('m_\\nu^2 (eV^2)'));
ylabel(sprintf('\\chi ^2 / 23 dof'));
%title('Katrin First Neutrino Mass Science Run');
grid on
PrettyFigureFormat
subplot(1,2,2)
p=plot(myscan(:,1),smooth(myscan(:,2)-min(myscan(:,2))),'Color',rgb('DarkBlue'),'LineWidth',4);
hold on
e=myscan(:,1);
ep=e(e>=-1);
em=e(e<=-1);
p1=plot(ep,(ep+1).^2./(0.89)^2,'Color',rgb('DarkGreen'),'LineWidth',2,'LineStyle','--');
p2=plot(em,(em+1).^2./(1.06)^2,'Color',rgb('IndianRed'),'LineWidth',2,'LineStyle','--');
hold off
xlabel(sprintf('m_\\nu^2 (eV^2)'));
ylabel(sprintf('\\Delta \\chi ^2 '));
%title('Katrin First Neutrino Mass Science Run');
legend([p p1 p2],'True profile','Parabolic \sigma=0.9 eV^2','Parabolic \sigma=1.1 eV^2','Location','NorthWest');
legend('boxoff');
grid on
PrettyFigureFormat
export_fig(1,'KatrinFirstNeutrinoMassScienceRun_Chi2Scan.png','-png','-r300');

% Table 
table(myscan(:,1),round(myscan(:,2),2),'VariableNames',{'m2','chisquare'})

% 