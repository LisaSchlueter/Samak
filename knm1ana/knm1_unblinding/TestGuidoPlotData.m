%% settings
RunList               = 'KNM1';
exclDataStart         = 14; % 27 subruns

%% Init Model Object and covariance matrix object
tic;
Real = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Real',...
    'exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11',...
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1.064,...
    'minuitOpt','min ; minos',...
    'FSDFlag','SibilleFull',...
    'ELossFlag','KatrinT2',...
    'SysBudget',22);
Real.Fit;
toc;
%%
Real.chi2 = 'chi2CMShape';
Real.ComputeCM;
Real.Fit;
TBDIS_mNuSqfit = Real.ModelObj.TBDIS(Real.exclDataStart:end);
Data = Real.RunData.TBDIS(Real.exclDataStart:end);
E0fit = Real.FitResult.par(2)+Real.ModelObj.Q_i;
%%
Real.fixPar = '1 5 6 7 8 9 10 11';
Real.ModelObj.mnuSq_i = 0;
%Real.ModelObj.Q_i = E0fit;
Real.Fit;
TBDIS0free = Real.ModelObj.TBDIS(Real.exclDataStart:end);

Real.fixPar = '1 2 5 6 7 8 9 10 11';
Real.ModelObj.mnuSq_i = 0;
Real.ModelObj.Q_i = E0fit;
Real.Fit;
TBDIS0 = Real.ModelObj.TBDIS(Real.exclDataStart:end);

%%

fig5 = figure('Renderer','openGL');
set(fig5, 'Units', 'normalized', 'Position', [0.001, 0.001,0.5, 0.8]);

subplot(2,1,1);
qU = Real.ModelObj.qU(Real.exclDataStart:end)-18574;
%ResidualsErr = (  (Data.^(0.5)./TBDIS0).^2 + (Data.*TBDIS0.^(0.5)./TBDIS0.^2).^2  ).^(0.5);
ResidualsErr = Data.^(0.5)./TBDIS0free; % assuming only data has an uncertainty
e1 = errorbar(qU,(Data-TBDIS0free)./TBDIS0free,ResidualsErr,... 
    'o','Color',rgb('DarkSlateGray'),'MarkerFaceColor',rgb('DimGray'),'MarkerSize',5); % MC data
hold on;
p1 = plot(qU,(TBDIS_mNuSqfit-TBDIS0free)./TBDIS0free,'-','Color',rgb('IndianRed'),'LineWidth',2); % fit
p2 = plot(qU,(TBDIS0free-TBDIS0free)./TBDIS0free,'-','Color',rgb('Silver'),'LineWidth',2); % expectation at zero nu-mass
hold off;
PRLFormat;
xlabel('retarding energy - 18574 (eV)');
ylabel(sprintf('rel. difference'))
title(sprintf('Data - rel. difference of integral spectra with respect to m_\\beta = 0 eV^2'));
xlim([-42,49]);
leg = legend([e1,p1,p2],'Data','fit',sprintf('Fit with fixed m_\\beta = 0, E_0 free'));%E_0 = best fit')); %E_0 free'));
legend boxoff;

subplot(2,1,2);
qU = Real.ModelObj.qU(Real.exclDataStart:end)-18574;

%ResidualsErr = (  (Data.^(0.5)./TBDIS0).^2 + (Data.*TBDIS0.^(0.5)./TBDIS0.^2).^2  ).^(0.5);
ResidualsErr = Data.^(0.5)./TBDIS0; % assuming only data has an uncertainty
e1 = errorbar(qU,(Data-TBDIS0)./TBDIS0,ResidualsErr,... 
    'o','Color',rgb('DarkSlateGray'),'MarkerFaceColor',rgb('DimGray'),'MarkerSize',5); % MC data
hold on;
p1 = plot(qU,(TBDIS_mNuSqfit-TBDIS0)./TBDIS0,'-','Color',rgb('IndianRed'),'LineWidth',2); % fit
p2 = plot(qU,(TBDIS0-TBDIS0)./TBDIS0,'-','Color',rgb('Silver'),'LineWidth',2); % expectation at zero nu-mass
hold off;
PRLFormat;
xlabel('retarding energy - 18574 (eV)');
ylabel(sprintf('rel. difference'))
xlim([-42,49]);
leg = legend([e1,p1,p2],'Data','fit',sprintf('Fit with fixed m_\\beta = 0, E_0 = best fit')); %E_0 free'));
legend boxoff;


print(gcf,'GuidosPlotSamak_E0_Data.png','-dpng','-r100')