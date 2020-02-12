% plot for KNM1: final fit result
% test FC 
D = MultiRunAnalysis('RunList','KNM1','FSDFlag','Sibille0p5eV','fitter',...
    'minuit','minuitOpt','min;migrad','exclDataStart',14,'chi2','chi2CMShape',...
    'TwinBias_mNuSq',0,'DataType','Twin');
%%
nSamples = 1000;
TBDIS = mvnrnd(D.RunData.TBDIS',D.FitCMShape,nSamples)';

for i=1:nSamples
    progressbar(i/nSamples)
    
    D.RunData.TBDIS = TBDIS(:,i);
    D.RunData.TBDISE = sqrt(TBDIS(:,i));
    
    D.Fit;
    FitPar  = D.FitResult.par;
    FitErr   = D.FitResult.err;
    FitChi2min = D.FitResult.chi2min;
    dof           = D.FitResult.dof;
    
    if FitPar(1)<=-0.9
        mNuIndex = i;
        return
    end
end
%%

TBDIS_mNuSqfit = D.ModelObj.TBDIS;
%%
D.fixPar = '1 2 5 6 7 8 9 10 11';
D.ModelObj.mnuSq_i = 0;
D.ModelObj.Q_i = D.ModelObj.Q_i + FitPar(2);
D.Fit;
TBDIS_mNuSq0eV = D.ModelObj.TBDIS;
TBDIS0 = TBDIS_mNuSq0eV(D.exclDataStart:end);

D.fixPar = '1 5 6 7 8 9 10 11';
D.ModelObj.mnuSq_i = 0;
D.Fit;
TBDIS0free = D.ModelObj.TBDIS(D.exclDataStart:end);

%%
 fig5 = figure('Renderer','openGL');
 set(fig5, 'Units', 'normalized', 'Position', [0.001, 0.001,0.5, 0.8]);
 
 subplot(2,1,2); 
qU = D.ModelObj.qU(D.exclDataStart:end)-18574;
Data = TBDIS(D.exclDataStart:end,mNuIndex);
ResidualsErr = Data.^(0.5)./TBDIS0;
e1 = errorbar(qU,(Data-TBDIS0)./TBDIS0,ResidualsErr,... 
    'o','Color',rgb('DarkSlateGray'),'MarkerFaceColor',rgb('DimGray'),'MarkerSize',5); % MC data
hold on;
p1 = plot(qU,(TBDIS_mNuSqfit(D.exclDataStart:end)-TBDIS0)./TBDIS0,'-','Color',rgb('DodgerBlue'),'LineWidth',2); % fit
p2 = plot(qU,(TBDIS_mNuSq0eV(D.exclDataStart:end)-TBDIS0)./TBDIS0,'-','Color',rgb('Silver'),'LineWidth',2); % expectation at zero nu-mass
hold off;
PRLFormat;
xlabel('retarding energy - 18574 (eV)');
ylabel(sprintf('rel. difference'))
title(sprintf('MC - rel. difference of integral spectra with respect to m_\\beta = 0 eV^2'));
xlim([-42,49]);
leg = legend([e1,p1,p2],'MC data','fit',sprintf('fit with fixed m_\\beta = 0, E_0 = best fit'));
legend boxoff;

 subplot(2,1,1); 
ResidualsErr = Data.^(0.5)./TBDIS0free;
e1 = errorbar(qU,(Data-TBDIS0free)./TBDIS0free,ResidualsErr,... 
    'o','Color',rgb('DarkSlateGray'),'MarkerFaceColor',rgb('DimGray'),'MarkerSize',5); % MC data
hold on;
p1 = plot(qU,(TBDIS_mNuSqfit(D.exclDataStart:end)-TBDIS0free)./TBDIS0free,'-','Color',rgb('DodgerBlue'),'LineWidth',2); % fit
p2 = plot(qU,(TBDIS0free-TBDIS0free)./TBDIS0free,'-','Color',rgb('Silver'),'LineWidth',2); % expectation at zero nu-mass
hold off;
PRLFormat;
xlabel('retarding energy - 18574 (eV)');
ylabel(sprintf('rel. difference'))
title(sprintf('MC - rel. difference of integral spectra with respect to m_\\beta = 0 eV^2'));
xlim([-42,49]);
leg = legend([e1,p1,p2],'MC data','fit',sprintf('fit with fixed m_\\beta = 0, E_0 free'));%E_0 = best fit'));
legend boxoff;
print(gcf,'GuidosPlotSamak_E0.png','-dpng','-r100')