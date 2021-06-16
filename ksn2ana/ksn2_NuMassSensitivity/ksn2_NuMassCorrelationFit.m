% estimate correlation m4-m_nu
% method:
% 1. select MC truth for m_4,Ue4
% 2. vary m_4 slightly around MC truth & fit m_nu
DataType = 'Twin';
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/CorrFits/'];
MakeDir(savedir);
nStepsAll = 25;
nSteps = 5;


savefileCombi = sprintf('%sksn2_NuMassCorrelationFit_%s_%.0fGrid_%.0ffits.mat',savedir,DataType,nStepsAll,nSteps);

if exist(savefileCombi,'file') 
    load(savefileCombi)
else
    % init variables
    mNu4Sq_allgrid = repmat(logspace(-1,log10(40^2),nStepsAll),nStepsAll,1);
    mNu4Sq_all = reshape(mNu4Sq_allgrid,nStepsAll^2,1);
    sin2T4_allgrid = repmat(logspace(-3,log10(0.5),nStepsAll)',1,nStepsAll);
    sin2T4_all = reshape(sin2T4_allgrid,nStepsAll^2,1);
    CorrCoeff = zeros(nStepsAll^2,1);
    Covar = zeros(nStepsAll^2,1);
    SlopeFit = zeros(nStepsAll^2,1);
    OffsetFit = zeros(nStepsAll^2,1);
    mNuSq_fit = zeros(nStepsAll^2,nSteps);
    mNu4Sq_fit = zeros(nStepsAll^2,nSteps);
    
    for j=1:numel(mNu4Sq_all)
        progressbar(j./numel(mNu4Sq_all));
        mNu4Sq = mNu4Sq_all(j);
        sin2T4 = sin2T4_all(j);
        
        savefile = sprintf('%sksn2_NuMassCorrelationFit_%s_MCtruth_m4Sq%.3geV2_sin2T4%.3g_nSteps%.0f.mat',savedir,DataType,mNu4Sq,sin2T4,nSteps);
        if exist(savefile,'file')
            %  load(savefile);
            d = importdata(savefile);
            CorrCoeff(j) = d.CorrMat(2);
            Covar(j)     = d.CovMat(2);
            mNuSq_fit(j,:) = d.mNuSq;
            mNu4Sq_fit(j,:) = d.mNu4Sq_test;
            p=polyfit(mNu4Sq_fit(j,:), mNuSq_fit(j,:),1);
            SlopeFit(j) = p(1); 
            OffsetFit(j) = p(2);
        else
            %% settings that might change
            chi2 = 'chi2CMShape';
            nGridSteps = 50;
            range = 40;
            
            %% configure RunAnalysis object
            if strcmp(chi2,'chi2Stat')
                NonPoissonScaleFactor = 1;
            elseif  strcmp(chi2,'chi2CMShape')
                NonPoissonScaleFactor = 1.112;
            end
            RunAnaArg = {'RunList','KNM2_Prompt',...
                'DataType',DataType,...
                'fixPar','mNu E0 Norm Bkg',...%free par
                'SysBudget',40,...
                'fitter','minuit',...
                'minuitOpt','min;migrad',...
                'RadiativeFlag','ON',...
                'FSDFlag','KNM2_0p1eV',...
                'ELossFlag','KatrinT2A20',...
                'AnaFlag','StackPixel',...
                'chi2',chi2,...
                'NonPoissonScaleFactor',NonPoissonScaleFactor,...
                'FSD_Sigma',sqrt(0.0124+0.0025),...
                'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
                'TwinBias_Q',18573.7,...
                'PullFlag',99,...;%99 = no pull
                'BKG_PtSlope',3*1e-06,...
                'TwinBias_BKG_PtSlope',3*1e-06,...
                'DopplerEffectFlag','FSD'};
            A = MultiRunAnalysis(RunAnaArg{:});
            A.InitModelObj_Norm_BKG('RecomputeFlag','ON');
            %%
            
            A.ModelObj.SetFitBiasSterile(mNu4Sq,sin2T4);
            A.ModelObj.ComputeTBDDS(...
                'E0_bias',A.FitResult.par(2),...
                'B_bias',A.FitResult.par(3),...
                'N_bias',A.FitResult.par(4));
            A.ModelObj.ComputeTBDIS;
            TBDIS_i = A.ModelObj.TBDIS;
            A.RunData.TBDIS = TBDIS_i;
            
            
            mNu4Sq_test = mNu4Sq+linspace(-1,1,nSteps);
            FitResult = cell(nSteps,1);
            mNuSq = zeros(nSteps,1);
            for i=1:nSteps
                progressbar((i-1)./nSteps)
                A.ModelObj.SetFitBias(0);
                A.ModelObj.SetFitBiasSterile(mNu4Sq_test(i),sin2T4);
                A.Fit;
                FitResult{i} = A.FitResult;
                mNuSq(i) = A.FitResult.par(1);
            end
            
            CorrMat = corrcoef(mNu4Sq_test,mNuSq);
            CovMat  = cov(mNu4Sq_test,mNuSq);
            
            save(savefile,'mNuSq','mNu4Sq_test','mNu4Sq','sin2T4','FitResult','CorrMat','CovMat');
        end
    end
    %%
    CorrCoeff  = reshape(CorrCoeff,nStepsAll,nStepsAll);
    Covar      = reshape(Covar,nStepsAll,nStepsAll);
    SlopeFit   = reshape(SlopeFit,nStepsAll,nStepsAll);
    OffsetFit  = reshape(OffsetFit,nStepsAll,nStepsAll);
    mNuSq_fit  =  reshape(mNuSq_fit,nStepsAll,nStepsAll,nSteps);
    mNu4Sq_fit = reshape(mNu4Sq_fit,nStepsAll,nStepsAll,nSteps);
    
    %%
    save(savefileCombi,'sin2T4_allgrid','mNu4Sq_allgrid',...
        'CorrCoeff','Covar','mNuSq_fit','mNu4Sq_fit','SlopeFit','OffsetFit');
end

%% plotting
plotdir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/plots/'];
MakeDir(plotdir);
%% plot slope from linear fit
GetFigure
 %surf(sin2T4_allgrid,mNu4Sq_allgrid,SlopeFit,'EdgeColor','none');
% hold on;
% [~,ct] = contour3(sin2T4_allgrid,mNu4Sq_allgrid,SlopeFit,[-0.1,0,0.1],'Color',rgb('Black'),...
%     'ShowText','on','LineWidth',2.5);
ContourVec = [-2,-1,-0.5,-0.25,-0.1,-0.05,-0.01,0,0.01,0.05,0.1,0.2];
[M,ct] = contourf(sin2T4_allgrid,mNu4Sq_allgrid,SlopeFit,ContourVec,...
    'ShowText','on','LineWidth',1.5,'LabelSpacing',180,'LineColor',rgb('Black'));
cl = clabel(M,ct,'FontSize',12);
set(gca,'XScale','log');
set(gca,'YScale','log');
view(2);
c = colorbar; 
c.Ticks = (-1:0.2:0.1);
c.Label.String = sprintf('Slope {\\ita}  ({\\itm}_\\nu^2 = {\\ita} \\cdot {\\itm}_4^2 + offset)');
grid off; 
PrettyFigureFormat('FontSize',22);
c.Label.FontSize = 18;
ylim([0.1 1600]);
xlim([1e-03 0.5])
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
plotname = sprintf('%sksn2_NuMassCorrelationFit_%s_%.0fGrid_%.0ffits_LinSlope.png',plotdir,DataType,nStepsAll,nSteps);
print(gcf,plotname,'-dpng','-r350');

%% plot correlation
GetFigure
surf(sin2T4_allgrid,mNu4Sq_allgrid,CorrCoeff,'EdgeColor','none');
set(gca,'XScale','log');
set(gca,'YScale','log');
view(2);
c = colorbar;
c.Ticks = (-1:0.2:1);
c.Label.FontSize = 18;
c.Label.String = sprintf('Correlation coefficient');
grid off;
PrettyFigureFormat('FontSize',22);

ylim([0.1 1600]);
xlim([1e-03 0.5])
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
plotname = sprintf('%sksn2_NuMassCorrelationFit_%s_%.0fGrid_%.0ffits_CorrCoeff.png',plotdir,DataType,nStepsAll,nSteps);
print(gcf,plotname,'-dpng','-r350');

%% plot covariance
GetFigure
surf(sin2T4_allgrid,mNu4Sq_allgrid,Covar,'EdgeColor','none');
set(gca,'XScale','log');
set(gca,'YScale','log');
view(2);
c = colorbar;
c.Label.FontSize = 18;
c.Label.String = sprintf('Covariance (eV^4)');
grid off;
PrettyFigureFormat('FontSize',22);

ylim([0.1 1600]);
xlim([1e-03 0.5])
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
plotname = sprintf('%sksn2_NuMassCorrelationFit_%s_%.0fGrid_%.0ffits_Covar.png',plotdir,DataType,nStepsAll,nSteps);
print(gcf,plotname,'-dpng','-r350');

%% 1 example fit
GetFigure
IdxmNu = 14;
Idxsin = 1;
x = linspace(min(squeeze(mNu4Sq_fit(Idxsin,IdxmNu,:)))-0.5,max(squeeze(mNu4Sq_fit(Idxsin,IdxmNu,:)))+0.5,1e2);
pfit = plot(x,SlopeFit(Idxsin,IdxmNu).*x+OffsetFit(Idxsin,IdxmNu),'-','LineWidth',2,'Color',rgb('DodgerBlue'));
hold on;
px = plot(squeeze(mNu4Sq_fit(Idxsin,IdxmNu,:)),squeeze(mNuSq_fit(Idxsin,IdxmNu,:)),'k.','LineWidth',2,'MarkerSize',18);
PrettyFigureFormat('FontSize',20);
xlabel(sprintf('Fix: {\\itm}_4^2 (eV^2)'));
ylabel(sprintf('Fit: {\\itm}_\\nu^2 (eV^2)'));
leg = legend([px,pfit],sprintf('|{\\itU}_{e4}|^2 = %.3g',sin2T4_allgrid(Idxsin,IdxmNu)),...
    sprintf('Linear fit slope = %.2g ',SlopeFit(Idxsin,IdxmNu)),...
    'FontWeight','normal','FontSize',get(gca,'FontSize'));
xlim([min(x),max(x)])
legend boxoff
plotname = sprintf('%sksn2_NuMassCorrelationFit_%s_%.0fGrid_%.0ffits_mNu4Sq%.1f_sin2T4%.3g_LinFit.png',...
    plotdir,DataType,nStepsAll,nSteps,mNu4Sq_allgrid(Idxsin,IdxmNu),sin2T4_allgrid(Idxsin,IdxmNu));
print(gcf,plotname,'-dpng','-r350');
fprintf('save plot to %s \n',plotname)