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
close all
GetFigure
 %surf(sin2T4_allgrid,mNu4Sq_allgrid,SlopeFit,'EdgeColor','none');
% hold on;
% [~,ct] = contour3(sin2T4_allgrid,mNu4Sq_allgrid,SlopeFit,[-0.1,0,0.1],'Color',rgb('Black'),...
%     'ShowText','on','LineWidth',2.5);
ContourVec = [-2,-1,-0.5,-0.25,-0.1,-0.05,-0.01,0,0.01,0.05,0.1,0.2];
surf(sin2T4_allgrid,mNu4Sq_allgrid,SlopeFit-1e-02,'EdgeColor','interp','FaceColor','interp');
hold on;
[M,ct] = contour3(sin2T4_allgrid,mNu4Sq_allgrid,SlopeFit,ContourVec,...
    'ShowText','off','LineWidth',2.5,'LabelSpacing',180,'LineColor',rgb('Black'));
%clabel(M,ct,'FontSize',26,'FontName','Times New Roman','Color',ct.LineColor);
set(gca,'XScale','log');
set(gca,'YScale','log');
view(2);
c = colorbar; 
c.Ticks = (-1:0.2:0.1);
c.Label.String = sprintf('\\alpha_{slope}');%sprintf('Slope {\\ita}  ({\\itm}_\\nu^2 = {\\ita} \\cdot {\\itm}_4^2 + offset)');
grid off; 
PRLFormat;
set(gca,'FontSize',30);
ax = gca;
c.Label.FontSize = ax.XLabel.FontSize;
yticks([[0.1,1,10, 1e2 ,1e3]])
ylim([0.1 1600]);
xlim([1e-03 0.5])
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
ax.Position(3) = 0.65;
ax.Position(1) = 0.165;

% by hand: text boxes:
% t.delete;
% t2.delete;
% t3.delete;
% t4.delete;
% t5.delete;
% t6.delete;
% t7.delete;
% t8.delete;
% t9.delete;
% t10.delete;
t = text(0.0086,0.15,0.17,'-0.01','FontSize',28,'FontName','Times New Roman','Rotation',90) ;
t2 = text(0.041,0.15,0.17,'-0.05','FontSize',28,'FontName','Times New Roman','Rotation',90) ;
t3 = text(0.078,0.15,0.17,'-0.1','FontSize',28,'FontName','Times New Roman','Rotation',90) ;
t4 = text(0.168,0.15,0.17,'-0.25','FontSize',28,'FontName','Times New Roman','Rotation',90) ;
t5 = text(0.275,0.15,0.17,'-0.5','FontSize',28,'FontName','Times New Roman','Rotation',90) ;
t6 = text(0.415,0.15,0.17,'-1','FontSize',28,'FontName','Times New Roman','Rotation',90) ;
t7 = text(0.0015,47,0.17,'0','FontSize',28,'FontName','Times New Roman','Rotation',0) ;
t8 = text(0.025,190,0.17,'0.01','FontSize',28,'FontName','Times New Roman','Rotation',30) ;
t9 = text(0.115,160,0.17,'0.05','FontSize',28,'FontName','Times New Roman','Rotation',40) ;
t10 = text(0.225,150,0.17,'0.1','FontSize',28,'FontName','Times New Roman','Rotation',50) ;
ax = gca;
ax.SortMethod='ChildOrder';
%
plotname = sprintf('%sksn2_NuMassCorrelationFit_%s_%.0fGrid_%.0ffits_LinSlope.pdf',plotdir,DataType,nStepsAll,nSteps);
export_fig(plotname);

%print(gcf,plotname,'-dpng','-r350');


return
%% plot correlation
GetFigure
surf(sin2T4_allgrid,mNu4Sq_allgrid,CorrCoeff,'EdgeColor','none');
% [M,ct] = contourf(sin2T4_allgrid,mNu4Sq_allgrid,CorrCoeff,[-1,-0.9999,0,0.999,1],...
%     'ShowText','on','LineWidth',1.5,'LabelSpacing',180,'LineColor',rgb('Black'));
%cl = clabel(M,ct,'FontSize',12);
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
%surf(sin2T4_allgrid,mNu4Sq_allgrid,Covar,'EdgeColor','none');
[M,ct] = contourf(sin2T4_allgrid,mNu4Sq_allgrid,Covar,[-0.7:0.1:0.2,-0.01,0.01,0.05],...
    'ShowText','on','LineWidth',1.5,'LabelSpacing',250,'LineColor',rgb('Black'));
cl = clabel(M,ct,'FontSize',14);
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