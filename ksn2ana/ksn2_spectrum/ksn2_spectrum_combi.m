% ksn2 calculate chi2 grid search
%% settings that might change

range = 40;
LocalFontSize = 18;
LocalLineWidth = 1.5;
 sin2T4 = 0.01;    % works for all 
 mNu4Sq_plt = 10^2; % works from (1:40).^2;
 
 SavePlot = 'ON';
 
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_spectrum/results/'];
savefile = sprintf('%sksn2_spectrum_combi.mat',savedir);

if exist(savefile,'file')
    load(savefile)
else
%% configure RunAnalysis object
chi2 = 'chi2CMShape';
DataType = 'Real';

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
    'FSDFlag','KNM2_0p5eV',...
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
A.exclDataStart = A.GetexclDataStart(range);
A.InitModelObj_Norm_BKG;
A.Fit;
FitResultsNH = A.FitResult;
%%
qU = A.ModelObj.qU(A.exclDataStart:end);
Time = A.ModelObj.qUfrac(A.exclDataStart:end).*A.ModelObj.TimeSec;
TBDIS_NH = A.ModelObj.TBDIS(A.exclDataStart:end);
Rate_NH = TBDIS_NH./Time;
Data = A.RunData.TBDIS(A.exclDataStart:end);
DataRate = Data./Time;
DataRateErr = sqrt(Data)./Time;
Q_i   = A.ModelObj.Q_i;
BKG_i = A.ModelObj.BKG_RateSec_i;
BKG_RateSec_tot = A.ModelObj.BKG_RateSec+0.5.*A.ModelObj.BKG_PtSlope.*Time./A.nRuns;
          
%% get sterile hypothesis
A.ModelObj.mnuSq_i = A.FitResult.par(1);
A.ModelObj.Q_i = Q_i+A.FitResult.par(2);
A.ModelObj.BKG_RateSec_i = BKG_i + A.FitResult.par(3);
A.ModelObj.normFit_i = A.FitResult.par(4);

mNu4Sq = (1:40).^2;
H1 = zeros(numel(mNu4Sq),28);
H1Err = zeros(numel(mNu4Sq),28);
for i=1:numel(mNu4Sq)
    A.ModelObj.SetFitBiasSterile(mNu4Sq(i),1);
    A.ModelObj.ComputeTBDDS;
    A.ModelObj.ComputeTBDIS;
    H1(i,:) = A.ModelObj.TBDIS(A.exclDataStart:end)./Time;
    H1Err(i,:) =  sqrt(A.ModelObj.TBDIS(A.exclDataStart:end))./Time;
end

MakeDir(savedir);
save(savefile,'H1','H1Err','mNu4Sq','qU','Time','Rate_NH','TBDIS_NH','DataRate',...
    'DataRateErr','Q_i','BKG_i','BKG_RateSec_tot','FitResultsNH');
end
%% Plot Spectrum with NH model
% Spectrum + Fit with Residuals
fig5 = figure('Renderer','painters');
set(fig5, 'Units', 'normalized', 'Position', [0.001, 0.001,0.45, 0.7]);
s1= subplot(4,1,[1 2]);

ystr = 'Count rate (cps)';
qUDisp = 'Rel';
ErrorBarScaling = 50;

if strcmp(qUDisp,'Rel')
    x =  qU - 18574;
    xstr = sprintf('Retarding energy - 18574 (eV)');
    textx = -38;
elseif strcmp(qUDisp,'Abs')
    x = qU ;
    xstr = sprintf('Retarding energy (eV)');
    myxticks = (round(min(qU),0):20:max(qU));
    textx =min(qU)+0.5;
end

% subplot 1
pfit = plot(x,Rate_NH,'Color',rgb('DodgerBlue'),'LineWidth',LocalLineWidth);
hold on;
pdata = errorbar(x,DataRate,ErrorBarScaling .*DataRateErr,'k.','CapSize',0,'LineWidth',LocalLineWidth,...
    'MarkerSize',10);
PRLFormat;
set(gca,'FontSize',LocalFontSize);

% subplot1: legend
if ErrorBarScaling==1
    datalabel = sprintf(' KATRIN data with %.0f\\sigma error bars',ErrorBarScaling);
else
    datalabel = sprintf(' KATRIN data with 1 \\sigma error bars \\times %.0f',ErrorBarScaling);
end
leg = legend([pdata,pfit],datalabel,sprintf('3-\\nu model'));
PrettyLegendFormat(leg);
leg.FontSize = LocalFontSize;

ax = gca;
mypos = ax.Position;
ax.Position = [mypos(1)+0.05 mypos(2)+0.03 mypos(3) mypos(4)+0.02];
mylim = ylim;
text(textx(1)+1,max(mylim)*1.1,'a)','FontSize',LocalFontSize+2,'FontName',get(gca,'FontName'));
ylabel('Count rate (cps)','FontSize',LocalFontSize+4);
set(gca,'yscale','log');


if strcmp(qUDisp,'Abs')
    xticks(myxticks);
    ax = gca;
    ax.XAxis.Exponent = 0;
end
             
             
 %% subplot 2
 s2= subplot(4,1,3);
 RateH1_plt = sin2T4.*H1(find(mNu4Sq==mNu4Sq_plt),:)'+ (1-sin2T4).*Rate_NH;
 RateH1err_plt = sqrt(sin2T4.*Time.*H1(find(mNu4Sq==mNu4Sq_plt),:)'+ (1-sin2T4).*Rate_NH.*Time)./Time;
 Ratio = RateH1_plt./Rate_NH;

 pH1  = plot(x,Ratio,'LineWidth',LocalLineWidth,'Color',rgb('LightCoral'));
 hold on;
 pH0 = plot(x,Rate_NH./Rate_NH,'--','Color',rgb('DodgerBlue'),'LineWidth',LocalLineWidth);
 
 fitH1 = errorbar(x,Ratio, RateH1err_plt./Rate_NH,'k.',...
     'CapSize',0,'LineWidth',LocalLineWidth, 'MarkerSize',10);
 pnone = plot(qU,zeros(numel(qU),1),'LineStyle','none');
 hold off
 if strcmp(qUDisp,'Abs')
     xticks(myxticks);
     ax = gca;
     ax.XAxis.Exponent = 0;
 end
 ylim([0.98 1.02]);
   PRLFormat;
 set(gca,'FontSize',LocalFontSize);
 ylabel('Ratio','FontSize',LocalFontSize+4);

 
 katrinsim   = sprintf('3+1 simulation {\\itm}_{4} = %.1f eV   |{\\itU}_{e4}|^2 = %.2f',sqrt( mNu4Sq_plt),sin2T4);
 sterilemod  = sprintf('3+1 model');
 hl=legend([pnone,fitH1, pH1,pH0 ],{'',katrinsim,'3-\nu model',sterilemod},'Location','southeast','box','off');
 ylim([0.978 1.012]);
 hl.FontSize = LocalFontSize;
 hl.NumColumns=2;
 ax2 = gca;


mypos2 = ax2.Position;
ax2.Position = [ax.Position(1) mypos2(2)+0.005 ax.Position(3) mypos2(4)+0.035];
text(ax2.YLabel.Position(1)+18,1.007,'b)','FontSize',LocalFontSize+2,'FontName',get(gca,'FontName'));
%hl.Position(2) = 0.333;
%% 
%% 
s3 = subplot(4,1,4);
bT = bar(x,Time./(60*60),0.5,'FaceColor',rgb('DodgerBlue'),'EdgeColor','none');
bT.BarWidth = 0.7;
% apperance: legend, labels etc.
xlabel(xstr);
ylh = ylabel('Time (h)');
ylh.Position(1) = ax2.YLabel.Position(1);%
yl1.Position(1) = ax2.YLabel.Position(1);
ylh.Position(2) =  35;
ylim([0 70])
yticks([0 35 70])
xlim([-range-2 50])
PRLFormat;
set(gca,'FontSize',LocalFontSize);
set(get(gca,'XLabel'),'FontSize',LocalFontSize+2);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+2);
text(ax2.YLabel.Position(1)+18,60,'c)','FontSize',LocalFontSize+2,'FontName',get(gca,'FontName'));

linkaxes([s1,s2,s3],'x');

ax3 = gca;
mypos2 = ax3.Position;
ax3.Position = [ax.Position(1) mypos2(2)-0.02 ax.Position(3) mypos2(4)+0.035];

if strcmp(SavePlot,'ON')
    plotdir  = strrep(savedir,'results','plots');
    MakeDir(plotdir);
    plotname = sprintf('%sksn2_spectrum_combi_mNu4Sq%.3geV2_sin2T4%.3g.pdf',plotdir,mNu4Sq_plt,sin2T4);
end
export_fig(plotname);
fprintf('save plot to %s \n',plotname);