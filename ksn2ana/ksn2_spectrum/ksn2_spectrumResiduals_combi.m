% KSN2 paper plot
%% settings that might change

range = 40;
LocalFontSize = 20;
LocalLineWidth = 2.5;

 SavePlot = 'ON';
 
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_spectrum/results/'];
savefile = sprintf('%sksn2_spectrumBF_combi.mat',savedir);

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

%% fit with 3nu model
A.InitModelObj_Norm_BKG;
A.Fit;
FitResultsNH = A.FitResult;

% some info
qU = A.ModelObj.qU(A.exclDataStart:end);
Time = A.ModelObj.qUfrac(A.exclDataStart:end).*A.ModelObj.TimeSec;
Q_i   = A.ModelObj.Q_i;
BKG_i = A.ModelObj.BKG_RateSec_i;
BKG_RateSec_tot = A.ModelObj.BKG_RateSec+0.5.*A.ModelObj.BKG_PtSlope.*Time./A.nRuns;

% data
DataCounts = A.RunData.TBDIS(A.exclDataStart:end);
DataRate = DataCounts./Time;
DataRateErr = sqrt(DataCounts)./Time;
 
% null hypothesis best fit with 3nu model
NullBfCounts = A.ModelObj.TBDIS(A.exclDataStart:end);
NullBfRate   = NullBfCounts./Time;
NullBfRateErr = sqrt(NullBfCounts)./Time;
 

%% get sterile hypothesis best fit parameter
savedir = sprintf('%sksn2ana/ksn2_BestFit/results/',getenv('SamakPath'));
savefile2 = sprintf('%sksn2_BestFitPar_mNuE0NormBkg_%s_%s.mat',...
    savedir,DataType,chi2);
if exist(savefile2,'file')
    dbf = importdata(savefile2);
end
A.ModelObj.mnuSq_i = dbf.FitResult.par(1);
A.ModelObj.Q_i = Q_i+dbf.FitResult.par(2);
A.ModelObj.BKG_RateSec_i = BKG_i + dbf.FitResult.par(3);
A.ModelObj.normFit_i = dbf.FitResult.par(4);

A.ModelObj.SetFitBiasSterile(dbf.mNu4Sq_bf,dbf.sin2T4_bf);
A.ModelObj.ComputeTBDDS;
A.ModelObj.ComputeTBDIS;
H1BfRate    = A.ModelObj.TBDIS(A.exclDataStart:end)./Time;
H1BfRateErr =  sqrt(A.ModelObj.TBDIS(A.exclDataStart:end))./Time;

% divide active and sterile branch in 3nu+1 fit
% sterile
A.ModelObj.SetFitBiasSterile(dbf.mNu4Sq_bf,1); % only sterile
A.ModelObj.ComputeTBDDS;
A.ModelObj.TBDDS = A.ModelObj.TBDDS.*dbf.sin2T4_bf;
A.ModelObj.ComputeTBDIS;
RateSterileB    = A.ModelObj.TBDIS(A.exclDataStart:end)./Time-BKG_RateSec_tot;
RateSterileB(RateSterileB<0)=0;

% active
%A.ModelObj.normFit_i = dbf.FitResult.par(4).*(1-dbf.sin2T4_bf);
A.ModelObj.SetFitBiasSterile(0,0); % no sterile
A.ModelObj.ComputeTBDDS;
A.ModelObj.TBDDS = A.ModelObj.TBDDS.*(1-dbf.sin2T4_bf);
A.ModelObj.ComputeTBDIS;
RateActiveB    = A.ModelObj.TBDIS(A.exclDataStart:end)./Time-BKG_RateSec_tot;
RateActiveB(RateActiveB<0) = 0; % numerical noise

%% some simulation for ratio plot. Use null hypothesis best fit parameter + some sterile hypothesis
A.ModelObj.mnuSq_i = FitResultsNH.par(1);
A.ModelObj.Q_i = Q_i+FitResultsNH.par(2);
A.ModelObj.BKG_RateSec_i = BKG_i + FitResultsNH.par(3);
A.ModelObj.normFit_i = FitResultsNH.par(4);

A.ModelObj.SetFitBiasSterile(dbf.mNu4Sq_bf,dbf.sin2T4_bf);
A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
SimRateBf =  A.ModelObj.TBDIS(A.exclDataStart:end)./Time;

SimmNu4Sq = [10,10,20,20,30,30].^2;
Simsin2T4 = [0.01,0.1,0.01,0.1,0.01,0.1];
SimRate = zeros(numel(SimmNu4Sq),numel(SimRateBf));

for i=1:numel(SimmNu4Sq)
    A.ModelObj.SetFitBiasSterile(SimmNu4Sq(i),Simsin2T4(i));
    A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
    SimRate(i,:) =  A.ModelObj.TBDIS(A.exclDataStart:end)./Time;
end

%% save
MakeDir(savedir);
save(savefile,...
     'Q_i','BKG_i','BKG_RateSec_tot','qU','Time',...
     'DataCounts','DataRate','DataRateErr',...
    'FitResultsNH','NullBfCounts','NullBfRate','NullBfRateErr',...
     'dbf','H1BfRate','H1BfRateErr','RateActiveB','RateSterileB',...
   'SimRateBf','SimRate','SimmNu4Sq','Simsin2T4');
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
pfit = plot(x,H1BfRate,'Color',rgb('DodgerBlue'),'LineWidth',LocalLineWidth);
hold on;
pactive = plot(x,RateActiveB,'-.','LineWidth',LocalLineWidth,'Color',rgb('Orange'));
%hold on;
psterile = plot(x,RateSterileB,':','LineWidth',LocalLineWidth,'Color',rgb('FireBrick'));
%psterile = plot(x,RateSterileB+BKG_RateSec_tot-H1,'-','LineWidth',LocalLineWidth,'Color',rgb('Black'));
pdata = errorbar(x,DataRate,ErrorBarScaling .*DataRateErr,'k.','CapSize',0,'LineWidth',LocalLineWidth-0.5,...
    'MarkerSize',10);
PRLFormat;
set(gca,'FontSize',LocalFontSize);

% subplot1: legend
if ErrorBarScaling==1
    datalabel = sprintf(' KATRIN data with %.0f\\sigma error bars',ErrorBarScaling);
else
    datalabel = sprintf(' KATRIN data with 1 \\sigma error bars \\times %.0f',ErrorBarScaling);
end
leg = legend([pdata,pfit,pactive,psterile],datalabel,...
    sprintf('3\\nu+1 best fit model'),...
    sprintf('Active branch'),...%: {\\itm}_\\nu^2 = %.1f eV^2',dbf.FitResult.par(1)),...
    sprintf('Sterile branch: {\\itm}_\\nu^2 = %.1f eV^2, |{\\itU}_{e4}|^2 = %.3f',dbf.mNu4Sq_bf,dbf.sin2T4_bf));
PrettyLegendFormat(leg);
leg.FontSize = LocalFontSize-0.5;
set(gca,'yscale','log');
ylabel('Count rate (cps)','FontSize',LocalFontSize+6);

% axis
ax = gca;
xlim([-43 140]);
ylim([0.1, 100]);
mypos = ax.Position;
ax.Position = [mypos(1)+0.05 mypos(2)+0.035 mypos(3) mypos(4)+0.02];
mylim = ylim;
text(textx(1)+1,65,'a)','FontSize',LocalFontSize+2,'FontName',get(gca,'FontName'));


if strcmp(qUDisp,'Abs')
    xticks(myxticks);
    ax = gca;
    ax.XAxis.Exponent = 0;
end
             
           
%% subplot 2
Mode = 'Ratio';

s2= subplot(4,1,3);
%RateH1_plt = sin2T4.*H1(find(mNu4Sq==mNu4Sq_plt),:);
%RateH1err_plt = sqrt(sin2T4.*Time.*H1(find(mNu4Sq==mNu4Sq_plt),:)'+ (1-sin2T4).*Rate_NH.*Time)./Time;
switch Mode
    case 'Residuals'
        y1 = (DataRate-H1BfRate)./DataRateErr;
        y0 = (DataRate-NullBfRate)./DataRateErr;
        
        yref = 0;
        yStr = sprintf('Residuals (\\sigma)');
        MarkerStyle = '.';
    case 'Ratio'  
        y1 = SimRateBf./NullBfRate;
        y0 = NullBfRate./NullBfRate; 
      
        yref = 1;
        yStr = sprintf('Ratio');
        MarkerStyle = '-';
end

%pref = plot(linspace(-50,150,10),yref.*ones(10,1),'-','Color',rgb('Black'),'LineWidth',1.5);
pnone = plot(NaN,NaN,'.','Color',rgb('White'));
hold on;
if strcmp(Mode,'Residuals')
    plot(linspace(-50,150,10),yref.*ones(10,1)-1,'-','Color',rgb('Gray'),'LineWidth',1.5);
    plot(linspace(-50,150,10),yref.*ones(10,1)+1,'-','Color',rgb('Gray'),'LineWidth',1.5);
end

pH0  = plot(linspace(-50,150,numel(y1)),y0,':','MarkerSize',10,'LineWidth',LocalLineWidth-0.5,'Color',rgb('Black'));
pH1  = plot(x,y1,MarkerStyle,'MarkerSize',10,'LineWidth',LocalLineWidth,'Color',rgb('DodgerBlue'));

%eSimData = errorbar(x,SimRateBf./NullBfRate,(sqrt(SimRateBf.*Time)./Time)./NullBfRate,...
%    '.','Color',rgb('Black'),'CapSize',0,'LineWidth',LocalLineWidth-0.5, 'MarkerSize',10);
hold off
if strcmp(qUDisp,'Abs')
    xticks(myxticks);
    ax = gca;
    ax.XAxis.Exponent = 0;
end

PRLFormat;
set(gca,'FontSize',LocalFontSize);
ylabel(yStr,'FontSize',LocalFontSize+6);

katrinsim   = sprintf('3\\nu+1 simulation {\\itm}_4^2 = %.1f eV^2, |{\\itU}_{e4}|^2 = %.3f',dbf.mNu4Sq_bf,dbf.sin2T4_bf);

%ylim([0.972 1.013]);
ylim([0.977 1.015]);
xlim([-43 140]);

ax2 = gca;
%axis
mypos2 = ax2.Position;
ax2.Position = [ax.Position(1) mypos2(2)+0.012 ax.Position(3) mypos2(4)+0.035];
text(textx(1)+1,1.009,'b)','FontSize',LocalFontSize+2,'FontName',get(gca,'FontName'));


% legend
% hl=legend([eSimData,pH1],...
%     {katrinsim,'3\nu+1 model'},'Location','southeast','box','off');
hl=legend(pH1,...
   katrinsim,'Location','northeast','box','off');
hl.FontSize = LocalFontSize-0.5;

% axleg=axes('Position',get(gca,'Position'),'Visible','Off');
% hl2 = legend(axleg,pH0,'3\nu model','box','off','Location','southeast');
% hl2.FontSize = LocalFontSize-0.5;
PRLFormat;

%% 
s3 = subplot(4,1,4);
bT = bar(x,Time./(60*60),0.5,'FaceColor',rgb('Black'),'EdgeColor','none');
bT.BarWidth = 0.7;

% apperance: legend, labels etc.
xlabel(xstr);
ylh = ylabel('Time (h)');
ylh.Position(1) = ax2.YLabel.Position(1);%
yl1.Position(1) = ax2.YLabel.Position(1);
ylh.Position(2) =  35;
ylim([0 70])
yticks([0 35 70])
PRLFormat;
set(gca,'FontSize',LocalFontSize);
set(get(gca,'XLabel'),'FontSize',LocalFontSize+6);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+6);
text(textx(1)+1,60,'c)','FontSize',LocalFontSize+2,'FontName',get(gca,'FontName'));

linkaxes([s1,s2,s3],'x');
xlim([-43 140]);
ax3 = gca;
mypos2 = ax3.Position;
ax3.Position = [ax.Position(1) mypos2(2)-0.01 ax.Position(3) mypos2(4)+0.035];
% hl2.Position(2) = hl.Position(2);
 
if strcmp(SavePlot,'ON')
    plotdir  = strrep(savedir,'results','plots');
    MakeDir(plotdir);
    plotname = sprintf('%sksn2_spectrum.pdf',plotdir);
end
export_fig(plotname);
fprintf('save plot to %s \n',plotname);