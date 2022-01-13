% KSN2 paper plot, but with different middle panel
% instead of simluation: data!
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
    fprintf('run ksn2_spectrumPaperPlot.m to generate file %s \n',savefile);
    return
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

%% subplot 1
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
    sprintf('3\\nu+1 best-fit model'),...
    sprintf('Active branch: {\\itm}_\\nu^2 = %.1f eV^2',dbf.FitResult.par(1)),...
    sprintf('Sterile branch: {\\itm}_4^2 = %.1f eV^2, |{\\itU}_{e4}|^2 = %.3f',dbf.mNu4Sq_bf,dbf.sin2T4_bf));
PrettyLegendFormat(leg);
leg.FontSize = LocalFontSize-0.5;
set(gca,'yscale','log');
yl1 = ylabel('Count rate (cps)','FontSize',LocalFontSize+6);

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
Mode = 'Ratio2';

s2= subplot(4,1,3);

MarkerStyle0 = {'d','MarkerSize',5,'LineWidth',1.5,'Color',rgb('LimeGreen')};
MarkerStyle1 = {'.','MarkerSize',12,'LineWidth',1.5,'Color',rgb('DodgerBlue')};

if strcmp(Mode,'Residuals')
    plot(linspace(-50,150,10),zeros(10,1),':','Color',rgb('Black'),'LineWidth',2);
    hold on;
    pH0 = plot(x,(DataRate-NullBfRate)./DataRateErr,MarkerStyle0{:});
    pH1 = plot(x,(DataRate-H1BfRate)./DataRateErr,MarkerStyle1{:});
    
     yStr = sprintf('Residuals (\\sigma)');
elseif strcmp(Mode,'Ratio')
    plot(linspace(-50,150,10),ones(10,1),':','Color',rgb('Black'),'LineWidth',2); 
    hold on;
    pH0  = errorbar(x,DataRate./NullBfRate,DataRateErr./NullBfRate,MarkerStyle0{:},'CapSize',0);
    pH1  = errorbar(x,DataRate./H1BfRate,DataRateErr./H1BfRate,MarkerStyle1{:},'CapSize',0);    
     yStr = sprintf('Ratio');
elseif strcmp(Mode,'Ratio2')
     plot(linspace(-50,150,10),ones(10,1),':','Color',rgb('Black'),'LineWidth',2); 
    hold on;
    pH0  = plot(x,NullBfRate./NullBfRate,'LineWidth',1.5,'Color',rgb('LimeGreen'));
    pH1  = plot(x,H1BfRate./NullBfRate,'LineWidth',1.5,'Color',rgb('DodgerBlue'));  
    pD   = errorbar(x,DataRate./NullBfRate,DataRateErr./NullBfRate,'.k','CapSize',0,'MarkerSize',12,'LineWidth',1.5);
     yStr = sprintf('Ratio');
end

hold off
if strcmp(qUDisp,'Abs')
    xticks(myxticks);
    ax = gca;
    ax.XAxis.Exponent = 0;
end

PRLFormat;
set(gca,'FontSize',LocalFontSize);
yl2 = ylabel(yStr,'FontSize',LocalFontSize+6);


if strcmp(Mode,'Residuals')
ylim([-3,3]);
yb = -2.1;
elseif contains(Mode,'Ratio')
    ylim([0.98 1.02]);
    yb = 0.9865;
 
end
xlim([-43 140]);

ax2 = gca;
%axis
mypos2 = ax2.Position;
ax2.Position = [ax.Position(1) mypos2(2)+0.012 ax.Position(3) mypos2(4)+0.035];
text(textx(1)+1,yb,'b)','FontSize',LocalFontSize+2,'FontName',get(gca,'FontName'));


hl=legend([pH0,pH1],...
    sprintf('3\\nu best-fit model'),...
    sprintf('3\\nu+1 best-fit model'),...
   'Location','southeast','box','off');
hl.FontSize = LocalFontSize-0.5;
%hl.NumColumns = 2;

%% mtd
s3 = subplot(4,1,4);
bT = bar(x,Time./(60*60),0.5,'FaceColor',rgb('Black'),'EdgeColor','none');
bT.BarWidth = 0.7;

% apperance: legend, labels etc.
xlabel(xstr);
ylim([0 70])
yticks([0 35 70])

yl3 = ylabel('Time (h)');
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

% align y labels
yl2.Position(1) = yl1.Position(1);
yl3.Position(1) = yl1.Position(1);
yl3.Position(2) =  35;



if strcmp(SavePlot,'ON')
    plotdir  = strrep(savedir,'results','plots');
    MakeDir(plotdir);
    plotname = sprintf('%sksn2_spectrum_Data%s.pdf',plotdir,Mode);
    
    export_fig(plotname);
    fprintf('save plot to %s \n',plotname);
end