% KNM2 - -300V Rate Monitor FPD
%  
% Stakced-pixel Evolution 
%
% Last Modified: 04/02/2020
% T. Lasserre
% 

%% Read Data
DataType  = 'Real';
RunList   = 'KNM2_RW1';
FSDFlag   = 'BlindingKNM2';
ELossFlag = 'KatrinT2';
AnaFlag   = 'StackPixel'; % uniform FPD
RunAnaArg = {'RunList',RunList,'DataType',DataType,...
            'FSDFlag',FSDFlag,'ELossFlag',ELossFlag,'AnaFlag',AnaFlag,'RingMerge','Full'};
MR        = MultiRunAnalysis(RunAnaArg{:});
A         = RingAnalysis('RunAnaObj',MR,'RingList',1:4);

%%
R         = A.MultiObj(1);

%% Slow Control Data
p1 =(R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')./mean((R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')).*R.SingleRunData.WGTS_CD_MolPerCm2'./mean(R.SingleRunData.WGTS_CD_MolPerCm2');
p2 = mean(R.SingleRunData.qU_RM,1); p2=p2-mean(p2);           

%% Time in days
StartTimeStampDays = days(R.SingleRunData.StartTimeStamp-R.SingleRunData.StartTimeStamp(1));

%% Stacked Pixel Data for each patch
count = zeros(A.nRings,numel(A.RunAnaObj.RunList));
rate  = zeros(A.nRings,numel(A.RunAnaObj.RunList));
rateE = zeros(A.nRings,numel(A.RunAnaObj.RunList));
cf    = zeros(A.nRings,numel(A.RunAnaObj.RunList));
for i=1:A.nRings
R           = A.MultiObj(i);
count(i,:)  = R.SingleRunData.TBDIS_RM;
sstime      = mean(R.SingleRunData.qUfrac_RM,1).*R.SingleRunData.TimeSec;
rate(i,:)   = count(i,:)./sstime;
rateE(i,:)  = sqrt(count(i,:))./sstime;
cf(i,:)     = R.RMRateErosCorrectionqUActivity;
end

%% Rate Evolution --> mV equivalent
myMainTitle = sprintf('KATRIN - %s - FPD Rate Evolution 300eV below Endpoint',RunList);
maintitle   = myMainTitle;
savefile    = sprintf('plots/test.png');
fig2      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
    'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

for i=1:A.nRings

rateEquivalentmV      =   - (rate(i,:) - mean(rate(i,:))) ./737.9 *1e3;
rateEquivalentmV_E    =   rateE(i,:).*cf(i,:) ./737.9 *1e3;
CrateEquivalentmV     =   - (rate(i,:).*cf(i,:) - mean(rate(i,:).*cf(i,:))) ./737.9 *1e3;
subplot(A.nRings,1,i)
hc=plot(StartTimeStampDays,CrateEquivalentmV,'s--','Color',rgb('IndianRed'),'LineWidth',1,'MarkerSize',12,'markerfacecolor',rgb('IndianRed'));

% Fit: Constant 
fitType = fittype('a'); % The equation for your fit goes here
[f,gofCons] = fit(StartTimeStampDays',CrateEquivalentmV',fitType,...
    'Weights',(1./rateEquivalentmV_E).^2,...
    'StartPoint',[0],...
    'Robust','Bisquare'); 
ci = confint(f,0.68); uncertaintyCons = (ci(2,1)-ci(1,1))/2;
hold on
hc = plot(StartTimeStampDays,a.*ones(1,numel(StartTimeStampDays)),'--','Color',rgb('Black'),'LineWidth',2);
hold off

% Fit: Linear 
fitType = fittype('a*x + b'); % The equation for your fit goes here
[f,gofLin] = fit(StartTimeStampDays',CrateEquivalentmV',fitType,...
    'Weights',(1./rateEquivalentmV_E).^2,...
    'StartPoint',[0 CrateEquivalentmV(1)],...
    'Robust','Bisquare'); 
ci = confint(f,0.68); uncertaintyLin = (ci(2,1)-ci(1,1))/2;
hold on
hl = plot(StartTimeStampDays,f.a*StartTimeStampDays+f.b,'-','Color',rgb('DarkBlue'),'LineWidth',2);
hold off

% Fit: Quadratic 
fitType = fittype('a*x^2 + b*x + c'); % The equation for your fit goes here
[f,gofQuad] = fit(StartTimeStampDays',CrateEquivalentmV',fitType,...
    'Weights',(1./rateEquivalentmV_E).^2,...
    'StartPoint',[0 0 CrateEquivalentmV(1)],...
    'Robust','Bisquare'); 
ci = confint(f,0.68); uncertaintyQuad = (ci(2,2)-ci(1,2))/2;
hold on
hq = plot(StartTimeStampDays,f.a*StartTimeStampDays.^2 + f.b*StartTimeStampDays+f.c,'-','Color',rgb('DarkGreen'),'LineWidth',2);
hold off


ylabel('\Delta U (meV)');
xlabel('Days since last RW setting');
leg=legend([hc hl hq],...
    sprintf('constant fit - pvalue = %g',chi2pvalue(gofCons.sse,gofConss.dfe)),...
    sprintf('linear fit = %.1f +- %.1f mV/day - pvalue = %g',f.a,uncertaintyLin,chi2pvalue(gofLin.sse,gofLin.dfe)),...
    sprintf('quadratic fit - pvalue = %g',chi2pvalue(gofQuad.sse,gofQuad.dfe)),...
    'location','best','FontSize',8);
leg.Color = 'none'; legend boxoff;
PrettyFigureFormat
end

%% Rate Evolution
myMainTitle = sprintf('KATRIN - %s - FPD Rate Evolution 300eV below Endpoint',RunList);
maintitle   = myMainTitle;
savefile    = sprintf('plots/test.png');
fig2      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
    'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

for i=1:A.nRings

subplot(A.nRings,1,i)
hnc=plot(R.SingleRunData.StartTimeStamp,rate(i,:),'s--','Color',rgb('Green'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('DarkGreen'));
hold on
hc=plot(R.SingleRunData.StartTimeStamp,rate(i,:).*cf(i,:),'s--','Color',rgb('Red'),'LineWidth',1,'MarkerSize',12,'markerfacecolor',rgb('IndianRed'));
hold off
ylabel('cps');
xlabel('Scan Start Time');
leg=legend([hnc hc],...
    sprintf('before correction'),...
    sprintf('after correction'),...
    'location','northeast');
leg.Color = 'none'; legend boxoff;
PrettyFigureFormat
end

