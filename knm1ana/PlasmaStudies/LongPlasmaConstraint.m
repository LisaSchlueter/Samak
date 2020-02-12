%% Attempt to constrain plasma bias with FPD data 
%% Proxy = rate @-200eV

%% Definition of a MultiRunAnalysis
option = {...
    'DataType','Real',...
    'RunList','KNM1_m149mvRW',...
    'exclDataStart',1,...
    'fixPar','1 5 6 7 8 9 10',...
    'chi2','chi2Stat',...
    'ELossFlag','KatrinD2',...
    'StackTolerance',0.2,...
    'NonPoissonScaleFactor',1.5,...
    'Debug','ON',...
    'AnaFlag','StackPixel',...
    'PixList',[]...
    };
PC_NC = MultiRunAnalysis(option{:},'TestRateCorr','OFF');

%% Extract Reference Rate @-200V 
%% Apply EROS Rate Correction
rateNC = zeros(numel(PC_NC.RunData.qU),numel(PC_NC.RunList));
for i=1:1:numel(PC_NC.RunData.qU)
    rateNC(i,:) = PC_NC.SingleRunData.TBDIS(i,:);
end
RaterefNC       = rateNC(1,:)./PC_NC.SingleRunData.TimeperSubRun(1,:);
RateErrrefNC    = sqrt(rateNC(1,:))./PC_NC.SingleRunData.TimeperSubRun(1,:);
RaterefC        = RaterefNC.*PC_NC.R200RateErosCorrection;
RateErrrefC     = sqrt(rateNC(1,:))./PC_NC.SingleRunData.TimeperSubRun(1,:);
qUref           = PC_NC.SingleRunData.qU(1,:);

%% Plot Rates Before / After Correction
fig2000 = figure('Renderer','opengl');
set(fig2000,'units','normalized','pos',[0.1, 0.1,0.9,0.5]);
pc=plot(PC_NC.SingleRunData.StartTimeStamp,RaterefC,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('IndianRed'),'MarkerSize',8);
hold on
pnc=plot(PC_NC.SingleRunData.StartTimeStamp,RaterefNC,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('SlateGray'),'MarkerSize',8);
hold off
legend([pnc,pc],'Raw Data', 'Samak Corrected Data');legend boxoff;
ylabel('Reference Rate (cps)');
xlabel('Start Time Stamp (s)');
xlim([min(PC_NC.SingleRunData.StartTimeStamp) max(PC_NC.SingleRunData.StartTimeStamp)]);
PrettyFigureFormat;set(gca,'FontSize',20);

%% Plot Rates Before / After Correction with errorbars
fig2000 = figure('Renderer','opengl');
set(fig2000,'units','normalized','pos',[0.1, 0.1,0.9,0.7]);
 pc=errorbar(PC_NC.SingleRunData.StartTimeStamp,RaterefC,RateErrrefC,'s',...
                 'MarkerSize',5,'MarkerEdgeColor', rgb('IndianRed'),'MarkerFaceColor', rgb('IndianRed'),'Color',rgb('IndianRed'));
pc.CapSize = .2;
% pc=errorbar(PC_NC.SingleRunData.StartTimeStamp,RaterefC,RateErrrefC,'s',...
%                 'MarkerSize',3,'MarkerEdgeColor', rgb('White'),'MarkerFaceColor', rgb('White'),'Color',rgb('White'));
hold on
pnc=errorbar(PC_NC.SingleRunData.StartTimeStamp,RaterefNC,RateErrrefNC,'o',...
                'MarkerSize',5,'MarkerEdgeColor', rgb('SlateGray'),'MarkerFaceColor', rgb('SlateGray'),'Color',rgb('SlateGray'));
pnc.CapSize = .2;
hold off
legend([pnc,pc],'Raw Data', 'Samak Corrected Rates');legend boxoff;
ylabel('Reference Rate (cps)');
xlabel('Start Time Stamp (s)');
xlim([min(PC_NC.SingleRunData.StartTimeStamp) max(PC_NC.SingleRunData.StartTimeStamp)]);
PrettyFigureFormat;set(gca,'FontSize',20);

%% Rates Histograms
% Plot Rate / Corrected Rate
fig2000 = figure('Renderer','opengl');
set(fig2000,'units','normalized','pos',[0.1, 0.1,0.8,0.8]);
% not corrected
h1=histfit((RaterefNC)-mean((RaterefNC)));%,'facecolor',rgb('IndianRed'),'facealpha',.6,'edgecolor','none');
set(h1(1),'facecolor',rgb('SlateGrey'),'facealpha',.4,'edgecolor','none'); set(h1(2),'color',rgb('SlateGrey'));
hold on 
% corrected
h2=histfit((RaterefC)-mean((RaterefC)));%,'facecolor',rgb('IndianRed'),'facealpha',.6,'edgecolor','none');
set(h2(1),'facecolor',rgb('IndianRed'),'facealpha',.8,'edgecolor','none'); set(h2(2),'color',rgb('IndianRed'));
pd = fitdist(((RaterefC)-((RaterefC)))','Normal');
hold off
xlabel('Reference Rate at E_0-200V (cps)');
grid on
box off
nc=sprintf('Uncorrected Rates: std = %.1f cps',std((RaterefNC)));
c=sprintf('Corrected Rates: std = %.1f cps',std((RaterefC)));
legend([h1(2), h2(2)],nc,c,'location','northwest');
legend boxoff
PrettyFigureFormat;set(gca,'FontSize',20);

%% Compute a function relating qU-offset and Reference Rate
%% Similar to Pro-KATRIN
PC_NC.fixPar              = '1 5 6 7 8 ; fix 9 10 11';
PC_NC.exclDataStart       = 2;
PC_NC.ModelObj.qUOffset_i = 0;
PC_NC.ModelObj.qUOffset   = 0;
PC_NC.ModelObj.ComputeTBDDS('qUOffset_bias',0);
PC_NC.ModelObj.ComputeTBDIS;
PC_NC.Fit; % PC_NC.PlotFit('Mode','Rate','saveplot','png')
RefRateZeroOffset = PC_NC.ModelObj.TBDIS(1)./PC_NC.ModelObj.TimeSec./PC_NC.ModelObj.qUfrac(1);
fprintf('qU=%.3f - qUOffset = %.1f V (%f / %f) --> %.5f\n',...
    PC_NC.ModelObj.qU(1),...
    PC_NC.ModelObj.qUOffset,...
    PC_NC.ModelObj.qUOffset_i,...
    PC_NC.ModelObj.qUOffset_bias,...
    PC_NC.ModelObj.TBDIS(1)./PC_NC.ModelObj.TimeSec./PC_NC.ModelObj.qUfrac(1));
PC_NC.ModelObj.qUOffset_bias = 0;

OffSets         = -3:0.1:3;
ModelRateOffset = zeros(1,numel(OffSets));
counter=0;
for i=OffSets
counter=counter+1;
PC_NC.ModelObj.qUOffset_i = 0;
PC_NC.ModelObj.qUOffset   = 0;
PC_NC.ModelObj.ComputeTBDDS('qUOffset_bias',i);
PC_NC.ModelObj.ComputeTBDIS;
ModelRateOffset(counter) = (PC_NC.ModelObj.TBDIS(1)./PC_NC.ModelObj.TimeSec./PC_NC.ModelObj.qUfrac(1)*(PC_NC.FitResult.par(4)+1));
fprintf('qU=%.3f - qUOffset = %.1f V (%f / %f) --> %.5f\n',...
    PC_NC.ModelObj.qU(1),...
    PC_NC.ModelObj.qUOffset,...
    PC_NC.ModelObj.qUOffset_i,...
    PC_NC.ModelObj.qUOffset_bias,...
    ModelRateOffset(counter));
PC_NC.ModelObj.qUOffset_bias = 0;
end

%% Plot
fig2000 = figure('Renderer','opengl');
set(fig2000,'units','normalized','pos',[0.1, 0.1,0.5,0.4]);
p=plot(OffSets,ModelRateOffset,'d','MarkerSize',10);
PrettyFigureFormat
xlabel('qU-offset (V)');
ylabel('Rate @-200V (cps)');
P = polyfit(OffSets,ModelRateOffset,1);
yfit = P(1)*OffSets+P(2);
hold on;
f=plot(OffSets,yfit,'r--','LineWidth',4,'Color',rgb('Amethyst'));
legend([p,f],'Samak Simulation', 'Linear Model');legend boxoff;
hold off
PrettyFigureFormat;set(gca,'FontSize',20);

%% Interpolation Rate/qUoffset (pro-KATRIN)
offset = @(rate) interp1(yfit,(yfit - P(2))./P(1),rate);
% Plot Rate / Corrected Rate
fig2000 = figure('Renderer','opengl');
set(fig2000,'units','normalized','pos',[0.1, 0.1,0.8,0.8]);
% not corrected
h1=histfit(offset(RaterefNC)-mean(offset(RaterefNC)));%,'facecolor',rgb('IndianRed'),'facealpha',.6,'edgecolor','none');
set(h1(1),'facecolor',rgb('SlateGrey'),'facealpha',.4,'edgecolor','none'); set(h1(2),'color',rgb('SlateGrey'));
hold on 
% corrected
h2=histfit(offset(RaterefC)-mean(offset(RaterefC)));%,'facecolor',rgb('IndianRed'),'facealpha',.6,'edgecolor','none');
set(h2(1),'facecolor',rgb('IndianRed'),'facealpha',.8,'edgecolor','none'); set(h2(2),'color',rgb('IndianRed'));
pd = fitdist((offset(RaterefC)-mean(offset(RaterefC)))','Normal');
% Simulation of the expected fluctuation
SimRefRate = PC_NC.ModelObj.TBDIS(1) + sqrt(PC_NC.ModelObj.TBDIS(1)).*randn(1,numel(RaterefC));
SimRefRate = SimRefRate./PC_NC.ModelObj.TimeSec./PC_NC.ModelObj.qUfrac(1)*(PC_NC.FitResult.par(4)+1);
h3=histfit(offset(SimRefRate)-mean(offset(SimRefRate)));%,'facecolor',rgb('IndianRed'),'facealpha',.6,'edgecolor','none');
set(h3(1),'facecolor',rgb('CadetBlue'),'facealpha',.6,'edgecolor','none'); set(h3(2),'color',rgb('CadetBlue'));
hold off
xlabel('Plasma Bias Voltage (V)');
grid on
box off
nc=sprintf('qU-offset (Uncorrected Rates) std = %.1f mV',std(offset(RaterefNC))*1000);
c=sprintf('qU-offset (Corrected Rates) std = %.1f mV',std(offset(RaterefC))*1000);
s=sprintf('Stat. Expected: std = %.1f mV',std(offset(SimRefRate))*1000);
legend([h1(2), h2(2), h3(2)],nc,c,s,'location','northwest');
%xlim([-0.3 0.3]);
legend boxoff
PrettyFigureFormat;set(gca,'FontSize',20);


%% Periodicity Check
[ff,PP,prob] = lomb(double(PC_NC.SingleRunData.StartTimeStamp)',RaterefC',2,.5);
figure(123456)
plot(1./ff./86400,PP,...
    'ks-.','MarkerSize',12,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
xlabel('Period (d^{-1})');
ylabel('Strength');
[Pmax,jmax] = max(PP);
[Psort Pindex] = sort(PP);
for i=1:numel(Pindex)
    if prob(Pindex(i))<1e-4
        disp(['Most significant period is ',num2str(1/ff(Pindex(i))/86400),...
            ' days with FAP of ',num2str(prob(Pindex(i)))]);
    end
end
PrettyFigureFormat ; set(gca,'FontSize',16);

X = double(PC_NC.SingleRunData.StartTimeStamp);
X = X';
Y = (RaterefC-mean(RaterefC))';
Y = (mean(RaterefC) .* (1+10*sin(1e-5*X')))';
myeq = '(1+a*sin(e*x+b))*c';
%myeq='sin1';
figure(555)
foo = fit(X,Y,myeq)
scatter(X,Y);
hold on
plot(foo)