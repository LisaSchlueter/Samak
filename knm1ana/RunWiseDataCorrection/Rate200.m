option = {...
    'DataType','Real',...
    'RunList','KNM1_m149mvRW',...
    'exclDataStart',2,...
    'fixPar','1 5 6 7 8 9 10',...
    'chi2','chi2Stat',...
    'ELossFlag','KatrinD2',...
    'StackTolerance',0.2,...
    'NonPoissonScaleFactor',1.5,...
    'Debug','ON',...
    'AnaFlag','StackPixel',...
    'PixList',[]...
    };

SPR=MultiRunAnalysis(option{:});

%% Extract Rates
rate = zeros(numel(SPR.RunData.qU),numel(SPR.RunList));
for i=1:1:numel(SPR.RunData.qU)
        rate(i,:)=SPR.SingleRunData.TBDIS(i,:);%./mean(SPR.SingleRunData.qUfrac(i,:,SPR.PixList),3);
end
Rateref = rate(2,:)./SPR.SingleRunData.TimeperSubRun(2,:);
qU200   = SPR.SingleRunData.qU(2,:);

%% Plot Rate / RhoD / TT Verus Times
fig1000 = figure('Renderer','opengl');
set(fig1000,'units','normalized','pos',[0.1, 0.1,0.9,0.8]);
s1=subplot(3,1,1)
plot(SPR.SingleRunData.StartTimeStamp,Rateref,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('DarkGreen'),'MarkerSize',8);
ylabel('Rate (cps)');
xlabel('Start Time Stamp (s)');
PrettyFigureFormat;set(gca,'FontSize',20);
s2=subplot(3,1,2)
plot(SPR.SingleRunData.StartTimeStamp,SPR.SingleRunData.WGTS_CD_MolPerCm2,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('DarkBlue'),'MarkerSize',8);
ylabel('\rho d (mol/cm^2)');
xlabel('Start Time Stamp (s)');
PrettyFigureFormat;set(gca,'FontSize',20);
s3=subplot(3,1,3)
plot(SPR.SingleRunData.StartTimeStamp,SPR.SingleRunData.WGTS_MolFrac_TT,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('DarkRed'),'MarkerSize',8);
ylabel('[T-T]');
xlabel('Start Time Stamp (s)');
PrettyFigureFormat;set(gca,'FontSize',20);
            linkaxes([s1,s2,s3],'x');

            
%% Plot Rate Verus RhoD
fig1001 = figure('Renderer','opengl');
set(fig1001,'units','normalized','pos',[0.1, 0.1,0.9,0.5]);
plot(SPR.SingleRunData.WGTS_CD_MolPerCm2,Rateref,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('DarkBlue'),'MarkerSize',8);
ylabel('rate (cps)');
xlabel('\rho d (mol/cm^2)');
PrettyFigureFormat
set(gca,'FontSize',24);
P = polyfit(SPR.SingleRunData.WGTS_CD_MolPerCm2,Rateref,1);
yfit = P(1)*SPR.SingleRunData.WGTS_CD_MolPerCm2+P(2);
hold on;
plot(SPR.SingleRunData.WGTS_CD_MolPerCm2,yfit,'r-.','LineWidth',5,'Color',rgb('SteelBlue'));
hold off

%% Plot Rate Verus qU
fig1001 = figure('Renderer','opengl');
set(fig1001,'units','normalized','pos',[0.1, 0.1,0.9,0.5]);
plot(qU200,Rateref,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('DarkGrey'),'MarkerSize',8);
ylabel('rate (cps)');
xlabel('qU200');
PrettyFigureFormat
set(gca,'FontSize',24);
P = polyfit(qU200,Rateref,1);
yfit = P(1)*qU200+P(2);
hold on;
plot(qU200,yfit,'r-.','LineWidth',5,'Color',rgb('SteelBlue'));
hold off

%% Plot Rate Verus TT
fig1001 = figure('Renderer','opengl');
set(fig1001,'units','normalized','pos',[0.1, 0.1,0.9,0.5]);
plot(SPR.SingleRunData.WGTS_MolFrac_TT,Rateref,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('DarkRed'),'MarkerSize',8);
ylabel('rate (cps)');
xlabel('[TT]');
PrettyFigureFormat
set(gca,'FontSize',24);
P = polyfit(SPR.SingleRunData.WGTS_MolFrac_TT,Rateref,1);
yfit = P(1)*SPR.SingleRunData.WGTS_MolFrac_TT+P(2);
hold on;
plot(SPR.SingleRunData.WGTS_MolFrac_TT,yfit,'r-.','LineWidth',5,'Color',rgb('SteelBlue'));
hold off

%% Plot Rate Verus HT
fig1001 = figure('Renderer','opengl');
set(fig1001,'units','normalized','pos',[0.1, 0.1,0.9,0.5]);
plot(SPR.SingleRunData.WGTS_MolFrac_HT,Rateref,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('Amethyst'),'MarkerSize',8);
ylabel('rate (cps)');
xlabel('[HT]');
PrettyFigureFormat
set(gca,'FontSize',24);
P = polyfit(SPR.SingleRunData.WGTS_MolFrac_HT,Rateref,1);
yfit = P(1)*SPR.SingleRunData.WGTS_MolFrac_HT+P(2);
hold on;
plot(SPR.SingleRunData.WGTS_MolFrac_HT,yfit,'r-.','LineWidth',5,'Color',rgb('SteelBlue'));
hold off

%% Plot Rate Verus DT
fig1001 = figure('Renderer','opengl');
set(fig1001,'units','normalized','pos',[0.1, 0.1,0.9,0.5]);
plot(SPR.SingleRunData.WGTS_MolFrac_DT,Rateref,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('IndianRed'),'MarkerSize',8);
ylabel('rate (cps)');
xlabel('[HT]');
PrettyFigureFormat
set(gca,'FontSize',24);
P = polyfit(SPR.SingleRunData.WGTS_MolFrac_DT,Rateref,1);
yfit = P(1)*SPR.SingleRunData.WGTS_MolFrac_DT+P(2);
hold on;
plot(SPR.SingleRunData.WGTS_MolFrac_DT,yfit,'r-.','LineWidth',5,'Color',rgb('SteelBlue'));
hold off

%% Plot TT Verus RhoD
fig1001 = figure('Renderer','opengl');
set(fig1001,'units','normalized','pos',[0.1, 0.1,0.9,0.5]);
plot(SPR.SingleRunData.WGTS_CD_MolPerCm2,SPR.SingleRunData.WGTS_MolFrac_TT,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('DarkBlue'),'MarkerSize',8);
ylabel('[TT]');
xlabel('\rho d (mol/cm^2)');
PrettyFigureFormat
set(gca,'FontSize',24);
P = polyfit(SPR.SingleRunData.WGTS_CD_MolPerCm2,SPR.SingleRunData.WGTS_MolFrac_TT,1);
yfit = P(1)*SPR.SingleRunData.WGTS_CD_MolPerCm2+P(2);
hold on;
plot(SPR.SingleRunData.WGTS_CD_MolPerCm2,yfit,'r-.','LineWidth',5,'Color',rgb('SteelBlue'));
hold off


%% Compute Mean Coefficients
MeanRhoD  = mean(SPR.SingleRunData.WGTS_CD_MolPerCm2);

MeanTT    = mean(SPR.SingleRunData.WGTS_MolFrac_TT);

MeanHT    = mean(SPR.SingleRunData.WGTS_MolFrac_HT);

MeanDT    = mean(SPR.SingleRunData.WGTS_MolFrac_DT);

MeanqU200 = mean(qU200);

%% Compute Standard Deviation Coefficients
StdRate  = std(Rateref);

StdRhoD  = std(SPR.SingleRunData.WGTS_CD_MolPerCm2);

StdTT    = std(SPR.SingleRunData.WGTS_MolFrac_TT);

StdHT    = std(SPR.SingleRunData.WGTS_MolFrac_HT);

StdDT    = std(SPR.SingleRunData.WGTS_MolFrac_DT);

StdqU200 = std(qU200);

%% Compute Correlation Coefficients

% Rate VS RhoD
CorrRateRhoD = corrcoef(Rateref,SPR.SingleRunData.WGTS_CD_MolPerCm2);
CorrRateRhoD = CorrRateRhoD(1,2);

% Rate VS qU200
CorrRateqU200 = corrcoef(Rateref,qU200);
CorrRateqU200 = CorrRateqU200(1,2);

% Rate VS TT
CorrRateTT = corrcoef(Rateref,SPR.SingleRunData.WGTS_MolFrac_TT);
CorrRateTT = CorrRateTT(1,2);

% Rate VS HT
CorrRateHT = corrcoef(Rateref,SPR.SingleRunData.WGTS_MolFrac_HT);
CorrRateHT = CorrRateHT(1,2);

% Rate VS DT
CorrRateDT = corrcoef(Rateref,SPR.SingleRunData.WGTS_MolFrac_DT);
CorrRateDT = CorrRateDT(1,2);

% RhoD Vs TT
CorrRhoDTT = corrcoef(SPR.SingleRunData.WGTS_CD_MolPerCm2,SPR.SingleRunData.WGTS_MolFrac_TT);
CorrRhoDTT = CorrRhoDTT(1,2);

% RhoD Vs HT
CorrRhoDHT = corrcoef(SPR.SingleRunData.WGTS_CD_MolPerCm2,SPR.SingleRunData.WGTS_MolFrac_HT);
CorrRhoDHT = CorrRhoDHT(1,2);

% RhoD Vs DT
CorrRhoDDT = corrcoef(SPR.SingleRunData.WGTS_CD_MolPerCm2,SPR.SingleRunData.WGTS_MolFrac_DT);
CorrRhoDDT = CorrRhoDDT(1,2);

% RhoD Vs qU200
CorrRhoDqU200 = corrcoef(SPR.SingleRunData.WGTS_CD_MolPerCm2,qU200);
CorrRhoDqU200 = CorrRhoDqU200(1,2);

% Rate VS qU200
CorrTTqU200 = corrcoef(SPR.SingleRunData.WGTS_MolFrac_TT,qU200);
CorrTTqU200 = CorrTTqU200(1,2);

CorrHTqU200 = corrcoef(SPR.SingleRunData.WGTS_MolFrac_HT,qU200);
CorrHTqU200 = CorrHTqU200(1,2);

CorrDTqU200 = corrcoef(SPR.SingleRunData.WGTS_MolFrac_DT,qU200);
CorrDTqU200 = CorrDTqU200(1,2);

% TT VS HT
CorrTTHT = corrcoef(SPR.SingleRunData.WGTS_MolFrac_TT,SPR.SingleRunData.WGTS_MolFrac_HT);
CorrTTHT = CorrTTHT(1,2);

% DT VS HT
CorrHTDT = corrcoef(SPR.SingleRunData.WGTS_MolFrac_HT,SPR.SingleRunData.WGTS_MolFrac_DT);
CorrHTDT = CorrHTDT(1,2);

% TT VS DT
CorrTTDT = corrcoef(SPR.SingleRunData.WGTS_MolFrac_TT,SPR.SingleRunData.WGTS_MolFrac_DT);
CorrTTDT = CorrTTDT(1,2);

%% Definition of Parameters in use
p1=SPR.SingleRunData.WGTS_CD_MolPerCm2;
p2=SPR.SingleRunData.WGTS_MolFrac_HT;
p3=SPR.SingleRunData.WGTS_MolFrac_TT;
p4=qU200;
p5=SPR.SingleRunData.WGTS_MolFrac_DT;

m1=MeanRhoD;
m2=MeanHT;
m3=MeanTT;
m4=MeanqU200;
m5=MeanDT;

std1=StdRhoD;
std2=StdHT;
std3=StdTT;
std4=StdqU200;
std5=StdDT;

rhoP1=CorrRateRhoD;
rhoP2=CorrRateHT;
rhoP3=CorrRateTT;
rhoP4=CorrRateqU200;
rhoP5=CorrRateDT;

rho12=CorrRhoDHT;
rho13=CorrRhoDTT;
rho23=CorrTTHT;
rho14=CorrRhoDqU200;
rho24=CorrHTqU200;
rho34=CorrTTqU200;
rho15=CorrRhoDDT;
rho25=CorrHTDT;
rho35=CorrTTDT;
rho45=CorrDTqU200;

%% Correction Coefficientss
TM = [...
      std1         rho12*std2     rho13*std3    rho14*std4  rho15*std5      ; ...
      rho12*std1   std2           rho23*std3    rho24*std4  rho25*std5      ; ...
      rho13*std1   rho23*std2     std3          rho34*std4  rho35*std5      ; ...
      rho14*std1   rho24*std2     rho34*std3    std4        rho45*std5      ; ...
      rho15*std1   rho25*std2     rho35*std3    rho45*std4  std5            ; ...
    ];
TX = StdRate * [rhoP1 rhoP2 rhoP3 rhoP4 rhoP5]';

TS = TM\TX;

RaterefCorrected = Rateref  - TS(1) * (p1-m1) - TS(2) * (p2-m2) - TS(3) * (p3-m3) - TS(4) * (p4-m4) - TS(5) * (p5-m5);

%% Rate / Corrected Rate
fig2000 = figure('Renderer','opengl');
set(fig2000,'units','normalized','pos',[0.1, 0.1,0.9,0.8]);
subplot(3,1,1)
plot(SPR.SingleRunData.StartTimeStamp,Rateref,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('Black'),'MarkerSize',8);
ylabel('Rate (cps)');
xlabel('Start Time Stamp (s)');
PrettyFigureFormat;set(gca,'FontSize',20);
subplot(3,1,2)
plot(SPR.SingleRunData.StartTimeStamp,RaterefCorrected,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('Red'),'MarkerSize',8);
ylabel('Corrected Rate (cps)');
xlabel('Start Time Stamp (s)');
PrettyFigureFormat;set(gca,'FontSize',20);
subplot(3,1,3)
CorrectionFactor = RaterefCorrected./Rateref;
plot(SPR.SingleRunData.StartTimeStamp,CorrectionFactor,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('Blue'),'MarkerSize',8);
ylabel('Correction');
xlabel('Start Time Stamp (s)');
PrettyFigureFormat;set(gca,'FontSize',20);

%% Rate / Corrected Rate
subrun=1;
r=rate(subrun,:)./SPR.SingleRunData.TimeperSubRun(subrun,:);
rc=rate(subrun,:)./SPR.SingleRunData.TimeperSubRun(subrun,:).*CorrectionFactor;
           
fig2000 = figure('Renderer','opengl');
set(fig2000,'units','normalized','pos',[0.1, 0.1,0.8,0.8]);
%histogram(r,'facecolor',rgb('CadetBlue'),'facealpha',.6,'edgecolor','none');
h=histfit(r);%,'facecolor',rgb('IndianRed'),'facealpha',.6,'edgecolor','none');
set(h(1),'facecolor',rgb('CadetBlue'),'facealpha',.6,'edgecolor','none'); set(h(2),'color',rgb('CadetBlue'));
hold on 
h=histfit(rc);%,'facecolor',rgb('IndianRed'),'facealpha',.6,'edgecolor','none');
set(h(1),'facecolor',rgb('IndianRed'),'facealpha',.6,'edgecolor','none'); set(h(2),'color',rgb('IndianRed'));
pd = fitdist(rc','Normal');
hold off
xlabel('Rate (cps)');
grid on
box off
nc=sprintf('Not Corrected: std/mean=%.3f %%',std(r)./mean(r)*100);
c=sprintf('Corrected: std/mean=%.3f %% \n Poisson Expected =%.3f %%',std(rc)./mean(r)*100,1/sqrt(pd.mu*mean(SPR.SingleRunData.TimeperSubRun(subrun,:)))*100);
legend(nc,c,'location','northwest')
legend boxoff
PrettyFigureFormat;set(gca,'FontSize',20);


%% Rate / Corrected Rate
fig2000 = figure('Renderer','opengl');
set(fig2000,'units','normalized','pos',[0.1, 0.1,0.9,0.5]);
plot(SPR.SingleRunData.StartTimeStamp,Rateref,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('Black'),'MarkerSize',8);
hold on
plot(SPR.SingleRunData.StartTimeStamp,RaterefCorrected,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('Red'),'MarkerSize',8);
hold off
ylabel('Corrected Rate (cps)');
xlabel('Start Time Stamp (s)');
PrettyFigureFormat;set(gca,'FontSize',20);

%% Correction Factors
fig2000 = figure('Renderer','opengl');
set(fig2000,'units','normalized','pos',[0.1, 0.1,0.8,0.9]);
s1=subplot(6,1,1)
f1=plot(SPR.SingleRunData.StartTimeStamp,...
    Rateref  - TS(1) * (p1-m1),'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('Red'),'MarkerSize',6);
l=legend(f1,'\rhod corrected','Location','best');
PrettyFigureFormat;set(gca,'FontSize',12);
ylabel('Rate (cps)');
s2=subplot(6,1,2)
f2=plot(SPR.SingleRunData.StartTimeStamp,...
    Rateref  - TS(1) * (p1-m1) - TS(2) * (p2-m2) ,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('Red'),'MarkerSize',6);
l=legend(f2,'\rhod + HT corrected','Location','best');
PrettyFigureFormat;set(gca,'FontSize',12);
ylabel('Rate (cps)');
s3=subplot(6,1,3)
f3=plot(SPR.SingleRunData.StartTimeStamp,...
    Rateref  - TS(1) * (p1-m1) - TS(3) * (p3-m3) ,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('Red'),'MarkerSize',6);
l=legend(f3,'\rhod + TT corrected','Location','best');
PrettyFigureFormat;set(gca,'FontSize',12);
ylabel('Rate (cps)');
s4=subplot(6,1,4)
f4=plot(SPR.SingleRunData.StartTimeStamp,...
    Rateref  - TS(1) * (p1-m1) - TS(5) * (p5-m5) ,'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('Red'),'MarkerSize',6);
l=legend(f4,'\rhod + DT corrected','Location','best');
PrettyFigureFormat;set(gca,'FontSize',12);
ylabel('Rate (cps)');
s5=subplot(6,1,5)
f5=plot(SPR.SingleRunData.StartTimeStamp,...
    Rateref - TS(1) * (p1-m1) - TS(4) * (p4-m4),'s',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('Red'),'MarkerSize',6);
l=legend(f5,'\rhod + qU corrected','Location','best');
ylabel('Rate (cps)');
xlabel('Start Time Stamp (s)');
PrettyFigureFormat;set(gca,'FontSize',12);

return;

%% Test on Stacked Spectrum
option = {...
    'DataType','Real',...
    'RunList','KNM1_m149mvRW',...
    'exclDataStart',2,...
    'fixPar','1 5 6 7 8 9 10',...
    'chi2','chi2Stat',...
    'ELossFlag','KatrinD2',...
    'StackTolerance',0.2,...
    'NonPoissonScaleFactor',1.5,...
    'Debug','ON',...
    'AnaFlag','StackPixel',...
    'PixList',[]...
    };
SPR_OFF=MultiRunAnalysis(option{:},'TestRateCorr','OFF');
SPR_ActivityScaling=MultiRunAnalysis(option{:},'TestRateCorr','ActivityScaling');
SPR_EROS=MultiRunAnalysis(option{:},'TestRateCorr','EROS');

%%
for i=1:1:15
    SPR_OFF.exclDataStart=i;
    SPR_EROS.exclDataStart=i;
    SPR_OFF.Fit; offChi2(i) = SPR_OFF.FitResult.chi2min;
    SPR_EROS.Fit;  erosChi2(i) = SPR_EROS.FitResult.chi2min;
end

%%
figure(1)
h2=stairs(SPR_ActivityScaling.RunData.TBDIS(2:end)./SPR_OFF.RunData.TBDIS(2:end));
hold on
h3=stairs(SPR_EROS.RunData.TBDIS(2:end)./SPR_OFF.RunData.TBDIS(2:end));
stairs( 1+1./sqrt(1*SPR_EROS.RunData.TBDIS(2:end)))
stairs( 1-1./sqrt(1*SPR_EROS.RunData.TBDIS(2:end)))
legend([h2,h3],'ActivityScaling','EROS')
hold off
PrettyFigureFormat

%%
figure(111)
subplot(6,1,1)
stairs(SPR.SingleRunData.WGTS_MolFrac_TT);
subplot(6,1,2)
stairs(SPR.SingleRunData.WGTS_MolFrac_HT);
subplot(6,1,3)
stairs(SPR.SingleRunData.WGTS_MolFrac_DT);
subplot(6,1,4)
stairs(SPR.SingleRunData.WGTS_MolFrac_DT+SPR.SingleRunData.WGTS_MolFrac_HT+SPR.SingleRunData.WGTS_MolFrac_TT);
subplot(6,1,5)
stairs(SPR.SingleRunData.WGTS_MolFrac_DT/2+SPR.SingleRunData.WGTS_MolFrac_HT/2+SPR.SingleRunData.WGTS_MolFrac_TT);
subplot(6,1,6)
stairs((SPR.SingleRunData.WGTS_MolFrac_DT/2+SPR.SingleRunData.WGTS_MolFrac_HT/2+SPR.SingleRunData.WGTS_MolFrac_TT).*SPR.SingleRunData.WGTS_CD_MolPerCm2);

%% 
figure(112)
rhodActivity = (SPR.SingleRunData.WGTS_MolFrac_DT/2+SPR.SingleRunData.WGTS_MolFrac_HT/2+SPR.SingleRunData.WGTS_MolFrac_TT).*SPR.SingleRunData.WGTS_CD_MolPerCm2;
stairs(rhodActivity./mean(rhodActivity));
hold on
hold off