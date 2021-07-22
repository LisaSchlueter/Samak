% plot signal to background as a function of qU
savedir = sprintf('%sknm12Combi/knm12_General/results/',getenv('SamakPath'));
savename = sprintf('%sknm12_DataModel.mat',savedir);

if ~exist(savename,'file') 
    M2 = MultiRunAnalysis('RunList','KNM2_Prompt',...
        'NonPoissonScaleFactor',1.112,...
        'chi2','chi2Stat',...
        'BKG_PtSlope',3*1e-06,...
        'FSDFlag','KNM2_0p1eV',...
        'DopplerEffectFlag','FSD',...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'ELossFlag','KatrinT2A20',...
        'RadiativeFlag','ON',...
        'AngularTFFlag','ON');%,...
       % 'fixPar','mNu E0 Norm Bkg');
    M2.exclDataStart=  M2.GetexclDataStart(40);
   % M2.Fit;
    M2.InitModelObj_Norm_BKG('RecomputeFlag','ON');
    
    M1 = MultiRunAnalysis('RunList','KNM1',...
        'NonPoissonScaleFactor',1.064,...
        'chi2','chi2Stat',...
        'BKG_PtSlope',-2.2*1e-06,...
        'FSDFlag','KNM2_0p1eV',...
        'DopplerEffectFlag','FSD',...
        'FSD_Sigma',0,...
        'ELossFlag','KatrinT2A20',...
        'RadiativeFlag','ON',...
        'AngularTFFlag','ON');
    M1.exclDataStart=  M1.GetexclDataStart(40);
    M1.InitModelObj_Norm_BKG('RecomputeFlag','ON');
    
    StackData1 = M1.RunData;
    Data1      = M1.SingleRunData;
    
    StackData2 = M2.RunData;
    Data2      = M2.SingleRunData;
    
    Model1 = M1.ModelObj;
    Model2 = M2.ModelObj;
    
    % Get background only (model)
    TimeTotSubrun2    = Model2.TimeSec.*Model2.qUfrac;
    TimeAvSubrun2     = Model2.TimeSec.*Model2.qUfrac./Model2.nRuns;
    BkgRate_PngSlope2 = 0.5.*Model2.BKG_PtSlope.*TimeAvSubrun2;
    Bkg_PtSlope2      = BkgRate_PngSlope2.*TimeTotSubrun2;
    BkgCounts2        =  Model2.BKG_RateSec.*TimeTotSubrun2 + Bkg_PtSlope2 ;
    BkgRate2          = Model2.BKG_RateSec+BkgRate_PngSlope2;
    SignalRate2       = Model2.TBDIS./TimeTotSubrun2-BkgRate2;
    SignalCount2      = Model2.TBDIS-BkgCounts2;
    
    TimeTotSubrun1    = Model1.TimeSec.*Model1.qUfrac;
    TimeAvSubrun1     = Model1.TimeSec.*Model1.qUfrac./Model1.nRuns;
    BkgRate_PngSlope1 = 0.5.*Model1.BKG_PtSlope.*TimeAvSubrun1;
    Bkg_PtSlope1      = BkgRate_PngSlope1.*TimeTotSubrun1;
    BkgCounts1        =  Model1.BKG_RateSec.*TimeTotSubrun1 + Bkg_PtSlope1 ;
    BkgRate1          = Model1.BKG_RateSec+BkgRate_PngSlope1;
    SignalRate1       = Model1.TBDIS./TimeTotSubrun1-BkgRate1;
    SignalCount1      = Model1.TBDIS-BkgCounts1;
    
    Sig2Bkg1 = SignalCount1./BkgCounts1; Sig2Bkg1(Sig2Bkg1<0)=0;
    Sig2Bkg2 = SignalCount2./BkgCounts2; Sig2Bkg2(Sig2Bkg2<0)=0;
    
    qU1 = Model1.qU;
    qU2 = Model2.qU;
    
    MakeDir(savedir);
    save(savename,'Data1','Data2','StackData1','StackData2',...
                   'Model1','Model2',...
                   'SignalRate2','SignalRate2','SignalCount1','SignalCount2',...
                   'BkgRate1','BkgRate2','BkgCounts1','BkgCounts2',...
                   'TimeTotSubrun1','TimeTotSubrun2',...
                   'Sig2Bkg1','Sig2Bkg2','qU1','qU2');
else
    load(savename);
end
%% select range [E0-40,40] eV
StartIdx1 = 13;
StartIdx2 = 11;
StopIdx1  = 34;
StopIdx2  = 33;
%% absolute (sum) model
SigSum2 = sum(SignalCount2(StartIdx2:StopIdx2));
BkgSum2 = sum(BkgCounts2(StartIdx2:StopIdx2));
Sig2Bkgtot2 = SigSum2/BkgSum2;

SigSum1 = sum(SignalCount1(StartIdx1:StopIdx1));
BkgSum1 = sum(BkgCounts1(StartIdx1:StopIdx1));
Sig2Bkgtot1 = SigSum1/BkgSum1;

%% something weird with data a la christian
% CountsAll2      = sum(StackData2.TBDIS(StartIdx2-1:StopIdx2));
% MeanBkgRateData = mean(StackData2.TBDIS(StopIdx2+1:end)./TimeTotSubrun2(StopIdx2+1:end));
% BkgCounts2      = sum(0.22.*TimeTotSubrun2(StartIdx2:StopIdx2)); %MeanBkgRateData
% (CountsAll2-BkgCounts2)/BkgCounts2

%% 

%% plt
GetFigure;
plot(qU1-18574,ones(39,1),':','LineWidth',2,'Color',rgb('Black'));
hold on
plot(qU1-18574,Sig2Bkgtot1.*ones(39,1),':','LineWidth',2,'Color',rgb('Orange'));
plot(qU2-18574,Sig2Bkgtot2.*ones(38,1),':','LineWidth',2,'Color',rgb('ForestGreen'));

p1 = plot(qU1-18574,Sig2Bkg1,'.-','LineWidth',3,'Color',rgb('Orange'),'MarkerSize',18);

p2 = plot(qU2-18574,Sig2Bkg2,'.-','LineWidth',3,'Color',rgb('ForestGreen'),'MarkerSize',18);

xlabel('Retarding energy - 18574 (eV)');
ylabel('Signal to background');
PrettyFigureFormat('FontSize',22);
xlim([-42 5]);
ylim([1e-03 350])
set(gca,'YScale','log');
leg = legend([p1,p2],...
    'KNM1',...
    'KNM2',...
    'Location','northeast');
PrettyLegendFormat(leg);
t1 = text(-9,Sig2Bkgtot1+1.5,...
    sprintf('\\Sigma sig / \\Sigma bkg = %.1f',Sig2Bkgtot1),'Color',rgb('Orange'),...
    'FontSize',get(gca,'FontSize')-1,'FontName',get(gca,'FontName'));
t2 = text(-5,1.5,...
    sprintf('sig = bkg'),'Color',rgb('Black'),...
    'FontSize',get(gca,'FontSize')-1,'FontName',get(gca,'FontName'));
t0 = text(-9,Sig2Bkgtot2+4,...
    sprintf('\\Sigma sig / \\Sigma bkg = %.1f',Sig2Bkgtot2),'Color',rgb('ForestGreen'),...
    'FontSize',get(gca,'FontSize')-1,'FontName',get(gca,'FontName'));

%% interpolate for KSN2:
S2B = interp1(qU2-18574,Sig2Bkg2,[-40,-30,-20,-10,0],'spline');
fprintf('Signal to background ratio of 2 at %.1f eV \n',interp1(Sig2Bkg2(1:end-5),qU2(1:end-5)-18574,1,'spline'))

