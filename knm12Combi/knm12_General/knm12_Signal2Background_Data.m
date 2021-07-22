% plot signal to background as a function of qU
savedir = sprintf('%sknm12Combi/knm12_General/results/',getenv('SamakPath'));
savename = sprintf('%sknm12_Data.mat',savedir);

if ~exist(savename,'file')
    M2 = MultiRunAnalysis('RunList','KNM2_Prompt','NonPoissonScaleFactor',1);
    M1 = MultiRunAnalysis('RunList','KNM1','NonPoissonScaleFactor',1);
    
    StackData1 = M1.RunData;
    Data1      = M1.SingleRunData;
    
    StackData2 = M2.RunData;
    Data2      = M2.SingleRunData;
    
    MakeDir(savedir);
    save(savename,'Data1','Data2','StackData1','StackData2');
else
    load(savename);
end
%% select range [E0-40,40] eV
StartIdx1 = 13;
StartIdx2 = 11;
StopIdx1  = 34;
StopIdx2  = 33;

Time1 = StackData1.TimeSec.*StackData1.qUfrac;
Time2 = StackData2.TimeSec.*StackData2.qUfrac;

CountsSum2      = sum(StackData2.TBDIS(StartIdx2:StopIdx2)); % [only E0-40,E0]
CountsSum1      = sum(StackData1.TBDIS(StartIdx1:StopIdx1)); % [only E0-40,E0]
Counts2      = StackData2.TBDIS; % entire range
Counts1      = StackData1.TBDIS; % entire range
  
CountsSum2all   = sum(StackData2.TBDIS(StartIdx2:end));  % all 28 subruns
CountsSum1all   = sum(StackData1.TBDIS(StartIdx1:end)); % all 27 subruns

%Estimate Background: use mean and scale to all sub-runs
BkgCounts2_OnlyBkg      = sum(StackData2.TBDIS(StopIdx2+1:end));
MeanBkgRate2            = mean(StackData2.TBDIS(StopIdx2+1:end)./Time2(StopIdx2+1:end));
Bkg2 = MeanBkgRate2.*Time2;
BkgSum2 = sum(MeanBkgRate2.*Time2(StartIdx2:StopIdx2));

BkgCounts1_OnlyBkg      = sum(StackData1.TBDIS(StopIdx1+1:end));
MeanBkgRate1            = mean(StackData1.TBDIS(StopIdx1+1:end)./Time1(StopIdx1+1:end));
Bkg1 = MeanBkgRate1.*Time1;
BkgSum1 = sum(MeanBkgRate1.*Time1(StartIdx1:StopIdx1));

% subtract background to get signal- only
SigSum2 = CountsSum2-BkgSum2;
SigSum1 = CountsSum1-BkgSum1;
Sig1 = Counts1-Bkg1;
Sig2 = Counts2-Bkg2;

Sig2BkgSum2 = SigSum2/BkgSum2;
Sig2BkgSum1 = SigSum1/BkgSum1;
Sig2Bkg2 = Sig2./Bkg2; Sig2Bkg2(Sig2Bkg2<0)=0;
Sig2Bkg1 = Sig1./Bkg1;Sig2Bkg1(Sig2Bkg1<0)=0;


fprintf('KNM1 data in range [E0-40eV , E0] \n');
fprintf('Counts [E0-40eV , E0+47]    %.2f millon \n',CountsSum1all./1e6);
fprintf('Counts [E0-40eV , E0]       %.2f millon \n',CountsSum1./1e6);
fprintf('Counts signal:              %.2f millon \n',SigSum1./1e6);
fprintf('Counts bg (only bg points): %.2f millon \n',sum(StackData1.TBDIS(StopIdx1+1:end))./1e6);
fprintf('Counts background:          %.2f millon \n',BkgSum1./1e6);
fprintf('S/B:                        %.1f  \n',Sig2BkgSum1);

fprintf('KNM2 data in range [E0-40eV , E0] \n');
fprintf('Counts [E0-40eV , E0+135]   %.2f millon \n',CountsSum2all./1e6);
fprintf('Counts [E0-40eV , E0]       %.2f millon \n',CountsSum2./1e6);
fprintf('Counts signal:              %.2f millon \n',SigSum2./1e6);
fprintf('Counts bg (only bg points): %.2f millon \n',sum(StackData2.TBDIS(StopIdx2+1:end))./1e6);
fprintf('Counts background:          %.2f millon \n',BkgSum2./1e6);
fprintf('S/B:                        %.1f  \n',Sig2BkgSum2);

%% plt
GetFigure;
%plot(qU1-18574,Sig2Bkgtot2.*ones(39,1),':','LineWidth',2,'Color',rgb('ForestGreen'));

p1 = plot(StackData1.qU-18574,Sig2Bkg1,'-','LineWidth',3);
hold on;
p2 = plot(StackData2.qU-18574,Sig2Bkg2,'-.','LineWidth',3);

xlabel('Retarding energy - 18574 (eV)');
ylabel('Signal to background');
PrettyFigureFormat('FontSize',22);
xlim([-42 5]);
ylim([1e-03 300])
set(gca,'YScale','log');
leg = legend([p1,p2],...
    'KNM1',...
    'KNM2',...
    'Location','northeast');
PrettyLegendFormat(leg);




