KSN1 = MultiRunAnalysis('RunList','KNM1',...
    'chi2','chi2Stat',...
    'DataType','Real',...
    'fixPar','',...
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1.064,...
    'minuitOpt','min ; migrad',...
    'FSDFlag','SibilleFull',...
    'ELossFlag','KatrinT2',...
    'SysBudget',22,...
    'exclDataStart',1); % 7 = --> 65 eV  /  13 = --> 40 eV / 

%KSN1.Fit;
%KSN1.PlotDataModel_KSN1;

%% Statistics
KSN1.PrintDataStatistics

%% Display Counting Figures
timeSec          = sum(KSN1.RunData.qUfrac(KSN1.exclDataStart:end)*KSN1.RunData.TimeSec); 
timeBelowE0Sec   = sum(KSN1.RunData.qUfrac(KSN1.exclDataStart:end-5).*KSN1.RunData.TimeSec); 
splusb           = sum(KSN1.RunData.TBDIS(KSN1.exclDataStart:end)); 
brate            = 0.293; %cps
b                = timeSec * brate;
s                = splusb - b;
bbelowE0         = timeBelowE0Sec * 0.293;

fprintf(2,'\n');

fprintf(2,'KSN1 - ExclDataStart=%.0f - subscan=%0.f -  qU=%.1f V - qU-E0=%.1f eV\n',...
    KSN1.exclDataStart,numel(KSN1.exclDataStart:KSN1.exclDataStop), KSN1.RunData.qU(KSN1.exclDataStart),KSN1.RunData.qU(KSN1.exclDataStart)-18573.73);
fprintf(2,'KSN1 - subscan below E0 = %0.f -  subscan above E0 = %0.f\n',...
    numel(KSN1.exclDataStart:KSN1.exclDataStop-5),numel(KSN1.exclDataStop-4:KSN1.exclDataStop));
fprintf(2,'KSN1 - Total Measurement Time: %.0f seconds - %.1f hours - %.2f days \n',timeSec,timeSec/3600,timeSec/86400);
fprintf(2,'KSN1 - Measurement Time Below E0: %.0f seconds - %.1f hours - %.2f days \n',timeBelowE0Sec,timeBelowE0Sec/3600,timeBelowE0Sec/86400);
fprintf(2,'KSN1 - Tritium Beta Decay + Background Electrons: %.3g counts \n',splusb);
fprintf(2,'KSN1 - Tritium Beta Decay Electrons: %.3g counts \n',s);
fprintf(2,'KSN1 - Total Background Electrons: %.3g counts \n',b);
fprintf(2,'KSN1 - Total Background Below E0: %.3g counts \n',bbelowE0);
fprintf(2,'KSN1 - Max S/B: %.4g  \n',(KSN1.RunData.TBDIS(KSN1.exclDataStart)./KSN1.RunData.qUfrac(KSN1.exclDataStart)./KSN1.RunData.TimeSec - brate)./brate);
fprintf(2,'\n');

%% Signal / Background
hv               = KSN1.RunData.qU(KSN1.exclDataStart:end-5)-18573.73;
sbratio          = ((KSN1.RunData.TBDIS(KSN1.exclDataStart:end-5)./KSN1.RunData.qUfrac(KSN1.exclDataStart:end-5)./KSN1.RunData.TimeSec) - brate) ./ brate;
fig = figure('Renderer','painters');
set(fig, 'Units', 'normalized', 'Position', [0.001, 0.001,0.8, 0.45]);
semilogy(hv,((sbratio)),'MarkerSize',2,'MarkerFaceColor',rgb('Black'),'MarkerEdgeColor',rgb('Black'),'LineWidth',3);
xlabel('Retarding energy - 18574 (eV)');
ylabel('log(Tritium / Background e^- Rate)');
ylim([0.1 2000])
grid on
PRLFormat     

%% Subscan Average Times
KSN1.SingleRunData.qUfrac(1).*(KSN1.SingleRunData.TimeSec);