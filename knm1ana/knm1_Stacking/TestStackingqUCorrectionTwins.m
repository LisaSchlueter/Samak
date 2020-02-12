CommongArg = {'RunList','KNM1','DataType','Twin','fixPar','5 6 7 8 9 10 11'};

if ~exist('Moff','var')
Moff = MultiRunAnalysis(CommongArg{:},'StackqUCorr','OFF');
StackqU_i = Moff.RunData.qU;
Moff.Fit;
end

if ~exist('Mon','var')
Mon = MultiRunAnalysis(CommongArg{:},'StackqUCorr','ON','Debug','ON');
%StackqU_corr = Mon.RunData.qU;
Mon.Fit;
end

%[SingleRunDataRate,SingleRunqU] = Mon.ComputeStackqU;
%% plot old and new <qU> for 1 subrun
qu = 30;
%SingleRunDataRate = Moff.SingleRunData.TBDIS(qu,:)./(Moff.SingleRunData.qUfrac(qu,:).*M.SingleRunData.TimeSec);
StackDataRate = Moff.RunData.TBDIS(qu)./(Moff.RunData.qUfrac(qu).*Moff.ModelObj.TimeSec);

psingle = plot(SingleRunDataRate(qu,:),SingleRunqU(qu,:)-18575,'x','MarkerSize',9,'Color',rgb('CadetBlue'));
hold on;
pStackoff = plot(StackDataRate,Moff.RunData.qU(qu)-18575,...
    'o','MarkerFaceColor',rgb('SlateGray'),'MarkerSize',9,'Color',rgb('SlateGray'));
pStackon = plot(StackDataRate,Mon.RunData.qU(qu)-18575,...
    'o','MarkerFaceColor',rgb('FireBrick'),'MarkerSize',9,'Color',rgb('FireBrick'));

PrettyFigureFormat;
ylabel('qU below E0 (eV)'); xlabel('count rate (cps)');
legend('single runs', 'stack wmean',sprintf('stack interpol')); legend boxoff
grid on;
set(gca,'FontSize',18); set(gca,'YScale','log'); title('ComputeStackqU');
hold off;