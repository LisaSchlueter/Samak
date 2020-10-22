
RunDataR = GetRunData('Real');
RunDataT =  GetRunData('Twin');

range = 40;
Start =  find((RunDataR.qU)>=18573.7-range,1);
TimeSubRun = RunDataR.qUfrac(Start:end)*RunDataR.TimeSec;

f4 = figure('Units','normalized','Position',[0.1,0.1,0.45,0.5]);
plot(RunDataR.qU(Start:end)-18574,RunDataR.TBDIS(Start:end)./TimeSubRun,'-','LineWidth',3,'MarkerSize',25,'Color',rgb('DodgerBlue'));
hold on;
plot(RunDataT.qU(Start:end)-18574,RunDataT.TBDIS(Start:end)./TimeSubRun,'--','LineWidth',3,'MarkerSize',25,'Color',rgb('GoldenRod'));
xlim([-40 5])
set(gca,'YScale','log');
PrettyFigureFormat('FontSize',22);
ylabel('Count rate (cps)');
xlabel('Retarding potential - 18574 (eV)');
leg = legend(sprintf('\\Sigma data = %.3e counts',sum(RunDataR.TBDIS(Start:end))),...
             sprintf('\\Sigma twin^{ } = %.3e counts',sum(RunDataT.TBDIS(Start:end))),'EdgeColor',rgb('Silver'));
hold off
savename = [getenv('SamakPath'),'knm2ana/knm2_unblindingFinal/plots/knm2ub2_OverlayTwinDataTBDIS.png'];
print(savename,'-dpng','-r350');
%% function
function RunData = GetRunData(DataType)
savedir = [getenv('SamakPath'),'knm2ana/knm2_unblindingFinal/results/'];
savename = sprintf('%sknm2ub1_OverlayTwinDataTBDIS_%s.mat',...
    savedir,DataType);

if exist(savename,'file')
    RunData = importdata(savename);
    
else
    SigmaSq =  0.0124+0.0025;
    
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'TwinBias_Q',18573.7,...
        'FSDFlag','Sibille0p5eV'};
    A = MultiRunAnalysis(RunAnaArg{:});
    RunData = A.RunData;
    save(savename,'RunData');
end
end