% response function combi plot
savedir = [getenv('SamakPath'),'knm2ana/knm2_ResponseFunction/results/'];
savename = sprintf('%sResponseFunction.mat',savedir);

if exist(savename,'file')
    load(savename);
else
    DataType = 'Real';
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag','KNM2',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',.....
        'AngularTFFlag','ON',...
        'SynchrotronFlag','ON'....
        'RadiativeFlag','ON'};
    
    A2 = MultiRunAnalysis(RunAnaArg{:});
    Te_2 = A2.ModelObj.Te;
    qU_2 = A2.ModelObj.qU;
    RF_2 = A2.ModelObj.RF;
    DataRate_2 = A2.RunData.TBDIS./(A2.RunData.qUfrac.*A2.RunData.TimeSec);
    RunData2 = A2.RunData;
    A2.InitModelObj_Norm_BKG;
   ModelRate_2 = A2.ModelObj.TBDIS./(A2.RunData.qUfrac.*A2.RunData.TimeSec);
   
    A1 = MultiRunAnalysis('RunList','KNM1',...
        'DataType',DataType,...
        'minuitOpt','min ; minos',...
        'FSDFlag','Sibille0p5eV',...
        'ELossFlag','KatrinT2A20',...
        'AngularTFFlag','ON',...
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',1);
    Te_1 = A1.ModelObj.Te;
    qU_1 = A1.ModelObj.qU;
    RF_1 = A1.ModelObj.RF;
    DataRate_1 = A1.RunData.TBDIS./(A1.RunData.qUfrac.*A1.RunData.TimeSec);
    RunData1 = A1.RunData;
    A1.InitModelObj_Norm_BKG;
    ModelRate_1 = A1.ModelObj.TBDIS./(A1.RunData.qUfrac.*A1.RunData.TimeSec);
    
    save(savename,'Te_1','Te_2','qU_1','qU_2','RF_1','RF_2',...
        'DataRate_1','DataRate_2','RunData1','RunData2',...
        'ModelRate_1','ModelRate_2');
end

%% plot response function

GetFigure;
p1 = plot(Te_1-qU_1(30),RF_1(:,30),'-','LineWidth',2.5,'Color',rgb('DodgerBlue'));
hold on;
p2 = plot(Te_2-qU_2(30),RF_2(:,30),'-.','LineWidth',2.5,'Color',rgb('Orange'));
PrettyFigureFormat('FontSize',28);
xlabel(sprintf('Energy {\\itE} - {\\it qU} (eV)'));
ylabel('Transmission probability');
leg = legend([p1,p2],'KNM-1','KNM-2','Location','northwest');
PrettyLegendFormat(leg);
xlim([-5 40])

pltdir = strrep(savedir,'results','plots');
pltname = sprintf('%sResponseFunction.png',pltdir);
print(gcf,pltname,'-dpng','-r350');

%%
%% plot response function

GetFigure;
p1 = plot(qU_1-18574,ModelRate_1,'-','LineWidth',2.5,'Color',rgb('DodgerBlue'));
hold on;
p2 = plot(qU_2-18574,ModelRate_2,'-.','LineWidth',2.5,'Color',rgb('Orange'));
PrettyFigureFormat('FontSize',28);
xlabel(sprintf('Retarding energy {\\it qU} - {\\it E}_0 (eV)'));
ylabel('Rate (cps)');
leg = legend([p1,p2],'KNM-1','KNM-2','Location','northeast');
PrettyLegendFormat(leg);
xlim([-42 10])
set(gca,'YScale','log');
ylim([0.1 200])
pltdir = strrep(savedir,'results','plots');
pltname = sprintf('%sIntSpec.png',pltdir);
print(gcf,pltname,'-dpng','-r350');

