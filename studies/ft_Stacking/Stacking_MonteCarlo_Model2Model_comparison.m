%------------------------------------------------------------------------%
% Little Study to compare Modelling of Stacking to Modelling of Single Runs
% for first Tritium
% Fit to Normalization and Background, Rest fixed
% L.Schlueter June/2018
%------------------------------------------------------------------------%
close all;
%Stacked Model
RunList = [40538:40543,40603,40604,40610:40613,40667:40693];%[40538:40543];%, 40603,40604,40610:40613,40667:40693];
A = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','ringCutFlag','ex2b','fixPar','1 2');

%Single Runs
SingleRunObj = cell(length(A.StackedRuns),4);
TBDIS_Stack = struct('AllAverage', zeros(A.ModelObj.nqU,1),...
    'SingleRun_qU',zeros(A.ModelObj.nqU,1),...
    'SingleRun_Rho_DTTTHT',zeros(A.ModelObj.nqU,1),...
    'AllSingleRun',zeros(A.ModelObj.nqU,1));

for r=1:length(A.StackedRuns)
 % Single Runs with Average qU, column density, DT, HT, TT - fractions
 SingleRunObj{r,1} = ref_RunSummaries_StackPix(...
     A.StackedRuns(r),A.ringCutFlag,...
     'TD',['Run',A.StackFileName,A.ringCutFlag],...
     'TimeSec',(A.ModelObj.TimeSec/length(A.StackedRuns)),...
     'ISCS','Theory','recomputeRF','OFF',...
     'WGTS_CD_MolPerCm2', A.ModelObj.WGTS_CD_MolPerCm2,...
     'WGTS_MolFrac_TT', A.ModelObj.WGTS_MolFrac_TT,...
     'WGTS_MolFrac_HT', A.ModelObj.WGTS_MolFrac_HT,...
     'WGTS_MolFrac_DT', A.ModelObj.WGTS_MolFrac_DT);
    SingleRunObj{r,1}.ComputeTBDDS;
    SingleRunObj{r,1}.ComputeTBDIS;
    TBDIS_Stack.AllAverage = TBDIS_Stack.AllAverage + SingleRunObj{r,1}.TBDIS;
    
 % Single Runs with average column density, DT, HT, TT - fractions
 % and individual qU
    SingleRunObj{r,2} = ref_RunSummaries_StackPix(A.StackedRuns(r),A.ringCutFlag,...
        'ISCS','Theory','recomputeRF','OFF',...
        'WGTS_CD_MolPerCm2', A.ModelObj.WGTS_CD_MolPerCm2,...
        'WGTS_MolFrac_TT', A.ModelObj.WGTS_MolFrac_TT,...
        'WGTS_MolFrac_HT', A.ModelObj.WGTS_MolFrac_HT,...
        'WGTS_MolFrac_DT', A.ModelObj.WGTS_MolFrac_DT);
    
    SingleRunObj{r,2}.ComputeTBDDS;
    SingleRunObj{r,2}.ComputeTBDIS;
    TBDIS_Stack.SingleRun_qU = TBDIS_Stack.SingleRun_qU + SingleRunObj{r,2}.TBDIS;
    
 % Single Runs with average qU
 % and individual DT,HT,TT, column density
    SingleRunObj{r,3} = ref_RunSummaries_StackPix(A.StackedRuns(r),A.ringCutFlag,...
        'ISCS','Theory','recomputeRF','OFF',...
        'TD',['Run',A.StackFileName,A.ringCutFlag] ,...
        'TimeSec',A.ModelObj.TimeSec/length(A.StackedRuns));   
    SingleRunObj{r,3}.ComputeTBDDS;
    SingleRunObj{r,3}.ComputeTBDIS;
    TBDIS_Stack.SingleRun_Rho_DTTTHT = TBDIS_Stack.SingleRun_Rho_DTTTHT + SingleRunObj{r,3}.TBDIS;
    
  % Single Runs everything individual 
    SingleRunObj{r,4} = ref_RunSummaries_StackPix(A.StackedRuns(r),A.ringCutFlag,...
        'ISCS','Theory','recomputeRF','OFF');  
    SingleRunObj{r,4}.ComputeTBDDS;
    SingleRunObj{r,4}.ComputeTBDIS;
    TBDIS_Stack.AllSingleRun = TBDIS_Stack.AllSingleRun + SingleRunObj{r,4}.TBDIS;
       
end

%%
mycolor ={rgb('CadetBlue'); rgb('Orange'); rgb('DarkRed'); rgb('Green')};
fieldnames = fields(TBDIS_Stack);

for i=1:numel(fields(TBDIS_Stack))
 A.RunData.qU = A.ModelObj.qU;
 A.RunData.TBDIS = TBDIS_Stack.(fieldnames{i});
 A.fixPar = '1 2';
 A.Fit;
%A.InitModelObj_Norm_BKG;
%Plot
f8 = figure(i);
set(f8, 'Units', 'normalized', 'Position', [0.9, 0.9, 1, 0.9]);

subplot(2,1,1)
plot(A.ModelObj.qU-A.ModelObj.Q_i,TBDIS_Stack.(fieldnames{i}),'--',...
     A.ModelObj.qU-A.ModelObj.Q_i,A.ModelObj.TBDIS,'-','Color',mycolor{i});
 set(gca,'FontSize',18);
grid on;
legend('\Sigma Single Run Models', 'Stack Model');
xlabel('retarding potential - 18575 (eV)');
ylabel('counts');
xlim([A.RunData.qU(A.exclDataStart)-A.ModelObj.Q_i, A.RunData.qU(end)-A.ModelObj.Q_i]);
mytitle = sprintf('\n%s',A.GetRunTitle);
title(['\Sigma Single Run Simulation vs Stacked Run Simulation',mytitle])

subplot(2,1,2);
plot(A.ModelObj.qU-A.ModelObj.Q_i,...
    (TBDIS_Stack.(fieldnames{i})-A.ModelObj.TBDIS)./sqrt(TBDIS_Stack.(fieldnames{i})),...
    's-','Color',mycolor{i},'LineWidth',3,'MarkerSize',6);
hold on;   
plot(A.ModelObj.qU-A.ModelObj.Q_i, zeros(A.ModelObj.nqU,1),'k--');
%PrettyFigureFormat;
xlabel('retarding potential - 18575 (eV)');
ylabel('norm. residuals');
xlim([A.RunData.qU(A.exclDataStart)-A.ModelObj.Q_i, A.RunData.qU(end)-A.ModelObj.Q_i]);
chi2leg = '.';%sprintf('%.4f (N, B fit)\n',A.FitResult.chi2min);
legend(['\chi^2=',chi2leg,strrep(fieldnames{i},'_','-')],'Location','northwest');
set(gca,'FontSize',18);
grid on;
hold off;
publish_figure(f8,['./MonteCarlo-Test/',fieldnames{i},'.eps']);
end


% Display
 %A.PlotStackingModel;
 fprintf(2,'------------------------------------------------------------------------\n')
 fprintf(2,'TD StackedRuns =  %s\n', A.ModelObj.TD)
 fprintf(2,'TD SingleRuns  =  %s\n' , SingleRunObj{1,1}.TD)
 fprintf(2,' Difference WGTS_CD_molcm2 = %g \n', A.ModelObj.WGTS_CD_MolPerCm2-SingleRunObj{1,1}.WGTS_CD_MolPerCm2)
 fprintf(2,' Difference DT-fraction = %.5f \n', A.ModelObj.WGTS_MolFrac_DT-SingleRunObj{1,1}.WGTS_MolFrac_DT)
% fprintf(2,' Difference Time = %0.5f \n', A.ModelObj.TimeSec-total_time)
 fprintf(2,'------------------------------------------------------------------------\n')