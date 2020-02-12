M = MultiRunAnalysis('RunList','KNM1','DataType','Real');
WGTS_CD_MolPerCm2_i = M.ModelObj.WGTS_CD_MolPerCm2;

WGTS_CD_MolPerCm2_Rel = [0.95:0.01:1.05];
Rate200     = zeros(numel(WGTS_CD_MolPerCm2_Rel),1);
RateErr200 = zeros(numel(WGTS_CD_MolPerCm2_Rel),1);

for i=1:numel(WGTS_CD_MolPerCm2_Rel)
M.SimulateStackRuns('WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2_Rel(i)*WGTS_CD_MolPerCm2_i);
Rate200(i) = M.ModelObj.TBDIS(1)./(M.ModelObj.qUfrac(1).*M.ModelObj.TimeSec);
RateErr200(i) = M.ModelObj.TBDISE(1)./(M.ModelObj.qUfrac(1).*M.ModelObj.TimeSec);
end

Rate200Data    =  M.SingleRunData.TBDIS(1,:)./(M.SingleRunData.qUfrac(1,:).*M.SingleRunData.TimeSec);
RateErr200Data = M.SingleRunData.TBDISE(1,:)./(M.SingleRunData.qUfrac(1,:).*M.SingleRunData.TimeSec);
DataTpurity =  M.SingleRunData.WGTS_MolFrac_TT + 0.5*M.SingleRunData.WGTS_MolFrac_DT + ...
                0.5*M.SingleRunData.WGTS_MolFrac_HT;

%%
%M.LoadSingleRunObj;
%%
SingleObjActivity = cell2mat(cellfun(@(x) x.WGTS_epsT,M.SingleRunObj,'UniformOutput',0)).*...
                         M.SingleRunData.WGTS_CD_MolPerCm2';
SingleObjRate = cell2mat(cellfun(@(x) x.TBDIS(1)./(x.qUfrac(1)*x.TimeSec),M.SingleRunObj,...
    'UniformOutput',0));
%% fit
[a,b,c,d] = linFit(((M.ModelObj.WGTS_epsT.*WGTS_CD_MolPerCm2_Rel.*WGTS_CD_MolPerCm2_i)./mean((M.ModelObj.WGTS_epsT.*WGTS_CD_MolPerCm2_Rel.*WGTS_CD_MolPerCm2_i)))',Rate200./mean(Rate200),...
    RateErr200./mean(Rate200))
    
%
plot(M.ModelObj.WGTS_epsT.*WGTS_CD_MolPerCm2_Rel.*WGTS_CD_MolPerCm2_i,Rate200./mean(Rate200),...
    '--x','LineWidth',2);
hold on;
errorbar(DataTpurity.*M.SingleRunData.WGTS_CD_MolPerCm2,...
    Rate200Data./mean(Rate200Data),RateErr200Data./mean(Rate200Data),'o');
%plot(SingleObjActivity,SingleObjRate./mean(SingleObjRate),'s');
  hold off;  
PrettyFigureFormat;
xlabel('column density * tritium purity');
ylabel('-200eV FPD rate')