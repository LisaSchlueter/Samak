M = MultiRunAnalysis('RunList','KNM1','DataType','Twin','fixPar','5 6 7 8 9 10 11',...
    'exclDataStart',14,'SysBudget',22);
S = RunSensitivity('RunAnaObj',M);
M.chi2 = 'chi2CMShape';
S.PlotSysBreakdownBars('Ranges',[2,10,11,12,14,17]);
% try
%     TestStackingqUCorrectionTwins;
% catch
% end
% try
%     M = MultiRunAnalysis('RunList','KNM1','DataType','Twin','fixPar','5 6 7 8 9 10 11','exclDataStart',14,'SysBudget',1);
%     S = RunSensitivity('RunAnaObj',M);
%     S.RecomputeFlag='ON';
%     S.PlotSysBreakdownBars('Ranges',[2,10,11,12,14,17]);
% catch
% end
% 
% try
%     M = MultiRunAnalysis('RunList','KNM1','DataType','Twin','fixPar','5 6 7 8 9 10 11','exclDataStart',14,'SysBudget',2);
%     S = RunSensitivity('RunAnaObj',M);
%     S.RecomputeFlag='ON';
%     S.PlotSysBreakdownBars('Ranges',[2,10,11,12,14,17]);
% catch
% end
% 
% try
%     ColumnDensityScan;
% catch
% end