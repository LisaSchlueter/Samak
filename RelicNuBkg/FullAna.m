DataSet='KNM1';

if strcmp(DataSet,'KNM1')
    SB=24;
elseif strcmp(DataSet,'KNM2')
    SB=40;
elseif strcmp(DataSet,'TDR')
    SB=67;
end
A=RelicNuAnalysis('Params',DataSet);

UpperLimitvsmnuSq('Params',DataSet);                                                                        %Sensitivity limit of twins

A.EtaFit('Mode','SinglemnuSq','SysBudget',SB,'fitPar','E0 Norm Bkg eta');                                   %Raster scan twins

A.Chi2Scan_2D('RunList',DataSet,'SystBudget',SB);                                                           %2D scan twins

A.SystBreakdown('TwinBias_mnuSq',1);                                                                        %Systematics breakdown twins

PlotBestFit('RunList',DataSet,'DataType','Twin','Nfit',1000);                                               %p-value

if ~strcmp(DataSet,'TDR')
    PlotBestFit('RunList',DataSet,'DataType','Real','saveplot','ON');                                           %Best fit data

    A.EtaFit('Mode','SinglemnuSq','SysBudget',SB,'DataType','Real','fitPar','E0 Norm Bkg eta','DeltaChi2',1);   %Raster Scan data

    A.EtaFit('Mode','SinglemnuSq','SysBudget',SB,'DataType','Real','fitPar','E0 Norm Bkg eta','DeltaChi2',4);

    A.EtaFit('Mode','SinglemnuSq','SysBudget',SB,'DataType','Real','fitPar','E0 Norm Bkg eta','DeltaChi2',9);

    A.Chi2Scan_2D('RunList',DataSet,'SystBudget',SB,'DataType','Real','Nmnubins',30,'Netabins',30);                                         %2D scan data

    A.SystBreakdown('DataType','Real');                                                                         %Systematics breakdown data
end