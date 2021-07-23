DataSet    = 'KNM2_Prompt';
Recompute  = 'ON';
CorrectErr = 'OFF';
fitter     = 'minuit';

if strcmp(DataSet,'KNM1')
    SB=24;
elseif strcmp(DataSet,'KNM2_Prompt')
    SB=40;
elseif strcmp(DataSet,'TDR')
    SB=67;
end
A=RelicNuAnalysis('Params',DataSet);

%UpperLimitvsmnuSq('Params',DataSet,'Recompute',Recompute);                                                              %Sensitivity limit of twins

%A.EtaFit('Mode','SinglemnuSq','SysBudget',SB,'fitPar','E0 Norm Bkg eta','CorrectErr',CorrectErr,'Recompute',Recompute);  %Raster scan twins

%A.Chi2Scan_2D('RunList',DataSet,'SystBudget',SB,'Recompute',Recompute);                                                 %2D scan twins

%A.SystBreakdown('TwinBias_mnuSq',1,'CorrectErr',CorrectErr,'Recompute',Recompute);                                       %Systematics breakdown twins

%PlotBestFit('RunList',DataSet,'DataType','Twin','Nfit',1000);                                                           %p-value

%ProfileLikelihood('RunList',DataSet,'SysBudget',SB);

if ~strcmp(DataSet,'TDR')
    %PlotBestFit('RunList',DataSet,'DataType','Real','Plot','ON','saveplot','ON');                                                                           %Best fit data

    %A.EtaFit('Mode','SinglemnuSq','SysBudget',SB,'DataType','Real','fitPar','E0 Norm Bkg eta','DeltaChi2',1,'CorrectErr',CorrectErr,'Recompute',Recompute); %Raster Scan data

    %A.EtaFit('Mode','SinglemnuSq','SysBudget',SB,'DataType','Real','fitPar','E0 Norm Bkg eta','DeltaChi2',4,'CorrectErr',CorrectErr,'Recompute',Recompute);

    %A.EtaFit('Mode','SinglemnuSq','SysBudget',SB,'DataType','Real','fitPar','E0 Norm Bkg eta','DeltaChi2',9,'CorrectErr',CorrectErr,'Recompute',Recompute);

    A.Chi2Scan_2D('RunList',DataSet,'SystBudget',SB,'DataType','Real','Nmnubins',50,'Netabins',51,'Recompute',Recompute,'fitter',fitter,'pullFlag',31);     %2D scan data
    B=RelicNuAnalysis('Params','KNM1'); B.Chi2Scan_2D('RunList','KNM1','SystBudget',24,'DataType','Real','Nmnubins',50,'Netabins',51,'Recompute',Recompute,'fitter',fitter,'pullFlag',31);

    A.SystBreakdown('DataType','Real','CorrectErr',CorrectErr,'Recompute',Recompute);                                                         %Systematics breakdown data
end