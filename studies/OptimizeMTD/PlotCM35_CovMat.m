%% RF
rf=load('WGTSMACE_CovMat_MTDcreator60_CD5e+17molPercm2-0.005err_xsection0.005err_Bs-2.52T_Bmax-4.20T_BF0.005err_Ba7G-0.0054err_1000Trials.mat');
rf.obj.StudyObject.TimeSec = 86400*365;
rf.obj.PlotCM('ConvergenceTest','OFF');
publish_figurePDF(1,'tmpRF.pdf');

%% FSD
fsd=load('FSD_DT-HT-TT-CovMat_1000Trials_MTDcreator60_0.01NormErr_0.04GS_0.18ES_ShapeErr.mat');
fsd.obj.StudyObject.TimeSec = 86400*365;
fsd.obj.PlotCM('ConvergenceTest','OFF');
publish_figurePDF(1,'tmpFSD.pdf');

%% TC
tc=load('TCoff_CovMat_MTDcreator60.mat');
tc.obj.StudyObject.TimeSec = 86400*365;
tc.obj.PlotCM('ConvergenceTest','OFF');
publish_figurePDF(1,'tmpTC.pdf');

% Combi

%% TASR
ta=load('TCoff_CovMat_RunStack_538_540_541_542_543_603_604_611_613_667_668_669_670_671_672_673_674_675_676_677_678_679_680_681_682_683_684_685_686_687_688_689_690_691_692_693ex2b_1-Runs.mat');
%ta.obj.StudyObject.TimeSec = 86400;
ta.obj.PlotCM('ConvergenceTest','OFF');
publish_figurePDF(1,'tmpTASR.pdf');

%% Combi
c=load('CombiCM_MTDcreator60_1000-Trials_RF_0.002CDErr-FSD_0.01NormErr0.04ShapeGSErr0.18ShapeESErr-TCoff_RAD-TCoff_OTHER.mat');
c.obj.StudyObject.TimeSec = 86400*365;
c.obj.PlotCM('ConvergenceTest','OFF');


