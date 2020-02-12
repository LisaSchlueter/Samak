addpath(genpath('../../../Samak2.0'));
resultsCM = load('FitResultsCM_hvdata_all_17-Mar-2018_newB-FieldCorr40Pixels.mat');
results = load('FitResults_hvdata_all_16-Mar-2018_newB-FieldCorr40Pixels.mat');

FitE_CM = resultsCM.LineE_i + resultsCM.par(1);
FitEerr_CM = resultsCM.err(1);
FitW_CM = resultsCM.LineW_i + resultsCM.par(2);
FitWerr_CM = resultsCM.err(2);
FitPhi0_CM = resultsCM.LinePhi0_i(1:40) + resultsCM.par(3:42);
FitPhi0err_CM= resultsCM.err(3:42);
FitOffset_CM = resultsCM.LineBKG_i(1:40) + resultsCM.par(43:82);
FitOffseterr_CM = resultsCM.err(43:82);
FitOffsetWMean_CM = wmean(FitOffset_CM,FitOffseterr_CM);
FitPhi0WMean_CM = wmean(FitPhi0_CM,FitPhi0err_CM);


FitE = results.LineE_i + results.par(1);
FitEerr = results.err(1);
FitW = results.LineW_i + results.par(2);
FitWerr = results.err(2);
FitPhi0 = results.LinePhi0_i(1:40) + results.par(3:42);
FitPhi0err= results.err(3:42);
FitOffset = results.LineBKG_i(1:40) + results.par(43:82);
FitOffseterr = results.err(43:82);
FitOffsetWMean = wmean(FitOffset,FitOffseterr);
FitPhi0WMean = wmean(FitPhi0,FitPhi0err);