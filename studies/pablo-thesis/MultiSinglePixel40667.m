%addpath(genpath('../../../Samak2.0'));

fitter = 'minuit';

chi2name = 'chi2P';

dataStart = 9; 

R40667Model = ref_RunAnalysis('40677','','mpix',...
    'FPD_Segmentation','SINGLEPIXEL','FPD_Pixel',1:124,...
    'nTeBinningFactor',5,'recomputeRF','OFF','useParallelRF','OFF');
R40667Model.ComputeTBDDS();
R40667Model.ComputeTBDIS();

R40667Data = load('40667mpix.mat');

Data = [R40667Data.qU,R40667Data.TBDIS,sqrt(R40667Data.TBDIS)];

% Single Pixel

for p = 1:1

F = FITC('SO',R40667Model,'DATA',Data,'fitter',fitter,...
    'chi2name',chi2name,...
    'fixPar','5 6',...
    'exclDataStart',dataStart);
end