format long; 
addpath(genpath('../../../Samak2.0'));
A= InitKrKATRIN_krl332();
%MyCovMat = cell(30,1)
tic
progressbar('Building CovMat')
for i= 1:15
  progressbar(i/15)
  A.FPD_Pixel = i;
  MyCovMat{i,1} = BuildCovMat_HVRipple(A, 1e4);
  mystr = sprintf('./CovMat_uniform/CovMatrix_KrL3LineFit_Pixel%u.mat', A.FPD_Pixel);
  save(mystr,'MyCovMat','-mat'); 
end
toc
