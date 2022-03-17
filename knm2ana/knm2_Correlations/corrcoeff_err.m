function [MeanCorrMat,StdCorrMat,CorrMat_i] = corrcoeff_err(FitPar)
CorrMat_i = corrcoef(FitPar'); % correlation matrix

nReSample = 1e4;
RandIdx = randi(size(FitPar,2),[size(FitPar,2),nReSample]);

CorrMat_sample = zeros(4,4,nReSample);
for i=1:nReSample
    CorrMat_sample(:,:,i) = corrcoef(FitPar(:,RandIdx(:,i))');
end
MeanCorrMat = mean(CorrMat_sample,3);
StdCorrMat = std(CorrMat_sample,0,3);


fprintf('Correlation(m2, E0): %.2f +- %.3f \n',MeanCorrMat(1,2),StdCorrMat(1,2));
fprintf('Correlation(m2, Bkg): %.2f +- %.3f \n',MeanCorrMat(1,3),StdCorrMat(1,3));
fprintf('Correlation(m2, Norm): %.2f +- %.3f \n',MeanCorrMat(1,4),StdCorrMat(1,4));
fprintf('Correlation(E0, Bkg): %.2f +- %.3f \n',MeanCorrMat(2,3),StdCorrMat(2,3));
fprintf('Correlation(E0, Norm): %.2f +- %.3f \n',MeanCorrMat(2,4),StdCorrMat(2,4));
fprintf('Correlation(Bkg, Norm): %.2f +- %.3f \n',MeanCorrMat(3,4),StdCorrMat(3,4));
end
