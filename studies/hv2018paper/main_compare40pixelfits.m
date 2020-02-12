load('./results/FitResults_hvdata_all40Pixels_4March_nonvectorized.mat')
fprintf('non-vectorized: W = %.3f\n', (LineW_i+par(2)))
FitPhi0 = LinePhi0_i(1:40) + par(3:42);
FitOffset = LineBKG_i(1:40) + par(43:82);

load('./results/FitResultsCM_hvdata_all40Pixels_vectorizedv2.mat')
fprintf('vectorized: W = %.3f\n', (LineW_i+par(2)))
FitPhi0_vec = LinePhi0_i(1:40) + par(3:42);
FitOffset_vec = LineBKG_i(1:40) + par(43:82);

figure(5)
subplot(2,1,1);
Phi0Diff = FitPhi0-FitPhi0_vec;
nhist(Phi0Diff);
xlabel('Phi0 Diff');
title('Non Vec - Vec');
subplot(2,1,2);
OffsetDiff = FitOffset-FitOffset_vec;
nhist(OffsetDiff);
xlabel('Offset Diff');

%FPDViewer(Phi0Diff);