R = RunAnalysis('RunNr',51043);
exclDataStart = 2:30;
E0 = zeros(numel(exclDataStart),1);
qUStart = zeros(1,numel(exclDataStart));
errorStart = zeros(numel(exclDataStart),1);

for i = 1:numel(exclDataStart);
    R.exclDataStart = exclDataStart(i);
    R.Fit;
    E0(i) = R.FitResult.par(2) + 18575;
    errorStart(i) = R.FitResult.err(2);
    qUStart(i) = R.RunData.qU(i);
end;
figure(1);
errorbar(qUStart,E0,errorStart)
hold on
M = MultiRunAnalysis('RunList', 'KNM1_m149mvRW');
MultiexclDataStart = 2:30;
MultiE0 = zeros(numel(MultiexclDataStart),1);
MultiqUStart = zeros(1,numel(MultiexclDataStart));
MultierrorStart = zeros(numel(MultiexclDataStart),1);
for i = 1:numel(MultiexclDataStart);
    M.exclDataStart = MultiexclDataStart(i);
    M.Fit;
    MultiE0(i) = M.FitResult.par(2) + 18575;
    MultierrorStart(i) = M.FitResult.err(2);
    MultiqUStart(i) = M.RunData.qU(i);
end;
figure(1);
errorbar(MultiqUStart,MultiE0,MultierrorStart)