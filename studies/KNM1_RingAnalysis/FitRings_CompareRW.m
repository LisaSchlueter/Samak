% Compare rear wall (RW) voltages: 300mV and 175mV
Rings = 1:12;
nRings =numel(Rings);
E0300 = zeros(nRings,1);
E0Err300 = zeros(nRings,1);
chi2min300 = zeros(nRings,1);
filepath = sprintf('../KNM1_RingAnalysis/results/');
if ~exist(filepath,'dir')
    system(['mkdir ',filepath]);
end
filename = sprintf('FitRingList.mat');

if exist([filepath,filename],'file')
    load([filepath,filename]);
else
%%
for i=1:nRings
    M = MultiRunAnalysis('RunList','KNM1_300mvRW','DataType','Real','exclDataStart',1,...
        'AnaFlag','Ring','RingList',Rings(i));
   M.Fit;
   E0300(i) = M.FitResult.par(2);
   E0Err300(i) = M.FitResult.err(2);
   chi2min300(i) = M.FitResult.chi2min;
end
%%
E0175 = zeros(nRings,1);
E0Err175 = zeros(nRings,1);
chi2min175 = zeros(nRings,1);
for i=1:nRings
    M = MultiRunAnalysis('RunList','KNM1_175mvRW','DataType','Real','exclDataStart',1,...
        'AnaFlag','Ring','RingList',Rings(i));
   M.Fit;
   E0175(i) = M.FitResult.par(2);
   E0Err175(i) = M.FitResult.err(2);
   chi2min175(i) = M.FitResult.chi2min;
end
save([filepath,filename],'Rings','nRings','E0300','E0Err300','chi2min300','E0175','E0Err175','chi2min175');
end
%%
e300 = errorbar(Rings,E0300,E0Err300,'LineWidth',3,'LineStyle','--');
hold on;
e175 = errorbar(Rings,E0175,E0Err175,'LineWidth',3,'LineStyle','--');
PrettyFigureFormat;
xlabel('Ring');
ylabel('E0 -18575.0 (eV)');
hold off;
leg = legend([e300,e175],'300 mV','175 mV');
leg.Title.String = 'real wall offset';
set(gca,'FontSize',20);

