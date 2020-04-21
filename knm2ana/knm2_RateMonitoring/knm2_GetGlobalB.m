
Mode = 'Ring';
Weights_v = [171,93,97]./361; % approx. weights -> numbers of scans
switch Mode
    case 'Ring'
        nRings = 4;
        Sigma_M = (1e-03.*[18.7, 35.7, 29.0, 16.4;...
                           22.9, 14.1, 8.9, 10.5;...
                           26.8, 26.5, 32.9, 30.9]);    
        Shift_M = 1e-03.*[73.5, 84.6, 90.5, 120.3;...
                           0 , 43.2, 60.2, 106.8;...
                          120.5, 148.9, 161.9, 197.7];
    case 'Uniform'
        nRings = 1;
        Sigma_M = 1e-03.*[31.9,  20.3 , 28.25]';
        Shift_M = [];
end
SigmaGlobal = zeros(nRings,1);
plotdir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/plots/'];
MakeDir(plotdir);

for i=1:nRings
    
     plotname = sprintf('%sknm2_GetGlobalSigma_%s.png',plotdir,Mode);
    if nRings>1
     plotname = strrep(plotname,'.png',sprintf('_Ring%.0f.png',i));
    end
    
    SigmaGlobal(i)  =knm2_ConvertShiftDrift2GlobalSigma('shifts_v',Shift_M(:,i),...
        'sigmas_v',Sigma_M(:,i),...
        'weights_v',Weights_v,...
        'SanityPlot','OFF','PlotName',plotname);
  
end