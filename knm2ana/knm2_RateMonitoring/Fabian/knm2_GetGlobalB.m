% convert period-wise broadening + drift into 1 global broadening
% method: MC simulation. samples according to plasma model (3 gaussians)
% corresponding uncertainties are calculating in script "knm2_GetGlobalBerr"
% Lisa, April 2020

%% setttings
SanityPlot    = 'ON'; % plot and save plot
Mode          = 'Ring';       % Ring or Uniform
RecomputeFlag = 'ON';

savedir  = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/'];
MakeDir(savedir);
savename = sprintf('%sknm2_GetGlobalB_%s.mat',savedir,Mode);

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
else
    %% inputs from Fabian: broadenings + shifts
   savenameF = sprintf('%sknm2_RManalysis_%s.mat',savedir,Mode);
   load(savenameF,'Sigma','Shift');
 
    Weights_v = [171,93,97]./361; % approx. weights -> numbers of scans
    switch Mode
        case 'Ring'
            nRings = 4;     
        case 'Uniform'
            nRings = 1;  
    end
    
    %% calculate global broadening for uniform or pseudo-ring wise
    SigmaGlobal = zeros(nRings,1);
    plotdir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/plots/'];
    MakeDir(plotdir);
    for i=1:nRings
        
        plotname = sprintf('%sknm2_GetGlobalSigma_%s.pdf',plotdir,Mode);
        if nRings>1
            plotname = strrep(plotname,'.pdf',sprintf('_%.0f.pdf',i));
        end
        
        SigmaGlobal(i)  = knm2_ConvertShiftDrift2GlobalSigma('shifts_v',Shift(:,i),...
            'sigmas_v',Sigma(:,i),...
            'weights_v',Weights_v,...
            'SanityPlot',SanityPlot,'PlotName',plotname);
    end
    
    save(savename,'SigmaGlobal','Shift','Sigma','Weights_v');
end
%% display result
fprintf('Global broadening = %.1f meV \n',1e3.*SigmaGlobal);
fprintf('Global broadening = %.3g eV^2 \n',SigmaGlobal.^2);