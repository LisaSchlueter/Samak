% error propagation of global sigma (script knm2_GetGlobalB"
% method: MC simulation. draw parameters within uncertainties 
% and sample according to plasma model (3 gaussians)
%
% Lisa, April 2020

%% settings
nSamples = 1e3;
Mode = 'Uniform';

RecomputeFlag = 'ON';

savedir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/'];
MakeDir(savedir);
savename = sprintf('%sknm2_GetGlobalBerr_%s.mat',savedir,Mode);

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
else
%% inputs from Fabian: broadenings + shifts + uncertainties
Weights_v = [171,93,97]./361; % approx. weights -> numbers of scans

 savenameF = sprintf('%sknm2_RManalysis_%s.mat',savedir,Mode);
   load(savenameF,'Sigma','Shift','SigmaErr','ShiftErr');
switch Mode
    case 'Ring'
        nRings = 4;  
    case 'Uniform'
        nRings = 1;
end

%% calculate uncertainty
SigmaGlobal = zeros(nRings,nSamples);
for r=1:nRings
    
    % draw random samples
    Sigma_samples = Sigma(:,r)+SigmaErr(:,r).*randn(size(SigmaErr,1),nSamples);
    Shift_samples = Shift(:,r)+ShiftErr(:,r).*randn(size(ShiftErr,1),nSamples);
    
    for i=1:nSamples
        progressbar(i/nSamples);
        SigmaGlobal(r,i) = knm2_ConvertShiftDrift2GlobalSigma(...
            'shifts_v',Shift_samples(:,i),...
            'sigmas_v',Sigma_samples(:,i),...
            'weights_v',Weights_v,...
            'SanityPlot','OFF');
    end
end
% uncertainty on global sigma
SigmaGlobalErr = std(SigmaGlobal,0,2);

save(savename,'SigmaGlobal','SigmaGlobalErr',...
     'Shift','Sigma','Weights_v','ShiftErr','SigmaErr');
end
%% display result
fprintf('Global broadening uncertainty %s \n',Mode)
fprintf('%.1f meV \n',1e3.*SigmaGlobalErr);


