% error propagation of global sigma (script knm2_GetGlobalB"
% method: MC simulation. draw parameters within uncertainties 
% and sample according to plasma model (3 gaussians)
%
% Lisa, April 2020

%% settings
nSamples = 1e3;
Mode = 'Uniform';

RecomputeFlag = 'OFF';

savedir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/'];
MakeDir(savedir);
savename = sprintf('%sknm2_GetGlobalBerr_%s.mat',savedir,Mode);

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
else
%% inputs from Fabian: broadenings + shifts + uncertainties
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
        
        Sigma_M_err = (1e-03.*[7.2, 4.1, 3.6, 2.4;...
            5.1, 1.8, 3.1, 0.7;...
            0.9, 0.6, 0.1, 0.1]);
        Shift_M_err = 1e-03.*[4.9, 4.3, 4.5, 6.0; ...
            0, 5.4, 5.5, 7.3; ...
            6.1, 5.4, 5.5 ,7.3];
    case 'Uniform'
        nRings = 1;
        Sigma_M     = 1e-03.*[31.9,20.3,28.25]';
        Shift_M     = 1e-03.*[41.8,0,105.8]';
        Sigma_M_err = 1e-03.*[8.1,4.7,0]';
        Shift_M_err = 1e-03.*[1.8,0,2.1]';
     
end

%% calculate uncertainty
SigmaGlobal = zeros(nRings,nSamples);
for r=1:nRings
    
    % draw random samples
    Sigma_M_samples = Sigma_M(:,r)+Sigma_M_err(:,r).*randn(size(Sigma_M_err,1),nSamples);
    Shift_M_samples = Shift_M(:,r)+Shift_M_err(:,r).*randn(size(Shift_M_err,1),nSamples);
    
    for i=1:nSamples
        progressbar(i/nSamples);
        SigmaGlobal(r,i) = knm2_ConvertShiftDrift2GlobalSigma(...
            'shifts_v',Shift_M_samples(:,i),...
            'sigmas_v',Sigma_M_samples(:,i),...
            'weights_v',Weights_v,...
            'SanityPlot','OFF');
    end
end
% uncertainty on global sigma
SigmaGlobalErr = std(SigmaGlobal,0,2);

save(savename,'SigmaGlobal','SigmaGlobalErr',...
     'Shift_M','Sigma_M','Weights_v','Shift_M_err','Sigma_M_err');
end
%% display result
fprintf('Global broadening uncertainty %s \n',Mode)
fprintf('%.1f meV \n',1e3.*SigmaGlobalErr);


