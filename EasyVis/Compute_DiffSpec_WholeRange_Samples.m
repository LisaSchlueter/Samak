% calculate spectrum (samples) for trititum spectrum
% save for later use
% nominal KATRIN
% based ony script from (Lisa - June 2019)
% script to calculate whole differential spectrum from 0-endpoint
%-> takes a lot of time for this range)
% may want to comment "InitalizeRF" in TBD.m constructor
savedir = [getenv('SamakPath'),'EasyVis/results/'];
savename = sprintf('%sDiffSpecWholeRange.mat',savedir);

if exist(savename,'file')
else
    TD = 'Flat18575';
    qU = linspace(10,18576,10000)';
    qUfrac = (1/numel(qU)).*ones(numel(qU),1);
    save([getenv('SamakPath'),'simulation/katrinsetup/TD_DataBank/Flat18575.mat'],...
        'TD','qU','qUfrac');
    
    A = ref_TBD_NominalKATRIN('TD',TD);
    %% Draw random energies
    nSamples = 1e6; % number of entries in final histogram
    
    % diff. spec
    A.ComputeTBDDS;
    TBDDS = A.TBDDS(1:end-1);
    TBDDS_cdf = cumsum(TBDDS);
    TBDDS_cdf = TBDDS_cdf./TBDDS_cdf(end);
    
    % unique for interpolation
    [TBDDS_cdf,Idx_unique,~] = unique(TBDDS_cdf);
    Te = A.Te(1:end-1);
    Te = Te(Idx_unique);
    TBDDS = TBDDS(Idx_unique);
    
    % uniform sampling over inverse function
    % --> draw random energie values following shape of diff. spec
    Energy_samples = interp1(TBDDS_cdf,Te,rand(nSamples,1),'spline');
    
    %% save
    
    save(savename,'Te','TBDDS','TBDDS_cdf','Energy_samples');
    fprintf('save %s \n',savename);
end


