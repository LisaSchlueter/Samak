% calculate spectrum (samples) for trititum spectrum
% save for later use
% nominal KATRIN
% based ony script from (Lisa - June 2019)
% script to calculate differential spectrum - samples


savedir = [getenv('SamakPath'),'EasyVis/results/'];
savename = sprintf('%sDiffSpec40eVRange.mat',savedir);

if exist(savename,'file')
else
    TD = 'Flat40l';
    qU = linspace(18575-40,18575+10,10000)';
    qUfrac = (1/numel(qU)).*ones(numel(qU),1);
    save([getenv('SamakPath'),'simulation/katrinsetup/TD_DataBank/Flat40l.mat'],...
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

return

%% display
f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,07]);
PrettyFigureFormat;
xlabel('Energie der Elektronen (eV)');
ylabel('Anzahl Elektronen');
set(gca,'FontSize',32);
ylim([0 5500]);
xlim([0 18575]);

filename = 'DifferentialSpectrum.gif';
nIteration = 100;
gif(filename,'frame',gcf,'DelayTime',0.2);%,'LoopCount',1)

samples = Energy_samples(1:nSamples/nIteration);

for i=1:nIteration
    samples = [samples,Energy_samples(((i-1).*nSamples/nIteration+1):(i.*nSamples/nIteration))];
    histogram(samples,'BinWidth',500,'FaceColor',rgb('SteelBlue'),'EdgeColor',rgb('DarkSlateGray'));
    xlabel('Energie des Elektrons (eV)');
    ylabel('Anzahl Elektronen');
    PrettyFigureFormat;
    set(gca,'FontSize',32);
    grid on;
    ylim([0 5500])
    xlim([0 18575])
    gif
end


%% second gif 20 electrons one after the other
f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,07]);
PrettyFigureFormat;
xlabel('Energie der Elektronen (eV)');
ylabel('Anzahl Elektronen');
set(gca,'FontSize',32);
ylim([0 2]);
xlim([0 18575]);
filename = 'DifferentialSpectrum_20slow.gif';
nSamples = 10;
nIteration = nSamples; % one after the other
gif(filename,'frame',gcf,'DelayTime',1.5);%,'LoopCount',1)

samples = Energy_samples(1);
for i=1:nIteration
    histogram(samples,'BinWidth',500,'FaceColor',rgb('SteelBlue'),'EdgeColor',rgb('DarkSlateGray'));
    xlabel('Energie des Elektrons (eV)');
    ylabel('Anzahl Elektronen');
    PrettyFigureFormat;
    set(gca,'FontSize',32);
    grid on;
    xlim([0 18575])
    ylim([0 2]);
    gif
    samples = [samples,Energy_samples(i)];
end



