% Lisa - June 2019
% script to plot whole differential spectrum from 0-endpoint
% attion: comment InitializeRF in TBD (otherwise it will compute the response function 
% -> takes a lot of time for this range)

TD = 'Flat18575';
qU = linspace(10,18576,10000)';
qUfrac = (1/numel(qU)).*ones(numel(qU),1);
save([getenv('SamakPath'),'simulation/katrinsetup/TD_DataBank/Flat18575.mat'],...
       'TD','qU','qUfrac');

A = ref_TBD_NominalKATRIN('TD',TD);
%% Draw random energies
A.ComputeTBDDS;

cdfTBDDS = cumsum(A.TBDDS(1:end-1));
cdfTBDDS = cdfTBDDS./cdfTBDDS(end);
nSamples = 100000;

[cdfTBDDS_u,Idx_u,~] = unique(cdfTBDDS);
Te_u = A.Te(1:end-1);
Te_u = Te_u(Idx_u);
cdfInvTBDDS = interp1(cdfTBDDS_u,Te_u,rand(nSamples,1),'spline');

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

samples = cdfInvTBDDS(1:nSamples/nIteration);

for i=1:nIteration
    samples = [samples,cdfInvTBDDS(((i-1).*nSamples/nIteration+1):(i.*nSamples/nIteration))];
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

samples = cdfInvTBDDS(1);
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
    samples = [samples,cdfInvTBDDS(i)];
end



