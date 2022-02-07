% plot/movie of diff. spec (nominal KATRIN)
% (based on Lisa - June 2019)

% load samples from "Compute_SpectrumWholeRange.m"
savedir = [getenv('SamakPath'),'EasyVis/results/'];
savename = sprintf('%sDiffSpecWholeRange.mat',savedir);
load(savename);
%% generate more samples if needed
%nSamples = 1e6;
%Energy_samples = interp1(TBDDS_cdf,Te,rand(nSamples,1),'spline');
   
%% gif1:  plot
plotdir = [getenv('SamakPath'),'EasyVis/plots/'];

HistStyle = {'FaceColor',rgb('DeepSkyBlue'),'FaceAlpha',1,...
    'EdgeColor','none'};
f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,07]);

MakePretty1;

filename = [plotdir,'DifferentialSpectrum.gif'];
system(sprintf('rm %s',filename));

nIteration = 10;
nSamples = numel(Energy_samples); % number of entries in final histogram
RandEnergies = 1e-03.*reshape(Energy_samples,nIteration,nSamples/nIteration);

gif(filename,'frame',gcf,'DelayTime',0.2);%,'LoopCount',1)

%CumIdx = 1;
for i=1:nIteration
    histogram(RandEnergies(1:i,:),'BinWidth',0.5.*1e-03,...
       HistStyle{:});
    MakePretty1;

    gif
    
   % CumIdx = i+1;
end

%%
% %% second gif 20 electrons one after the other
% f3 = figure('Renderer','opengl');
% set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,07]);
% PrettyFigureFormat;
% xlabel('Energie der Elektronen (eV)');
% ylabel('Anzahl Elektronen');
% set(gca,'FontSize',32);
% ylim([0 2]);
% xlim([0 18575]);
% filename = 'DifferentialSpectrum_20slow.gif';
% nSamples = 10;
% nIteration = nSamples; % one after the other
% gif(filename,'frame',gcf,'DelayTime',1.5);%,'LoopCount',1)
% 
% samples = Energy_samples(1);
% for i=1:nIteration
%     histogram(samples,'BinWidth',500,'FaceColor',rgb('SteelBlue'),'EdgeColor',rgb('DarkSlateGray'));
%     xlabel('Energie des Elektrons (eV)');
%     ylabel('Anzahl Elektronen');
%     PrettyFigureFormat;
%     set(gca,'FontSize',32);
%     grid on;
%     xlim([0 18575])
%     ylim([0 2]);
%     gif
%     samples = [samples,Energy_samples(i)];
% end

%% make pretty
function MakePretty1()
xlabel('Energie des Elektrons (keV)');
ylabel('Anzahl Elektronen');
PrettyFigureFormat;
set(gca,'FontSize',32);
grid off;
ylim([0 4.2e4]);
xlim([18575-40 18575].*1e-03)
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
ax = gca;
end


