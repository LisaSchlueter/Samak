% plot/movie of diff. spec (nominal KATRIN)
% (based on Lisa - June 2019)

% load samples from "Compute_SpectrumWholeRange.m"
savedir = [getenv('SamakPath'),'EasyVis/results/'];
savename = sprintf('%sDiffSpecWholeRange.mat',savedir);
load(savename);
%% generate more samples if needed
%nSamples = 1e6;
%Energy_samples = interp1(TBDDS_cdf,Te,rand(nSamples,1),'spline');
   
%% gif1:  histogram diff spec all 1e6 samples
plotdir = [getenv('SamakPath'),'EasyVis/plots/'];

HistStyle = {'FaceColor',rgb('DeepSkyBlue'),'FaceAlpha',1,...
    'EdgeColor','none'};
f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,07]);

MakePretty1;

filename = [plotdir,'DifferentialSpectrum_1e6.gif'];
system(sprintf('rm %s',filename));

nIteration = 100;
nSamples = numel(Energy_samples); % number of entries in final histogram
RandEnergies = 1e-03.*reshape(Energy_samples,nIteration,nSamples/nIteration);

gif(filename,'frame',gcf,'DelayTime',0.2);%,'LoopCount',1)

for i=1:nIteration
    histogram(RandEnergies(1:i,:),'BinWidth',500.*1e-03,...
   HistStyle{:});
    MakePretty1;
    gif
end


%% plot 2: only 20 electrons one after the other
f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,07]);
MakePretty2;


filename2 = [plotdir,'DifferentialSpectrum_20.gif'];
system(sprintf('rm %s',filename2));

gif(filename2,'frame',gcf,'DelayTime',1.5);%,'LoopCount',1)


for i=1:20
    histogram(1e-03.*Energy_samples(1:i),'BinWidth',500.*1e-03,HistStyle{:});
    MakePretty2
    gif
end

%% make pretty
function MakePretty1()
xlabel('Energie (keV)');
ylabel('Anzahl Elektronen');
PrettyFigureFormat;
set(gca,'FontSize',32);
grid off;
ylim([0 5.4e4]);%5500
xlim([0 18575].*1e-03)
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
end


function MakePretty2()
xlabel('Energie (keV)');
ylabel('Anzahl Elektronen');
PrettyFigureFormat;
set(gca,'FontSize',32);
grid off;
ylim([0 2]);%
xlim([0 18575].*1e-03)
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
end


