function FT_SimMTD_multiruns(varargin)
%
% Read a Given TD from Databank
% Create N TD's, for N runs, with fluctuations of the qU's
% Assuming an instrumental dispersion of the actual qU values
%
% Thierry Lasserre
%  Lisa Schlueter
% Last updated: May 24 2018
%

% Miscellaneous
close all;
format long;

p = inputParser;
p.addParameter('TDi','FT-TL3');
p.addParameter('nRuns',44,@(x)isfloat(x) && x>0);
p.addParameter('qUoffset_l',-0.075,@(x)isfloat(x)); % fluctuation, eV
p.addParameter('qUoffset_r',+0.075,@(x)isfloat(x)); % fluctuation, eV
p.addParameter('ClearTDruns','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('Display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('SanityPlots','OFF',@(x)ismember(x,{'ON','OFF'}));
p.parse(varargin{:});
TDi                = p.Results.TDi;
nRuns              = p.Results.nRuns;
qUoffset_l         = p.Results.qUoffset_l;
qUoffset_r         = p.Results.qUoffset_r;
ClearTDruns        = p.Results.ClearTDruns;
Display            = p.Results.Display;
SanityPlots        = p.Results.SanityPlots;

% Clear TD Files
folder = '../../simulation/katrinsetup/TD_DataBank/';
switch ClearTDruns
    case 'ON'
        for i=1:1:nRuns
            TD=sprintf('%s_run%g',TDi,i);
            filename = sprintf('%s.mat',TD);
            if exist(filename, 'file') == 2
                command = sprintf('rm %s%s.mat',folder,TD);
                system(command);
                fprintf('createTDruns: Delete %s%s\n', folder,filename);
            end
        end
        return;
end

fprintf('-------------------START: createTDruns ---------------------- \n')

% Read TD Databank
tdfile = sprintf('%s%s.mat',folder,TDi);
td = load(tdfile);
qUfrac = td.qUfrac;
%init
qUrun = zeros(numel(td.qU),nRuns);
shift = zeros(numel(td.qU),nRuns); % (nTrial, qU)

for i=1:1:nRuns % Loop on runs
    for qui =1:numel(td.qU) % Loop on qU
    shift(qui,i) = unifrnd(qUoffset_l(qui),qUoffset_r(qui),1);
    end 
    qUrun(:,i) = td.qU + shift(:,i);
    
    if strcmp(SanityPlots,'ON')
    fig1Scatter = figure(3);
    set(fig1Scatter, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.8, 1.2]);
    c = repmat(colormap(hsv),ceil(nRuns/64),1);
    scatter(td.qU-18575,shift(:,i),80,c(i,:),'o','filled','MarkerEdgeColor',rgb('DimGray'));
    PrettyFigureFormat;
    hold on;
    set(gca,'FontSize',16);
    xlabel('<qU> - 18575  (eV)','FontSize',20);
    ylabel('qU - <qU> (eV)','FontSize',20);
    title(sprintf('Stacking Runs \n qU-distribution for %u Sample Runs',nRuns),'FontSize',15);
    end
      
    % Write new TD file, for run i, in TD Databank
    qU = qUrun(:,i);
    TD=sprintf('%s_run%g',TDi,i);
    tdname = sprintf('%s%s.mat',folder,TD);
    save(tdname,'TD','qU','qUfrac');
    fprintf('createTDruns: creating and saving %s\n', tdname);
end
plot([0.9*min(td.qU-18575) 1.1*max(td.qU-18575)],[0 0],'k--');
xlim([0.9*min(td.qU-18575) 1.1*max(td.qU-18575)]);
grid on;
hold off;
export_fig('./plots/CMStacking_Sample-qUDistribution.eps');
 
switch Display
    case 'ON'
        qUrunMean = zeros(numel(td.qU),1);
        qUrunStd = zeros(numel(td.qU),1);
        index = 1:1:numel(td.qU);
        for h=index
            qUrunMean(h) = mean(qUrun(h,:));
            qUrunStd(h)  = std(qUrun(h,:));
        end
        figure(999)
        errorbar(index,td.qU-qUrunMean,qUrunStd,'d','MarkerSize',3,'MarkerFaceColor',.5*[1 1 1],'LineWidth',1);
        plt = Plot();
        plt.LineWidth = 1;
        plt.LineStyle = '--';
        plt.Markers = {'s'};
        pltstr     = sprintf('Stacking of %g runs',nRuns);
        plt.Title  = pltstr; % plot title
        plt.XLabel = 'qU index'; % xlabel
        plt.YLabel = 'qU_{ref} - mean(qU_{sim})'; %ylabel
        plt.YScale = 'lin'; % 'linear' or 'log'
        %plt.YLim       = [0 1.1];
        plt.XScale = 'lin'; % 'linear' or 'log'
        plt.FontSize = 16;
        plttitle   = sprintf('StackRun%g-%s-Mean.png',nRuns,TD);
        plt.export(plttitle);
        
        figure(9999)
        hold on
        for i=1:1:nRuns
        errorbar(index,td.qU-qUrun(:,i),index*0);
        end
        hold off
        clear plt2; plt2 = Plot();
        plt2str     = sprintf('Stacking of %g runs',nRuns);
        plt2.LineStyle  = 'none';
        plt2.Title  = plt2str; % plot title
        plt2.XLabel = 'qU index'; % xlabel
        plt2.YLabel = 'qU_{ref} - qU_{sim}'; %ylabel
        plt2.YScale = 'lin'; % 'linear' or 'log'
        plt2.XScale = 'lin'; % 'linear' or 'log'
        plt2.FontSize = 16;
        plt2title   = sprintf('StackRun%g-%ss-Dist.png',nRuns,TD);
        plt2.export(plt2title);
end

fprintf('-------------------END:  createTDruns ---------------------- \n')
end
