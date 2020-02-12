function createTDruns(varargin)
%
% Read a Given TD from Databank
% Create N TD's, for N runs, with fluctuations of the qU's
% Assuming an instrumental dispersion of the actual qU values
%
% Thierry Lasserre
% Last updated: April 12 2018
%

% Miscellaneous
format long;

p = inputParser;
p.addParameter('TDi','FT-TL4');
p.addParameter('nRuns',10,@(x)isfloat(x) && x>0);
p.addParameter('qUoffset',1,@(x)isfloat(x) && x>=0); % fluctuation, eV
p.addParameter('ClearTDruns','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('Display','ON',@(x)ismember(x,{'ON','OFF'}));
p.parse(varargin{:});
TDi                = p.Results.TDi;
nRuns              = p.Results.nRuns;
qUoffset           = p.Results.qUoffset;
ClearTDruns        = p.Results.ClearTDruns;
Display            = p.Results.Display;

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

qUrun = zeros(numel(td.qU),nRuns);
% Loop on runs
for i=1:1:nRuns
    % Flat fluctuation, within qU-fluct - qU+fluct - asymetric
    qUrun(:,i) = td.qU + qUoffset .* ( -0.5 + 2.*rand(numel(td.qU),1));
    % Write new TD file, for run i, in TD Databank
    qU = qUrun(:,i);
    TD=sprintf('%s_run%g',TDi,i);
    tdname = sprintf('%s%s.mat',folder,TD);
    save(tdname,'TD','qU','qUfrac');
    fprintf('createTDruns: creating and saving %s\n', tdname);
end

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
        strtitle = sprintf('Stack of %g runs - %g V spread (uniform)',nRuns,qUoffset);
        title(strtitle);       
        subplot(2,1,1)
        errorbar(index,td.qU,qUrunStd,'d','MarkerSize',3,'MarkerFaceColor',.5*[1 1 1],'LineWidth',1);
        grid on
        xlabel('qU index','FontSize',14);
        ylabel('qU (V)','FontSize',14);
        PrettyFigureFormat
        title(strtitle);       
        subplot(2,1,2)
        errorbar(index,td.qU-qUrunMean,qUrunStd,'d','MarkerSize',3,'MarkerFaceColor',.5*[1 1 1],'LineWidth',1);
        grid on
        xlabel('qU index','FontSize',14);
        ylabel('qU - mean(qU) ','FontSize',14);
        PrettyFigureFormat
end

fprintf('-------------------END:  createTDruns ---------------------- \n')
end
