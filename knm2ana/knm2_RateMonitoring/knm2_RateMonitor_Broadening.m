% 300 eV analysis
% calculate period-wise broadenings

savedir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/'];
savename = sprintf('%sknm2_RateMonitoring_GetCounts_%s.mat',savedir,'StackPixel');
load(savename);

%
SigmaSqNorm  = zeros(3,1);
SigmaSqPoiss = zeros(3,1);
SigmaSqBroad = zeros(3,1);
ax = cell(3,1);
%%
f2 = figure('Units','normalized','Position',[0.8,0.8,1,0.4]);% breite, höhe

for i=1:3
    if i==1
        Idx_Rw = Idx_Rw1;
    elseif i==2
        Idx_Rw = Idx_Rw2;
    elseif i==3
        Idx_Rw = Idx_Rw3;
    end
    subplot(1,3,i)
    
    tmpname = sprintf('%sknm2_RateMonitor_Broadening_RW%.0f.mat',savedir,i);
    if exist(tmpname,'file')
        d = importdata(tmpname);
        x1 = d.x1; plotx = d.plotx;
        pdfPois = d.pdfPois; pdfNorm = d.pdfNorm;
        SigmaSqNorm(i) = d.SigmaSqNorm_t;
        SigmaSqPoiss(i) = d.SigmaSqPoiss_t;
        fprintf('load %s \n',tmpname);
    else
        [x1,plotx,pdfNorm,pdfPois,SigmaSqNorm(i),SigmaSqPoiss(i)] = Fit(Idx_Rw,RatesCorr,TimeSec_ScanStep);
        SigmaSqNorm_t = SigmaSqNorm(i);
        SigmaSqPoiss_t = SigmaSqPoiss(i);
        save(tmpname,'x1','plotx','pdfNorm','pdfPois','SigmaSqPoiss_t','SigmaSqNorm_t');
    end
    h1 = histogram(x1,'FaceColor',rgb('Silver'),'FaceAlpha',1,'EdgeColor',rgb('DarkGray'));
    hold on;
    lN = plot(plotx,pdfNorm.*h1.BinWidth.*numel(x1),'LineWidth',3,'Color',rgb('Crimson'));
    lP = plot(plotx,pdfPois.*h1.BinWidth.*numel(x1),'-.','LineWidth',3,'Color',[0.929,0.694,0.125]);
    xlim([min(x1),max(x1)]);
    PrettyFigureFormat('FontSize',22);
    xlabel(sprintf('Counts in %.0f s',mean(TimeSec_ScanStep)));
    ylabel('Occurence');
    
    leg = legend([lN,lP],....
        sprintf('\\sigma_G^2 = %.3f \\times 10^7 counts^2',1e-07.*SigmaSqNorm(i)),...
        sprintf('\\sigma_P^2 = %.3f \\times 10^7 counts^2',1e-07.*SigmaSqPoiss(i)),...
        'Location','north');
    PrettyLegendFormat(leg);
    leg.Title.String = sprintf('Period %.0f',i);
    leg.Title.FontWeight = 'normal';
    leg.FontSize = get(gca,'FontSize');
    
    if i==1
        ylim([0 126])
    elseif i==2
         ylim([0 72])
    elseif i==3
          ylim([0 69])
    end
    ax{i} = gca;
end


pltdir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/plots/'];
MakeDir(pltdir);
pltname = sprintf('%sknm2_RateMonitoring_BroadeningHist_%s.pdf',pltdir,'StackPixel');
export_fig(pltname);
fprintf('save plot to %s \n',pltname);


%%

function [x1,plotx,pdfNorm,pdfPois,SigmaSqNorm,SigmaSqPoiss] = Fit(Idx,RatesCorr,TimeSec_ScanStep)
x1 = (RatesCorr(Idx).*mean(TimeSec_ScanStep))';
x1(abs(x1-mean(x1))>5*std(x1)) = [];

dPois = fitdist(x1,'Poisson');
dNorm = fitdist(x1,'Normal');

plotx = linspace(min(x1),max(x1),1e2);
pdfNorm = pdf(dNorm,plotx);
pdfPois= interp1(0:1:round(max(plotx)),poisspdf(0:1:round(max(plotx)),dPois.lambda),plotx,'spline');

SigmaSqNorm  = dNorm.sigma.^2;
SigmaSqPoiss = dPois.lambda;
end

