% illustrate chi2 minimization
% based on knm2 data and model

savedir = [getenv('SamakPath'),'EasyVis/results/'];
savename = sprintf('%sComputeIntSpec_FromRealChi2Profile.mat',savedir);

if exist(savename,'file')
    load(savename);
     fprintf('load %s\n',savename)
else
    fprintf(2,'file not found %s \n Run Compute_IntSpec_FitPar.m\n',savename)
    return
end

%% gif 1: spectrum with (random) fit parameter values and matching chi^2
pause('on')
plotdir = [getenv('SamakPath'),'EasyVis/plots/'];

f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.9 ,0.7]);
subplot(1,4,1:3)
MakePretty1;

nIteration = size(TBDIS,2);

filename = [plotdir,'FitPar_Chi2_RealProfile.gif'];
system(sprintf('rm %s',filename));
gif(filename,'frame',gcf,'DelayTime',0.3,'LoopCount',1)

%Chi2Idx = flipud(find(chi2min<1e6));
Chi2Idx = 1:numel(chi2min);

for i=1:numel(Chi2Idx)%nIteration:-2:1
    subplot(1,4,1:3)
   
    p = plot(1e-03.*qU,TBDIS(:,Chi2Idx(i))./Time,'-','LineWidth',4,'Color',rgb('Orange'));
     hold on;
    for j = 1:28
        if TBDIS_data(j)<TBDIS(j,Chi2Idx(i))
            plot(1e-03.*qU(j).*ones(1,100),linspace(TBDIS_data(j)./Time(j),TBDIS(j,Chi2Idx(i))./Time(j),100),':','LineWidth',4,'Color',rgb('Red'));
        else
            plot(1e-03.*qU(j).*ones(1,100),linspace(TBDIS(j,Chi2Idx(i))./Time(j),TBDIS_data(j)./Time(j),100),':','LineWidth',4,'Color',rgb('Red'));
            
        end
    end
    
      plot(1e-03.*qU,TBDIS_data./Time,'.','MarkerSize',15,'LineWidth',4,'Color',rgb('Black'));
  
    hold off
    MakePretty1;
    
    subplot(1,4,4)
    bar(chi2min(Chi2Idx(i)),'FaceAlpha',1,'FaceColor',rgb('Red'),'EdgeColor','none','BarWidth',0.2);
    PrettyFigureFormat('FontSize',28)
    xticks([]);
    yticks([]);
     ylim([0.1,max(chi2min(Chi2Idx))]);
   % ylim([0.1,max(chi2min)]);
  %  set(gca,'YScale','log');
    xlabel('Unterschied');
    ax = gca;
    ax.YAxis.Color = 'w';
   % ax.XAxis.Color = 'w';
    box off
    gif
end

%% best fit: do a couple of twice to last longer
  subplot(1,4,1:3)
   
    p = plot(1e-03.*qU,TBDIS_bf./Time,'-','LineWidth',4,'Color',rgb('DarkOrange'));
     MakePretty1;
     hold on;
    for j = 1:28
        if TBDIS_data(j)<TBDIS_bf(j)
            plot(1e-03.*qU(j).*ones(1,100),linspace(TBDIS_data(j)./Time(j),TBDIS_bf(j)./Time(j),100),':','LineWidth',4,'Color',rgb('Red'));
        else
            plot(1e-03.*qU(j).*ones(1,100),linspace(TBDIS_bf(j)./Time(j),TBDIS_data(j)./Time(j),100),':','LineWidth',4,'Color',rgb('Red'));
            
        end
    end
    
      plot(1e-03.*qU,TBDIS_data./Time,'.','MarkerSize',15,'LineWidth',4,'Color',rgb('Black'));
  text(18.57,20,sprintf('Bester Fit'),'FontSize',get(gca,'FontSize')+4,'Color',p.Color,'HorizontalAlignment','center');
  
    hold off
   
    
    subplot(1,4,4)
    bar(chi2_bf,'FaceAlpha',1,'FaceColor',rgb('Red'),'EdgeColor','none','BarWidth',0.2);
    PrettyFigureFormat('FontSize',26)
    xticks([]);
    yticks([]);
    ylim([0.1,max(chi2min(Chi2Idx))]);
  %  set(gca,'YScale','log');
    xlabel('Unterschied');
    ax = gca;
    ax.YAxis.Color = 'w';
   % ax.XAxis.Color = 'w';
    box off
    
    for k=1:20 % write frame a couple of time, to pause before repeating gif
    gif
    end
%%
function MakePretty1()
xlabel('Energie (keV)');
ylabel('Elektronen pro Sekunde');
PrettyFigureFormat('FontSize',26);
set(gca,'YScale','log');
%set(gca,'FontSize',28);
grid off;
ylim([0.1 100]);%5500
xlim(1e-03.*[18574-39.8 18575+10])
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
xticks([18.53:0.01:18.58]);
yticklabels({'0.1','1','10','100'});

end
