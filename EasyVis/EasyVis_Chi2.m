% illustrate chi2 minimization
% based on knm2 data and model

savedir = [getenv('SamakPath'),'EasyVis/results/'];
savename = sprintf('%sComputeIntSpec_Chi2.mat',savedir);

if exist(savename,'file')
    load(savename);
     fprintf('load %s\n',savename)
else
    fprintf(2,'file not found %s \n Run Compute_IntSpec_FitPar.m\n',savename)
    return
end

Format = 'video';
Language = 'eng';
color_background = 'w';

if strcmp(color_background,'k')
    color_text = 'w';
else
    color_text = 'k';
end
%% gif 1: spectrum with (random) fit parameter values and matching chi^2
pause('on')
plotdir = [getenv('SamakPath'),'EasyVis/plots/'];

f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8 ,0.7]);
subplot(1,4,1:3)
MakePretty1(color_background,Language);

nIteration = size(TBDIS,2);

pltname = [plotdir,'FitPar_Chi2_',color_background];
if strcmp(Language,'eng')
    pltname = [pltname,'_eng'];
end

if strcmp(Format,'gif')
    filename = [pltname,'.gif'];
    system(sprintf('rm %s',filename));
    gif(filename,'frame',gcf,'DelayTime',0.3,'LoopCount',1)
    
elseif strcmp(Format,'video')
    video = VideoWriter(pltname,'MPEG-4');
    video.FrameRate = 3;
    open(video);
end

Chi2Idx = flipud(find(chi2min<1e6));

for i=1:numel(Chi2Idx)%nIteration:-2:1
    subplot(1,4,1:3)
   
    p = plot(1e-03.*qU,TBDIS(:,Chi2Idx(i))./Time,'-','LineWidth',6,'Color',rgb('SlateGray'));
     hold on;
    for j = 1:28
        if TBDIS_data(j)<TBDIS(j,Chi2Idx(i))
            plot(1e-03.*qU(j).*ones(1,100),linspace(TBDIS_data(j)./Time(j),TBDIS(j,Chi2Idx(i))./Time(j),100),':','LineWidth',4,'Color',rgb('Red'));
        else
            plot(1e-03.*qU(j).*ones(1,100),linspace(TBDIS(j,Chi2Idx(i))./Time(j),TBDIS_data(j)./Time(j),100),':','LineWidth',4,'Color',rgb('Red'));
            
        end
    end
    
    plot(1e-03.*qU,TBDIS_data./Time,'.','MarkerSize',25,'Color',color_text);
     hold off
    MakePretty1(color_background,Language);
    ax1 = gca;
 
    subplot(1,4,4)
    bar(chi2min(Chi2Idx(i)),'FaceAlpha',1,'FaceColor',rgb('OrangeRed'),'EdgeColor','none','BarWidth',0.2);
    PrettyFigureFormat('FontSize',28)
     set(gca,'Color',color_background);
    set(gcf,'Color',color_background);
    xticks([]);
    yticks([]);
     ylim([0.1,max(chi2min(Chi2Idx))]);
   % ylim([0.1,max(chi2min)]);
  %  set(gca,'YScale','log');
  if strcmp(Language,'eng')
      xlabel('Difference');
  else
      xlabel('Unterschied');
  end
  set(gca,'LineWidth',3);
  ax2 = gca;
  ax2.YAxis.Color = color_background;
  ax2.XAxis.Color = color_text;
  box off
    

    %%
    
%      if i==1
%         print(gcf,[plotdir,'EasyVis_Chi2_1.png'],'-dpng','-r500');
%      end
    if strcmp(Format,'gif')
        gif
    elseif strcmp(Format,'video')
      % NFrame = ceil((numel(Chi2Idx)-i+1)/10);
        F = getframe(gcf);
      %  for j=1:NFrame
            writeVideo(video,F);
      %  end
    end
end

%% best fit: do a couple of twice to last longer
  subplot(1,4,1:3)
   
    p = plot(1e-03.*qU,TBDIS_bf./Time,'-','LineWidth',6,'Color',rgb('SlateGray'));
     MakePretty1(color_background,Language);
     hold on;
    for j = 1:28
        if TBDIS_data(j)<TBDIS_bf(j)
            plot(1e-03.*qU(j).*ones(1,100),linspace(TBDIS_data(j)./Time(j),TBDIS_bf(j)./Time(j),100),':','LineWidth',5,'Color',rgb('Red'));
        else
            plot(1e-03.*qU(j).*ones(1,100),linspace(TBDIS_bf(j)./Time(j),TBDIS_data(j)./Time(j),100),':','LineWidth',5,'Color',rgb('Red'));
            
        end
    end
    
    plot(1e-03.*qU,TBDIS_data./Time,'.','MarkerSize',25,'LineWidth',4,'Color',color_text);
    if ~strcmp(Language,'eng')
        t =text(18.57,12,sprintf('Beste Übereinstimmung \nmit {\\itm}_\\nu^2 = 0.26 eV^2'),'FontSize',get(gca,'FontSize')+4,'Color',p.Color,'HorizontalAlignment','center','FontWeight','bold');
    end
    ax1 = gca;
    hold off
    
    
    subplot(1,4,4)
    bar(chi2_bf,'FaceAlpha',1,'FaceColor',rgb('Red'),'EdgeColor','none','BarWidth',0.2);
    PrettyFigureFormat('FontSize',26)
    set(gca,'Color',color_background);
    set(gcf,'Color',color_background);
    xticks([]);
    yticks([]);
    ylim([0.1,max(chi2min(Chi2Idx))]);
    %  set(gca,'YScale','log');
    if strcmp(Language,'eng')
        xlabel('Difference');
    else
        xlabel('Unterschied');
    end
    set(gca,'LineWidth',3);
    ax2 = gca;
    ax2.YAxis.Color = color_background;
    ax2.XAxis.Color = color_text;
    box off
    
    
    if strcmp(Format,'gif')
        for k=1:20 % write frame a couple of time, to pause before repeating gif
            gif
        end
    elseif strcmp(Format,'video')
        for k=1:20
            F = getframe(gcf);
            writeVideo(video,F);
        end
    end
    
     % save
  if strcmp(Format,'video')
      close(video);
  end
%%
function MakePretty1(color,Language)
if strcmp(Language,'eng')
    xlabel('Energy (keV)');
    ylabel('Electrons per second');
else
    xlabel('Energie (keV)');
    ylabel('Elektronen pro Sekunde');
end
PrettyFigureFormat('FontSize',28);
set(gca,'YScale','log');
%set(gca,'FontSize',28);
grid off;
ylim([0.1 100]);%5500
xlim(1e-03.*[18574-40.5 18575+10])
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
yticks([1 10 100])
yticklabels({'1','10','100'});
set(gca,'LineWidth',4);
xticks([18.54 18.56 18.58])
box off

set(gcf,'Color',color);
set(gca,'Color',color);
ax = gca;
if strcmp(color,'k') ||  strcmp(color,'none')
    ax.YAxis.Color = 'w';
    ax.XAxis.Color = 'w';
end
end

function DrawCircle
    R = 1. ;
    % center of cricle
    C = [0. 0.] ;
    % angle
    th = 90:-1:0 ;
    % points of cricle
    xc = C(1)+R*cosd(th) ;
    yc = C(2)+R*sind(th) ;
    plot(xc,yc)
    pmNu =area(-xc,yc,'FaceColor',rgb('Red'),'EdgeColor','none');
   
    hold on;
    pE0 =area(xc,yc,'FaceColor',rgb('Orange'),'EdgeColor','none');
    pN =area(-xc,-yc,'FaceColor',rgb('LimeGreen'),'EdgeColor','none');
    pN =area(xc,-yc,'FaceColor',rgb('DeepSkyBlue'),'EdgeColor','none');
    
    text(-0.5,0.5,sprintf('%.0f',randi(100)-1),'FontSize',20,'HorizontalAlignment','center');
    text(-0.5,-0.5,sprintf('%.0f',randi(100)-1),'FontSize',20,'HorizontalAlignment','center');
    text(0.5,-0.5,sprintf('%.0f',randi(100)-1),'FontSize',20,'HorizontalAlignment','center');
    text(0.5,0.5,sprintf('%.0f',randi(100)-1),'FontSize',20,'HorizontalAlignment','center');
    
    ax3 = gca;
    ax3.YAxis.Color = 'w';
    ax3.XAxis.Color = 'w';
    ax3.Position = [0.6,0.7,0.09,0.16];
    
end