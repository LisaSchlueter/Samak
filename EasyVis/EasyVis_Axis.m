% plot beta-spectrum axis
savedir = [getenv('SamakPath'),'EasyVis/results/'];
savename = sprintf('%sComputeIntSpec_FitPar_100.mat',savedir);

color_background = 'w';

if strcmp(color_background,'k')
    color_text = 'w';
else
    color_text = 'k';
end

if exist(savename,'file')
    load(savename);
     fprintf('load %s\n',savename)
else
    fprintf(2,'file not found %s \n Run Compute_IntSpec_FitPar.m\n',savename)
    return
end


%%
Language = 'eng';
Format = 'video';%

%% gif1:  signal normalization 
plotdir = [getenv('SamakPath'),'EasyVis/plots/'];

f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,07]);

if strcmp(Format,'gif')
    filename = [plotdir,'Axis.gif'];
    system(sprintf('rm %s',filename));
    gif(filename,'frame',gcf,'DelayTime',1e-10);%,'LoopCount',1)
elseif strcmp(Format,'video')
    if strcmp(Language,'eng')
        video = VideoWriter([plotdir,'Axis_',color_background,'_eng'],'MPEG-4');
    else
        video = VideoWriter([plotdir,'Axis_',color_background],'MPEG-4');
    end
    video.FrameRate = 10;
    open(video);
end

set(gcf,'Color',color_background);

WriteVideo(Format,video,1);
MakePretty1(color_background);

if strcmp(Language,'eng')
     ylabel('Electrons per second');
else
    ylabel('Elektronen pro Sekunde');
end
ax = gca;
ax.XAxis.Color = color_background;

if strcmp(color_background,'k')
    ax.YAxis.Color = color_text;
end

WriteVideo(Format,video,20);
ax.XAxis.Color = color_text;
if strcmp(Language,'eng')
    xlabel('Energy');
else
    xlabel('Energie');
end

WriteVideo(Format,video,20);

for i=1:28
    plot(1e-03.*qU(1:i),TBDIS_data(1:i)./Time(1:i),'x','Color',color_text,'LineWidth',7,'MarkerSize',20);
    MakePretty1(color_background);
    ax.YAxis.Color = color_text;
    ax.XAxis.Color = color_text;
    
    if strcmp(Language,'eng')
        xlabel('Energy');
        ylabel('Electrons per second');
    else
        xlabel('Energie');
        ylabel('Elektronen pro Sekunde');
    end
    
    WriteVideo(Format,video,1);
end

% save several seconds at the end
WriteVideo(Format,video,120);

% save
  if strcmp(Format,'video')
      close(video);
  end
  

%%
function WriteVideo(Format,video,nFrame)
if strcmp(Format,'gif')
    gif
elseif strcmp(Format,'video')
    for i=1:nFrame
        F = getframe(gcf);
        writeVideo(video,F);
    end
end

end

function MakePretty1(color)
PrettyFigureFormat('FontSize',40);
set(gca,'YScale','log');
%set(gca,'FontSize',28);
grid off;
ylim([0.1 100]);%5500
xlim(1e-03.*[18574-41 18575+10])
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
yticks([])
%yticklabels({'1','10','100'});
set(gca,'LineWidth',8);
xticks([])
box off
set(gca,'Color',color);
set(gcf,'Color',color);
end
