% plot beta-spectrum model slowly

savedir = [getenv('SamakPath'),'EasyVis/results/'];
savename = sprintf('%sComputeIntSpec_Chi2.mat',savedir);

if exist(savename,'file')
    load(savename);
     fprintf('load %s\n',savename)
else
    fprintf(2,'file not found %s \n Run Compute_IntSpec_FitPar.m\n',savename)
    return
end


%%
Language = 'eng';
bf = 'ON';
Format = 'video';%
Chi2Idx = flipud(find(chi2min<1e6));
TBDIS = TBDIS(:,Chi2Idx);
%% gif1:  signal normalization 
plotdir = [getenv('SamakPath'),'EasyVis/plots/'];

f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,07]);

filename = 'Model';
if strcmp(bf,'ON')
    filename = [filename,'_bf'];
end

if strcmp(Language,'eng')
    filename = [filename,'_eng'];
end

if strcmp(Format,'gif')
    filename = [plotdir,filename,'.gif'];
    system(sprintf('rm %s',filename));
    gif(filename,'frame',gcf,'DelayTime',1e-10);%,'LoopCount',1)
elseif strcmp(Format,'video')
    video = VideoWriter([plotdir,filename],'MPEG-4');
    video.FrameRate = 10;
    open(video);
end
WriteVideo(Format,video,1);

plot(1e-03.*qU,TBDIS_data./Time,'kx','LineWidth',7,'MarkerSize',20);
MakePretty1(Language);
WriteVideo(Format,video,1);

qU_inter = linspace(qU(1),qU(end-3),50);
if strcmp(bf,'ON')
    TBDIS_inter = interp1(qU(1:end-3),TBDIS_bf(1:end-3)./Time(1:end-3),qU_inter,'spline');
else
    TBDIS_inter = interp1(qU(1:end-3),TBDIS(1:end-3,end-20)./Time(1:end-3),qU_inter,'spline');
end

for i=1:numel(qU_inter)
    plot(1e-03.*qU,TBDIS_data./Time,'kx','LineWidth',7,'MarkerSize',20);
    hold on;
    plot(1e-03.*qU_inter(1:i),TBDIS_inter(1:i),'-','Color',rgb('Red'),'LineWidth',7,'MarkerSize',20);
    hold off;
    MakePretty1(Language);
    WriteVideo(Format,video,1);
end

  WriteVideo(Format,video,100);
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

function MakePretty1(language)
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
if strcmp(language,'eng')
    ylabel('Electrons per second');
    xlabel('Energy');
else
    ylabel('Elektronen pro Sekunde');
    xlabel('Energie');
end
end
