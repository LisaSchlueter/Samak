% illustrate influence of fit parameter on integral spectrum
% based on knm2 data

savedir = [getenv('SamakPath'),'EasyVis/results/'];
savename = sprintf('%sComputeIntSpec_FitPar_100.mat',savedir);

if exist(savename,'file')
    load(savename);
     fprintf('load %s\n',savename)
else
    fprintf(2,'file not found %s \n Run Compute_IntSpec_FitPar.m\n',savename)
    return
end

% re-sort: 
nIteration = size(TBDIS_N,2);
TBDIS_N = TBDIS_N(:,[nIteration/2:nIteration,nIteration/2-1:-1:1]);
TBDIS_B = TBDIS_B(:,[nIteration/2:nIteration,nIteration/2-1:-1:1]);
TBDIS_E0 = TBDIS_E0(:,[nIteration/2:nIteration,nIteration/2-1:-1:1]);
TBDIS_mNuSq = TBDIS_mNuSq(:,[nIteration/2:nIteration,nIteration/2-1:-1:1]);
%%
Format = 'video';%

Mode = '';%half';%

color_background = 'w';

if strcmp(color_background,'k')
    color_text = 'w';
else
    color_text = 'k';
end


%% gif1:  signal normalization 
plotdir = [getenv('SamakPath'),'EasyVis/plots/'];

f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,07]);

MakePretty1(color_background);

filename = [plotdir,'FitPar_N_',color_background];
if strcmp(Mode,'half')
  filename = [filename,'_half'];
  nLoop = nIteration/2;
else
     nLoop = nIteration;
end

if strcmp(Format,'gif')
    filename = [filename,'.gif'];
    system(sprintf('rm %s',filename));
    gif(filename,'frame',gcf,'DelayTime',1e-10);%,'LoopCount',1)
elseif strcmp(Format,'video')
    video = VideoWriter(filename,'MPEG-4');
    video.FrameRate = 20;
    open(video);
end

% first frame
p = plot(1e-03.*qU,TBDIS_bf./Time,'-','LineWidth',6,'Color',rgb('Black'));
MakePretty1(color_background);
ax = gca;
if strcmp(color_background,'k')
    ax.YAxis.Color = color_text;
    ax.XAxis.Color = color_text;
end

if strcmp(Format,'gif')
    gif
elseif strcmp(Format,'video')
    F = getframe(gcf);
    writeVideo(video,F);
end

% moving frames
for i=1:nLoop
    plot(1e-03.*qU,TBDIS_bf./Time,'-','LineWidth',6,'Color',color_text);
    hold on
    p = plot(1e-03.*qU,TBDIS_N(:,i)./Time,'-.','LineWidth',6,'Color',rgb('Orange'));
    
    hold off
    MakePretty1(color_background);
    
    FontArg = {'FontSize',50,'Color',p.Color,'FontWeight','bold'};
    if i>nIteration/2+1
        Scale = 0.74;
        t =text(1e-03.*qU(5),TBDIS_N(5,i)./Time(5).*Scale,sprintf('\\downarrow'),FontArg{:});
        text(1e-03.*qU(3),TBDIS_N(3,i)./Time(3).*Scale,sprintf('\\downarrow'),FontArg{:});
        text(1e-03.*qU(7),TBDIS_N(7,i)./Time(7).*Scale,sprintf('\\downarrow'),FontArg{:});
        text(1e-03.*qU(9),TBDIS_N(9,i)./Time(9).*(Scale-0.01),sprintf('\\downarrow'),FontArg{:});
        text(1e-03.*qU(11)-0.0004,TBDIS_N(11,i)./Time(11).*(Scale+0.02),sprintf('\\downarrow'),FontArg{:});
        
        FontSize = 55-0.7*(i-50);% i = [50,100] von groß nach klein:60->20
        text(18.572,15,sprintf(' Signalstärke \n nimmt ab'),'FontSize',FontSize,'Color',p.Color,'HorizontalAlignment','center');
  
    else
        Scale = 1.3;
        t =text(1e-03.*qU(5),TBDIS_N(5,i)./Time(5).*Scale,sprintf('\\uparrow'),FontArg{:});
        text(1e-03.*qU(3),TBDIS_N(3,i)./Time(3).*Scale,sprintf('\\uparrow'),FontArg{:});
        text(1e-03.*qU(7),TBDIS_N(7,i)./Time(7).*(Scale-0.02),sprintf('\\uparrow'),FontArg{:});
        text(1e-03.*qU(9),TBDIS_N(9,i)./Time(9).*(Scale-0.05),sprintf('\\uparrow'),FontArg{:});
        text(1e-03.*qU(11)-0.00075,TBDIS_N(11,i)./Time(11).*(Scale+0.11),sprintf('\\uparrow'),FontArg{:});
        FontSize = 20+0.7*(i);% i = [1,50] von klein nach groß:20->60
        text(18.572,15,sprintf(' Signalstärke \n nimmt zu'),'FontSize',FontSize,'Color',p.Color,'HorizontalAlignment','center');
        
    end
    
    if strcmp(color_background,'k')
        ax.YAxis.Color = color_text;
        ax.XAxis.Color = color_text;
    end

    if i==nIteration/2 || i==nIteration
        nFrames = 40;
    else
        nFrames = 1;
    end
    
    for f=1:nFrames
        if strcmp(Format,'gif')
            gif
        elseif strcmp(Format,'video')
            F = getframe(gcf);
            writeVideo(video,F);
        end
    end
end

if  ~strcmp(Mode,'half')  
    %% last frame
    Scale = 0.74;
    plot(1e-03.*qU,TBDIS_bf./Time,'-','LineWidth',6,'Color',rgb('Black'));
    MakePretty1(color_background);
    text(1e-03.*qU(3),1.22.*TBDIS_bf(3)./Time(3).*Scale,sprintf('$$\\big\\updownarrow$$'),'Interpreter','latex',FontArg{:});
    text(1e-03.*qU(5),1.22.*TBDIS_bf(5)./Time(5).*Scale,sprintf('$$\\big\\updownarrow$$'),'Interpreter','latex',FontArg{:});
    text(1e-03.*qU(7),1.2.*TBDIS_bf(7)./Time(7).*Scale,sprintf('$$\\big\\updownarrow$$'),'Interpreter','latex',FontArg{:});
    text(1e-03.*qU(9),1.18.*TBDIS_bf(9)./Time(9).*Scale,sprintf('$$\\big\\updownarrow$$'),'Interpreter','latex',FontArg{:});
    text(1e-03.*qU(11),1.15.*TBDIS_bf(11)./Time(11).*Scale,sprintf('$$\\big\\updownarrow$$'),'Interpreter','latex',FontArg{:});
    text(18.572,15,sprintf(' Signalstärke'),'FontSize',40,'Color',rgb('Orange'),'HorizontalAlignment','center');
end

 if strcmp(Format,'gif')
            gif
        elseif strcmp(Format,'video')
            F = getframe(gcf);
            writeVideo(video,F);
 end
        
 % save
  if strcmp(Format,'video')
      close(video);
  end
%% gif 2 background

f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,07]);

MakePretty1(color_background);
filename = [plotdir,'FitPar_B_',color_background];
if strcmp(Mode,'half')
  filename = [filename,'_half'];
  nLoop = nIteration/2;
else
     nLoop = nIteration;
end

if strcmp(Format,'gif')
    filename = [filename,'.gif'];
    system(sprintf('rm %s',filename));
    gif(filename,'frame',gcf,'DelayTime',1e-10);%,'LoopCount',1)
elseif strcmp(Format,'video')
    video = VideoWriter(filename,'MPEG-4');
    video.FrameRate = 20;
    open(video);
end

% first frame
p = plot(1e-03.*qU,TBDIS_bf./Time,'-','LineWidth',6,'Color',rgb('Black'));
MakePretty1(color_background);
ax = gca;
if strcmp(color_background,'k')
    ax.YAxis.Color = color_text;
    ax.XAxis.Color = color_text;
end

if strcmp(Format,'gif')
    gif
elseif strcmp(Format,'video')
    F = getframe(gcf);
    writeVideo(video,F);
end

for i=1:nLoop
    plot(1e-03.*qU,TBDIS_bf./Time,'-','LineWidth',6,'Color',color_text);
    hold on;
    p = plot(1e-03.*qU,TBDIS_B(:,i)./Time,'-.','LineWidth',6,'Color',rgb('LimeGreen'));
    hold off
    MakePretty1(color_background);
    
    FontArg = {'FontSize',50,'Color',p.Color,'FontWeight','bold'};
    if i>nIteration/2+1
        Scale = 0.8;
        t =text(1e-03.*18571,TBDIS_B(end-6,i)./Time(end-6).*Scale,sprintf('\\downarrow'),FontArg{:});
        text(1e-03.*18576,TBDIS_B(end-5,i)./Time(end-5).*Scale,sprintf('\\downarrow'),FontArg{:});
        text(1e-03.*18581,TBDIS_B(end-4,i)./Time(end-4).*Scale,sprintf('\\downarrow'),FontArg{:});
        FontSize = 55-0.7*(i-50);% i = [50,100] von groß nach klein:60->20
        
        text(18.552,0.4,sprintf(' Untergrund \n nimmt ab'),'FontSize',FontSize,'Color',p.Color,'HorizontalAlignment','center');
        
    else
        Scale = 1.45;
        t =text(1e-03.*18571,TBDIS_B(end-6,i)./Time(end-6).*Scale,sprintf('\\uparrow'),FontArg{:});
        text(1e-03.*18576,TBDIS_B(end-5,i)./Time(end-5).*Scale,sprintf('\\uparrow'),FontArg{:});
        text(1e-03.*18581,TBDIS_B(end-4,i)./Time(end-4).*Scale,sprintf('\\uparrow'),FontArg{:});
        FontSize = 20+0.7*(i);% i = [1,50] von klein nach groß:20->60
        text(18.552,0.4,sprintf(' Untergrund \n nimmt zu'),'FontSize',FontSize,'Color',p.Color,'HorizontalAlignment','center');
    end
    
    if strcmp(color_background,'k')
        ax.YAxis.Color = color_text;
        ax.XAxis.Color = color_text;
    end
    
    if i==nIteration/2 || i==nIteration
        nFrames = 40;
    else
        nFrames = 1;
    end
    
    for f=1:nFrames
        if strcmp(Format,'gif')
            gif
        elseif strcmp(Format,'video')
            F = getframe(gcf);
            writeVideo(video,F);
        end
    end
end

if  ~strcmp(Mode,'half')
    % last frame
    plot(1e-03.*qU,TBDIS_bf./Time,'-','LineWidth',6,'Color',color_text);
    MakePretty1(color_background);
    text(1e-03.*18571,1.22.*TBDIS_bf(end-6)./Time(end-6).*Scale,sprintf('$$\\big\\updownarrow$$'),'Interpreter','latex',FontArg{:});
    text(1e-03.*18576,1.22.*TBDIS_bf(end-5)./Time(end-5).*Scale,sprintf('$$\\big\\updownarrow$$'),'Interpreter','latex',FontArg{:});
    text(1e-03.*18581,1.2.*TBDIS_bf(end-4)./Time(end-4).*Scale,sprintf('$$\\big\\updownarrow$$'),'Interpreter','latex',FontArg{:});
    text(18.552,0.4,sprintf('Untergrund'),'FontSize',40,'Color',rgb('LimeGreen'),'HorizontalAlignment','center');
end

if strcmp(Format,'gif')
    gif
elseif strcmp(Format,'video')
    F = getframe(gcf);
    writeVideo(video,F);
end

% save
if strcmp(Format,'video')
    close(video);
end

%% gif 3 E0 

f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,07]);

filename = [plotdir,'FitPar_E0_',color_background];
if strcmp(Mode,'half')
  filename = [filename,'_half'];
  nLoop = nIteration/2;
else
     nLoop = nIteration;
end

if strcmp(Format,'gif')
    filename = [filename,'_E0.gif'];
    system(sprintf('rm %s',filename));
    gif(filename,'frame',gcf,'DelayTime',1e-10);%,'LoopCount',1)
elseif strcmp(Format,'video')
    video = VideoWriter(filename,'MPEG-4');
    video.FrameRate = 20;
    open(video);
end

% first frame
p = plot(1e-03.*qU,TBDIS_bf./Time,'-','LineWidth',6,'Color',color_text);
MakePretty1(color_background);
ax = gca;
if strcmp(color_background,'k')
    ax.YAxis.Color = color_text;
    ax.XAxis.Color = color_text;
end

if strcmp(Format,'gif')
    gif
elseif strcmp(Format,'video')
    F = getframe(gcf);
    writeVideo(video,F);
end


for i=1:nLoop
     plot(1e-03.*qU,TBDIS_bf./Time,'-','LineWidth',6,'Color',color_text);
     hold on
    p = plot(1e-03.*qU,TBDIS_E0(:,i)./Time,'-.','LineWidth',6,'Color',rgb('Red'));
   % hold on
   
    hold off
    MakePretty1(color_background);
    
    FontArg = {'FontSize',50,'Color',p.Color,'FontWeight','bold'};
    Idx = [10,13,16];
    
    CountsArrow = logspace(log10(0.4),log10(5),4);%[0.4,,1.8,5,10];
    qUArrow = interp1(TBDIS_E0(:,i)./Time,qU,CountsArrow,'lin');
    
    if i>nIteration/2+1
        for  j=1:numel(CountsArrow)
             t =text(1e-03.*qUArrow(j)-0.0031,CountsArrow(j),sprintf('\\leftarrow'),FontArg{:});
        end
        FontSize = 55-0.7*(i-50);% i = [50,100] von groß nach klein:60->20
       
         text(18.572,17,sprintf(' Endpunkt \n nimmt ab'),'FontSize',FontSize,'Color',p.Color,'HorizontalAlignment','center');
    
    else
        for  j=1:numel(CountsArrow)
            t =text(1e-03.*qUArrow(j)+0.0001,CountsArrow(j),sprintf('\\rightarrow'),FontArg{:});
        end
          FontSize = 20+0.7*(i);% i = [1,50] von klein nach groß:20->60
          text(18.572,17,sprintf(' Endpunkt \n nimmt zu'),'FontSize',FontSize,'Color',p.Color,'HorizontalAlignment','center');
    end
    
    if strcmp(color_background,'k')
        ax.YAxis.Color = color_text;
        ax.XAxis.Color = color_text;
    end
    
    if i==nIteration/2 || i==nIteration
        nFrames = 40;
    else
        nFrames = 1;
    end
    
    for f=1:nFrames
        if strcmp(Format,'gif')
            gif
        elseif strcmp(Format,'video')
            F = getframe(gcf);
            writeVideo(video,F);
        end
    end
    
end

if ~strcmp(Mode,'half')
    % last frame
    plot(1e-03.*qU,TBDIS_bf./Time,'-','LineWidth',6,'Color',color_text);
    MakePretty1(color_background);
    qUArrow = interp1(TBDIS_bf./Time,qU,CountsArrow,'lin');
    
    for  j=1:numel(CountsArrow)
        t =text(1e-03.*qUArrow(j)-0.003,CountsArrow(j),sprintf('$$\\longleftrightarrow$$'),'Interpreter','latex',FontArg{:});
    end
    text(18.572,17,sprintf('Endpunkt'),'FontSize',40,'Color',rgb('Red'),'HorizontalAlignment','center');
end

if strcmp(Format,'gif')
    gif
elseif strcmp(Format,'video')
    F = getframe(gcf);
    writeVideo(video,F);
end

% save
if strcmp(Format,'video')
    close(video);
end
%% gif 4 numass

f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,07]);

filename = [plotdir,'FitPar_mNuSq_',color_background];
if strcmp(Mode,'half')
  filename = [filename,'_half'];
  nLoop = nIteration/2;
else
     nLoop = nIteration;
end

if strcmp(Format,'gif')
    filename = [filename,'_mNuSq.gif'];
    system(sprintf('rm %s',filename));
    gif(filename,'frame',gcf,'DelayTime',1e-10);%,'LoopCount',1)
elseif strcmp(Format,'video')
    video = VideoWriter(filename,'MPEG-4');
    video.FrameRate = 20;
    open(video);
end

% first frame
p = plot(1e-03.*qU,TBDIS_bf./Time,'-','LineWidth',6,'Color',color_text);
MakePretty1(color_background);
ax = gca;
if strcmp(color_background,'k')
    ax.YAxis.Color = color_text;
    ax.XAxis.Color = color_text;
end

if strcmp(Format,'gif')
    gif
elseif strcmp(Format,'video')
    F = getframe(gcf);
    writeVideo(video,F);
end


for i=1:nLoop
    plot(1e-03.*qU,TBDIS_bf./Time,'-','LineWidth',6,'Color',color_text);
    hold on
    p = plot(1e-03.*qU,TBDIS_mNuSq(:,i)./Time,'-.','LineWidth',6,'Color',rgb('DeepSkyBlue'));
    hold off
    MakePretty1(color_background);
    
    FontArg = {'FontSize',50,'Color',p.Color,'FontWeight','bold'};
    Idx = [10,13,16];
    
    CountsArrow = logspace(log10(0.4),log10(5),4);%[0.4,,1.8,5,10];
    qUArrow = interp1(TBDIS_mNuSq(:,i)./Time,qU,CountsArrow,'lin');
    
    if i>nIteration/2+1
        for  j=1:numel(CountsArrow)
            t =text(1e-03.*qUArrow(j)+0.0001,(1.11-0.0006.*i).*CountsArrow(j),sprintf('$$\\nearrow$$'),'Interpreter','latex',FontArg{:});
           % if CountsArrow(j)<5
                t.Rotation = 0.3.*i;
           % end
        end
         FontSize = 55-0.7*(i-50);% i = [50,100] von groß nach klein:60->20
         text(18.573,4,sprintf(' Neutrinomasse \n nimmt ab'),'FontSize',FontSize,'Color',p.Color,'HorizontalAlignment','center');
  
    else
        for  j=1:numel(CountsArrow) 
            if CountsArrow(j)>4% from (y-axis)  bottom to top
                t =text(1e-03.*qUArrow(j)-(0.0033-i*0.00002),0.77.*CountsArrow(j),sprintf('$$\\swarrow$$'),'Interpreter','latex',FontArg{:});
            else
                t =text(1e-03.*qUArrow(j)-(0.0032-i.*0.00001),0.77.*CountsArrow(j),sprintf('$$\\swarrow$$'),'Interpreter','latex',FontArg{:}); 
            end
             t.Rotation = -0.5.*i;
        end
        FontSize = 20+0.7*(i);% i = [1,50] von klein nach groß:20->60
        text(18.573,4,sprintf(' Neutrinomasse \n nimmt zu'),'FontSize',FontSize,'Color',p.Color,'HorizontalAlignment','center');
        
    end
    
    if strcmp(color_background,'k')
        ax.YAxis.Color = color_text;
        ax.XAxis.Color = color_text;
    end
    
    if i==nIteration/2 || i==nIteration
        nFrames = 40;
    else
        nFrames = 1;
    end
    
    for f=1:nFrames
        if strcmp(Format,'gif')
            gif
        elseif strcmp(Format,'video')
            F = getframe(gcf);
            writeVideo(video,F);
        end
    end

end

if  ~strcmp(Mode,'half')
    % last frame
    plot(1e-03.*qU,TBDIS_bf./Time,'-','LineWidth',6,'Color',color_text);
       if strcmp(color_background,'k')
        ax.YAxis.Color = color_text;
        ax.XAxis.Color = color_text;
    end
    MakePretty1(color_background);
    qUArrow = interp1(TBDIS_bf./Time,qU,CountsArrow,'lin');
    
    for  j=1:numel(CountsArrow)
        t =text(1e-03.*qUArrow(j)-0.0015,CountsArrow(j),sprintf('$$\\nearrow$$'),'Interpreter','latex',FontArg{:});
        t =text(1e-03.*qUArrow(j)-0.0015,CountsArrow(j),sprintf('$$\\swarrow$$'),'Interpreter','latex',FontArg{:});
    end
    text(18.573,4,sprintf('Neutrinomasse'),'FontSize',40,'Color',rgb('DeepSkyBlue'),'HorizontalAlignment','center');
end
% save
for f=1:80
    if strcmp(Format,'gif')
        gif
    elseif strcmp(Format,'video')
        
        
        F = getframe(gcf);
        writeVideo(video,F);
    end
end
% save
if strcmp(Format,'video')
    close(video);
end

%%
function MakePretty1(color)
xlabel('Energie (keV)');
ylabel('Elektronen pro Sekunde');
PrettyFigureFormat('FontSize',30);
set(gca,'YScale','log');
%set(gca,'FontSize',28);
grid off;
ylim([0.1 100]);%5500
xlim(1e-03.*[18574-39.8 18575+10])
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
yticks([1 10 100])
yticklabels({'1','10','100'});
set(gca,'LineWidth',3);
xticks([18.54 18.56 18.58])
box off
set(gcf,'Color',color);
set(gca,'Color',color);
end
