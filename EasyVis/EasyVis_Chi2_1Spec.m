% illustrate chi2 minimization
% based on knm2 data and model
% 1 example spectrum - slow
savedir = [getenv('SamakPath'),'EasyVis/results/'];
savename = sprintf('%sComputeIntSpec_Chi2.mat',savedir);

if exist(savename,'file')
    load(savename);
    fprintf('load %s\n',savename)
else
    fprintf(2,'file not found %s \n Run Compute_IntSpec_FitPar.m\n',savename)
    return
end

Language = 'eng';

i=1;
Chi2Idx = flipud(find(chi2min<1e6));
Chi2SubRun = (TBDIS_data-TBDIS(:,Chi2Idx(i))).^2./TBDIS_data;

color_background = 'w';

if strcmp(color_background,'k')
    color_text = 'w';
else
    color_text = 'k';
end
%% gif 1: spectrum with (random) fit parameter values and matching chi^2
plotdir = [getenv('SamakPath'),'EasyVis/plots/'];
if strcmp(Language,'eng')
    video = VideoWriter([plotdir,'FitPar_Chi2_1Spec_',color_background,'_eng'],'MPEG-4');
else
video = VideoWriter([plotdir,'FitPar_Chi2_1Spec_',color_background],'MPEG-4');
end
video.FrameRate = 2;
open(video);

f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8 ,0.7]);
subplot(1,4,1:3)
plot(1e-03.*qU,TBDIS_data./Time,'.','MarkerSize',25,'Color',color_text);
MakePretty1(color_background,Language);
hold on,
p = plot(1e-03.*qU,TBDIS(:,Chi2Idx(i))./Time,'-','LineWidth',6,'Color',rgb('SlateGray'));
WriteVideo(video,5);

for j = 1:28
    subplot(1,4,1:3)

    k=j;
   % for k=1:j
     plot(1e-03.*qU,TBDIS_data./Time,'.','MarkerSize',25,'Color',color_text);
       hold on;
        if TBDIS_data(j)<TBDIS(k,Chi2Idx(i))
            plot(1e-03.*qU(k).*ones(1,100),linspace(TBDIS_data(k)./Time(k),TBDIS(k,Chi2Idx(i))./Time(k),100),':','LineWidth',5,'Color',rgb('Red'));
        else
            plot(1e-03.*qU(k).*ones(1,100),linspace(TBDIS(k,Chi2Idx(i))./Time(k),TBDIS_data(k)./Time(k),100),':','LineWidth',5,'Color',rgb('Red'));
        end
        %end
      %  hold on;
       % p = plot(1e-03.*qU,TBDIS(:,Chi2Idx(i))./Time,'-','LineWidth',6,'Color',rgb('SlateGray'));
        
        % hold off
        MakePretty1(color_background,Language);
    
    subplot(1,4,4)
    bar(sum(Chi2SubRun(1:j)),'FaceAlpha',1,'FaceColor',rgb('OrangeRed'),'EdgeColor','none','BarWidth',0.2);
    PrettyFigureFormat('FontSize',28)
    set(gca,'Color',color_background);
    set(gcf,'Color',color_background);
     ax2 = gca;
     xticks([]);
     yticks([]);
     ylim([0.1,sum(Chi2SubRun)]);
     if strcmp(Language,'eng')
         xlabel('Difference');
     else
         xlabel('Unterschied');
     end
    set(gca,'LineWidth',3);
   
    ax2.YAxis.Color = color_background;
     ax2.XAxis.Color = color_text;
    box off
    WriteVideo(video,1);
end

% save
 WriteVideo(video,20);
close(video);

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
if strcmp(color,'k')
    ax.YAxis.Color = 'w';
    ax.XAxis.Color = 'w';
end
end

function WriteVideo(video,nFrame)
    for i=1:nFrame
        F = getframe(gcf);
        writeVideo(video,F);
    end
end
