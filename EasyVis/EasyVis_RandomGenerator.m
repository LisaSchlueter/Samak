savedir = [getenv('SamakPath'),'EasyVis/results/'];
savename = sprintf('%sComputeIntSpec_Chi2.mat',savedir);

if exist(savename,'file')
    load(savename);
     fprintf('load %s\n',savename)
else
    fprintf(2,'file not found %s \n Run Compute_IntSpec_FitPar.m\n',savename)
    return
end


Chi2Idx = flipud(find(chi2min<1e6));


%% 
close all
FontSize = 100;

filename = [plotdir,'FitPar_Chi2.gif'];
system(sprintf('rm %s',filename));
figure('Units','normalized','Position',[0.1,0.1,0.9,0.4]);
MakePrety1(1);


gif(filename,'frame',gcf,'DelayTime',0.3,'LoopCount',1)


for i=1:10

subplot(1,4,1);
r = rectangle('Curvature',0.1,'LineWidth',6,'EdgeColor',rgb('Red'),'FaceColor','none');
MakePrety1(1)
t(1) = text(0.5,0.5,sprintf('%.1f',mNuSq(Chi2Idx(i))),'HorizontalAlignment','center','FontSize',FontSize,'FontWeight','bold','Color',r.EdgeColor);

subplot(1,4,2);
r = rectangle('Curvature',0.1,'LineWidth',6,'EdgeColor',rgb('Gold'),'FaceColor','none');
MakePrety1(2)
t(2) =text(0.5,0.5,sprintf('%.2f',1e-03.*(A.ModelObj.Q_i+E0(Chi2Idx(i)))),'HorizontalAlignment','center','FontSize',FontSize,'FontWeight','bold','Color',r.EdgeColor);

subplot(1,4,3);
r = rectangle('Curvature',0.1,'LineWidth',6,'EdgeColor',rgb('LimeGreen'),'FaceColor','none');
MakePrety1(3)
t(3) =text(0.5,0.5,sprintf('%.0f',1e3.*(B(Chi2Idx(i))+A.ModelObj.BKG_RateSec_i)),'HorizontalAlignment','center','FontSize',FontSize,'FontWeight','bold','Color',r.EdgeColor);


subplot(1,4,4);
r = rectangle('Curvature',0.1,'LineWidth',6,'EdgeColor',rgb('DeepSkyBlue'),'FaceColor','none');
MakePrety1(4)
t(4) = text(0.5,0.5,sprintf('%.1f',1+N(Chi2Idx(i))),'HorizontalAlignment','center','FontSize',FontSize,'FontWeight','bold','Color',r.EdgeColor);

gif

t.delete;
end

%%
function MakePrety1(sb)
xticks([]);
yticks([]);
ax  = gca;
ax.YAxis.Color = 'w';
ax.XAxis.Color = 'w';

if sb==2 %subplot number
    ax.Position(1)  = 0.3;
elseif sb==3
    ax.Position(1) = 0.47;
elseif sb==4
    ax.Position(1) = 0.64;
end
end
