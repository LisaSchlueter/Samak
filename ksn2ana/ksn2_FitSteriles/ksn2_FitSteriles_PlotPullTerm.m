% illustrate pull-term
pullFlag = 27;
range = 40;
mNu4Sq = linspace(-5,range^2+50,1e3);
sin2t4 = linspace(-0.5,1.5,1e3);

if pullFlag == 27
    sinU = 0.5;     % sin2t4 upper bound
    sinD = 0;       % sin2t4 lower bound
    m4U   = range^2; % mSq4 upper bound
    m4D   = 0;       % mSq4 lower bound
elseif pullFlag == 28
    sinU = 1;     % sin2t4 upper bound
    sinD = 0;       % sin2t4 lower bound
    m4U   = range^2; % mSq4 upper bound
    m4D   = 0;       % mSq4 lower bound
end

PullTermSin =  exp(1e3*(sin2t4-sinU-0.01)) + ...    % sin2t4 upper bound
               exp(-1e3*(sin2t4-sinD+0.01)) ;     % sin2t4 lower bound
    
PullTermM4 =  exp(-1e3*(mNu4Sq-m4D+0.01)) + ... % mSq4 lower bound
              exp(1e3*(mNu4Sq-m4U-0.01));       % mSq4 upper bound
          
t = tiledlayout(1,1);
ax1 = axes(t);
plot(ax1,mNu4Sq,PullTermM4,'LineWidth',2,'Color',rgb('FireBrick'));
ax1.XColor = rgb('FireBrick');
 ax1.YLim = [-0.2 5];
 xlim([-100 range^2+100]);
PrettyFigureFormat('FontSize',22)
 ax1.Box = 'off';
xlabel(sprintf('{\\itm}_4^2 (eV^2)'));
ylabel(sprintf('pull term (\\chi^2)'));
ax2 = axes(t);
plot(ax2,sin2t4,PullTermSin,'-.','LineWidth',2,'Color',rgb('DodgerBlue'))
PrettyFigureFormat('FontSize',22)
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax2.Box = 'off';
ax2.YTickLabel = '';
ax2.XColor = rgb('DodgerBlue');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
 ax2.YLim = [-0.2 5];
xlim([-0.1 sinU+0.1]);
pltdir = sprintf('%sksn2ana/ksn2_FitSteriles/plots/',getenv('SamakPath'));
MakeDir(pltdir);
pltname = sprintf('%sksn2_FitSteriles_PlotPullTerm.png',pltdir);
print(gcf,pltname,'-dpng','-r350');

%%