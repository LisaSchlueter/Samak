% plot: points of original grid  and interpolation grid
nGridSteps = 30;
range = 40;
ExtmNu4Sq = 'log';
mNu4SqTestGrid = 'OFF';%5;
nInter = 1e3; 
SavePlt = 'ON';

%% define grid
if strcmp(ExtmNu4Sq,'ON')
    mnu4Sq_ex = [0.1;0.35;0.7];
    nGridSteps_i = nGridSteps;
    nGridSteps = nGridSteps_i-3;
end

if mNu4SqTestGrid==1
    mnu4SqSmall  = logspace(0,log10((range-11)^2),nGridSteps-5)';
    mnu4SqLarge  = logspace(log10((range-10)^2),log10((range)^2),5)';
    mNu4Sq       = sort([mnu4SqSmall;mnu4SqLarge]);
elseif mNu4SqTestGrid==2
    mnu4SqSmall  = logspace(0,log10((range-11)^2),nGridSteps-7)';
    mnu4SqLarge  = logspace(log10((range-10)^2),log10((range)^2),7)';
    mNu4Sq       = sort([mnu4SqSmall;mnu4SqLarge]);
elseif mNu4SqTestGrid==3
    mnu4SqSmall  = logspace(0,log10((range-11)^2),nGridSteps-10)';
    mnu4SqLarge  = logspace(log10((range-10)^2),log10((range)^2),10)';
    mNu4Sq       = sort([mnu4SqSmall;mnu4SqLarge]);
elseif mNu4SqTestGrid==4
    mnu4SqSmall  = logspace(0,log10((range-11)^2),nGridSteps-5)';
    mnu4SqLarge  = logspace(log10((range-10)^2),log10((range)^2),5)';
    mNu4Sq       = sort([mnu4SqSmall;mnu4SqLarge]);
elseif any(mNu4SqTestGrid==5) || any(mNu4SqTestGrid==5.5)
    mnu4SqSmall  = logspace(0,log10((range-11)^2),nGridSteps-15)';
    mnu4SqLarge  = logspace(log10((range-10)^2),log10((range)^2),15)';
    mNu4Sq       = sort([mnu4SqSmall;mnu4SqLarge]);
elseif strcmp(ExtmNu4Sq,'log')
     mNu4Sq      = logspace(-1,log10((range)^2),nGridSteps)';
else
    mNu4Sq      = logspace(0,log10((range)^2),nGridSteps)';
end

if ismember(ExtmNu4Sq,{'ON','0.01'})
    mNu4Sq = [mnu4Sq_ex;mNu4Sq];%;logspace(0,log10((obj.range)^2),obj.nGridSteps-3)'];
    nGridSteps = nGridSteps_i;
end
sin2T4      = logspace(-3,log10(0.5),nGridSteps);
mNu4Sq      = repmat(mNu4Sq,1,nGridSteps);

% interpolation grid
mNu4_tmp = logspace(log10(min(mNu4Sq(:,1))),log10(max(mNu4Sq(:,1))),nInter);
mNu4Sq_inter = repmat(mNu4_tmp,nInter,1);
sin2T4_inter = logspace(log10(min(min(sin2T4))),log10(max(max(sin2T4))),nInter);

%% plot
GetFigure;
pInter = plot(sin2T4_inter,mNu4Sq_inter,'.','Color',rgb('Silver'),'MarkerSize',1);
hold on;
pOrig = plot(sin2T4,mNu4Sq','k.');
set(gca,'XScale','log');
set(gca,'YScale','log');
xlim([7e-04 0.7]);
ylim([0.05,100^2]);
hold off;
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'))
PrettyFigureFormat;

leg = legend([pOrig(1),pInter(1)],sprintf('Original grid:  %.0f \\times %.0f',nGridSteps,nGridSteps),...
    sprintf('Interpolation grid: %.0f \\times %.0f',nInter,nInter),'Location','northwest' );
PrettyLegendFormat(leg);

% save plot
if strcmp(SavePlt,'ON')
    pltdir = sprintf('%sksn2ana/ksn2_GridSize/plots/',getenv('SamakPath'));
    MakeDir(pltdir);
    pltname = [pltdir,'ksn2_GridPoints.png'];
    print(gcf,pltname,'-dpng','-r350');
    fprintf('save plot to %s \n',pltname);
end