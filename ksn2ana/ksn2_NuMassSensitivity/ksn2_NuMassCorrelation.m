% investigate the correlation between nu-mass and m4
DataType = 'Real';
InterpMode = 'spline'; 
%% load/compute file
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];
savefile = [savedir,'ksn2_NuMassCorrelation',DataType,'_',InterpMode,'.mat'];

if exist(savefile,'file') 
    load(savefile);
else
    %% settings that might change
    chi2 = 'chi2CMShape'; 
    nGridSteps = 40;
    range = 40;
    
    %% configure RunAnalysis object
    if strcmp(chi2,'chi2Stat')
        NonPoissonScaleFactor = 1;
    elseif  strcmp(chi2,'chi2CMShape')
        NonPoissonScaleFactor = 1.112;
    end
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar','mNu E0 Norm Bkg',...%free par
        'SysBudget',40,...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'RadiativeFlag','ON',...
        'FSDFlag','KNM2_0p5eV',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'chi2',chi2,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'TwinBias_Q',18573.7,...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',3*1e-06,...
        'TwinBias_BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
    A = MultiRunAnalysis(RunAnaArg{:});
    %% configure Sterile analysis object
    SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'nGridSteps',nGridSteps,...
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',...
        'range',range,...
        'LoadGridArg',{'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'}};
    %%
    S = SterileAnalysis(SterileArg{:});
    S.InterpMode = InterpMode;
    S.LoadGridFile(S.LoadGridArg{:});
     mNu4Sq_i = S.mNu4Sq;
    sin2T4_i = S.sin2T4;
    chi2_i   = S.chi2;
    mNuSq_i = S.mNuSq;
    
    S.Interp1Grid;
    S.ContourPlot('BestFit','ON'); close;
    %%
    mNu4Sq = S.mNu4Sq;
    sin2T4 = S.sin2T4;
    chi2   = S.chi2;
    mNuSq = S.mNuSq;
    
    sin2T4_contour  =  S.sin2T4_contour;
    mNu4Sq_contour = S.mNu4Sq_contour;
    sin2T4_contour_bf  =  S.sin2T4_bf;
    mNu4Sq_contour_bf = S.mNu4Sq_bf;
    save(savefile,'mNu4Sq','sin2T4','chi2','mNuSq',...
        'mNu4Sq_i','sin2T4_i','chi2_i','mNuSq_i',...
        'sin2T4_contour','mNu4Sq_contour','sin2T4_bf','mNu4Sq_contour_bf');
end
% surf(sin2T4,mNu4Sq,mNuSq,'EdgeColor','none');
% set(gca,'XScale','log');
% set(gca,'YScale','log');
% view(2)


%%
myColors = jet(7);
legStr  = cell(size(myColors,1),1);
pHandle = cell(size(myColors,1),1);
InterMode = 'Inter';
switch InterMode
    case 'Orig'
        x = mNuSq_i;
        y = mNu4Sq_i(:,1);
        z = sin2T4_i(1,:);
    case 'Inter'
        x = mNuSq;
        y = mNu4Sq(1,:);
        z = sin2T4(:,1);
end
PltMode = 'Gradient';

GetFigure;
for i=1:numel(z)
    if strcmp(PltMode,'Absolut')
        p = plot(x(i,:),y,'LineWidth',2);
    elseif strcmp(PltMode,'Gradient')
        p = plot(diffxy(y,x(i,:)),y,'LineWidth',2);
    end
    
    hold on;
    if z(i)<=1e-02 % 0.001 to 0.01
        p.Color = myColors(1,:);
        legStr{1} = sprintf('0.001 \\leq |{\\itU}_{e4}|^2 \\leq 0.01');
        pHandle{1} = p;
    elseif z(i)<=5e-02 %  0.01 to 0.05
        p.Color = myColors(2,:);
        legStr{2} = sprintf('0.01  < |{\\itU}_{e4}|^2 \\leq 0.05');
        pHandle{2} = p;
    elseif z(i)<=0.1 % 0.05 to 0.1
        p.Color = myColors(3,:);
        legStr{3} = sprintf('0.05  < |{\\itU}_{e4}|^2 \\leq 0.1');
        pHandle{3} = p;
    elseif z(i)<=0.2 % 0.1 to 0.25
        p.Color = myColors(4,:);
        legStr{4} = sprintf('0.1   < |{\\itU}_{e4}|^2 \\leq 0.2');
        pHandle{4} = p;
    elseif z(i)<=0.3
        p.Color = myColors(5,:);
        legStr{5} = sprintf('0.2   < |{\\itU}_{e4}|^2 \\leq 0.3');
        pHandle{5} = p;
    elseif z(i)<=0.4
        p.Color = myColors(6,:);
        legStr{6} = sprintf('0.3   < |{\\itU}_{e4}|^2 \\leq 0.4');
        pHandle{6} = p;
    else
        p.Color = myColors(7,:);
        legStr{7} = sprintf('0.4   < |{\\itU}_{e4}|^2 \\leq 0.5');
        pHandle{7} = p;
    end
end
set(gca,'YScale','log');

if strcmp(PltMode,'Absolut')
    xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)')); 
elseif strcmp(PltMode,'Gradient')
    xlabel(sprintf('\\delta{\\itm}_\\nu^2/\\delta{\\itm}_4^2'));
    xlim([-3 3])
end

ylim([0.1,1600]) 
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));

leg = legend([pHandle{:}],legStr,'Location','southeast');
PrettyLegendFormat(leg);
leg.Title.String = 'Active-to-sterile mixing';
leg.Title.FontWeight = 'normal';
PrettyFigureFormat('FontSize',20);

leg.Title.FontSize = get(gca,'FontSize');

pltdir = strrep(savedir,'results','plots');
pltname = [pltdir,sprintf('ksn2_NuMassCorrelation%s_%s.png',PltMode,DataType)];
MakeDir(pltdir)
print(gcf,pltname,'-dpng','-r350');

%% 
       
Corr = zeros(1e3,1e3);
for i=1:numel(z)
    if strcmp(PltMode,'Absolut')    
    elseif strcmp(PltMode,'Gradient')
        Corr(i,:) = diffxy(mNu4Sq(1,:),mNuSq(i,:));
    end
end
%%
GetFigure;
%Corr(abs(Corr)<0.1) = NaN;
s1 = surf(sin2T4,mNu4Sq,Corr,'EdgeColor','none');
view(2);
set(gca,'YScale','log')
set(gca,'XScale','log')
c = colorbar;
xlim([1e-03 0.5]);
ylim([0.1 1600]);
PrettyFigureFormat('FontSize',22);
grid off;
xlabel(sprintf('|{\\itU}_{e4}|^2'))
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
%%
% Idx = 500;
% GetFigure
% s1 =  subplot(2,1,1);
% plot(mNu4Sq(Idx,:),mNuSq(Idx,:),'LineWidth',2);
% PrettyFigureFormat('FontSize',22);
% ylabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
% xlabel(sprintf('{\\itm}_4^2 (eV^2)'));
% set(gca,'XScale','log');
% 
% s2 = subplot(2,1,2);
% plot(mNu4Sq(1,:),diffxy(mNu4Sq(Idx,:),mNuSq(Idx,:)),'LineWidth',2);
% PrettyFigureFormat('FontSize',22);
% xlabel(sprintf('{\\itm}_4^2 (eV^2)'));
% ylabel(sprintf('\\delta{\\itm}_\\nu^2/\\delta{\\itm}_4^2'));
% 
% linkaxes([s1,s2],'x');
% set(gca,'XScale','log');
% xlim([0.1 1600])
% sin2T4(Idx,1);