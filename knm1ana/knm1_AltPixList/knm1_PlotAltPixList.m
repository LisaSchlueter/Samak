
% plot alternative pixel lists

% get standard uniform results
chi2 = 'chi2Stat';
NP = 1.064;

fUniform = [getenv('SamakPath'),'knm1ana/knm1_AlternativeRunLists/results/',...
    sprintf('knm1_AltRunList_%s_NP%2g_KNM1.mat',chi2,NP)];

d = importdata(fUniform);
nPix_u = 117;
mNuSq_u = d.FitResult.par(1);
mNuSqErr_u =  d.FitResult.err(1);
chi2min_u =  d.FitResult.chi2min;
p_u = 1-chi2cdf(d.FitResult.chi2min,d.FitResult.dof);

AltPixLists = {'Half','AziHalfNS','AziHalfEW'};
PixList_Labels = {'Golden','Inner','Outer','North','South','East','West'};

nPix = zeros(numel(AltPixLists),2);
mNuSq = zeros(numel(AltPixLists),2);
mNuSqErr = zeros(numel(AltPixLists),2);
chi2min = zeros(numel(AltPixLists),2);
p     = zeros(numel(AltPixLists),2);
PixNumList = cell(2*numel(AltPixLists),1);
%% fit run lists
LocalFontSize = 18;

for i=1:numel(AltPixLists)
    AltPixList =AltPixLists{i};  %'Half';%'AziHalfEW';  % defines alternative pixel list
    
    % label
    savedir = [getenv('SamakPath'),'knm1ana/knm1_AltPixList/results/'];
    savename = sprintf('%sknm1_PixListAlt_%s_%s_%s_%.0feV_%s_%2g.mat',...
        savedir,AltPixList,'Real','mNuE0BkgNorm',40,chi2,NP);
    
    d = importdata(savename);
    nPix(i,:) = cellfun(@(x) numel(x),d.PixList,'UniformOutput',1);
    PixNumList(2*(i-1)+1:2*i) = d.PixList;
    mNuSq(i,:) = d.FitResult.par(:,1);
   % mNuSqErr(i,:) =  d.FitResult.err(:,1);
    mNuSqErr(i,:) =  0.5.*(d.FitResult.errPos(:,1)-d.FitResult.errNeg(:,1));
    chi2min(i,:) =  d.FitResult.chi2min;
    p(i,:) = 1-chi2cdf(d.FitResult.chi2min,d.FitResult.dof);
end

%% plot

RandList = 'ON';
nPix = [nPix_u;reshape(nPix',[2*numel(AltPixLists),1])];
x = [mNuSq_u; reshape(mNuSq',[2*numel(AltPixLists),1])];
xerr = [mNuSqErr_u;reshape(mNuSqErr',[2*numel(AltPixLists),1])];
pp =[p_u;reshape(p',[2*numel(AltPixLists),1])];

figureHandle = figure('Units','normalized','Position',[0.1,0.1,0.6,0.5]);
s1 = subplot(1,4,1:3);

if strcmp(RandList,'ON')
    nPixList = 2000;
    savenameRand = [savedir,sprintf('RandomHalfPixList_Unblinded_%s_NP%2g_%.0ffits_merged.mat',chi2,NP,nPixList)];
    frand = importdata(savenameRand);
    mNuSq_rand =  cell2mat(cellfun(@(x) x.par(1),frand.FitResults,'UniformOutput',0));
    [l,a_rand] = boundedline(mean(mNuSq_rand).*ones(10,1),linspace(-1,3.*numel(AltPixLists),10),std(mNuSq_rand).*ones(10,1),'orientation','horiz');
    hold on
    l.LineStyle = 'none'; l.Color = rgb('DimGray'); l.LineWidth = 2;
    a_rand.FaceColor = rgb('Silver');
    a_rand.FaceAlpha = 0.5;   
end

pref = plot(mNuSq_u.*ones(10,1),linspace(0,numel(x)+1,10),':','LineWidth',2,'Color',rgb('LimeGreen'));
hold on;
e_fits= errorbar(x,1:numel(x),xerr,'horizontal',...
    '.','MarkerSize',20,'LineStyle','none','Color',rgb('DodgerBlue'),'CapSize',0,'LineWidth',2);
e_bf = errorbar(x(1),1,xerr(1),'horizontal',...
    '.','MarkerSize',20,'LineStyle','none','Color',pref.Color,'CapSize',0,'LineWidth',2);


xlabel(sprintf('{\\itm}_\\nu^{ 2} (eV^{ 2})'));
set(gca,'YMinorTick','off');
yticks(1:numel(x));
yticklabels(PixList_Labels)
ax1 = gca;
xlim([-5.5 2.7]);
ylabel('Pixel selection');
PrettyFigureFormat('FontSize',LocalFontSize);

if strcmp(RandList,'ON')
    leg = legend(a_rand,sprintf('Random half pixels: 1\\sigma band (%.0f samples)',nPixList));
    PrettyLegendFormat(leg,'alpha',0.8); leg.EdgeColor = 'none';
    leg.Location= 'northwest';
    leg.FontSize = LocalFontSize-2;
end


s2 = subplot(1,4,4);
area(linspace(0,0.05,10),numel(x)+1.*ones(10,1),'FaceColor',rgb('Red'));%(0.05.*ones(10,1),linspace(0,numel(RunLists)+1,10),'-','MarkerSize',20,'LineWidth',2,'Color',rgb('Red'));
hold on;
plot(pp,1:numel(x),'.','MarkerSize',20,'LineWidth',2,'Color',rgb('DodgerBlue'));
plot(pp(1),1,'.','MarkerSize',20,'LineWidth',2,'Color',pref.Color);

yticklabels('');
xlabel(sprintf('{\\it p}'));
ax2 = gca;
PrettyFigureFormat('FontSize',LocalFontSize);

linkaxes([s1,s2],'y');

ax1.Position(2) = 0.17;
ax2.Position(2) = ax1.Position(2);
ax2.Position(4) = ax1.Position(4);
ax1.Position(1) = 0.2;
ax2.Position(1) = 0.79;

if strcmp(RandList,'ON')
    ylim([0.5 numel(x)+1])
else
    ylim([0.5 numel(x)+0.5]);
end
%%
pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);
pltname = [pltdir,'knm1_AltPixLists.pdf'];
export_fig(gcf,pltname);


%%

FPDView = 'OFF';
ViewList = 'AziHalfEW';%{'Half','AziHalfNS','AziHalfEW'};
if strcmp(FPDView,'ON')
    close all
    Pixels = NaN.*zeros(148,1);
    if strcmp(ViewList,'Half')
        Pixels(PixNumList{1}) = 1;
        Pixels(PixNumList{2}) = 0;
    elseif strcmp(ViewList,'AziHalfNS')
        Pixels(PixNumList{3}) = 1;
        Pixels(PixNumList{4}) = 0;
    elseif strcmp(ViewList,'AziHalfEW')
        Pixels(PixNumList{5}) = 1;
        Pixels(PixNumList{6}) = 0;
    end
    [plotHandle, cbHandle] = FPDViewer(Pixels,'Label','ON','ReDrawSkeleton','ON');
    cbHandle.delete;
    colormap(summer);
    colormap(spring);
     ax = gca;
     
    if strcmp(ViewList,'Half')
        xlabel('Inner/Outer pixels');      
       ax.XLabel.FontSize = 24;
        export_fig(gcf,[pltdir,'FPDView_InOut.pdf']);
    elseif strcmp(ViewList,'AziHalfNS') 
        xlabel('Northern/Southern pixels');
         ax.XLabel.FontSize = 24;
         export_fig(gcf,[pltdir,'FPDView_NS.pdf']);
    elseif strcmp(ViewList,'AziHalfEW')
        xlabel('Eastern/Western pixels');
         ax.XLabel.FontSize = 24;
         export_fig(gcf,[pltdir,'FPDView_EW.pdf']);
    end
    
   
end


%% signficiance
fprintf('Inner vs. Outer: %.2f sigma \n',abs(x(3)-x(2))./mean(xerr(2:3)));
fprintf('North vs. South: %.2f sigma \n',abs(x(5)-x(4))./mean(xerr(4:5)));
fprintf('East vs. West  : %.2f sigma \n',abs(x(7)-x(6))./mean(xerr(6:7)));

% sigma with respect to random half
Sigma = (x-mean(mNuSq_rand))./std(mNuSq_rand);
