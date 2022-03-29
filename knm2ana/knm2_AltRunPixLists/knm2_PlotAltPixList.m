% plot alternative pixel list results
nRandFits = 2e3;
freePar = 'mNu E0 Bkg Norm';
DataType = 'Real';
Parameter = 'mNuSq';
RunList = 'KNM2_Prompt';
range = 40;                % fit range in eV below endpoint
AltPixLists = {'Half','AziHalfNS','AziHalfEW'};
PixList_Labels = {'Golden','Inner','Outer','North','South','East','West','Bullseye','1','2','3','4'};

savedir = [getenv('SamakPath'),'knm2ana/knm2_AltRunPixLists/results/'];

%% get uniform result
savedir_u = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savename_u = sprintf('%sknm2ubfinal_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s.mat',...
    savedir_u,3,'Real',40,'mNuE0BkgNorm','chi2Stat+','StackPixel','KNM2');
d_u = importdata(savename_u);
switch Parameter
    case 'mNuSq'
        Par_u = d_u.FitResult.par(1);
        Err_u = 0.5*(d_u.FitResult.errPos(1)-d_u.FitResult.errNeg(1));
    case 'E0'
        Par_u = d_u.FitResult.par(2);
        Err_u = d_u.FitResult.err(2);
end
nPix_u = numel(d_u.A.PixList);
p_u = 1-chi2cdf(d_u.FitResult.chi2min,d_u.FitResult.dof);

%% load random half
savenameRand = sprintf('%sknm2_PixListRandHalf_%s_%s_%.0feV_%.0ffits_KNM2_0p1eV_BkgPT3mucpsPers.mat',...
    savedir,DataType,strrep(freePar,' ',''),range,nRandFits);
if exist(savenameRand,'file')
    dRand = importdata(savenameRand);
else
    fprintf(2,'Random half pix list not computed yet \n')
    return
end
switch Parameter
    case 'mNuSq'
        x_rand     = cellfun(@(x) x.par(1),dRand.FitResults);
        x_rand     = x_rand-(median(x_rand)-d_u.FitResult.par(1));
        x_rand_err = 0.5.*(cellfun(@(x) x.errPos(1),dRand.FitResults)-cellfun(@(x) x.errNeg(1),dRand.FitResults));
        x_rand(x_rand_err<0.1) = []; % sometime uncertainty estimation failed
        x_rand_err(x_rand_err<0.1) = [];
       % x_rand_err(x_rand_err<0.1) = median(x_rand_err); % sometime uncertainty estimation failed
        xStr = sprintf('{\\itm}_\\nu^2');
        Unit = sprintf('eV^2');
    case 'E0'
        x_rand    = cellfun(@(x) x.par(2),dRand.FitResults)+Q_i;
        x_rand_err = cellfun(@(x) x.err(2),dRand.FitResults);
        xStr = sprintf('{\\itE}_0^{fit}');
        Unit = sprintf('eV');
end



%% load Alt run lists
ParAlt = zeros(numel(AltPixLists),2);
ErrAlt = zeros(numel(AltPixLists),2);
nSeg   = zeros(numel(AltPixLists),1);
nPix   = zeros(numel(AltPixLists),2);
p      = zeros(numel(AltPixLists),2); % pvalue
PixNumList = cell(2*numel(AltPixLists),1);

for i=1:numel(AltPixLists)
    savename = sprintf('%sknm2_PixListAlt_%s_%s_%s_%.0feV_KNM2.mat',...
        savedir,AltPixLists{i},DataType,strrep(freePar,' ',''),range);
    
    if exist(savename,'file')
        d = importdata(savename);
        switch Parameter
            case 'mNuSq'
                ParAlt(i,:) = d.mNuSq;
                ErrAlt(i,:) = 0.5*(d.FitResult.errPos(:,1)-d.FitResult.errNeg(:,1)); % mean asymmetric err
            case 'E0'
                ParAlt(i,:) = d.E0;
                ErrAlt(i,:) = d.E0Err;
        end
        
        p(i,:) = 1-chi2cdf(d.FitResult.chi2min,d.FitResult.dof);
        nSeg(i) = numel(d.E0);
        nPix(i,:) = cellfun(@(x) numel(x),d.PixList)';   
          PixNumList(2*(i-1)+1:2*i) = d.PixList;
    else
        fprintf('alternative pix list %s not computed yet \n',AltPixLists{i})
    end
end

nPix = [nPix_u;reshape(nPix',[2*numel(AltPixLists),1])];
x = [Par_u; reshape(ParAlt',[2*numel(AltPixLists),1])];
xerr = [Err_u;reshape(ErrAlt',[2*numel(AltPixLists),1])];

pp =[p_u;reshape(p',[2*numel(AltPixLists),1])];

%% plot
LocalFontSize = 18;

figureHandle = figure('Units','normalized','Position',[0.1,0.1,0.6,0.5]);
s1 = subplot(1,4,1:3);
[l,a_rand] = boundedline(mean(x_rand).*ones(10,1),linspace(-1,3.*numel(AltPixLists),10),std(x_rand).*ones(10,1),'orientation','horiz');
hold on
l.LineStyle = 'none'; l.Color = rgb('DimGray'); l.LineWidth = 2;
a_rand.FaceColor = rgb('Silver');
a_rand.FaceAlpha = 0.5;
pref = plot(Par_u.*ones(10,1),linspace(0,numel(x)+1,10),':','LineWidth',2,'Color',rgb('LimeGreen'));
e_fits= errorbar(x,1:numel(x),xerr,'horizontal',...
    '.','MarkerSize',20,'LineStyle','none','Color',rgb('DodgerBlue'),'CapSize',0,'LineWidth',2);
e_bf = errorbar(x(1),1,xerr(1),'horizontal',...
    '.','MarkerSize',20,'LineStyle','none','Color',pref.Color,'CapSize',0,'LineWidth',2);


xlabel(sprintf('{\\itm}_\\nu^{ 2} (eV^{ 2})'));
set(gca,'YMinorTick','off');
yticks(1:numel(x));
yticklabels(PixList_Labels)
ax1 = gca;
xlim([-1.2 1.6]);
ylabel('Pixel selection');
PrettyFigureFormat('FontSize',LocalFontSize);


leg = legend(a_rand,sprintf('Random half pixels: 1\\sigma band (%.0f samples)',nRandFits));
PrettyLegendFormat(leg,'alpha',0.8); leg.EdgeColor = 'none';
leg.Location= 'northwest';
leg.FontSize = LocalFontSize-2;



s2 = subplot(1,4,4);
area(linspace(0,0.05,10),numel(x)+1.*ones(10,1),'FaceColor',rgb('Red'));%(0.05.*ones(10,1),linspace(0,numel(RunLists)+1,10),'-','MarkerSize',20,'LineWidth',2,'Color',rgb('Red'));
hold on;
plot(pp,1:numel(x),'.','MarkerSize',20,'LineWidth',2,'Color',rgb('DodgerBlue'));
plot(pp(1),1,'.','MarkerSize',20,'LineWidth',2,'Color',pref.Color);

yticklabels('');
xlabel(sprintf('{\\it p}'));
xlim([0 1]);
ax2 = gca;
PrettyFigureFormat('FontSize',LocalFontSize);

linkaxes([s1,s2],'y');

ax1.Position(2) = 0.17;
ax2.Position(2) = ax1.Position(2);
ax2.Position(4) = ax1.Position(4);
ax1.Position(1) = 0.2;
ax2.Position(1) = 0.79;

ylim([0.5 numel(x)+1])

pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);
pltname = [pltdir,'knm2_AltPixLists.pdf'];
export_fig(gcf,pltname);


%% plot FPD viewer with pixel selection
FPDView = 'ON';
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
        export_fig(gcf,[pltdir,'knm2_FPDView_InOut.pdf']);
    elseif strcmp(ViewList,'AziHalfNS') 
        xlabel('Northern/Southern pixels');
         ax.XLabel.FontSize = 24;
         export_fig(gcf,[pltdir,'knm2_FPDView_NS.pdf']);
    elseif strcmp(ViewList,'AziHalfEW')
        xlabel('Eastern/Western pixels');
         ax.XLabel.FontSize = 24;
         export_fig(gcf,[pltdir,'knm2_FPDView_EW.pdf']);
    end
    
   
end


%% compatibility of uncorrelated pixel selections 

% uncertainty on difference through error prob. 
Sigma1 =  abs(x(3)-x(2))./sqrt(xerr(2)^2+xerr(3)^2);
Sigma2 =  abs(x(4)-x(5))./sqrt(xerr(4)^2+xerr(5)^2);
Sigma3 =  abs(x(6)-x(7))./sqrt(xerr(6)^2+xerr(7)^2);

fprintf('Inner vs. Outer: %.2f sigma \n',Sigma1);%abs(x(3)-x(2))./mean(xerr(2:3)));
fprintf('North vs. South: %.2f sigma \n',Sigma2);%abs(x(5)-x(4))./mean(xerr(4:5)));
fprintf('East vs. West  : %.2f sigma \n',Sigma3);%abs(x(7)-x(6))./mean(xerr(6:7)));

%%
% sigma with respect to random half
Sigma = abs(x-median(x_rand))./std(x_rand);%sqrt(xerr.^2+var(x_rand));

for i=1:numel(Sigma)
    fprintf('%s within %.1f sigma random half uncertainty band \n',PixList_Labels{i},Sigma(i));
end

