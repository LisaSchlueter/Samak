% plt file calculated in ksn2_NuMassSensititiyFit_Grid.m
chi2 = 'chi2CMShape';
DataType = 'Real';%'Twin';
SavePlt = 'ON';
PullFlag = 27;
nGridSteps = 10;
if strcmp(DataType,'Real')
    mNuSq = -1:0.2:2;%.5;
else
    mNuSq = -1:0.1:1;
end

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];
savefile = sprintf('%sksn2_NuMassSensitivityFits_%s_Pull%.0f_%s_mNuSqMin%.2f_mNuSqMax%.2f_mNuSteps%.0f_Grid%.0f.mat',...
    savedir,DataType,PullFlag,chi2,min(mNuSq),max(mNuSq),numel(mNuSq),nGridSteps);
d = importdata(savefile);

chi2 = reshape(d.chi2min,numel(d.mNuSq),nGridSteps,nGridSteps);

pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);

%%
chi2min = min(d.chi2min');
mNuSq_inter = linspace(min(d.mNuSq),max(d.mNuSq),5e3);
chi2min_inter = interp1(d.mNuSq,chi2min,mNuSq_inter,'spline');

MinIdx = find(chi2min_inter==min(chi2min_inter));

mNuSqUp = interp1(chi2min_inter(mNuSq_inter>1.1),mNuSq_inter(mNuSq_inter>1.1),min(chi2min_inter)+1,'spline');
mNuSqDown = interp1(chi2min_inter(mNuSq_inter<1.1),mNuSq_inter(mNuSq_inter<1.1),min(chi2min_inter)+1,'spline');
mNuSqErrUp = mNuSqUp-mNuSq_inter(MinIdx);
mNuSqErrDown = mNuSq_inter(MinIdx)-mNuSqDown;
mNuSq_bf = mNuSq_inter(MinIdx);
MeanErr = 0.5.*(mNuSqErrDown+mNuSqErrUp);
GetFigure
pbf = plot(d.mNuSq,ones(numel(d.mNuSq),1).*(min(chi2min)+1),':k','LineWidth',1.5);
hold on;
plot(d.mNuSq,ones(numel(d.mNuSq),1).*min(chi2min),':k','LineWidth',1.5);
pchi2 = plot(mNuSq_inter,chi2min_inter,'LineWidth',2.5);
plot(d.mNuSq,chi2min,'.','LineWidth',2.5,'MarkerSize',12);
xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
ylabel(sprintf('\\chi^2'));
PrettyFigureFormat('FontSize',22);
ylim([min(chi2min)-0.5 28])
leg = legend(pchi2,sprintf('{\\itm}_\\nu^2 = %.1f \\pm %.2f (^{-%.2f}_{+%.2f}) eV^2',mNuSq_bf,MeanErr,mNuSqErrDown,mNuSqErrUp));
PrettyLegendFormat(leg);
pltnamechi2 = sprintf('%sksn2_NuMassSensitivity_Fit_chi2%s.png',pltdir,DataType);
print(pltnamechi2,'-dpng','-r350');
%%
GetFigure;
histogram(d.chi2min(5,:),'BinWidth',1);
xlim([20 40])
%%
s = cell(numel(mNuSq),1);
for i=1:numel(mNuSq)
    PltIdx = mod(i,4);
    
    if mod(i,4)==1
        f1 = figure('Units','normalized','Position',[0.1,0.1,1,1]);
        pltname = sprintf('%sNuMassSensitivityFit_Grid_mNu%.1feV2.png',pltdir,mNuSq(i));
    elseif mod(i,4)==0
        PltIdx = 4;
    end
    s{i} = subplot(2,2,PltIdx);
    chi2_tmp = squeeze(chi2(i,:,:));
    xThres = median(median(chi2_tmp)).*1.2;
    chi2_tmp(chi2_tmp>xThres) = NaN;
    surf(d.sin2T4_i,d.mNu4Sq_i,chi2_tmp,'EdgeColor','none');
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    view(2);
    xlim([1e-03 0.5])
    ylim([0.1 1600])
    c = colorbar;
    c.Label.String = sprintf('\\chi^2');
    grid off
    PrettyFigureFormat('FontSize',22);
    c.Label.FontSize = get(gca,'FontSize');
    c.Limits = [min(min(chi2_tmp)),xThres];
    zlim([min(min(chi2_tmp)),xThres])
    leg = legend(sprintf('{\\itm}_\\nu^2 = %.1f eV^2',mNuSq(i)),'Location','southwest');
    PrettyLegendFormat(leg); leg.FontSize = get(gca,'FontSize');
    xlabel(sprintf('Fit init.: |{\\itU}_{e4}|^2'));
    ylabel(sprintf('Fit init.: {\\itm}_{4}^2 (eV^2)'));
    if mod(i,4)==0
             print(gcf,pltname,'-dpng','-r350');
    end
end

%% best fit location

sin2T4_tmp = zeros(16,100);
mNu4Sq_tmp = zeros(16,100);
for i=1:size(d.chi2min,2)
tmp = d.FitResults{i};
sin2T4_tmp(:,i) = cellfun(@(x) x.par(16),tmp);
mNu4Sq_tmp(:,i) = cellfun(@(x) x.par(15),tmp);
end

BfIdx = zeros(numel(mNuSq),1);%(find(d.chi2min'== chi2min));
sin2T4_bf = zeros(16,1);
mNu4Sq_bf = zeros(16,1);
for i=1:numel(mNuSq)
BfIdx(i) = (find(d.chi2min(i,:)== chi2min(i)));
sin2T4_bf(i) = sin2T4_tmp(i,BfIdx(i));
mNu4Sq_bf(i) = mNu4Sq_tmp(i,BfIdx(i));
end
%%
GetFigure;
plot(sin2T4_bf,mNu4Sq_bf,'.','MarkerSize',15);
PrettyFigureFormat;
% set(gca,'XScale','log');
 set(gca,'YScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_{4}^2 (eV^2)'));