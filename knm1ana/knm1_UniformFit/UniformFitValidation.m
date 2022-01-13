% -------------------------------------------------------------------------
% script to validate uniform fit (== stacked FPD pixels) for KNM1
% method:
% fit stacked twin data with response functions for single pixels 
% compare to fit with response function for stacked pixels
% -------------------------------------------------------------------------

RunList = 'KNM1';
if ~exist('M','var')
    M = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','DataType','Twin','exclDataStart',13,'NonPoissonScaleFactor',1);
    M.exclDataStart = 14;
    M.fixPar = '5 6 7 8 9 10 11';
    M.Fit;
    FitParStack  = M.FitResult.par;
    FitErrStack  = M.FitResult.err;
    FitChi2Stack = M.FitResult.chi2min;
end

%% step 1: get pixelwise response functions
GoldenPixList = M.PixList;
savedir = [getenv('SamakPath'),'knm1ana/knm1_UniformFit/results/'];
if ~exist(savedir,'dir')
    system(['mkdir -p ',savedir]);
end
savename = [savedir,sprintf('ResponseFunction_Pixelwise%.0f_%s.mat',numel(GoldenPixList),RunList)];

if exist(savename,'file')
    load(savename);
else
    % Get Response Function for single pixels (data)
    RF = zeros([numel(GoldenPixList),size(M.ModelObj.RF)]);
    
    progressbar('compute different RF');
    for i=1:numel(GoldenPixList)
        progressbar(i/numel(GoldenPixList));
        T = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','DataType','Real','PixList',GoldenPixList(i));
        RF(i,:,:) = T.ModelObj.RF;
    end
    
    % Get Response Function for stacked pixels (data)
    T = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','DataType','Real');
    RFStack = T.ModelObj.RF;
    save(savename,'RF','RFStack');
end

%% step 1.2 get pixelwise qU
if ~exist('qUPix','varxc')
    T = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','DataType','Real');
    T.ReadSingleRunData;
    qUPix = T.StackWmean(T.SingleRunData.qU,T.SingleRunData.TimeperSubRunperPixel);
    save(savename,'qUPix','-append');
end
%% step 2: fit twins with these response functions
savename2 = [savedir,sprintf('FitWithResponseFunction_Pixelwise%.0f_%s.mat',numel(GoldenPixList),RunList)];

if exist(savename2,'file')
    load(savename2);
else
    FitPar = zeros(numel(GoldenPixList),M.nPar);
    FitErr = zeros(numel(GoldenPixList),M.nPar);
    FitChi2 =  zeros(numel(GoldenPixList),1);
    
    progressbar('Fit with different RF');
    for i=1:numel(GoldenPixList)
        progressbar(i/numel(GoldenPixList));
        M.SimulateStackRuns('qU',qUPix(:,GoldenPixList(i)));
        M.ModelObj.RF = squeeze(RF(i,:,:));
        M.Fit;
        FitPar(i,:) = M.FitResult.par;
        FitErr(i,:) = M.FitResult.err;
        FitChi2(i)  = M.FitResult.chi2min;
    end
    save(savename2,'FitPar','FitErr','FitChi2','FitParStack','FitErrStack','FitChi2Stack');   
end

%% plot fit results (default: neutrino mass)

plotdir = strrep(savedir,'results','plots');
if ~exist(plotdir,'dir')
    system(['mkdir -p ',plotdir]);
end

f55 = figure('Name','MultiBarPlot','Renderer','opengl');
set(f55, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
Parameter = 1;
hpixel = histogram(FitPar(:,Parameter)); hpixel.FaceColor = rgb('SteelBlue'); hpixel.FaceAlpha = 0.8; hpixel.EdgeColor = rgb('DarkSlateGray');
hpixel.BinWidth = 0.0007;
hold on;
pmean = plot(mean(FitPar(:,Parameter)).*ones(2,1),[0,30],'-.','Color',rgb('GoldenRod'),'LineWidth',4);
pstack = plot(FitParStack(Parameter).*ones(2,1),[0,30],'--','Color',rgb('FireBrick'),'LineWidth',4);
PrettyFigureFormat;
xlabel(sprintf('m_\\nu^2 (eV^2)'));
%xlabel(sprintf('E_0 (eV)'))
leg = legend([hpixel, pmean,pstack],'pixelwise RF','<pixelwise RF>','stacked pixel RF'); legend boxoff;
leg.Title.String = 'response function';
hold off
ylim([0 max(hpixel.BinCounts)+1]);
plotname = [plotdir,sprintf('FitResult_RFPixelwise%.0f_%s.png',numel(GoldenPixList),RunList)];
set(gca,'FontSize',24);
print(f55,plotname,'-dpng','-r450');
%% plot response function
f1 = figure('Name','MultiBarPlot','Renderer','opengl');
set(f1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
plotRF = zeros(numel(GoldenPixList),M.ModelObj.nTe);
for i=1:numel(GoldenPixList)
   plotRF(i,:) = interp1(M.ModelObj.Te-qUPix(27,GoldenPixList(i)),RF(i,:,27),M.ModelObj.Te-M.ModelObj.qU(27));
end
Residuals = (mean(plotRF));%-RFStack(:,27)')./(RFStack(:,27)');
Residuals(isnan(Residuals)) = 0;  Residuals(isinf(Residuals)) = 0;
Residuals(Residuals>1e5) = 0;
[l,a] = boundedline(M.ModelObj.Te-M.ModelObj.qU(27),Residuals,std(plotRF));
l.Color = rgb('SteelBlue'); l.LineStyle = '-.'; l.LineWidth = 2;
a.FaceAlpha = 0.3; a.FaceColor = rgb('SteelBlue');
hold on;
pstack2 = plot(M.ModelObj.Te-M.ModelObj.qU(27),RFStack(:,27),'--','Color',rgb('FireBrick'),'LineWidth',2);
leg2 = legend([a,l,pstack2],sprintf('1\\sigma error band pixelwise'),'< pixelwise >','stacked pixel'); legend boxoff;
leg2.Title.String = 'response function'; leg2.Location='northwest';
xlabel(sprintf('E_{kin} - qU  (eV)'));
ylabel('transmission probability');
PrettyFigureFormat;
xlim([-0.5,4]);
ylim([-0.02 0.82])
hold off;
set(gca,'FontSize',24);

plotname = [plotdir,sprintf('ResponseFunction_RFPixelwise%.0f_%s.png',numel(GoldenPixList),RunList)];
print(f1,plotname,'-dpng','-r450');
