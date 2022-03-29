% KNM1 Alternative run lists
% save or load results
% plot results

DataType    = 'Real';
freePar = 'mNu E0 Bkg Norm';
range       = 40;
FSDFlag     = 'KNM2_0p1eV';
BKG_PtSlope = 3e-06;

savedir = [getenv('SamakPath'),'knm2ana/knm2_AltRunPixLists/results/'];

% define run lists
AltRunLists = {'KNM2_Up',...
    'KNM2_Down',...
    'KNM2_RW1',...
    'KNM2_RW2',...
     'KNM2_RW3'};
RunLists_Labels = {'Golden scans','Up scans','Down scans',...
                 sprintf('{\\itU}_{RW} = - 49.6 mV'),...
                 sprintf('{\\itU}_{RW} = -7.7 mV'),...
                 sprintf('{\\itU}_{RW} = 193 mV')};
%%  prepare to load  results
nRuns_v = zeros(numel(AltRunLists)+1,1);
mNuSq_v = zeros(numel(AltRunLists)+1,1);
mNuSqErr_v = zeros(numel(AltRunLists)+1,1);
chi2min_v = zeros(numel(AltRunLists)+1,1);
p_v = zeros(numel(AltRunLists)+1,1);

%% uniform fit result
savedir_u = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savename_u = sprintf('%sknm2ubfinal_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s.mat',...
    savedir_u,BKG_PtSlope*1e6,DataType,range,strrep(freePar,' ',''),'chi2Stat+','StackPixel','KNM2');
d_u           = importdata(savename_u);
nRuns_v(1)    = numel(d_u.A.RunList);
mNuSq_v(1)    = d_u.FitResult.par(1);
mNuSqErr_v(1) = 0.5.*(d_u.FitResult.errPos(1)-d_u.FitResult.errNeg(1));
chi2min_v(1)  = d_u.FitResult.chi2min;
p_v(1)        = 1-chi2cdf(d_u.FitResult.chi2min,d_u.FitResult.dof);

%% load alternative run lists
for i=2:numel(AltRunLists)+1
    filename_tmp = sprintf('%sknm2_AltRunList_%s_%s_%s_%.0feV_%s_BkgPt%.2g.mat',...
    savedir,AltRunLists{i-1},DataType,strrep(freePar,' ',''),range,FSDFlag,BKG_PtSlope*1e6);

    if exist(filename_tmp,'file')
        fprintf('Load file %s\n',filename_tmp);
        ftmp = importdata(filename_tmp);
        nRuns_v(i) = numel(ftmp.RunList);
        mNuSq_v(i) = ftmp.FitResult.par(1);
        mNuSqErr_v(i) =  0.5.*(ftmp.FitResult.errPos(1)-ftmp.FitResult.errNeg(1));
        chi2min_v(i) =  ftmp.FitResult.chi2min;
        p_v(i) = 1-chi2cdf(ftmp.FitResult.chi2min,ftmp.FitResult.dof);  
    else
        fprintf(2,'file not found %s\n',filename_tmp);
        return
    end
end

%% print results on screen
for i=1:numel(AltRunLists)+1
    
    fprintf('m2 = %.2f +- %.2f eV^2 , chi2min = %.1f, p = %.2f, %.0f runs, %s \n',...
        mNuSq_v(i),mNuSqErr_v(i),chi2min_v(i),p_v(i),nRuns_v(i),RunLists_Labels{i});
end

%% plot
RandList = 'ON';
LocalFontSize = 18;
figureHandle = figure('Units','normalized','Position',[0.1,0.1,0.6,0.5]);
s1 = subplot(1,4,1:3);

if strcmp(RandList,'ON')
    nRunList = 2000;
    chi2 = 'chi2Stat';
    NP = 1.064;
    savenameRand = [savedir,'knm2_RunListRandHalf_Real_mNuE0BkgNorm_40eV_2000fits_KNM2_0p1eV_BkgPT3mucpsPers.mat'];    
   % sprintf('RandomHalfRunList_Unblinded_%s_NP%2g_%.0ffits_merged.mat',chi2,NP,nRunList)];
    frand = importdata(savenameRand);
    mNuSq_rand =  cell2mat(cellfun(@(x) x.par(1),frand.FitResults,'UniformOutput',0));
    [l,a_rand] = boundedline(mean(mNuSq_rand).*ones(10,1),linspace(-1,numel(AltRunLists)+3,10),std(mNuSq_rand).*ones(10,1),'orientation','horiz');
    hold on
    l.LineStyle = 'none'; l.Color = rgb('DimGray'); l.LineWidth = 2;
    a_rand.FaceColor = rgb('Silver');
    a_rand.FaceAlpha = 0.5;   
end

%
pref = plot(mNuSq_v(1).*ones(10,1),linspace(-1,numel(AltRunLists)+2,10),':','LineWidth',2.5,'Color',rgb('LimeGreen'));
hold on;
e_fits = errorbar(mNuSq_v,1:numel(AltRunLists)+1,mNuSqErr_v,'horizontal',...
    '.','MarkerSize',20,'LineStyle','none','Color',rgb('DodgerBlue'),'CapSize',0,'LineWidth',2);
e_bf = errorbar(mNuSq_v(1),1,mNuSqErr_v(1),'horizontal',...
    '.','MarkerSize',20,'LineStyle','none','Color',pref.Color,'CapSize',0,'LineWidth',2);
xlabel(sprintf('{\\itm}_\\nu^{ 2} (eV^{ 2})'));
ylabel('Scan selection');

PrettyFigureFormat('FontSize',LocalFontSize);
set(gca,'YMinorTick','off');
yticks(1:numel(AltRunLists)+1);
yticklabels(RunLists_Labels)
ax1 = gca;
xlim([-1 1.2]);

if strcmp(RandList,'ON')
    leg = legend([a_rand],...
       sprintf('Random half scans: 1\\sigma band (%.0f samples)',nRunList));
    PrettyLegendFormat(leg,'alpha',0.8); leg.EdgeColor = 'none';
    leg.Location= 'northwest';
    leg.FontSize = LocalFontSize-2;
end

s2 = subplot(1,4,4);
area(linspace(0,0.05,10),(numel(AltRunLists)+2).*ones(10,1),'FaceColor',rgb('Red'));%(0.05.*ones(10,1),linspace(0,numel(RunLists)+1,10),'-','MarkerSize',20,'LineWidth',2,'Color',rgb('Red'));
hold on;
plot(p_v,1:numel(AltRunLists)+1,'.','MarkerSize',20,'LineWidth',2,'Color',rgb('DodgerBlue'));
plot(p_v(1),1,'.','MarkerSize',20,'LineWidth',2,'Color',pref.Color);

yticklabels('');
xlabel(sprintf('{\\it p}'));
ax2 = gca;
xlim([0 1])
ylim([0.5 numel(AltRunLists)+1.5])
PrettyFigureFormat('FontSize',LocalFontSize);
linkaxes([s1,s2],'y');

ax1.Position(2) = 0.17;
ax2.Position(2) = ax1.Position(2);
ax2.Position(4) = ax1.Position(4);
ax1.Position(1) = 0.2;
ax2.Position(1) = 0.79;

if strcmp(RandList,'ON')
    ylim([0 numel(AltRunLists)+2.5])
else
    ylim([0.5 numel(AltRunLists)+1.5])
end

% save plot
pltname = [strrep(savedir,'results','plots'),'knm2_AltRunLists.pdf'];
 export_fig(gcf,pltname);

 %% how many results lie within the 1sigma band from the random half?
mNuSq_low =  mean(mNuSq_rand)-std(mNuSq_rand);
mNuSq_up =  mean(mNuSq_rand)+std(mNuSq_rand);

LogIdx = mNuSq_v>=mNuSq_low & mNuSq_v<=mNuSq_up;
Frac = sum(LogIdx)./numel(LogIdx);
fprintf('%.0f out of %.0f (%.1f%%) fit results are within 1sigma\n',sum(LogIdx),numel(LogIdx),1e2*Frac);
 
%% calculate sigmas
Sigmas = (mNuSq_v-mean(mNuSq_rand))./std(mNuSq_rand);


%% upward/downward scans
mNuSqDiff = abs(mNuSq_v(2)-mNuSq_v(3));
mNuSqDifferr = sqrt(mNuSqErr_v(2)^2+mNuSqErr_v(3)^2);
Sigma = mNuSqDiff./mNuSqDifferr;


