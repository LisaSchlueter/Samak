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
AltRunLists = {
    'KNM2_RW1',...
    'KNM2_RW2',...
     'KNM2_RW3'};
 
RunLists_Labels = {sprintf('{\\itU}_{RW} = - 49.6 mV'),...
                 sprintf('{\\itU}_{RW} = -7.7 mV'),...
                 sprintf('{\\itU}_{RW} = 193 mV')};
%%  prepare to load  results
nRuns_v = zeros(numel(AltRunLists),1);
mNuSq_v = zeros(numel(AltRunLists),1);
mNuSqErr_v = zeros(numel(AltRunLists),1);
chi2min_v = zeros(numel(AltRunLists),1);
p_v = zeros(numel(AltRunLists),1);
U_rw = zeros(numel(AltRunLists),1);

%% uniform fit result
savedir_u = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savename_u = sprintf('%sknm2ubfinal_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s.mat',...
    savedir_u,BKG_PtSlope*1e6,DataType,range,strrep(freePar,' ',''),'chi2Stat+','StackPixel','KNM2');
d_u           = importdata(savename_u);
nRuns_u    = numel(d_u.A.RunList);
mNuSq_u    = d_u.FitResult.par(1);
mNuSqErr_u = 0.5.*(d_u.FitResult.errPos(1)-d_u.FitResult.errNeg(1));
chi2min_u  = d_u.FitResult.chi2min;
p_u        = 1-chi2cdf(d_u.FitResult.chi2min,d_u.FitResult.dof);

%% load alternative run lists
for i=1:numel(AltRunLists)
    filename_tmp = sprintf('%sknm2_AltRunList_%s_%s_%s_%.0feV_%s_BkgPt%.2g.mat',...
    savedir,AltRunLists{i},DataType,strrep(freePar,' ',''),range,FSDFlag,BKG_PtSlope*1e6);
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
    if contains(AltRunLists{i},'RW1')
        U_rw(i) = -49.6; % mV
    elseif contains(AltRunLists{i},'RW2')
        U_rw(i) = -7.7;% mV
    elseif contains(AltRunLists{i},'RW3')
        U_rw(i) = 193;% mV
    end
end


%% print results on screen
for i=1:numel(AltRunLists)
    fprintf('m2 = %.2f +- %.2f eV^2 , chi2min = %.1f, p = %.2f, %.0f runs, %s \n',...
        mNuSq_v(i),mNuSqErr_v(i),chi2min_v(i),p_v(i),nRuns_v(i),RunLists_Labels{i});
end

%% linear Fit
 [linFitPar, linFitErr, linFitchi2min,linFitdof] = linFit(U_rw,mNuSq_v,mNuSqErr_v);
      xlin =   linspace(min(U_rw),max(U_rw),1e2);  
      ylin = xlin.*linFitPar(1)+linFitPar(2);
% plot
LocalFontSize = 18;
figureHandle = figure('Units','normalized','Position',[0.1,0.1,0.5,0.4]);


plin = plot(xlin,ylin,'-','LineWidth',2,'Color',rgb('DodgerBlue'));
hold on
e_fits = errorbar(U_rw,mNuSq_v,mNuSqErr_v, '.','MarkerSize',20,'LineStyle','none','Color',rgb('DodgerBlue'),'CapSize',0,'LineWidth',2);

ylabel(sprintf('{\\itm}_\\nu^{ 2} (eV^{ 2})'));
xlabel(sprintf('{\\itU}_{RW} (mV)'));
PrettyFigureFormat('FontSize',LocalFontSize);

ax1 = gca;
xlim([-65 210]);
ylim([-1.2 1.4]);

leg = legend([e_fits,plin],'Period-wise fits',sprintf('Slope = (%.1f \\pm %.1f)\\times10^{-3} eV^2/mV',...
    1e3.*linFitPar(1),1e3.*linFitErr(1)));
PrettyLegendFormat(leg);
leg.FontSize = get(gca,'FontSize')+2;

%save plot
pltname = [strrep(savedir,'results','plots'),'knm2_AltRunLists_RearWall.pdf'];
export_fig(gcf,pltname);
fprintf('save plot %s \n',pltname);
%%
fprintf('Slope significance %.1f sigma \n',abs(linFitPar(1)./linFitErr(1)));


