% combine knm1 and knm2
% by adding chi^2 curves
DataType = 'Real';
chi2 = 'chi2CMShape';
%% load chi2-profiles (pre-calculated)
% knm1
k1file = [getenv('SamakPath'),sprintf('tritium-data/fit/Knm1/Chi2Profile/Uniform/Chi2Profile_%s_UniformScan_mNu_Knm1_UniformFPD_%s_NP1.064_FitParE0BkgNorm_nFit20_min-2_max1.mat',DataType,chi2)];
d1 = importdata(k1file); fprintf('load knm1: %s \n',k1file);
ScanResults1 = d1.ScanResults;
mNuSq1 =[flipud(ScanResults1.ParScan(:,2));ScanResults1.ParScan(2:end,1)];
Chi21 = [flipud(ScanResults1.chi2min(:,2));ScanResults1.chi2min(2:end,1)];

% knm2
k2file = [getenv('SamakPath'),sprintf('tritium-data/fit/Knm2/Chi2Profile/Uniform/Chi2Profile_%s_UniformScan_mNu_Knm2_UniformFPD_%s_NP1.112_FitParE0BkgNorm_nFit20_min-2_max1.mat',DataType,chi2)];
d2 = importdata(k2file); fprintf('load knm2: %s \n',k2file)
ScanResults2 = d2.ScanResults;
mNuSq2 =[flipud(ScanResults2.ParScan(:,2));ScanResults2.ParScan(2:end,1)];
Chi22 = [flipud(ScanResults2.chi2min(:,2));ScanResults2.chi2min(2:end,1)];

% check if mNuSq binning is the same: 
if any(mNuSq2 ~= mNuSq1)
    fprintf('binning not the same - interpolation necessary! \n');
    return

end
%% interpolate
mNuSq        = linspace(min(mNuSq1),max(mNuSq1),5e3);
Chi21_plot   = interp1(mNuSq1,Chi21,mNuSq,'spline');
Chi22_plot   = interp1(mNuSq1,Chi22,mNuSq,'spline');
Chi2sum_plot = interp1(mNuSq1,Chi21+Chi22,mNuSq,'spline');
%% find minimum and uncertainty
%knm1
Chi21_min = min(Chi21_plot);
mNuSq1_bf = mNuSq(Chi21_plot==Chi21_min);
mNuSq1_errNeg = interp1(Chi21_plot(mNuSq<mNuSq1_bf),mNuSq(mNuSq<mNuSq1_bf),Chi21_min+1,'spline')-mNuSq1_bf;
mNuSq1_errPos = interp1(Chi21_plot(mNuSq>mNuSq1_bf),mNuSq(mNuSq>mNuSq1_bf),Chi21_min+1,'spline')-mNuSq1_bf;

%knm2
Chi22_min = min(Chi22_plot);
mNuSq2_bf = mNuSq(Chi22_plot==Chi22_min);
mNuSq2_errNeg = interp1(Chi22_plot(mNuSq<mNuSq2_bf),mNuSq(mNuSq<mNuSq2_bf),Chi22_min+1,'spline')-mNuSq2_bf;
mNuSq2_errPos = interp1(Chi22_plot(mNuSq>mNuSq2_bf),mNuSq(mNuSq>mNuSq2_bf),Chi22_min+1,'spline')-mNuSq2_bf;

%knm1+2
Chi2sum_min = min(Chi2sum_plot);
mNuSqsum_bf = mNuSq(Chi2sum_plot==Chi2sum_min);
mNuSqsum_errNeg = interp1(Chi2sum_plot(mNuSq<mNuSqsum_bf),mNuSq(mNuSq<mNuSqsum_bf),Chi2sum_min+1,'spline')-mNuSqsum_bf;
mNuSqsum_errPos = interp1(Chi2sum_plot(mNuSq>mNuSqsum_bf),mNuSq(mNuSq>mNuSqsum_bf),Chi2sum_min+1,'spline')-mNuSqsum_bf;
%%
GetFigure;
p1 = plot(mNuSq,Chi21_plot,':','LineWidth',2);
hold on
p1bf = errorbar(mNuSq1_bf,Chi21_min,0,0,mNuSq1_errNeg,mNuSq1_errPos,'.','MarkerSize',20,'Color',p1.Color,'LineWidth',p1.LineWidth,'CapSize',0);

p2 = plot(mNuSq1,Chi22,'-.','LineWidth',2);
p2bf = errorbar(mNuSq2_bf,Chi22_min,0,0,mNuSq2_errNeg,mNuSq2_errPos,'.','MarkerSize',20,'Color',p2.Color,'LineWidth',p2.LineWidth,'CapSize',0);

psum = plot(mNuSq1,Chi21+Chi22,'LineWidth',2);
psumbf = errorbar(mNuSqsum_bf,Chi2sum_min,0,0,mNuSqsum_errNeg,mNuSqsum_errPos,'.','MarkerSize',20,'Color',psum.Color,'LineWidth',psum.LineWidth,'CapSize',0);

PrettyFigureFormat('FontSize',22);
leg = legend([p1,p2,psum,p1bf,p2bf,psumbf],...
    'KNM-1','KNM-2','KNM-1 and KNM-2',...
    sprintf('best fit: {\\itm}_\\nu^2 = %.2f \\pm %.2f eV^2',mNuSq1_bf,0.5.*(mNuSq1_errPos-mNuSq1_errNeg)),...
   sprintf('best fit: {\\itm}_\\nu^2 = %.2f \\pm %.2f eV^2',mNuSq2_bf,0.5.*(-mNuSq2_errNeg+mNuSq2_errPos)),...
    sprintf('best fit: {\\itm}_\\nu^2 = %.2f \\pm %.2f eV^2',mNuSqsum_bf,0.5.*(-mNuSqsum_errNeg+mNuSqsum_errPos)));
leg.NumColumns = 2;
leg.FontSize = get(gca,'FontSize')-3;
xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
ylabel(sprintf('\\chi^2'));
PrettyLegendFormat(leg);
xlim([-2,1])

%% save plot
savedir = [getenv('SamakPath'),'knm2ana/knm2_Combination/plots/'];
MakeDir(savedir);
savename = sprintf('%sknm2_CombiChi2_%s_Uniform_%s.png',savedir,DataType,chi2);
print(gcf,savename,'-dpng','-r300');
fprintf('save plot to %s \n',savename);
