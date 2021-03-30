% combine knm1 and knm2
% by adding chi^2 curves
DataType = 'Real';
chi2 = 'chi2CMShape';
Knm2AnaFlag = 'Uniform';%MR-4';
nFit = 50;
%% load chi2-profiles (pre-calculated)
% knm1
if strcmp(chi2,'chi2CMShape')
    chi2Str1 = 'chi2CMShape_SysBudget24';
else
     chi2Str1 = chi2;
end

d1 = LoadChi2Profile('DataSet','Knm1','DataType',DataType,'chi2',chi2,'AnaStr','Uniform','nFit',nFit,'mNuSqMin',-2.6,'mNuSqMax',1);
ScanResults1 = d1.ScanResults;
mNuSq1 =[flipud(ScanResults1.ParScan(:,2));ScanResults1.ParScan(2:end,1)];
Chi21 = [flipud(ScanResults1.chi2min(:,2));ScanResults1.chi2min(2:end,1)];
mNuSq1_bf    = ScanResults1.BestFit.par;
mNuSq1_errNeg = ScanResults1.BestFit.errNeg;
mNuSq1_errPos = ScanResults1.BestFit.errPos;
Chi21_min     = ScanResults1.BestFit.chi2;

% knm2
d2 = LoadChi2Profile('DataSet','Knm2','DataType',DataType,'chi2',chi2,'AnaStr',Knm2AnaFlag,'nFit',nFit,'mNuSqMin',-2.6,'mNuSqMax',1);
ScanResults2 = d2.ScanResults;
mNuSq2 =[flipud(ScanResults2.ParScan(:,2));ScanResults2.ParScan(2:end,1)];
Chi22 = [flipud(ScanResults2.chi2min(:,2));ScanResults2.chi2min(2:end,1)];
mNuSq2_bf = ScanResults2.BestFit.par;
mNuSq2_errNeg = ScanResults2.BestFit.errNeg;
mNuSq2_errPos = ScanResults2.BestFit.errPos;
Chi22_min     = ScanResults2.BestFit.chi2;
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

%% find common minimum and uncertainty

ksumfile = [getenv('SamakPath'),sprintf('tritium-data/fit/Knm1/Chi2Profile/Uniform/Chi2ProfileCombi_%s_UniformScan_mNu_Knm1KNM2_UniformFPD_%s_FitParE0BkgNorm_nFit20_min-2.6_max1.mat',DataType,chi2)];
if exist(ksumfile,'file')
    dsum = importdata(ksumfile);
    Chi2sum_min = dsum.BestFit.chi2;
    mNuSqsum_bf = dsum.BestFit.mNuSq;
    mNuSqsum_errNeg = dsum.BestFit.mNuSqErrNeg;
    mNuSqsum_errPos = dsum.BestFit.mNuSqErrPos;
    mNuSqsum_err = dsum.BestFit.mNuSqErr;
else
    %knm1+2
    Chi2sum_min = min(Chi2sum_plot);
    mNuSqsum_bf = mNuSq(Chi2sum_plot==Chi2sum_min);
    mNuSqsum_errNeg = interp1(Chi2sum_plot(mNuSq<mNuSqsum_bf),mNuSq(mNuSq<mNuSqsum_bf),Chi2sum_min+1,'spline')-mNuSqsum_bf;
    mNuSqsum_errPos = interp1(Chi2sum_plot(mNuSq>mNuSqsum_bf),mNuSq(mNuSq>mNuSqsum_bf),Chi2sum_min+1,'spline')-mNuSqsum_bf;
    mNuSqsum_err = 0.5*(mNuSqsum_errPos-mNuSqsum_errNeg);
    
    BestFit = struct('chi2',Chi2sum_min,'mNuSq',mNuSqsum_bf,'mNuSqErrPos',mNuSqsum_errPos,...
        'mNuSqErrNeg',mNuSqsum_errNeg,'mNuSqErr',mNuSqsum_err);
    save(ksumfile,'mNuSq','Chi2sum_plot','BestFit');
end
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
    'KNM-1 (Uniform)',sprintf('KNM-2 (%s)',Knm2AnaFlag),'KNM-1 and KNM-2',...
    sprintf('best fit: {\\itm}_\\nu^2 = %.2f \\pm %.2f eV^2',mNuSq1_bf,ScanResults1.BestFit.errMean),...
   sprintf('best fit: {\\itm}_\\nu^2 = %.2f \\pm %.2f eV^2',mNuSq2_bf,ScanResults2.BestFit.errMean),...
    sprintf('best fit: {\\itm}_\\nu^2 = %.2f \\pm %.2f eV^2',mNuSqsum_bf,0.5.*(-mNuSqsum_errNeg+mNuSqsum_errPos)));
leg.NumColumns = 2;
leg.FontSize = get(gca,'FontSize')-3;
xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
ylabel(sprintf('\\chi^2'));
PrettyLegendFormat(leg);
xlim([-2.1,1])
if strcmp(Knm2AnaFlag,'MR-4')
    ylimMax = max(ylim)+33;
    ylim([min(ylim),ylimMax]);
end
%% save plot
savedir = [getenv('SamakPath'),'knm2ana/knm2_Combination/plots/'];
MakeDir(savedir);
savename = sprintf('%sknm2_CombiChi2_%s_Uniform_%s_KNM2%s.png',savedir,DataType,chi2,Knm2AnaFlag);
print(gcf,savename,'-dpng','-r300');
fprintf('save plot to %s \n',savename);


