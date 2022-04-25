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

if strcmp(DataType,'Twin')
    d1 = LoadChi2Profile('DataSet','Knm1','DataType',DataType,'chi2',chi2,'AnaStr','Uniform','nFit',100,'mNuSqMin',-5,'mNuSqMax',5);   
else
    d1 = LoadChi2Profile('DataSet','Knm1','DataType',DataType,'chi2',chi2,'AnaStr','Uniform','nFit',200,'mNuSqMin',-5,'mNuSqMax',5);
end
ScanResults1 = d1.ScanResults;
mNuSq1 =[flipud(ScanResults1.ParScan(:,2));ScanResults1.ParScan(2:end,1)];
Chi21 = [flipud(ScanResults1.chi2min(:,2));ScanResults1.chi2min(2:end,1)];
mNuSq1_bf    = ScanResults1.BestFit.par+0.002;
mNuSq1_errNeg = ScanResults1.BestFit.errNeg;
mNuSq1_errPos = ScanResults1.BestFit.errPos;
Chi21_min     = ScanResults1.BestFit.chi2;
dof1 = ScanResults1.dof(1)-1;
% knm2
if strcmp(DataType,'Twin')
   d2 = LoadChi2Profile('DataSet','Knm2','DataType',DataType,'chi2',chi2,'AnaStr',Knm2AnaFlag,'nFit',nFit,'mNuSqMin',-2.5,'mNuSqMax',2.5);  
else
    d2 = LoadChi2Profile('DataSet','Knm2','DataType',DataType,'chi2',chi2,'AnaStr',Knm2AnaFlag,'nFit',50,'mNuSqMin',-2.5,'mNuSqMax',2.5);
    
end
ScanResults2 = d2.ScanResults;
mNuSq2 =[flipud(ScanResults2.ParScan(:,2));ScanResults2.ParScan(2:end,1)];
Chi22 = [flipud(ScanResults2.chi2min(:,2));ScanResults2.chi2min(2:end,1)];
mNuSq2_bf = ScanResults2.BestFit.par;
mNuSq2_errNeg = ScanResults2.BestFit.errNeg;
mNuSq2_errPos = ScanResults2.BestFit.errPos;
Chi22_min     = ScanResults2.BestFit.chi2;
dof2 = ScanResults2.dof(1)-1;
% check if mNuSq binning is the same: 
% if any(mNuSq2 ~= mNuSq1)
%     fprintf('binning not the same - interpolation necessary! \n');
%     return
% end
%% interpolate
%mNuSq        = linspace(min(mNuSq1),max(mNuSq1),5e3);
 mNuSq        = linspace(-2.5,2.5,5e3);
Chi21_plot   = interp1(mNuSq1,Chi21,mNuSq,'spline');
Chi22_plot   = interp1(mNuSq2,Chi22,mNuSq,'spline');
Chi2sum_plot = Chi21_plot+Chi22_plot;%interp1(mNuSq1,Chi21+Chi22,mNuSq,'spline');

%% find common minimum and uncertainty

ksumfile = [getenv('SamakPath'),sprintf('tritium-data/fit/Knm1/Chi2Profile/Uniform/Chi2ProfileCombi_%s_UniformScan_mNu_Knm1KNM2_UniformFPD_%s_FitParE0BkgNorm_nFit%.0f_min-2.6_max1.mat',DataType,chi2,nFit)];
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
AreaFlag = 'ON'; % 1 sigma band
Colors = {rgb('DodgerBlue'),rgb('DarkOrange'),rgb('Crimson')};

f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.45]);
if strcmp(AreaFlag,'ON')
    [l,a1]= boundedline(mNuSq1_bf.*ones(10,1),linspace(0,100,10),...
        [ones(10,1).*mNuSq1_errNeg,ones(10,1).*mNuSq1_errPos],'orientation', 'horiz');
    hold on;
    l.delete; a1.FaceColor =Colors{1}; a1.FaceAlpha =0.2 ;
    
    [l,a2]= boundedline(mNuSq2_bf.*ones(10,1),linspace(0,100,10),...
        [ones(10,1).*mNuSq2_errNeg,ones(10,1).*mNuSq2_errPos],'orientation', 'horiz');
    l.delete;  a2.FaceColor = rgb('Orange'); a2.FaceAlpha = 0.9;
    
    [l,a3]= boundedline(mNuSqsum_bf.*ones(10,1),linspace(0,100,10),...
        [ones(10,1).*abs(mNuSqsum_errNeg),ones(10,1).*mNuSqsum_errPos],'orientation', 'horiz');
    l.delete;  a3.FaceColor = rgb('LightCoral'); a3.FaceAlpha = 1;
    
    if strcmp(DataType,'Real')
       a1.FaceAlpha = 0.2;
       a2.FaceAlpha = 0.35;
       a3.FaceAlpha = 0.4;
    end
end
p1 = plot(mNuSq,Chi21_plot,':','LineWidth',3,'Color',Colors{1});
hold on
p2 = plot(mNuSq,Chi22_plot,'-.','LineWidth',3,'Color',Colors{2});
psum = plot(mNuSq,Chi2sum_plot,'LineWidth',3,'Color',Colors{3});


if strcmp(DataType,'Real')
  %  p1bf = errorbar(mNuSq1_bf,Chi21_min,0,0,mNuSq1_errNeg,mNuSq1_errPos,'.','MarkerSize',20,'Color',p1.Color,'LineWidth',p1.LineWidth,'CapSize',10);   
  %  p2bf = errorbar(mNuSq2_bf,Chi22_min,0,0,mNuSq2_errNeg,mNuSq2_errPos,'.','MarkerSize',20,'Color',p2.Color,'LineWidth',p2.LineWidth,'CapSize',10);
  %  psumbf = errorbar(mNuSqsum_bf,Chi2sum_min,0,0,mNuSqsum_errNeg,mNuSqsum_errPos,'.','MarkerSize',20,'Color',psum.Color,'LineWidth',psum.LineWidth,'CapSize',10);
    p1bf = plot(mNuSq1_bf,Chi21_min,'.','MarkerSize',20,'Color',p1.Color,'LineWidth',p1.LineWidth);   
    p2bf = plot(mNuSq2_bf,Chi22_min,'.','MarkerSize',20,'Color',p2.Color,'LineWidth',p2.LineWidth);
    psumbf = plot(mNuSqsum_bf,Chi2sum_min,'.','MarkerSize',20,'Color',psum.Color,'LineWidth',psum.LineWidth);
   
    leg = legend([p1,p2,psum,p1bf,p2bf,psumbf],...
        sprintf('KNM1 (%.0f dof)',dof1),...
        sprintf('KNM2 (%.0f dof)',dof2),...
        sprintf('KNM1+2 (%.0f dof)',dof1+dof2+1),...
        sprintf('Best fit: {\\itm}_\\nu^2 = %.2f \\pm %.2f eV^2',mNuSq1_bf,ScanResults1.BestFit.errMean),...
        sprintf('Best fit: {\\itm}_\\nu^2 = %.2f \\pm %.2f eV^2',mNuSq2_bf,ScanResults2.BestFit.errMean),...
        sprintf('Best fit: {\\itm}_\\nu^2 = %.2f \\pm %.2f eV^2',mNuSqsum_bf,0.5.*(-mNuSqsum_errNeg+mNuSqsum_errPos)));
        leg.NumColumns = 2;
else
    
    leg = legend([p1,p2,psum],...
        sprintf('KNM1 (%.0f dof),     \\sigma({\\itm}_\\nu^2) = %.2f eV^2',dof1,ScanResults1.BestFit.errMean),...
        sprintf('KNM2 (%.0f dof),     \\sigma({\\itm}_\\nu^2) = %.2f eV^2',dof2,ScanResults2.BestFit.errMean),...
        sprintf('KNM1+2 (%.0f dof), \\sigma({\\itm}_\\nu^2) = %.2f eV^2',dof1+dof2+1,0.5.*(-mNuSqsum_errNeg+mNuSqsum_errPos)));
    
    leg.Title.String = 'Sensitivity at 68.3% C.L.';
    leg.Title.FontWeight = 'normal';
end

 leg.Location = 'north';
leg.FontSize = 15;
xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
ylabel(sprintf('\\chi^2'));

PrettyFigureFormat('FontSize',18);
PrettyLegendFormat(leg);
if strcmp(DataType,'Real')
    xlim([-2.5,2])
    ylim([16 99])
else
    xlim([-1.39,1.3])
    ylim([0 10])
end

if strcmp(Knm2AnaFlag,'MR-4')
    ylimMax = max(ylim)+33;
    ylim([min(ylim),ylimMax]);
end

%% save plot
savedir = [getenv('SamakPath'),'knm2ana/knm2_Combination/plots/'];
MakeDir(savedir);
savename = sprintf('%sknm2_CombiChi2_%s_Uniform_%s_KNM2%s.pdf',savedir,DataType,chi2,Knm2AnaFlag);
export_fig(savename);
%print(gcf,savename,'-dpng','-r300');
fprintf('save plot to %s \n',savename);

fprintf('mnu^2_bf common = %.4f eV2 \n',mNuSqsum_bf);
fprintf('chi2min common = %.2f \n',Chi2sum_min);

