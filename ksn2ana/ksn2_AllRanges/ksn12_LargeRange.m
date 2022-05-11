%combined ksn1+ksn2 large range constraints
freePar = 'E0 Norm Bkg';
savedir1 = [getenv('SamakPath'),'ksn1ana/ksn1_ReAna/results/'];
savename1 = sprintf('%sksn1_Object_%s.mat',savedir1,strrep(freePar,' ',''));
d1 = importdata(savename1);


savedir2 = [getenv('SamakPath'),'ksn2ana/ksn2_AllRanges/results/'];
savename2 = sprintf('%sksn2_Object_%s.mat',savedir2,strrep(freePar,' ',''));
d2 = importdata(savename2);


%%
d1.S.nGridSteps = 45;
d1.S.range = max(d1.ranges);
d1.S.LoadGridFile(d1.S.LoadGridArg{:});
d1.S.Interp1Grid('RecomputeFlag','ON','Minm4Sq',1,'Maxm4Sq',90^2);
d1.S.ContourPlot('BestFit','ON')

chi2min1 = d1.S.chi2_bf;
%%
d2.S.nGridSteps = 30;
d2.S.range = max(d2.ranges);
d2.S.LoadGridFile(d2.S.LoadGridArg{:});
d2.S.Interp1Grid('RecomputeFlag','ON','Minm4Sq',1,'Maxm4Sq',90^2);
d2.S.ContourPlot('BestFit','ON','HoldOn','ON','Color',rgb('Orange'));

chi2min2 = d2.S.chi2_bf;
%% check if they have same binning
 if sum(sum(d1.S.sin2T4~=d2.S.sin2T4))>0 || sum(sum(d1.S.mNu4Sq~=d2.S.mNu4Sq))>0
     fprintf(2, 'KSN-1 and KSN-2 do not have  the same binning - return \n');
     return;
 end
 
%% combi 
% sum

chi2_sum = d1.S.chi2+d2.S.chi2;
d1.S.chi2 = chi2_sum;
d1.S.chi2_ref = min(min(chi2_sum));%chi2ref_k1+chi2ref_k2;
[p12tot,sinMin,p12bf] = d1.S.ContourPlot('BestFit','ON','CL',95,'HoldOn',...
    'ON','Color',rgb('Crimson'),'LineStyle','-.','Marker','o');
BestFit = d1.S.FindBestFit;
mNu4Sq_contour = d1.S.mNu4Sq_contour;
sin2T4_contour = d1.S.sin2T4_contour;

chi2_Null = d1.S.chi2_Null+d2.S.chi2_Null;
nqU = d1.S.RunAnaObj.ModelObj.nqU+d2.S.RunAnaObj.ModelObj.nqU;
nFitPar = 2*3+2;
BestFit.dof = nqU-nFitPar;


savedir = [getenv('SamakPath'),'ksn2ana/ksn2_AllRanges/results/'];
savename = sprintf('%sksn12Combi_MaxRange_chi2CMShape_%s.mat',savedir2,strrep(freePar,' ',''));
save(savename,'mNu4Sq_contour','sin2T4_contour','BestFit','chi2_Null');




%% phat
chi2hat = BestFit.chi2min-(chi2min2+chi2min1);
phat = 1-chi2cdf(chi2hat,2);


