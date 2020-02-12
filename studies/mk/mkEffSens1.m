%% Variant 1) Nominal Analysis
%
% Sensitivity study concering active neutrino mass sensitivity at KATRIN
% and the impact of a sterile neutrino admixture. The code distinguishes three scenarios:
%
% 1) Nominal effective neutrino mass (meff) Analysis:
%    Minimization of the model with regards to the normalisation factor (N), the background rate
%    (Rbkg) and the spectrum energy endpoint (E0, ie the Q-value).
%
% 2) Maximal Estimation/extended marginalization:
%    Additional model minimization using a sterile neutrino mass (m4Sq) and the mixing (sinSqTh4).
%
% 3) Minimal Estimation:
%    Assuming the existence of sterile neutrino with given parameters (m4Sq,sinSqTh4), its parameter
%    space is scanned and for each value (m4Sq,sinSqTh4) the typical meff-analysis (as in 1))
%    repeated.
%
% Another characteristic of this class, is the option to switch between two fitters
% - Minuit (www.fis.unipr.it/~giuseppe.allodi/Fminuit/Fminuit_intro.html)
% - Matlab fitter provided within the curvefitting toolbox 
%   !fitting toolbox breaks down sometimes ... try if following works: 
%   restoredefaultpath; g = fittype( @(a, b, c, x) a*x.^2+b*x+c )
%
% May 2017 - M. Korzeczek - KIT, IEKP
%

clearvars; clear global; delete(findall(0,'Type','figure')); 
addpath('fminuit'); addpath('tools');  %Paths

%Spectrum options inherited from TBD()
CommonOpt = {'TD','DR30','TimeYear',3,'nTeBinningFactor',20, ...
           'BKG_RateAllFPDSec',100e-3,'TTFSD','DOSS'}; 
%DataOpt   = [CommonOpt,{'mnuSq_i',mtest_arr(lnm)^2}]; %options inherited from TBD()
TheoryOpt = [CommonOpt,{'mnuSq_i',0}]; %options inherited from TBD()


%% nominal analysis fitting loop
mtest_arr = linspace(0,.21,10); nmtest=length(mtest_arr); %Calculation masses for fit procedure
npar = 6;  par = ones(nmtest,npar);  err= zeros(nmtest,npar);   chi2 = ones(nmtest,1);

tic;
progressbar('effective mass fit');
for lnm=1:1:nmtest
    progressbar(lnm/nmtest);
    
    DataOpt   = [CommonOpt,{'mnuSq_i',mtest_arr(lnm)^2}];
    [chi2t,part,errt,TBDData]=mkChiFunc.minimizeTheoryToData(DataOpt, TheoryOpt, 'PlotsShow','ON','Fitter','Minuit');
    par(lnm,:)=part; err(lnm,:)=errt; chi2(lnm)=chi2t;
end
toc;

%% Standard plots
TBDData.PlotTD('fign',10);
%TBDData.PlotTBDDS('fign',11);
%TBDData.PlotTBDIS('fign',12);
%TBDData.PlotKTF('fign',13);
%TBDData.PlotFSD('fign',14,'TT','ON');
%TBDData.PlotPIS('fign',15);


%% chi2 values and corresponding mass sensitivities
Prob=[.90,.95,.99];  nScannedPar=1;

chiCL     = mkChiFunc.GetChiCL( Prob, nScannedPar);
mSensStat = mkChiFunc.GetMassStatCL(chi2,mtest_arr,Prob,nScannedPar);
mSensTot  = mkChiFunc.GetMassTotCL(chi2,mtest_arr,Prob,nScannedPar);

fprintf('\n|    CL    | mnu_stat |  mnu_tot |\n');
mkChiFunc.print([Prob; mSensStat; mSensTot]);

figure(2);
plot(mtest_arr, chi2); hold on;
plot(repmat([mtest_arr(1);mtest_arr(end)],1,3), [chiCL; chiCL]); hold off;
text(repmat(mtest_arr(end),1,3)/2.3,chiCL*1.1,strcat(num2str(Prob(:)*100),'\% CL'),'FontSize',12);
title('statistical sensitivity - nominal analysis');
xlabel('tested mass $m_\mathrm{eff}$ in eV'); ylabel('goodness of fit $\chi^2$');
legend(...
    sprintf('$%3.0f/%3.0f\\,$meV @$%2.0f$CL$^\\mathrm{stat/tot}$',...
    mSensStat(1)*1e3,mSensTot(1)*1e3,round(Prob(1)*100)),'Location','northwest');
FigureFormat;
