%% Variant 3) Minimal Estimation
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

%% Initialization
clearvars; clear global; delete(findall(0,'Type','figure')); 
addpath('fminuit'); addpath('tools');  %Paths

%Spectrum options inherited from TBD()
CommonOpt = {'TD','DR30','TimeYear',3,'nTeBinningFactor',20, ...
           'BKG_RateAllFPDSec',10e-3,'TTFSD','SAENZ'}; 

%% Fitting loop
mtest_arr     = [0,linspace(.15,0.23,10)]; nmtest   = length(mtest_arr);    
sinSq2Th4test = linspace(0,1,20);          nsintest = length(sinSq2Th4test); sinSqTh4test = (1-sqrt(1-sinSq2Th4test))/2; 
mnu4Sqtest    = [0,logspace(-2,2,19)];     nmnu4test= length(mnu4Sqtest);    mnu4test     = sqrt(mnu4Sqtest);            

npar = 6;  par  = ones(nsintest,nmnu4test,nmtest,npar);  err  = zeros(nsintest,nmnu4test,nmtest,npar);   
chi2 = repmat(linspace(1,nmtest,nmtest),nsintest,1,nmnu4test)*1e5;chi2=permute(chi2,[1,3,2]);

tic; cnt=0;
progressbar('minimial estimation');
for ls2=1:1:nsintest
    for lm4=1:1:nmnu4test
        for lnm=1:1:nmtest
            cnt=cnt+1; progressbar(cnt/nsintest/nmnu4test/nmtest);
            
            TheoryOpt = [CommonOpt,{'sin2T4_i',sinSqTh4test(ls2),'mnu4Sq_i',mnu4Sqtest(lm4)}]; 
            DataOpt   = [TheoryOpt,{'mnuSq_i',mtest_arr(lnm)^2}];
            
            [chi2t,part,errt,TBDData]=mkChiFunc.minimizeTheoryToData(DataOpt, TheoryOpt, 'PlotsShow','OFF');
            
            par(ls2,lm4,lnm,:)=part; chi2(ls2,lm4,lnm)=chi2t;  err(ls2,lm4,lnm,:)=errt;
                    
        end
    end
end
toc

%%
%save('matfiles/mkEffSens3b.mat');

%% Standard plots
TBDData.PlotTD('fign',10);
%TBDData.PlotTBDDS('fign',11);
%TBDData.PlotTBDIS('fign',12);
%TBDData.PlotKTF('fign',13);
%TBDData.PlotFSD('fign',14,'TT','ON');
%TBDData.PlotPIS('fign',15);


%% chi2 values and corresponding mass sensitivities for specific sterile neutrino parameters
Prob=[.90,.95,.99];  nParEstimated=1;  chitmp = squeeze( chi2(1,1,:) );

chiCL     = mkChiFunc.GetChiCL( Prob, nParEstimated);
mSensStat = mkChiFunc.GetMassStatCL(chitmp,mtest_arr,Prob,nParEstimated);
mSensTot  = mkChiFunc.GetMassTotCL(chitmp,mtest_arr,Prob,nParEstimated);

fprintf('\n|    CL    | mnu_stat |  mnu_tot |\n');
mkChiFunc.print([Prob; mSensStat; mSensTot]);

figure(2);
plot(mtest_arr, chitmp); hold on;
plot(repmat([mtest_arr(1);mtest_arr(end)],1,3), [chiCL; chiCL]); hold off;
text(repmat(mtest_arr(end),1,3)/2.3,chiCL*1.1,strcat(num2str(Prob(:)*100),'\% CL'),'FontSize',12);
xlabel('tested mass $m_\mathrm{eff}$ in eV'); ylabel('goodness of fit $\chi^2$');
FigureFormat;


%% Calculate the mass sensitivity from chi2 for all sterile neutrino parameters
Prob=0.9;  nParEstimated=1;

meffCL=zeros(nsintest,nmnu4test); cnt=0; progressbar('Calculate the mass sensitivity');
for ls2=1:1:nsintest
    for lm4=1:1:nmnu4test
        VecC2 = squeeze(chi2(ls2,lm4,:));              
        meffCL(ls2,lm4) = mkChiFunc.GetMassTotCL(VecC2, mtest_arr, Prob, nParEstimated);   
        
        cnt=cnt+1;  progressbar(cnt/nmnu4test/nsintest); 
    end 
end

%% Plot the
figure(3);
[C, h] = contour(sinSq2Th4test, mnu4Sqtest, meffCL'*1e3); colorbar('vert'); 
ax=gca;  set(ax,'YScale','log');  % set(ax,'XScale','log');   
xlabel('$\sin^2 (2\theta_{14})$'); ylabel('$\Delta m_{14}^2 (\mathrm{eV}^2)$');
set(h,'Fill','on','LineWidth',1,'LineColor',[.5 .5 .5],'LineStyle',':');
clabel(C,h,'LabelSpacing',200,'FontSize',11,'Color',[.3 .3 .3]);

cb=colorbar('vert'); zlab=get(cb,'Label');  
set(zlab,'String','$m_\mathrm{eff}@90CL (\mathrm{meV})$','FontSize',12);
EditColormap('Scale',[.1 1]); grid off;

FigureFormat();

