%% Variant 2) Maximal Estimation
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

%% general initialization

clearvars; clear global; delete(findall(0,'Type','figure'));
addpath('fminuit'); addpath('tools');  %Paths

%Spectrum options inherited from TBD()
CommonOpt = {'TD','DR30','TimeYear',3,'nTeBinningFactor',20, ...
    'BKG_RateAllFPDSec',10e-3,'TTFSD','SAENZ'};
%DataOpt   = [CommonOpt,{}]; %options inherited from TBD()
TheoryOpt = [CommonOpt,{}]; %options inherited from TBD()


%% Single marginalization with standard parameters
MinimizationArg={...
    'PlotsShow','ON','Fitter','Minuit',...
    'CeSoxType','OFF',... %'OFF','5','50','100','500'
    ... %         mSqB QB    RbkgB NormB m4SqB sinSqTh4B
    'FitParInit',[0    0     0     0     0     .3], ...
    'FitParFix', [true false false false false false] , ...
    ...
    'FitParUseLimits','ON', ...
    'FitParLower',[-10  -5     -10e-3 -1   -500    0] , ...
    'FitParUpper',[10   5      1      1    500   .5]};

mtest_arr = [0,linspace(0.15,1.1,2)];  nmtest = length(mtest_arr); %Calculation masses for fit procedure
npar = 6;  par = ones(nmtest,npar);    err= zeros(nmtest,npar);   chi2 = ones(nmtest,1);

progressbar('maximal estimation');
for lnm=1:1:nmtest
    progressbar(lnm/nmtest);
    DataOpt   = [CommonOpt,{'mnuSq_i',mtest_arr(lnm)^2}];
    [chi2t,part,errt,TBDData] = mkChiFunc.minimizeTheoryToData(DataOpt, TheoryOpt, MinimizationArg{:});
    par(lnm,:)=part; err(lnm,:)=errt; chi2(lnm)=chi2t;
end

%% -- chi2 values and corresponding mass sensitivities
Prob=[.90,.95,.99];
nParEstimated=1;

chiCL     = mkChiFunc.GetChiCL( Prob, nParEstimated);
mSensStat = mkChiFunc.GetMassStatCL(chi2,mtest_arr,Prob,nParEstimated);
mSensTot  = mkChiFunc.GetMassTotCL(chi2,mtest_arr,Prob,nParEstimated);

fprintf('\n|    CL    | mnu_stat |  mnu_tot |\n');
mkChiFunc.print([Prob; mSensStat; mSensTot]);


%% -- Standard plots
TBDData.PlotTD('fign',10);
%TBDData.PlotTBDDS('fign',11);
%TBDData.PlotTBDIS('fign',12);
%TBDData.PlotKTF('fign',13);
%TBDData.PlotFSD('fign',14,'TT','ON');
%TBDData.PlotPIS('fign',15);



%% Multiple Marginalizations (scan runtime & CeSOX)
mtest_arr = [0,linspace(0.15,1.1,2)];  nmtest = length(mtest_arr);
runtime=[.1,.2,.3,.4,.5,.6,.9, 1, 1.5, 2, 3]; nruntime=length(runtime);
Prob = .9; 

% - marginalization for different KATRIN & CeSOX runtimes
MyFitType = {'nominal','OFF','5','50','100','500'};  nMyFitType = length(MyFitType);
progressbar('Multiple Marginalizations'); cnt=0;
for typet = MyFitType
    type=typet{:};
    if(strcmp(type,'nominal')); name=type; MinimizationArg={'PlotsShow','ON','Fitter','Minuit'};
    else
        name=['CeSOX_',type,'_tot'];
        MinimizationArg={...
            'PlotsShow','ON','Fitter','Minuit',...
            'CeSoxType',type,... %'OFF','5','50','100','500'
            ... %         mSqB QB    RbkgB NormB m4SqB sinSqTh4B
            'FitParInit',[0    0     0     0     0     .3], ...
            'FitParFix', [true false false false false false] , ...
            ...
            'FitParUseLimits','ON', ...
            'FitParLower',[-10  -5     -10e-3 -1   -500    0] , ...
            'FitParUpper',[10   5      1      1    500   .5]};
    end
    
    tmpSens=[];
    %different runtime
    for time = runtime
        %marginalization
        for lnm=1:1:nmtest
            cnt=cnt+1; progressbar(cnt/nMyFitType/nruntime/nmtest);
            
            DataOpt   = [CommonOpt,{'TimeYear',time,'mnuSq_i',mtest_arr(lnm)^2}];
            TheoryOpt = [CommonOpt,{'TimeYear',time}];
            [chi2t,part,errt,TBDData] = mkChiFunc.minimizeTheoryToData(DataOpt, TheoryOpt,MinimizationArg{:});
            par(lnm,:)=part; err(lnm,:)=errt; chi2(lnm)=chi2t;
        end
        %mSensStat = mkChiFunc.GetMassStatCL(chi2,mtest_arr,Prob,nParEstimated); 
        mSensTot  = mkChiFunc.GetMassTotCL(chi2,mtest_arr,Prob,nParEstimated); 
        tmpSens =[tmpSens, mSensTot];
    end
    mEffSens.(name)  = tmpSens;
end

save('matfiles/mkEffSens2Mult.mat');


%% Plot result of multiple marginalization
figure(4)
plot(runtime, mEffSens.nominal_tot,'k-','LineWidth',3); hold on;
plot(runtime, mEffSens.CeSOX_OFF_tot,'c-','LineWidth',1.5);
plot(runtime, mEffSens.CeSOX_500_tot,'g--','LineWidth',1.5);
plot(runtime, mEffSens.CeSOX_100_tot,'r-.','LineWidth',1.5);
plot(runtime, mEffSens.CeSOX_50_tot,'b:','LineWidth',1.5);
plot(runtime, mEffSens.CeSOX_5_tot,'m--','LineWidth',1.5);
hold off;

legend('nominal','marginalization','w CeSOX 500d','w CeSOX 100d','w CeSOX 50d','w CeSOX 5d');
xlabel('measurement time $(\mathrm{y})$'); ylabel('$m_\mathrm{eff}@90CL (\mathrm{eV}^2)$');
%ax=gca;  set(ax,'YScale','log'); %set(ax,'YScale','log');
FigureFormat();