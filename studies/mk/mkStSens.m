%% Sensitivity study concering sterile neutrinos in the eV-mass range at KATRIN
%
% June 2017 - M. Korzeczek - KIT, IEKP
%

%% Initialization
clearvars; clear global; clear figure; %clear workspace
addpath('fminuit'); addpath('tools');  %Paths
%FigureFormat; %1st call sets latex interpreter

%global variables
clear Theory; global Theory;

%control variables
Fitter = 'Minuit';   %'Matlab','Minuit'
CovMat = 'OFF';      %'OFF','ON'  ... available if Fitter='Minuit'

%Calculation grid
sinSq2Th4test_arr = logspace(-2,0,21);  sinSqTh4test_arr = ( 1 - sqrt(1-sinSq2Th4test_arr) )/ 2;
mnu4Sqtest_arr    = logspace(-2,2,21);  mnu4test_arr     = sqrt(mnu4Sqtest_arr);

runtime = 5;  %in years
Rbkg = 10e-3; %in cps
SpecOpt = {...
    'TD','DR30','TimeYear',runtime,'nTeBinningFactor',20,...
    'BKG_RateAllFPDSec',Rbkg,'TTFSD','SAENZ','mnuSq_i',0};

%define theoretic spectrum
Theory=InitKATRIN(SpecOpt{:});
switch Fitter
    case 'Minuit'
        %Build covariance matrix
        switch CovMat
            case 'ON'
                nCMTrials = 100;
                fName = sprintf('./data/CovMat/WGTSMACE_CovMat_%gTrials_%s_Year_%g.mat',...
                    nCMTrials,Theory.TD,Theory.TimeYear);
                if exist(fName, 'file') == 2
                    CM=importdata(fName);
                else
                    Theory.ComputeTBDISCovarianceMatrix('nTrials',nCMTrials);
                    CM=importdata(fName);
                end
        end
    case 'Matlab'
        %define theory
        tmp = @(mSqb,qb,bb,nb,m4Sqb,s2T4b,qu) Theory.ComputeTBDIS4f(qu,mSqb,qb,bb,nb,0,0,m4Sqb,s2T4b);
        TBDISLIf = @(mSqb,qb,bb,nb,m4Sqb,s2T4b,e) interp1(Theory.qU,tmp(mSqb,qb,bb,nb,m4Sqb,s2T4b,Theory.qU),e);
        
        Theoryf = fittype(@(mSqf,qf,bf,nf,m4Sqf,sSqT4f,qu) TBDISLIf(mSqf,qf,bf,nf,m4Sqf,sSqT4f,qu),...
            'independent','qu','coefficients', {'mSqf','qf','bf','nf','m4Sqf','sSqT4f'});
end

%% Fitting loop
nsintest  = length(sinSqTh4test_arr);  nmnu4test = length(mnu4Sqtest_arr);  npar = 6;

par  = zeros(nsintest,nmnu4test,npar);  err  = zeros(nsintest,nmnu4test,npar);
chi2 = zeros(nsintest,nmnu4test);

tic; cnt=-1;
for ls2=1:1:nsintest
    for lm4=1:1:nmnu4test
        cnt=cnt+1;
        %define measured data
        TBDData=InitKATRIN(SpecOpt{:},...
            'sin2T4_i',sinSqTh4test_arr(ls2),'mnu4Sq_i',mnu4Sqtest_arr(lm4));
        Data = [TBDData.qU,TBDData.TBDIS,TBDData.TBDISE];
        % Set Fitter dependent fit procedure
        switch Fitter
            case 'Minuit'
                %Fit Options
                %                mSqB   QB    RbkgB  NormB  m4SqB  sinSqTh4B
                ParIni =   [      0      0      0      0      0        0     ];
                ParFix = 'fix     1                           5        6     ;';
                Args = {ParIni, Data, '-c',[ParFix,'set pri -10 ; set now; min']};
                %Start fit
                switch CovMat
                    case 'OFF'
                        [ fpar, ferr, chi2min, ~ ] = fminuit('mkChiFunc.calcChi2',Args{:});
                    case 'ON'
                        switch SystCorr
                            case 'OFF'
                                Theory.fitcov = diag(diag(CM.WGTSMACE_CovMat)) + diag(Theory.TritSpecData.TBDIS);
                            case 'ON'
                                Theory.fitcov = (CM.WGTSMACE_CovMat) + diag(Theory.TritSpecData.TBDIS);
                        end
                        [ fpar, ferr, chi2min, ~ ] = fminuit('mkChiFunc.calcChi2CovMat',Args{:});
                end
                par(ls2,lm4,:)=fpar; err(ls2,lm4,:)=ferr; chi2(ls2,lm4)=chi2min;
                
            case 'Matlab'
                %Fit options
                ParIni = {'StartPoint',[0 0 0 0 0 0]};
                ParLim = {'lower',[0 -5 -10e-3 -1 0 0],'upper',[0 5 1 1 0 0]};
                opts = fitoptions('Display','off','Robust','Off',ParIni{:},ParLim{:},...
                    'Method', 'NonlinearLeastSquares','Weights',1./Data(:,2));
                %Start fitting
                [fobj,gof,output] = fit(Data(:,1),Data(:,2),Theoryf,opts);
                par(ls2,lm4,:)=[fobj.mSqf fobj.qf fobj.bf fobj.nf fobj.m4Sqf fobj.sSqT4f];
                chi2(ls2,lm4)=gof.sse;  out(ls2,lm4)=output;
        end
        %print result
        fprintf(2,'--------------------------------------------------------------\n');
        fprintf(2,'  Processing = %g %% -- m4Sq=%g eV^2 \t sinSq2T=%g\n',...
            cnt/nmnu4test/nsintest*100,mnu4Sqtest_arr(lm4),sinSq2Th4test_arr(ls2));
        fprintf(2,'--------------------------------------------------------------\n');
        fprintf(2,'  dmnuSq   \t= %g ±\t%g\tmeV^2 \n', par(ls2,lm4,1)*1e6, err(ls2,lm4,1)*1e6);
        fprintf(2,'  dQ       \t= %g ±\t%g\tmeV \n', par(ls2,lm4,2)*1e3, err(ls2,lm4,2)*1e3);
        fprintf(2,'  dB       \t= %g ±\t%g\tmcps \n',par(ls2,lm4,3)*1e3, err(ls2,lm4,3)*1e3);
        fprintf(2,'  dN       \t= %g ±\t%g\n',       par(ls2,lm4,4), err(ls2,lm4,4));
        fprintf(2,'  dmnu4Sq  \t= %g ±\t%g\tmeV^2\n',  par(ls2,lm4,5)*1e6, err(ls2,lm4,5)*1e6);
        fprintf(2,'  dsinSqTh4\t= %g ±\t%g\n',       par(ls2,lm4,6), err(ls2,lm4,6));
        fprintf(2,'  Chi2 \t= %g / %g dof \n',       chi2(ls2,lm4),TBDData.nqU-4);
        fprintf(2,'--------------------------------------------------------------\n');
    end
end

toc

%% Standard plots
TBDData.PlotTD('fign',10);
%TBDData.PlotTBDDS('fign',11);
%TBDData.PlotTBDIS('fign',12);
%TBDData.PlotKTF('fign',13);
%TBDData.PlotFSD('fign',14,'TT','ON');
%TBDData.PlotPIS('fign',15);

%% Plot last fit vs data
figure(1);
x=Data(:,1)/1e3; norm=InitKATRIN(SpecOpt{:}); norm=norm.TBDIS;
errorbar(x,Data(:,2)./norm,Data(:,3)./norm,'k+'); hold on;
switch Fitter
    case 'Minuit'
        plot(x,mkChiFunc.calcModel(par(nsintest,nmnu4test,:))./norm); hold off;
    case 'Matlab'
        plot(x,Theoryf(par(nsintest,nmnu4test,:),Data(:,1))./norm); hold off;
end
title( num2str( par(nsintest,nmnu4test,:) ) );  legend('data/norm','theory/norm');
xlabel('electron kinetic energy in keV'); ylabel('ratio in a.u.');


%% Combined sensitivity contours KATRIN & CeSOX 
Prob=.9; CL=mkChiFunc.myInvChi2CdfRec(Prob,2); % Define contour lines %only single contour lines
%[X,Y] = meshgrid(sinSq2Th4test_arr,  mnu4Sqtest_arr);
 [X,Y] = meshgrid(logspace(-2,0,101), logspace(-1.3,2,101));  %Define fine grid for contour plot

%Organize chi2maps: KATRIN (from calculation above)          
%c2fKstat = @(xq,yq) interp2(sinSq2Th4test_arr,mnu4Sqtest_arr,chi2',xq,yq,'linear',0);
%zkstat = c2fKstat(X,Y);  zkstat = zkstat-min(min(zkstat));

mnu4SqTot = sqrt( mnu4Sqtest_arr.^2 + 0.017^2 *mkChiFunc.CL2Sig(Prob)^2 ) ;
c2fK  = @(xq,yq) interp2(sinSq2Th4test_arr,mnu4SqTot,chi2',xq,yq,'linear',0);
zk = c2fK(X,Y);  zk = zk-min(min(zk)); %added systematic uncertainty

%Organize chi2maps: CeSOX (from precalculated file by T. Lasserre)
t=load('Data/CeSOX/CeSOXsmall_100x100_chi2map500days.mat','dm2','sin22theta','chi2map');  
c2fCS = @(xq,yq) interp2(t.sin22theta,t.dm2,t.chi2map(:,:,1),xq,yq,'linear',0);
zc = c2fCS(X,Y);  zc = zc-min(min(zc));

%Plot contours for KATRIN, CeSOX and the combined case
figure(3); 
contour(X,Y, zk+zc,[CL CL],'k-','LineWidth',2); hold on;
contour(X,Y, zc   ,[CL CL],'g:','LineWidth',2);
contour(X,Y, zk   ,[CL CL],'b:','LineWidth',2); hold off;
ax=gca;  set(ax,'XScale','log'); set(ax,'YScale','log'); grid on; grid minor;

legend('Combined','CeSOX','KATRIN');
title(sprintf('%.2f y - %.2f mHz - %2.0f\\%%CL',runtime,Rbkg*1e3,Prob*100));
xlabel('$\sin^2 2\theta ()$');  ylabel('$\Delta m_\mathrm{new}^2 (\mathrm{eV}^2)$');


FigureFormat;

