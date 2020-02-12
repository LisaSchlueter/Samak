%
% Study of impact of stacking runs
%
% Thierry Lasserre
% Last updated: 13/04/2018
%
addpath(genpath('../../../Samak2.0'));

% Read baseline TD and Create TD for Nruns
%<<<<<<< HEAD
TD      = 'FT-TL4';
nRuns   = 4;
TimeAll = 21150*30;
% =======
% TD      = 'FT2';
% nRuns   = 350;
% 
% TimeAll = 17400*350;
% >>>>>>> 0f8f2b7f526ff1b28d2e3b676974242c7abe9b64
createTDruns('ClearTDruns','OFF','TDi',TD,'nRuns',nRuns,'qUoffset',0.15,'Display','ON');
%createTDruns('ClearTDruns','ON','TDi',TD,'nRuns',nRuns,'qUoffset',0.15,'Display','ON');return;


% Generate Reference TBD 
%if ~exist('Sref','var')
    Sref = ref_stackrun('TD',TD,'TimeSec',TimeAll);
    Sref.ComputeTBDDS;Sref.ComputeTBDIS
    fprintf(2,'Creating Reference TBD, Time=%g\n',TimeAll);
    disp(num2str(Sref.qU'));
    disp(num2str(Sref.qUfrac'));
    disp(num2str(Sref.TBDIS'/nRuns));
%end

%% Generate TBD for each runs, scale time by 1/nRuns
for i=1:1:nRuns
    tdname = sprintf('%s_run%g.mat',TD,i);
    %if ~exist('S','var') || numel(S) < i+1
        S(i)   = ref_stackrun('TD',tdname,'TimeSec',TimeAll/nRuns);
        fprintf(2,'Creating TDB object for run %g, TimeRun=%g \n',i,TimeAll/nRuns);
        S(i).ComputeTBDDS;S(i).ComputeTBDIS
        disp(num2str(S(i).qU'));
        disp(num2str(S(i).qUfrac'));
        disp(num2str(S(i).TBDIS'));
    %end
end

%% Computing and stacking Spectra for nRuns
clear StackTBDIS ; StackTBDIS = zeros(Sref.nqU,1);
for i=1:1:nRuns
StackTBDIS = StackTBDIS + S(i).TBDIS;
end

%% Plot Reference TBDIS and StackTBDIS
% Assuming the nominal qU - ignoring fluctuations
figure(100)
subplot(2,1,1)
Pref = errorbar(Sref.qU-Sref.Q,Sref.TBDIS,Sref.TBDISE,'d','MarkerSize',5,'MarkerFaceColor',.8*[1 1 1],'LineWidth',1);
hold on
Pstack = errorbar(Sref.qU-Sref.Q,StackTBDIS,Sref.TBDISE,'d','MarkerSize',5,'MarkerFaceColor',.5*[1 1 1],'LineWidth',1);
hold off
grid on
set(gca, 'YScale', 'log');
xlabel('E-E_0 (eV)','FontSize',14);
ylabel('counts','FontSize',14);
PrefStr   = sprintf('Reference Spectrum - %g seconds',TimeAll);
PstackStr = sprintf('Stack of %g runs x %g seconds',nRuns,TimeAll/nRuns);
a = legend([Pref,Pstack],PrefStr,PstackStr);
title('6h MTD - 150 mV dispersion (asymmetric and uniform)','FontSize',14)
set(gca,'FontSize',14);
PrettyFigureFormat;
% subplot(3,1,2)
% Pres = plot(Sref.qU-Sref.Q,((Sref.TBDIS./StackTBDIS)-1)*100,'ks','MarkerSize',7,'MarkerFaceColor',.8*[1 1 1],'LineWidth',1);
% xlabel('E-E_0 (eV)','FontSize',14);
% grid on
% ylabel('Rel. Difference (%)','FontSize',14);
% PrettyFigureFormat;
subplot(2,1,2)
Pres = plot(Sref.qU-Sref.Q,(Sref.TBDIS-StackTBDIS)./Sref.TBDISE,'ks','MarkerSize',5,'MarkerFaceColor',.5*[1 1 1],'LineWidth',1);
xlabel('E-E_0 (eV)','FontSize',14);
grid on
ylabel('Difference / \sigma_{stat}','FontSize',14);
PrettyFigureFormat;

%% Fit Data with Reference Spectrum
%  Data = StackTBDIS
%  Proxy = Endpoint - Marginalize over all other parameters
%Parameters
nfit = 1;
runtimes = TimeAll;
Ntimes = 1;
npar = 6;
fign = 1;
display = 'ON';
Fitter = 'Minuit';
time_dist = TD;
FitFlag         = 'ON';
savespectrum    = 'OFF';
plotresults     = 'ON';
plotdist        = 'OFF';
saveresults     = 'OFF';
savedist        = 'OFF';
%Initialize output matrices
fit_p = zeros(nfit,10);
muhat = zeros(Ntimes,1);
sigmahat = zeros(Ntimes,1);
fit_xi2min = zeros(nfit,1);
% Model
Sref = ref_stackrun('TD',TD,'TimeSec',TimeAll);
Sref.ComputeTBDDS;Sref.ComputeTBDIS
% Define Covariance Matrix;
clear FitCovErrMatrix ;  global FitCovErrMatrix;
clear SystCovErrMatrix ; global SystCovErrMatrix;
% Compute Data - No Statistical Fluctuations
clear DataRef;    DataRef   = [Sref.qU,Sref.TBDIS,sqrt(Sref.TBDIS)];
clear DataStack ; DataStack = [Sref.qU,StackTBDIS,sqrt(StackTBDIS)];
clear Data ; Data      = DataStack;
%SystCovErrMatrix=diag(sqrt(0.3)*Sref.TBDIS);
SystCovErrMatrix=diag((5e-4)^2*Sref.TBDIS.^2);
FitCovErrMatrix = diag(Sref.TBDIS);%+SystCovErrMatrix;

% Loop
progressbar('E0,m2,B,N Fits');
DoF = Sref.nqU - 4;
for kk = 1:nfit
    tic
    
    % Init Fit Parameters
    i_mnu      = 0; i_Q        = 0;
    i_B        = 0;%Sref.TBDIS(end)./(Sref.TimeSec*Sref.qUfrac(end));
    i_N        = 0;
    i_m4       = 0; i_s4       = 0;
    
    DataTBD = {Data,Sref};
    
    switch Fitter
        case 'Minuit'
            % Initializing Fit
            ParIni = [i_mnu i_Q i_B i_N i_m4 i_s4];
            parnames = ['mSq dQ B dN m4 s4'];
            tmparg = sprintf(['set pri -10 ;'...
                'fix 5 6 ;'...
                'set now; minos ; imp ; ' ]);
            Args = {ParIni, DataTBD, '-c',tmparg};
            [par, err, chi2min, errmat] = fminuit('NuMassChi2E0_Cov',Args{:});
            
        case 'Matlab'
            options = optimoptions('fminunc','Algorithm','quasi-newton',...
                'OptimalityTolerance',1e-7,'StepTolerance',1e-7,'FiniteDifferenceType','central');
            TBDfun = @(xx) NuMassChi2E0(xx,DataL);
            [par,chi2min,exitflag,output,grad,Hessian] = fminunc(TBDfun,ParIni,options);
            if exitflag ~= 1
                [par,chi2min,exitflag,output,grad,Hessian] = fminunc(TBDfun,par,options);
            end
            errmat = 0.5*Hessian;
            varcov = inv(errmat);
            err = sqrt(diag(varcov));
    end
    
    switch display
        case 'ON'
            fprintf(2,'----------------------------------------------\n');
            fprintf('  Processing \t= %g %% \n',kk/nfit*100);
            mnuSq_report = Sref.mnuSq_i+par(1);
            mnu_report = sqrt(abs(mnuSq_report));
            err_mnu = sqrt(abs((mnuSq_report + err(1)))) - mnu_report;
            BCK_fit = Sref.BKG_RateSec_i+par(3);
            fprintf('===============================================\n');
            fprintf('  m^2       = %g +- %g eV^2\n', mnuSq_report,err(1));
            fprintf('  m         = %g +- %g eV\n', mnu_report,err_mnu);
            fprintf('  dQ        = %g +- %g eV\n', par(2),err(2));
            fprintf('  B         = %g +- %g cps\n', BCK_fit,err(3));
            fprintf('  dN        = %g +- %g\n', par(4),err(4));
            fprintf('  Chi2/dof  = %g/%g\n', chi2min, Sref.nqU-4);
            fprintf('===============================================\n');
            
            figure(10)
            %fig = figure('Name','Fit and Residuals','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
            
            subplot(2,2,[1 2])
            hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),...
                'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
            hold on
            hdata2 = boundedline(Data(:,1),Data(:,2),sqrt(diag(FitCovErrMatrix)), 's','alpha');
            hfit1 = plot(Data(:,1),NuMassModel4par(par,Sref),...
                'Color','Black','LineWidth',1,'LineStyle','-');
            hfit3 = line([Data(1,1) Data(end,1)],[(Sref.BKG_RateSec_i+par(3))*Sref.TimeSec*Sref.qUfrac(end) (Sref.BKG_RateSec_i+par(3))*Sref.TimeSec*Sref.qUfrac(end)],'LineStyle','--','Color','Red');
            hold off
            grid on
            xlabel('qU (eV)','FontSize',14);
            ylabel('Counts/qU','FontSize',14);
            set(gca,'FontSize',14);
            set(gca,'yscale','lin');
            mydata = sprintf('Data: m=%.2f eV - E0=%.2f eV \n',...
                sqrt(Sref.mnuSq_i),18575);
            myfit = sprintf(' Fit: \\chi2 / dof=%.2f/%.0f \n m^2=%.3f \t\\pm %.3f eV \n E0=%.3f \t\\pm %.3f eV \n B=%.3f \t\\pm %.3f \n N=%.3f \t\\pm %.3f',...
                chi2min,Sref.nqU-4,((Sref.mnuSq_i+par(1))),err(1),(par(2)),err(2),...
                (Sref.BKG_RateSec_i+par(3))*1e3,err(3),par(4),err(4));
            mychi2 = sprintf('Data');
            legend([hdata  hfit1 hfit3],mychi2,myfit,'Offset','Location','NorthEast') ;
            axis([min(Sref.qU) max(Sref.qU)+1 0.*min(Data(:,2)) max(Data(:,2))*1.2])
            title(sprintf(' %s Data and Fit',...
                Sref.TD),'FontSize',14);
            
            subplot(2,2,[3 4])
            hdata = errorbar(Data(:,1),Data(:,2)-NuMassModel4par(par,Sref),Data(:,3),...
                'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
            hold on
            line([Data(1,1) Data(end,1)],[0 0 ],'LineStyle','--','Color','Blue');
            hdata2 = boundedline(Data(:,1),Data(:,2)-NuMassModel4par(par,Sref),sqrt(diag(FitCovErrMatrix)),'s','alpha');
            %set(hdata2, 'linestyle', ':', 'color', 'r', 'marker', '.');
            hold off
            grid on
            xlabel('qU (eV)','FontSize',14);
            ylabel('Residuals','FontSize',14);
            set(gca,'FontSize',14);
            set(gca,'yscale','lin');
            %axis([min(Sref.qU) max(Sref.qU)+1 -max(Data(:,3)) max(Data(:,3))]);
            %axis([min(Sref.qU) max(Sref.qU)+1 -2 2]);
            
    end
    toc
    fit_p(kk,:) = [par(1:4) err(1:4) chi2min DoF];  fit_xi2min(kk) = chi2min;
    disp(kk);
    progressbar(kk/nfit);
end



