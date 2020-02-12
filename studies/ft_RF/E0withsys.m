
%Initialization
addpath(genpath('../../../Samak2.0'));

%Parameters
nfit = 1000;
runtimes = 60*60*24;
Ntimes = 1;
npar = 6;
fign = 1;
display = 'ON';
Fitter = 'Minuit';
time_dist = 'Flat1000';

FitFlag         = 'ON';
savespectrum    = 'OFF';
plotresults     = 'ON';
plotdist        = 'OFF';
saveresults     = 'OFF';
savedist        = 'OFF';

A=ref_tfcm();
A.ComputeTBDDS(); A.ComputeTBDIS('IStype','SIMP');

%Initialize matrices
fit_p = zeros(nfit,10);
muhat = zeros(Ntimes,1);
sigmahat = zeros(Ntimes,1);
fit_xi2min = zeros(nfit,1);

% Covariance Matrix;
global FitCovErrMatrix;
cm=importdata('WGTSMACE_CovMat_200Trials_Simulation_Sec_86400.mat');
%FitCovErrMatrix = diag(A.TBDIS) + diag(diag(cm.WGTSMACE_CovMat));
FitCovErrMatrix = diag(A.TBDIS) + ((cm.WGTSMACE_CovMat));
%FitCovErrMatrix = diag(A.TBDIS);

% Loop
progressbar('Neutrino Mass Fits');    
    DoF = A.nqU - 4;
    
    for kk = 1:nfit
        tic
        A.ComputeTBDDS(); A.ComputeTBDIS('IStype','SIMP');
        A.AddStatFluctTBDIS();
        Data = [A.qU,A.TBDIS,A.TBDISE];
        
        % Init
        i_mnu      = 0; i_Q        = 0;
        i_B        = A.TBDIS(end)./(A.TimeSec*A.qUfrac(end));
        i_N        = 0;
        i_m4       = 0; i_s4       = 0;
        
        DataTBD = {Data,A};
        
        switch Fitter
            case 'Minuit'
                % Initializing Fit
                ParIni = [i_mnu i_Q i_B i_N i_m4 i_s4];
                parnames = ['mSq dQ B dN m4 s4'];
                % tmparg = sprintf(['set pri 1;'...
                % 'fix  5 6 ; min; minos; imp; scan 1 ']);
                tmparg = sprintf(['set pri -10 ;'...
                    'fix  5 6 ;'...
                    'set now; min ; ' ]);
                Args = {ParIni, DataTBD, '-c',tmparg};
                fprintf(2,'No Detector Systematics - Stat.\n');
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
                mnuSq_report = A.mnuSq_i+par(1);
                mnu_report = sqrt(abs(mnuSq_report));
                err_mnu = sqrt(abs((mnuSq_report + err(1)))) - mnu_report;
                BCK_fit = A.BKG_RateSec_i+par(3);
                fprintf('===============================================\n');
                fprintf('  m^2       = %g ± %g eV^2\n', mnuSq_report,err(1));
                fprintf('  m         = %g ± %g eV\n', mnu_report,err_mnu);
                fprintf('  dQ        = %g ± %g eV\n', par(2),err(2));
                fprintf('  B         = %g ± %g cps\n', BCK_fit,err(3));
                fprintf('  dN        = %g ± %g\n', par(4),err(4));
                fprintf('  Chi2/dof  = %g/%g\n', chi2min, A.nqU-4);
                fprintf('===============================================\n');
        end
        toc
        fit_p(kk,:) = [par(1:4) err(1:4) chi2min DoF];  fit_xi2min(kk) = chi2min;
        disp(kk);
        progressbar(kk/nfit);
    end
    
    %% plot
    ii=1
    fit_full(nfit*(ii-1)+1:nfit*ii,:) = fit_p;
    [muhat(ii),sigmahat(ii)] = normfit(fit_p(:,3));
    RMS = int32(rms(fit_p(:,3)));
    x = linspace(min(fit_p(:,3)),max(fit_p(:,3)),nfit);
    figure(ii+10)
    hold on
    [t,Nbins,Wbins] = nhist(fit_p(:,3));
    binWidth = Wbins(2) - Wbins(1);
    plot(x,nfit*binWidth*normpdf(x,muhat(ii),sigmahat(ii)));
    hold off
    title(['E_0 fit ',time_dist,...
        ' ',num2str(int32(runtimes/(60*60))),' h']);
    xlabel('E_0 [eV]')
    dim = [0.6 0.5 0.3 0.3];
    str = {['N_{fit} = ',num2str(nfit)],['\sigma = ',num2str(int32(sigmahat(ii)*1e3)),' meV'],...
        ['Time Dist. = ',time_dist],['Runtime = ',num2str(runtimes/(60*60)),'h'],['RMS = ',num2str(RMS*1e3),' [meV]']};
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
    fprintf(1,'--------------------------------------------------------------\n');
    fprintf(1,'sigma = %g [meV] \n', sigmahat(ii)*1e3);
    fprintf(1,'--------------------------------------------------------------\n');

    save(['results\fit_full_',time_dist,'_DPG.mat'],'fit_full');

figure(Ntimes+2);
scatter(runtimes'/(60*60),sigmahat*1e3,'.b');
title('Endpoint precision vs Time of measurment');
ylabel('Uncertainty in Endpoint [meV]');
xlabel('Time [hours]');

% 1 = norm; 3 = E0; 5 = mnu; 7 = Bck; 9 = Chi2;
