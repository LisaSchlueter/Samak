
%Initialization
addpath(genpath('../../../Samak2.0'));

%Parameters
nfit = 500;
runtimes = [1 2 3]*60*60*24;
Ntimes = length(runtimes);
npar = 6;
fign = 1;
display = 'ON';
nTeBinFac = 8;
time_dist = 'Flat30';
Fitter = 'Minuit';
FSD = 'DOSS'; FSDadd = ''; FSDplot = 'only T2';

FitFlag         = 'ON';
savespectrum    = 'OFF';
plotresults     = 'ON';
plotdist        = 'OFF';
saveresults     = 'OFF';
savedist        = 'OFF';



% Parametrization: True Value
mnuSq_t = (0.00)^2;

% Init
opt_bin = {...
    'nTeBinningFactor',nTeBinFac};

opt_wgts = {...
    'WGTS_CosMaxAAngle',0.6324,...
    'WGTS_Tp',1,...
    'WGTS_DTHTr',1,...
    'WGTS_FTR_cm',4.1,...
    'WGTS_CD_MolPerCm2',5e17,...
    'WGTS_B_T',3.6};

opt_mace = {...
    'MACE_Bmax_T',6,...
    'MACE_Ba_T',9e-4,...
    'MACE_R_eV',2.79};

opt_wgtsmace = {...
    'KTFFlag','SSCW_DCperqU'}; %MACE+WGTSIS %SSCW_DCOfficial %SSCW_DCperqU

opt_fpd = {...
    'FPD_Segmentation','OFF',...
    'FPD_Pixel',0,'nPixels',148};

opt_bkg = {...
    'BKG_Flag','ON',...
    'BKG_Type','FLAT',...
    'BKG_RateAllFPDSec',0.3};

opt_fsd= {...
    'TTFSD','DOSS',...
    'DTFSD','OFF',...
    'HTFSD','OFF'};

opt_theocorr = {...
    'ScreeningFlag','OFF',...
    'FiniteExtChargeFlag','OFF',...
    'EEexchangeFlag','OFF',...
    'RecoilCoulombFlag','OFF',...
    'RadiativeFlag','OFF',...
    'RecoilWmVmAFlag','OFF',...
    'WintFiniteSizeFlag','OFF'};

%Initialize matrices
fit_p = zeros(nfit,10);
muhat = zeros(Ntimes,1);
sigmahat = zeros(Ntimes,1);
fit_xi2min = zeros(nfit,1);


progressbar('Runtime Fits','Neutrino Mass Fits');
for ii = 1:Ntimes
    %% 
    opt_katrin = {...
        'TD',time_dist,...
        'TimeSec',runtimes(ii)};
    
    A = TBD('Mode', 'Sim',...
        opt_katrin{:},...
        opt_wgts{:},...
        opt_mace{:},...
        opt_wgtsmace{:},...
        opt_fsd{:},...
        opt_fpd{:},...
        opt_bkg{:},...
        opt_bin{:},...
        opt_theocorr{:},...
        'mnuSq_i',mnuSq_t);
    
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
                [par, err, chi2min, errmat] = fminuit('NuMassChi2E0',Args{:});
                
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
        disp(ii); disp(kk);
        progressbar(ii/Ntimes,kk/nfit);
    end
    
    %% plot
    
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
        ' ',num2str(int32(runtimes(ii)/(60*60))),' h']);
    xlabel('E_0 [eV]')
    dim = [0.6 0.5 0.3 0.3];
    str = {['N_{fit} = ',num2str(nfit)],['\sigma = ',num2str(int32(sigmahat(ii)*1e3)),' meV'],...
        ['Time Dist. = ',time_dist],['Runtime = ',num2str(runtimes(ii)/(60*60)),'h'],['RMS = ',num2str(RMS*1e3),' [meV]']};
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
    fprintf(1,'--------------------------------------------------------------\n');
    fprintf(1,'sigma = %g [meV] \n', sigmahat(ii)*1e3);
    fprintf(1,'--------------------------------------------------------------\n');
end

    save(['results\fit_full_',time_dist,'_DPG.mat'],'fit_full');

figure(Ntimes+2);
scatter(runtimes'/(60*60),sigmahat*1e3,'.b',);
title('Endpoint precision vs Time of measurment');
ylabel('Uncertainty in Endpoint [meV]');
xlabel('Time [hours]');

% 1 = norm; 3 = E0; 5 = mnu; 7 = Bck; 9 = Chi2;
