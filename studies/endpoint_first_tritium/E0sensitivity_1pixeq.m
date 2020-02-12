%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Endpoint sensitivity study
% Detector as 1 pixel equivalent
% Configuration for First Tritium May

% P. Morales 2018
% Last update 18/03/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initialization
%clear;
tstart = tic;
addpath(genpath('../../../Samak2.0'));
%Parameters
nfit            = 25;
runtimes        = 210.2189781*42;
%runtimes        = ([1 3 5 7 9 2 4 6 8 10])*60*60*24;
Ntimes          = length(runtimes);

nTeBinFac       = 5;
time_dist       = 'Flat30';
Fitter          = 'Minuit';

FSD             = 'DOSS';
FSDadd = ''; FSDplot = 'only DT';

BCK             = 'ON';

Mode            = 'Sim';

display         = 'ON';
FitFlag         = 'OFF';
savespectrum    = 'OFF';
plotresults     = 'OFF';
plotdist        = 'OFF';
saveresults     = 'OFF';
savedist        = 'OFF';


% Parametrization: True Value
mnuSq_t = (0.00)^2;

% Init
opt_bin = {...
    'nTeBinningFactor',nTeBinFac};

opt_wgtsmace = {...
    'KTFFlag','Compute'};

opt_bkg = {...
    'BKG_Flag',BCK,...
    'BKG_RateAllFPDSec',0.3};

opt_fsd= {'DTFSD',FSD};

%Initialize matrices
fit_p = zeros(nfit,10);
muhat = zeros(Ntimes,1);
sigmahat = zeros(Ntimes,1);
fit_xi2min = zeros(nfit,1);

%progressbar('Neutrino Mass Fits','Runtime Fits');
for ii = 1:Ntimes
    
    opt_katrin = {...
        'TD',time_dist,...
        'TimeSec',runtimes(ii),...
        'Mode',Mode};
    
    A = InitKATRINE0_1pixeq(...
        opt_katrin{:},...
        opt_wgtsmace{:},...
        opt_fsd{:},...
        opt_bkg{:},...
        opt_bin{:},...
        'mnuSq_i',mnuSq_t,...
        'qUmin',18575-200,'qUmax',18575+5,'qUStepSize',5);
    
    A.ComputeTBDDS();
    A.ComputeTBDIS();
    
    
    DoF = A.nqU - 4;
    if strcmp(FitFlag,'ON')
        for kk = 1:nfit
            tic
            % Create fake data
            A.ComputeTBDDS();
            A.ComputeTBDIS();
            A.AddStatFluctTBDIS();
            Data = [A.qU,A.TBDIS,A.TBDISE];
            
            % Create Model
            Afit = InitKATRINE0_1pixeq(...
                opt_katrin{:},...
                opt_wgtsmace{:},...
                opt_fsd{:},...
                'BKG_Flag','ON',...
                opt_bin{:},...
                'mnuSq_i',mnuSq_t,...
                'qUmin',18575-200,'qUmax',18575+5,'qUStepSize',5);
            Afit.ComputeTBDDS();
            Afit.ComputeTBDIS();
            
            % Init
            i_mnu      = 0;
            i_Q        = 0;
            i_B        = Data(end,2)/(A.qUfrac(1)*A.TimeSec);
            i_N        = 0;
            %i_N = Data(1,2)/(A.TBDIS(1)/(A.qUfrac(1)*A.TimeSec)+i_B) - 1;
            
            DataTBD = {Data,Afit};
            
            switch Fitter
                case 'Minuit'
                    % Initializing Fit
                    ParIni = [i_mnu i_Q i_B i_N];
                    parnames = ['mSq dQ B dN'];
                    % tmparg = sprintf(['set pri 1;'...
                    % 'fix  5 6 ; min; minos; imp; scan 1 ']);
                    tmparg = sprintf(['set pri -10;'...
                        'set now; min;']);
                    Args = {ParIni, DataTBD, '-c',tmparg};
                    [par, err, chi2min, errmat] = fminuit('Chi2E0_first',Args{:});
                    
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
                    mnuSq_report = Afit.mnuSq_i+par(1);
                    mnu_report = sqrt(abs(mnuSq_report));
                    err_mnu = sqrt(abs((mnuSq_report + err(1)))) - mnu_report;
                    BCK_fit = Afit.BKG_RateSec_i+par(3);
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
            progressbar(kk/nfit,ii/Ntimes);
        end
    end
    
    %% plot
    if strcmp(plotdist,'ON')
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
            ' ',num2str(int32(runtimes(ii)/(60*60*24))),' days']);
        xlabel('E_0 [eV]')
        dim = [0.6 0.5 0.3 0.3];
        str = {['N_{fit} = ',num2str(nfit)],['\sigma = ',num2str(int32(sigmahat(ii)*1e3)),' meV'],...
            ['Runtime = ',num2str(runtimes(ii)/(60*60)),'h']};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        
        fprintf(1,'--------------------------------------------------------------\n');
        fprintf(1,'sigma = %g [meV] \n', sigmahat(ii)*1e3);
        fprintf(1,'--------------------------------------------------------------\n');
        print(['results/E0',num2str(runtimes(ii)/(60*60*24)),'d'],'-dpng')
    end
end

if strcmp(saveresults,'ON'); save(['results\fit_full_',time_dist,'_DPG.mat'],'fit_full'); end

%% savage plots
if strcmp(plotresults,'ON')
    figure(Ntimes+101);
    scatter(runtimes'/(60*60*24),sigmahat*1e3,'k','filled');
    title('Endpoint precision vs Time of measurment');
    ylabel('Uncertainty in Endpoint [meV]');
    xlabel('Time [days]');
    ttotal = toc(tstart);
end
% 1 = norm; 3 = E0; 5 = mnu; 7 = Bck; 9 = Chi2;
