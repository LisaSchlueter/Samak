%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Endpoint sensitivity study
% Detector as 1 pixel equivalent
% Configuration for First Tritium May

% P. Morales 2018
% Last update 16/03/2018


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initialization
clear;
addpath(genpath('../../../Samak2.0'));

%Parameters
nfit            = 1;
%runtimes        = (1)*60*60*24;
runtimes        = 210.2189781*42;
%Ntimes          = length(runtimes);
Ntimes          = 1;

nTeBinFac       = 5;
time_dist       = 'Flat30';
Fitter          = 'Matlab';
nPixels         = 148;

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
    'KTFFlag','Compute'}; %MACE+WGTSIS %SSCW_DCOfficial %SSCW_DCperqU

opt_bkg = {...
    'BKG_Flag',BCK};

opt_fsd= {'DTFSD',FSD};

%Initialize matrices
muhat = zeros(Ntimes,1);
sigmahat = zeros(Ntimes,1);
fit_xi2min = zeros(nfit,1);

opt_katrin = {...
    'TD',time_dist,...
    'TimeSec',runtimes,...
    'Mode',Mode};

A = InitKATRINE0_allPixels(...
    opt_katrin{:},...
    opt_wgtsmace{:},...
    opt_fsd{:},...
    opt_bkg{:},...
    opt_bin{:},...
    'mnuSq_i',mnuSq_t,...
    'nPixels',nPixels,...
    'qUmin',18575-200,'qUmax',18575+5,'qUStepSize',5);

fit_p = zeros(nfit,4*A.nPixels+6);

A.ComputeTBDDSallPixels();
A.ComputeTBDISallPixels(0.3/148*ones(1,nPixels));

if strcmp(FitFlag,'ON')  
    for ii = 1:Ntimes
        
        tic
        A.ComputeTBDDSallPixels();
        A.ComputeTBDISallPixels(0.3/148*ones(1,nPixels));
        A.AddStatFluctTBDISallPixels();
        Data = [...
            reshape(repmat(A.qU,1,A.nPixels),1,A.nqU*A.nPixels)' , ...
            reshape((A.TBDISallPixels),1,A.nqU*A.nPixels)', ...
            reshape(sqrt(A.TBDISallPixels),1,A.nqU*A.nPixels)'];
        % Init
        i_mnu      = 0;
        i_Q        = 0;
        i_B        = A.TBDISallPixels(end,:)/(A.qUfrac(1)*A.TimeSec);
        i_N        = zeros(1,A.nPixels);
        
        DoF = A.nqU*A.nPixels - length(i_B) - length(i_N) - 2;
        
        Afit = InitKATRINE0_allPixels(...
            opt_katrin{:},...
            opt_wgtsmace{:},...
            opt_fsd{:},...
            opt_bkg{:},...
            opt_bin{:},...
            'mnuSq_i',mnuSq_t,...
            'nPixels',nPixels,...
            'qUmin',18575-200,'qUmax',18575+5,'qUStepSize',5);
        Afit.ComputeTBDDSallPixels();
        Afit.ComputeTBDISallPixels();
        
        DataTBD = {Data,Afit};
        ParIni = [i_mnu i_Q i_B i_N];
        
        switch Fitter
            case 'Minuit'
                %                 % Initializing Fit
                %
                %                 parnames = ['mSq dQ B dN'];
                %                 % tmparg = sprintf(['set pri 1;'...
                %                 % 'fix  5 6 ; min; minos; imp; scan 1 ']);
                %                 tmparg = sprintf(['set pri -10;'...
                %                     'set now; min;']);
                %                 Args = {ParIni, DataTBD, '-c',tmparg};
                %                 fprintf(2,'No Detector Systematics - Stat.\n');
                %                 [par, err, chi2min, errmat] = fminuit('Chi2E0_first',Args{:});
                
                fprintf('Not working yet. Use Matlab as Fitter. \n');
            case 'Matlab'
                options = optimoptions('fminunc','Algorithm','quasi-newton',...
                    'OptimalityTolerance',1e-7,'StepTolerance',1e-8,...
                    'FiniteDifferenceType','central',...
                    'MaxFunctionEvaluations',1e6);
                TBDfun = @(xx) Chi2E0_first_allPixels(xx,DataTBD);
                [par,chi2min,exitflag,output,grad,Hessian] = fminunc(TBDfun,ParIni,options);
                if exitflag ~= 1
                    [par,chi2min,exitflag,output,grad,Hessian] = fminunc(TBDfun,par,options);
                end
                errmat = 0.5*Hessian;
                varcov = inv(errmat);
                err = sqrt(diag(varcov))';
        end
        
        switch display
            case 'ON'
                mnuSq_report = Afit.mnuSq_i+par(1);
                mnu_report = sqrt(abs(mnuSq_report));
                err_mnu = sqrt(abs((mnuSq_report + err(1)))) - mnu_report;
                BCK_fit = sum(par(3:2+Afit.nPixels));
                BCK_err = sqrt(sumsqr(par(3:2+Afit.nPixels)));
                N_ave = mean(abs(par(3+Afit.nPixels:3+2*Afit.nPixels-1)));
                N_err = sqrt(sumsqr(par(3+Afit.nPixels:3+2*Afit.nPixels-1)))/Afit.nPixels;
                fprintf(2,'----------------------------------------------\n');
                fprintf('  Processing t= %g %% \n',ii/nfit*100);
                fprintf('===============================================\n');
                fprintf('  m^2       = %g ± %g eV^2\n', mnuSq_report,err(1));
                fprintf('  m         = %g ± %g eV\n', mnu_report,err_mnu);
                fprintf('  dQ        = %g ± %g eV\n', par(2),err(2));
                fprintf('  B total   = %g ± %g cps\n', BCK_fit,BCK_err);
                fprintf('  dN average= %g ± %g\n', N_ave,N_err);
                fprintf('  Chi2/dof  = %g/%g\n', chi2min, DoF);
                fprintf('===============================================\n');
        end
        toc
        fit_p(ii,:) = [par err chi2min DoF];
        disp(ii);
    end
end

%save(['results\fit_allPixels_',time_dist,'_DPG.mat'],'fit_p');

%% plot

%     fit_full(nfit*(ii-1)+1:nfit*ii,:) = fit_p;
%     [muhat(ii),sigmahat(ii)] = normfit(fit_p(:,3));
%     RMS = int32(rms(fit_p(:,3)));
%     x = linspace(min(fit_p(:,3)),max(fit_p(:,3)),nfit);
%     figure(ii+10)
%     hold on
%     [t,Nbins,Wbins] = nhist(fit_p(:,3));
%     binWidth = Wbins(2) - Wbins(1);
%     plot(x,nfit*binWidth*normpdf(x,muhat(ii),sigmahat(ii)));
%     hold off
%     title(['E_0 fit ',time_dist,...
%         ' ',num2str(int32(runtimes(ii)/(60*60))),' h']);
%     xlabel('E_0 [eV]')
%     dim = [0.6 0.5 0.3 0.3];
%     str = {['N_{fit} = ',num2str(nfit)],['\sigma = ',num2str(int32(sigmahat(ii)*1e3)),' meV'],...
%         ['Time Dist. = ',time_dist],['Runtime = ',num2str(runtimes(ii)/(60*60)),'h'],['RMS = ',num2str(RMS*1e3),' [meV]']};
%     annotation('textbox',dim,'String',str,'FitBoxToText','on');
%
%     fprintf(1,'--------------------------------------------------------------\n');
%     fprintf(1,'sigma = %g [meV] \n', sigmahat(ii)*1e3);
%     fprintf(1,'--------------------------------------------------------------\n');

% save(['results\fit_full_',time_dist,'_DPG.mat'],'fit_full');
%
% figure(Ntimes+2);
% scatter(runtimes'/(60*60),sigmahat*1e3,'.b');
% title('Endpoint precision vs Time of measurment');
% ylabel('Uncertainty in Endpoint [meV]');
% xlabel('Time [hours]');

% 1 = norm; 3 = E0; 5 = mnu; 7 = Bck; 9 = Chi2;
