%Initialization
addpath(genpath('../../../Samak2.0'));

%Parameters
nfit = 500;

runtimes = 60*60*24;
Ntimes = 1;
npar = 6;
fign = 1;
display = 'ON';
Fitter = 'Minuit';
time_dist = 'FT1';

FitFlag         = 'ON';
savespectrum    = 'OFF';
plotresults     = 'ON';
plotdist        = 'OFF';
saveresults     = 'OFF';
savedist        = 'OFF';

%
%1. Define your Study Obect (TBD Object)

A= ref_fttl2();
A.ComputeTBDDS(); A.ComputeTBDIS('IStype','SIMP');

%%2. Define your CM Class with your number of trials, choose SysEffect,...
C = CovarianceMatrix('StudyObject',A,'nTrials',1500,...
    'RecomputeFlag','ON',...
    'SysEffect','RF',...
    'MACE_Ba_T_RelErr',1e-2,...
    'WGTS_B_T_RelErr',1e-2,...
    'WGTS_CD_MolPerCm2_RelErr',0.1,...
    'ISXsection_RelErr',1e-2);

%%3.(Optional) For RF: Turn Contribution to CovMat ON/OFF (default: all ON)
 C.BFieldFlag = 'ON';
 C.ELossFlag  = 'ON';
 C.RhoXsectionFlag = 'ON';
 C.RecomputeFlag = 'OFF';
 
%%4.Compute Covariance Matrix for your SysEffect ; here RF (Response Function)
C.ComputeCM_RF;
C.PlotCM;        %Plot the CM
disp(C.CovMatFile)
C.ComputeFracCM; %Compute the fractional CM
C.ComputeCMRank; %Compute the Rank of the CM as a first sanity check

%%Initialize output matrices
%mu         -0.112479  0.0124407
%sigma       0.278182  0.0088101
nfit = 500;
fit_p = zeros(nfit,10);
muhat = zeros(Ntimes,1);
sigmahat = zeros(Ntimes,1);
fit_xi2min = zeros(nfit,1);

% Define Covariance Matrix;
clear FitCovErrMatrix ; global FitCovErrMatrix;
cm=load(C.CovMatFile);
%FitCovErrMatrix = diag(A.TBDIS) + diag(diag(cm.CovMat));
FitCovErrMatrix = diag(A.TBDIS) + ((cm.CovMat));
%FitCovErrMatrix = diag(A.TBDIS);


% Define Model, again...

A= ref_fttl2(...
    'WGTS_CD_MolPerCm2',2.5e17*(1+0.*randn(1)),...
    'WGTS_B_T',2.52*(1+0.0*randn(1)),...
    'MACE_Bmax_T',4.2*(1+0.0*randn(1)),...
    'MACE_Ba_T',6e-4*(1+0.0*randn(1)));
A.ComputeTBDDS(); A.ComputeTBDIS('IStype','SIMP');

% Compute Data - Fluctuate the Data - Syst + Stat !
D=ref_fttl2(...
    'WGTS_CD_MolPerCm2',2.5e17*(1.03),...
    'WGTS_B_T',2.52*(0.994),...
    'MACE_Bmax_T',4.2*(1.003),...
    'MACE_Ba_T',6e-4*(0.997));

% Loop
progressbar('E0,m2,B,N Fits');
DoF = A.nqU - 4;

%%
for kk = 1:nfit
    tic
    pause(5)
    % Compute Data - Fluctuate the Data - Syst + Stat !
    
%     D=ref_fttl2(...
%         'WGTS_CD_MolPerCm2',2.5e17*(1+0.1*randn(1)),...
%         'WGTS_B_T',2.52*(1+0.01*randn(1)),...
%         'MACE_Bmax_T',4.2*(1+0.01*randn(1)),...
%         'MACE_Ba_T',6e-4*(1+0.01*randn(1)));


    % Compute & Stat Fluctuate Data
    D.ComputeTBDDS(); D.ComputeTBDIS('IStype','SIMP');
    D.AddStatFluctTBDIS();
    Data = [D.qU,D.TBDIS,D.TBDISE];
    
    % Init Fit Parameters
    i_mnu      = 0; i_Q        = 0;
    i_B        = D.TBDIS(end)./(D.TimeSec*D.qUfrac(end));
    i_N        = 0;
    i_m4       = 0; i_s4       = 0;
    
    DataTBD = {Data,A};
    
    switch Fitter
        case 'Minuit'
            % Initializing Fit
            ParIni = [i_mnu i_Q i_B i_N i_m4 i_s4];
            parnames = ['mSq dQ B dN m4 s4'];
            tmparg = sprintf(['set pri -10 ;'...
                'fix  5 6 ;'...
                'set now; min ; ' ]);
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
            mnuSq_report = A.mnuSq_i+par(1);
            mnu_report = sqrt(abs(mnuSq_report));
            err_mnu = sqrt(abs((mnuSq_report + err(1)))) - mnu_report;
            BCK_fit = A.BKG_RateSec_i+par(3);
            fprintf('===============================================\n');
            fprintf('  m^2       = %g +- %g eV^2\n', mnuSq_report,err(1));
            fprintf('  m         = %g +- %g eV\n', mnu_report,err_mnu);
            fprintf('  dQ        = %g +- %g eV\n', par(2),err(2));
            fprintf('  B         = %g +- %g cps\n', BCK_fit,err(3));
            fprintf('  dN        = %g +- %g\n', par(4),err(4));
            fprintf('  Chi2/dof  = %g/%g\n', chi2min, A.nqU-4);
            fprintf('===============================================\n');
            
            figure(10)
            %fig = figure('Name','Fit and Residuals','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
            
            subplot(2,2,[1 2])
            hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),...
                'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
            hold on
            hdata2 = boundedline(Data(:,1),Data(:,2),sqrt(diag(FitCovErrMatrix)), 's','alpha');
            hfit1 = plot(Data(:,1),NuMassModel4par(par,A),...
                'Color','Black','LineWidth',1,'LineStyle','-');
            hfit3 = line([Data(1,1) Data(end,1)],[(A.BKG_RateSec_i+par(3))*A.TimeSec*A.qUfrac(end) (A.BKG_RateSec_i+par(3))*A.TimeSec*A.qUfrac(end)],'LineStyle','--','Color','Red');
            hold off
            grid on
            xlabel('qU (eV)','FontSize',14);
            ylabel('Counts/qU','FontSize',14);
            set(gca,'FontSize',14);
            set(gca,'yscale','lin');
            mydata = sprintf('Data: m=%.2f eV - E0=%.2f eV \n',...
                sqrt(A.mnuSq_i),18575);
            myfit = sprintf(' Fit: \\chi2 / dof=%.1f/%.0f \n m=%.3f \t\\pm %.3f eV \n E0=%.3f \t\\pm %.3f eV \n B=%.3f \t\\pm %.3f \n N=%.3f \t\\pm %.3f',...
                chi2min,A.nqU-4,sqrt(abs(A.mnuSq_i+par(1))),err(1),(par(2)),err(2),...
                (A.BKG_RateSec_i+par(3))*1e3,err(3),par(4),err(4));
            mychi2 = sprintf('Data');
            legend([hdata  hfit1 hfit3],mychi2,myfit,'Offset','Location','NorthEast') ;
            axis([min(A.qU) max(A.qU)+1 0.*min(Data(:,2)) max(Data(:,2))*1.2])
            title(sprintf(' %s Data and Fit',...
                A.TD),'FontSize',14);
            
            subplot(2,2,[3 4])
            hdata = errorbar(Data(:,1),Data(:,2)-NuMassModel4par(par,A),Data(:,3),...
                'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
            hold on
            line([Data(1,1) Data(end,1)],[0 0 ],'LineStyle','--','Color','Blue');
            hdata2 = boundedline(Data(:,1),Data(:,2)-NuMassModel4par(par,A),sqrt(diag(FitCovErrMatrix)),'s','alpha');
            %set(hdata2, 'linestyle', ':', 'color', 'r', 'marker', '.');
            hold off
            grid on
            xlabel('qU (eV)','FontSize',14);
            ylabel('Residuals','FontSize',14);
            set(gca,'FontSize',14);
            set(gca,'yscale','lin');
            axis([min(A.qU) max(A.qU)+1 min(Data(:,2)-NuMassModel4par(par,A))*5 max(Data(:,2)-NuMassModel4par(par,A))*5]);
            
    end
    toc
    fit_p(kk,:) = [par(1:4) err(1:4) chi2min DoF];  fit_xi2min(kk) = chi2min;
    disp(kk);
    progressbar(kk/nfit);
end

%% plot
ii=1;
fit_full(nfit*(ii-1)+1:nfit*ii,:) = fit_p;
[muhat(ii),sigmahat(ii)] = normfit(fit_p(:,2));
RMS = int32(rms(fit_p(:,2)));
x = linspace(min(fit_p(:,2)),max(fit_p(:,2)),nfit);
figure(ii+10)
hold on
[t,Nbins,Wbins] = nhist(fit_p(:,2));
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
