function fit_p = mynewfit(run,nfit)
%
% Study of impact of stacking runs
%
% Thierry Lasserre
% Last updated: 19/04/2018
%
addpath(genpath('../../../Samak2.0'));

% MTD
%run = 40257;

reffile = sprintf('Run%g.mat',run);
mtd = load(reffile);
TD      = mtd.TD;
TimeSec = mtd.RunTime;

%Parameters
%nfit = 1;
runtimes = TimeSec;
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

% Build Model
initrunfct = str2func(sprintf('ref_run%g',run));
Sref = initrunfct('TD',TD,'TimeSec',TimeSec,'WGTS_CD_MolPerCm2',5e17*1.1);
Sref.ComputeTBDDS;Sref.ComputeTBDIS

% Build Data VFT2
datafile = sprintf('data/spectrum_%g_outerRingsExcl_correctTrigger.txt',run);
runfile  = importdata(datafile);
if mod(run,2) == 0
    Count    = runfile(1:end,2);
else
    if run == 40257
        Count    = flip(runfile(1:end-2,2));
    else
        Count    = flip(runfile(1:end,2));
        
    end
end
Count    = Count.*Sref.TimeSec.*Sref.qUfrac; ErrCount = sqrt(Count);
DataRef  = [Sref.qU,Count,ErrCount];Data = DataRef;

% Define Covariance Matrix;
clear FitCovErrMatrix ; matrix = load('CM_VFT2ud_20percent2.mat');
%FitCovErrMatrix = ((Count .*  matrix.covmatfrac .* Count'))+diag(Count);
FitCovErrMatrix = diag(Count);

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
    
    % Lower  Bin 9  = -200 V
    FitLB = 1;
    % Higher Bin 26 = + 40 V
    FitHB = 26;
    % Data
    DataTBD = {DataRef,Sref,FitCovErrMatrix,FitLB,FitHB};
    
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
            hdata = errorbar(Data(FitLB:FitHB,1),Data(FitLB:FitHB,2),Data(FitLB:FitHB,3),...
                'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
            hold on
            hdata2 = boundedline(Data(FitLB:FitHB,1),Data(FitLB:FitHB,2),sqrt(diag(FitCovErrMatrix(FitLB:FitHB,FitLB:FitHB))), 's','alpha');
            m = NuMassModel4par(par,Sref);
            hfit1 = plot(Data(FitLB:FitHB,1),m(FitLB:FitHB),...
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
            axis([min(Sref.qU(FitLB:FitHB)) max(Sref.qU(FitLB:FitHB))+1 0.*min(Data(:,2)) max(Data(FitLB:FitHB,2))*1.2])
            title(sprintf(' %s Data and Fit',...
                Sref.TD),'FontSize',14);
            PrettyFigureFormat;
            
            subplot(2,2,[3 4])
            hdata = errorbar(Data(FitLB:FitHB,1),Data(FitLB:FitHB,2)-m(FitLB:FitHB),Data(FitLB:FitHB,3),...
                'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
            hold on
            line([Data(FitLB,1) Data(FitHB,1)],[0 0 ],'LineStyle','--','Color','Blue');
%            hdata2 = boundedline(Data(:,1),Data(:,2)-NuMassModel4par(par,Sref),sqrt(diag(FitCovErrMatrix)),'s','alpha');
            hdata2 = boundedline(Data(FitLB:FitHB,1),Data(FitLB:FitHB,1)./Data(FitLB:FitHB,1),sqrt(diag(FitCovErrMatrix(FitLB:FitHB,FitLB:FitHB))),'s','alpha');

            %set(hdata2, 'linestyle', ':', 'color', 'r', 'marker', '.');
            hold off
            grid on
            xlabel('qU (eV)','FontSize',14);
            ylabel('Residuals','FontSize',14);
            set(gca,'FontSize',14);
            set(gca,'yscale','lin');
            PrettyFigureFormat;
            
    end
    toc
    fit_p(kk,:) = [par(1:4) err(1:4) chi2min DoF];  fit_xi2min(kk) = chi2min;
    disp(kk);
    progressbar(kk/nfit);
end



