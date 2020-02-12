function fit_p = vftfit(run,nfit,sim)
%
% Study of impact of stacking runs
%
% Thierry Lasserre
% Last updated: 19/04/2018
%
addpath(genpath('../../../Samak2.0'));

  % Lower  Bin 9  = -200 V
    FitLB = 1;
    % Higher Bin 26 = + 40 V
    FitHB = 26;
    % Data
    
% Read datafile
datafile = sprintf('../../inputs/ftdata/%g.mat',run);
runfile  = load(datafile);

% MTD
reffile = sprintf('Run%g.mat',run);
mtd = load(reffile);
TD      = sprintf('Run%g',run);
TimeSec = runfile.TimeSec;

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
plotdist        = 'ON';
saveresults     = 'OFF';
savedist        = 'OFF';
%Initialize output matrices
fit_p = zeros(nfit,10);
muhat = zeros(Ntimes,1);
sigmahat = zeros(Ntimes,1);
fit_xi2min = zeros(nfit,1);

% Build Model
initrunfct = str2func(sprintf('ref_run%g',run));
Sref = ref_run40257('TD',TD,...
    'WGTS_CD_MolPerCm2',runfile.WGTS_CD_MolPerCm2,...
    'WGTS_MolFrac_DT',runfile.WGTS_MolFrac_DT,...
    'TimeSec',runfile.TimeSec);
Sref.ComputeTBDDS;Sref.ComputeTBDIS

if sim == 1
    Sref.AddStatFluctTBDIS;
    Count    = Sref.TBDIS; ErrCount = sqrt(Count); 
    DataRef  = [Sref.qU,Count,ErrCount];Data = DataRef;
else
    % Build Data VFT2
    if ((run == 40257) || (run == 40258))
        Count    = (runfile.TBDIS(3:end));
    else
        Count    = (runfile.TBDIS(1:end));
    end
    %Count    = Count.*Sref.TimeSec.*Sref.qUfrac; 
    ErrCount = sqrt(Count);
    ErrCount(4) = 99999;
    DataRef  = [Sref.qU,Count,ErrCount];Data = DataRef;
end


% Define Covariance Matrix - Signal
clear FitCovErrMatrixS ; matrixS = load('CM_VFT2ud_20percent2.mat');
FitCovErrMatrixS = ((Count .*  matrixS.covmatfrac .* Count'))+diag(ErrCount.^2);
%FitCovErrMatrixS = diag(ErrCount.^2);

% Define Covariance MAtrix - Background
clear FitCovErrMatrixB ; matrixB = load('covmatBfracRun4057.mat');
FitCovErrMatrixBfrac = matrixB.covmatBfrac;

% Sum 
ScaleB      = Sref.BKG_RateSec*Sref.qUfrac.*Sref.TimeSec;
%FitCovMatTB = FitCovErrMatrixS(FitLB:end,FitLB:end) + ( ScaleB(FitLB:end).*FitCovErrMatrixBfrac(FitLB:end,FitLB:end).*ScaleB(FitLB:end)');
%FitCovMatTB = FitCovErrMatrixS  + ( ScaleB .*FitCovErrMatrixBfrac.*ScaleB');
FitCovMatTB = FitCovErrMatrixS ;
%FitCovMatTB = diag(Count);


% Loop
progressbar('E0,m2,B,N Fits');
DoF = Sref.nqU - 4;
for kk = 1:nfit
    tic
    
    if sim == 1
        %Sref = initrunfct('TD',TD,'TimeSec',TimeSec,'WGTS_CD_MolPerCm2',5e17*0.9);
        SetFitBias(Sref,0); Sref.ComputeTBDDS; Sref.ComputeTBDIS;
        %Sref.AddStatFluctTBDIS; Count = Sref.TBDIS; ErrCount = sqrt(Count);
        Count = mvnrnd(Sref.TBDIS,FitCovErrMatrixS,1)';  ErrCount = sqrt(Count);
        DataRef  = [Sref.qU,Count,ErrCount];Data = DataRef;
    end
    
    % Init Fit Parameters
    i_mnu      = 0; i_Q        = 0;
    i_B        = 0;%Sref.TBDIS(end)./(Sref.TimeSec*Sref.qUfrac(end));
    i_N        = 0;
    i_DTGS     = 0; i_DTES       = 0;
    
    DataTBD = {DataRef,Sref,FitCovErrMatrixS,FitLB,FitHB,FitCovErrMatrixBfrac};
    
    switch Fitter
        case 'Minuit'
            % Initializing Fit
            ParIni = [i_mnu i_Q i_B i_N i_DTGS i_DTES];
            parnames = ['mSq dQ B dN i_DTGS i_DTES'];
            tmparg = sprintf(['set pri -10 ;'...
                'fix  5 6;'...
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
            fprintf('  DTGS      = %g +- %g cps\n', par(5),err(5));
            fprintf('  DTES      = %g +- %g\n', par(6),err(6));
            fprintf('  Chi2/dof  = %g/%g\n', chi2min, numel(Sref.qU(FitLB:FitHB))-6);
            fprintf('===============================================\n');
            
            figure(10)
            switch plotdist 
                case 'ON'
                fig = figure('Name','Fit and Residuals','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
                
                subplot(2,2,[1 2])
                hdata = errorbar(Data(FitLB:FitHB,1),Data(FitLB:FitHB,2),Data(FitLB:FitHB,3),...
                    'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
                hold on
                hdata2 = boundedline(Data(FitLB:FitHB,1),Data(FitLB:FitHB,2),sqrt(diag(FitCovMatTB(FitLB:FitHB,FitLB:FitHB))), 's','alpha');
                m = NuMassModel4par(par,Sref);
                hfit1 = plot(Data(FitLB:FitHB,1),m(FitLB:FitHB),...
                    'Color','Black','LineWidth',1,'LineStyle','-');
                hfit3 = line([Data(1,1) Data(end,1)],[(Sref.BKG_RateSec_i+par(3))*Sref.TimeSec*Sref.qUfrac(end) (Sref.BKG_RateSec_i+par(3))*Sref.TimeSec*Sref.qUfrac(end)],'LineStyle','--','Color','Red');
                hold off
                grid on
                xlabel('qU (eV)','FontSize',14);
                ylabel('Counts/qU','FontSize',14);
                set(gca,'FontSize',14);
                set(gca,'yscale','log');
                mydata = sprintf('Data: m=%.2f eV - E0=%.2f eV \n',...
                    sqrt(Sref.mnuSq_i),18575);
                myfit = sprintf(' Fit: \\chi2 / dof=%.2f/%.0f \n m^2=%.3f \t\\pm %.3f eV \n E0=%.3f \t\\pm %.3f eV \n B=%.3f \t\\pm %.3f \n N=%.3f \t\\pm %.3f',...
                    chi2min,numel(Sref.qU(FitLB:FitHB))-2,((Sref.mnuSq_i+par(1))),err(1),(par(2)),err(2),...
                    (Sref.BKG_RateSec_i+par(3))*1e3,err(3),par(4),err(4));
                mychi2 = sprintf('Data');
                legend([hdata  hfit1 hfit3],mychi2,myfit,'Offset','Location','NorthEast') ;
                axis([min(Sref.qU(FitLB:FitHB)) max(Sref.qU(FitLB:FitHB))+1 0.*min(Data(:,2)) max(Data(FitLB:FitHB,2))*1.2])
                title(sprintf(' %s Samak Data and Fit',...
                    Sref.TD),'FontSize',14);
                PrettyFigureFormat;
                
                subplot(2,2,[3 4])
                hdata = errorbar(Data(FitLB:FitHB,1),Data(FitLB:FitHB,2)-m(FitLB:FitHB),Data(FitLB:FitHB,3),...
                    'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
                hold on
                line([Data(FitLB,1) Data(FitHB,1)],[0 0 ],'LineStyle','--','Color','Blue');
                %            hdata2 = boundedline(Data(:,1),Data(:,2)-NuMassModel4par(par,Sref),sqrt(diag(FitCovErrMatrix)),'s','alpha');
                hdata2 = boundedline(Data(FitLB:FitHB,1),Data(FitLB:FitHB,1)-Data(FitLB:FitHB,1),sqrt(diag(FitCovMatTB(FitLB:FitHB,FitLB:FitHB))),'s','alpha');
                hdata = errorbar(Data(FitLB:FitHB,1),Data(FitLB:FitHB,2)-m(FitLB:FitHB),Data(FitLB:FitHB,3),...
                    'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
                %set(hdata2, 'linestyle', ':', 'color', 'r', 'marker', '.');
                hold off
                grid on
                xlabel('qU (eV)','FontSize',14);
                ylabel('Residuals','FontSize',14);
                set(gca,'FontSize',14);
                set(gca,'yscale','lin');
                PrettyFigureFormat;
            end
            
    end
    toc
    fit_p(kk,:) = [par(1:4) err(1:4) chi2min DoF];  fit_xi2min(kk) = chi2min;
    disp(kk);
    progressbar(kk/nfit);
end



