%Initialization
addpath(genpath('../../../Samak2.0'));
close all

% Run Parameters
runtimes  = 21150*4;    % Run Time
mtd       = 'FT-TL4';   % MTD

% Fits 
nfit   = 500;            % Per Trials 
npar   = 6;              % Number of Fit Parameters
Fitter = 'Minuit';
clear FitCovErrMatrix ; global FitCovErrMatrix;
    
% Tritium Column Density Case studied
WGTS_CD_MolPerCm2_RelErr = 0.1; % Fluctuation
WGTS_CD_MolPerCm2_local = 2.5e17*[0.25 0.5 1]'; %,2.5e17
Ntimes = numel(WGTS_CD_MolPerCm2_local);
% General Parameters
fign = 1;
display = 'ON';
FitFlag         = 'ON';
savespectrum    = 'ON';
plotresults     = 'ON';
plotdist        = 'OFF';
saveresults     = 'ON';
savedist        = 'ON';

% Initialize output matrices (E0 Fit Results / Column Density Trials)
muhat = zeros(Ntimes,1); sigmahat = zeros(Ntimes,1);

% Define  Study Object (TBD Object)

CovMat_rho = cell(Ntimes,1);

for ii =1:1:Ntimes 
    A=ref_fttl4('WGTS_CD_MolPerCm2', WGTS_CD_MolPerCm2_local(ii)); A.ComputeTBDDS(); A.ComputeTBDIS('IStype','SIMP');    
    fit_p = zeros(nfit,10);
    fit_xi2min = zeros(nfit,1);
    
    % Compute RF Covariance Matrix
    myEffects = struct(...
                 'RF_EL','ON',... % Response Function(RF) EnergyLoss
                 'RF_BF','ON',... % RF B-Fields
                 'RF_RX','ON',... % RF Column Density, Cross Section
                 'TASR','OFF',... % Tritium Column Density Stack Runs);
                  'FSD','OFF');
%    C = CovarianceMatrix('StudyObject',A, 'nTrials',1000,...
%       'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr,...%WGTS_CD_MolPerCm2_RelErr(ii),...
%        'RunFlag','OFF', 'RecomputeFlag','OFF');
%    C.ComputeCM_RF;
%    CovMat_rho{ii,1} = C.MultiCovMat.CM_RF;
    % Define Fit Covariance Matrix : Stat. + Syst. 
    FitCovErrMatrix = diag(A.TBDIS);% + ((C.MultiCovMat.CM_RF));
    
    % Define Model, again...
    A=ref_fttl4('WGTS_CD_MolPerCm2', WGTS_CD_MolPerCm2_local(ii)); A.ComputeTBDDS(); A.ComputeTBDIS('IStype','SIMP'); 

    % Loop of fits
    progressbar('E0,m2,B,N Fits');
    DoF = A.nqU - 4;  
    D = ref_fttl4('WGTS_CD_MolPerCm2', A.WGTS_CD_MolPerCm2*(1+WGTS_CD_MolPerCm2_RelErr));
    for kk = 1:nfit
        % Compute & Stat Fluctuate Data
     %  D = ref_ft2('WGTS_CD_MolPerCm2', A.WGTS_CD_MolPerCm2*(1+randn*WGTS_CD_MolPerCm2_RelErr));%(ii)));   % Compute Data - Fluctuate the Data: Sys     
        D.ComputeTBDDS(); D.ComputeTBDIS('IStype','SIMP');      
        D.AddStatFluctTBDIS;                                                                                            % Fluctuate the Data: stat
        Data = [D.qU,D.TBDIS,sqrt(D.TBDIS)];
       
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
                    'set now; minos ; ' ]);
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
                fprintf('  Chi2/DoF  = %g/%g\n', chi2min, A.nqU-4);
                fprintf('===============================================\n');
                
                switch plotdist
                    case 'ON'
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
                        mydata = sprintf('Data: m=%.2f eV - (E0)eff=%.2f eV \n',...
                            sqrt(A.mnuSq_i),18575);
                        myfit = sprintf(' Fit: \\chi2 / DoF=%.1f/%.0f \n m=%.3f \t\\pm %.3f eV \n (E0)eff=%.3f \t\\pm %.3f eV \n B=%.3f \t\\pm %.3f \n N=%.3f \t\\pm %.3f',...
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
                
        end
        fit_p(kk,:) = [par(1:4) err(1:4) chi2min DoF];  fit_xi2min(kk) = chi2min;
        progressbar(kk/nfit);
    end
    
    %% plot Results
    fit_full(nfit*(ii-1)+1:nfit*ii,:) = fit_p;
    eofit = fit_p(:,2);
    [muhat(ii),sigmahat(ii)] = normfit(eofit(fit_xi2min<3*DoF));
    RMS = round(rms(fit_p(:,2)),3);
    x = linspace(min(eofit(fit_xi2min<3*DoF)),max(eofit(fit_xi2min<3*DoF)),nfit);
    figure(ii+10)
    hold on
    [t,Nbins,Wbins] = nhist(eofit(fit_xi2min<3*DoF));
    binWidth = Wbins(2) - Wbins(1);
    plot(x,nfit*binWidth*normpdf(x,muhat(ii),sigmahat(ii)),'Color','Red','LineWidth',2,'LineStyle','-');
    hold off
    title(['(E0)_{eff} - ',mtd,...
        ' - ',num2str(int32(runtimes/(60*60))),' h - ',...
        ' Column Density =',num2str(100*WGTS_CD_MolPerCm2_local(ii)/(2.5e17)),'%']);
    xlabel('(E0)_{eff} - 18575 [eV]')
    ylabel('')
    dim = [0.7 0.6 0.3 0.3];
    str = {['N_{fit} = ',num2str(nfit)],['\mu =',num2str(int32(muhat(ii)*1e3))],['\sigma = ',num2str(int32(sigmahat(ii)*1e3)),' meV'],...
        ['Time Dist. = ',mtd],['Runtime = ',num2str(runtimes/(60*60)),'h'],['RMS = ',num2str(RMS*1e3),' [meV]']};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','red');
    PrettyFigureFormat
    
    figure(ii+20) %chi2 plot
    title(['(E0)_{eff} - ',mtd,...
        ' - ',num2str(int32(runtimes/(60*60))),' h - ',...
        ' Column Density = ',num2str(100*WGTS_CD_MolPerCm2_local(ii)/(2.5e17)), '%']);
    nhist(fit_xi2min(fit_xi2min<3*DoF));
    xlabel('$\mathbf{\chi^2}$/14 dof','Interpreter','latex','FontSize',18)
    ylabel('')
    dim = [0.7 0.6 0.3 0.3];
    str = {['N_{fit} = ',num2str(nfit)],['<\chi^2> = ',num2str(round(mean(fit_xi2min(fit_xi2min<3*DoF)),2))],...
        ['Time Dist. = ',mtd],['Runtime = ',num2str(runtimes/(60*60)),'h'],['RMS = ',num2str(rms(fit_xi2min(fit_xi2min<3*DoF)))]};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','red');
    PrettyFigureFormat
    
    
    fprintf(2,'--------------------------------------------------------------\n');
    fprintf(2,'sigma = %g [meV] \n', sigmahat(ii)*1e3);
    fprintf(2,'--------------------------------------------------------------\n');
    
    %save(['./results/fit_full_',mtd,'_WGTS.mat'],'fit_full','WGTS_CD_MolPerCm2_RelErr','muhat','sigmahat');
    fig1str = sprintf('./results/statE0eff_%s_%gs_%gfits_rho%g_%.2ferr_fit.eps',mtd,runtimes,nfit,WGTS_CD_MolPerCm2_local(ii),WGTS_CD_MolPerCm2_RelErr);
    fig11str = sprintf('./results/statE0eff_%s_%gs_%gfits_rho%g_%.2ferr_fit.fig',mtd,runtimes,nfit,WGTS_CD_MolPerCm2_local(ii),WGTS_CD_MolPerCm2_RelErr);
    fprintf('publishing %s ...\n',fig1str);
    publish_figure(ii+10,fig1str);
    export_fig(fig11str);
    fig2str = sprintf('./results/statE0eff_%s_%gs_%gfits_rho%g_%.2ferr_chi2.eps',mtd,runtimes,nfit,WGTS_CD_MolPerCm2_local(ii),WGTS_CD_MolPerCm2_RelErr);
    fig21str = sprintf('./results/statE0eff_%s_%gs_%gfits_rho%g_%.2ferr_chi2.fig',mtd,runtimes,nfit,WGTS_CD_MolPerCm2_local(ii),WGTS_CD_MolPerCm2_RelErr);
    publish_figure(ii+20,fig2str);
    export_fig(fig21str);
    fprintf('publishing %s ...\n',fig2str);

end

%% Results
fprintf(2,'--------------------------------------------------------------\n');
disp('               TA%                bias                sigma')
fprintf(2,'--------------------------------------------------------------\n');
disp([100*WGTS_CD_MolPerCm2_local/(2.5e17),muhat*1e3,sigmahat*1e3]);
fprintf(2,'--------------------------------------------------------------\n');

figure(999)
subplot(2,1,1)
strtitle = sprintf('(E0)eff - MTD=%s - RunTime=%.1f h - %g fits (No CM)',mtd,runtimes/(60*60),nfit);
errorbar(100*WGTS_CD_MolPerCm2_local/(2.5e17),muhat*1e3,sigmahat*1e3,'s','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','blue','LineWidth',2,'Color','Black');
grid on
xlabel('Nominal Column Density (%)')
ylabel('(E0)_{eff} bias (meV)')
grid on
set(gca,'FontSize',16);
title(strtitle);
PrettyFigureFormat

subplot(2,1,2)
strtitle = sprintf('(E0)eff - MTD=%s - RunTime=%.1f h - %g fits',mtd,runtimes/(60*60),nfit);
stairs(100*WGTS_CD_MolPerCm2_local/(2.5e17),sigmahat*1e3,'s','MarkerSize',1,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',2);
grid on
%xlabel('Column Density Fluctuation (%) ')
xlabel('Nominal Column Density (%)')
ylabel('(E0)_{eff}   uncertainty (meV)')
grid on
set(gca,'FontSize',16);
hold on
%xx = min(WGTS_CD_MolPerCm2_RelErr*100):0.1:max(WGTS_CD_MolPerCm2_RelErr*100);
%yy = spline(WGTS_CD_MolPerCm2_RelErr*100,sigmahat*1e3,xx);
xx = min(100*WGTS_CD_MolPerCm2_local/(2.5e17)):0.1:max(100*WGTS_CD_MolPerCm2_local/(2.5e17));
yy = spline(100*WGTS_CD_MolPerCm2_local/(2.5e17),sigmahat*1e3,xx);
plot(100*WGTS_CD_MolPerCm2_local/(2.5e17),sigmahat*1e3,'s',xx,yy,'LineWidth',2,'Color','black')
hold off
PrettyFigureFormat
fig999str = sprintf('./results/statE0eff_%s_%gs_%gfits_RFfluctuation_columndensity.eps',mtd,runtimes,nfit);
fig9991str = sprintf('./results/statE0eff_%s_%gs_%gfits_RFfluctuation_columndensity.fig',mtd,runtimes,nfit);
publish_figure(999,fig999str);
export_fig(fig9991str);
save('./results/statfite0_RFCovMat_columndensities_Workspace.mat'); %save whole workspace
