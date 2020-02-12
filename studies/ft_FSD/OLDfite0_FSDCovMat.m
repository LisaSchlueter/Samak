%Initialization
addpath(genpath('../../../Samak2.0'));
close all

savepath = sprintf('./FitResults/');
% Fits 
nfit   = 500;            % Per Trials 
npar   = 6;              % Number of Fit Parameters
Fitter = 'Minuit';
clear FitCovErrMatrix ; global FitCovErrMatrix;
    

% General Parameters
fign = 1;
display = 'ON';
FitFlag         = 'ON';
savespectrum    = 'ON';
plotresults     = 'ON';
plotdist        = 'OFF';
saveresults     = 'ON';
savedist        = 'ON';
Ntimes = 1;
% Initialize output matrices (E0 Fit Results / Column Density Trials)
muhat = zeros(Ntimes,1); sigmahat = zeros(Ntimes,1);

% Define  Study Object (TBD Object)
A=InitKatrin_ft_FSD(); A.ComputeTBDDS(); A.ComputeTBDIS('IStype','SIMP');
CovMat_FSD = cell(Ntimes,1);

    fit_p = zeros(nfit,10);
    fit_xi2min = zeros(nfit,1);
    
    % Compute RF Covariance Matrix
    myEffects = struct(...
                 'RF_EL','OFF',... % Response Function(RF) EnergyLoss
                 'RF_BF','OFF',... % RF B-Fields
                 'RF_RX','OFF',... % RF Column Density, Cross Section
                 'TASR','OFF',... % Tritium Column Density Stack Runs);
                 'FSD', 'OFF',...
                 'BM1S','ON',...
                 'TCoff_RAD','OFF',...  % Radiative corrections
                 'TCoff_OTHER','OFF');
      
    C = CovarianceMatrix('StudyObject',A, 'nTrials',1000, 'SysEffect',myEffects,...
       'FSDNorm_RelErr',0.01,'FSDShape_RelErr',0.05,...
        'RunFlag','OFF', 'RecomputeFlag','ON');
    C.ComputeCM_BM1S;
    C.PlotCM;
    C.DecomposeCM;
    
    C.ComputeCM_FSD;
    CovMat_FSD{1,1} = C.MultiCovMat.CM_FSD;
    
    C.ComputeCM_STAT;
    % Define Fit Covariance Matrix : Stat. + Syst.
    FitCovErrMatrix = (C.MultiCovMat.CM_STAT) + (C.MultiCovMat.CM_FSD);
    
    % Define Model, again...
    A=ref_fttl4(); A.ComputeTBDDS(); A.ComputeTBDIS('IStype','SIMP');
    
    % Run Parameters
    runtimes  = A.TimeSec;    % Run Time
    mtd       = A.TD;   % MTD
    % Loop of fits
    progressbar('E0,m2,B,N Fits');
    DoF = A.nqU - 4;
    %D = ref_fttl4();
    %NormBias = zeros(Ntimes,1);
    D = ref_fttl4();   % Compute Data - Fluctuate the Data: Sys
    NormBias = D.DTNormGS_i*C.FSDNorm_RelErr;
    DT_P = zeros(numel(D.DTexP),Ntimes);
    DT_P(:,1) = D.DTexP'.*(1+randn(numel(D.DTexP),1).*C.FSDShape_RelErr); 
    D.DTexP = DT_P(:,1)';
    for kk = 1:nfit
        % Compute & Stat Fluctuate Data               
        D.ComputeTBDDS('DTGS_bias',NormBias,'DTES_bias', -NormBias); 
        D.ComputeTBDIS('IStype','SIMP');      
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
    %fit_full(nfit*(ii-1)+1:nfit*ii,:) = fit_p;
    eofit = fit_p(:,2);
    [muhat,sigmahat] = normfit(eofit(fit_xi2min<3*DoF));
    RMS = round(rms(fit_p(:,2)),3);
    x = linspace(min(eofit(fit_xi2min<3*DoF)),max(eofit(fit_xi2min<3*DoF)),nfit);
    figure(10)
    hold on
    [t,Nbins,Wbins] = nhist(eofit(fit_xi2min<3*DoF));
    binWidth = Wbins(2) - Wbins(1);
    plot(x,nfit*binWidth*normpdf(x,muhat,sigmahat),'Color','Red','LineWidth',2,'LineStyle','-');
    hold off
    title(['(E0)_{eff} - ',mtd,...
        ' - ',num2str(int32(runtimes/(60*60))),' h - ',...
        'FSD','Normalization and Shape']);
    xlabel('(E0)_{eff} - 18575 [eV]')
    ylabel('')
    dim = [0.7 0.6 0.3 0.3];
    str = {['N_{fit} = ',num2str(nfit)],['\sigma = ',num2str(int32(sigmahat*1e3)),' meV'],...
        ['Time Dist. = ',mtd],['Runtime = ',num2str(runtimes/(60*60)),'h'],['RMS = ',num2str(RMS*1e3),' [meV]']};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','red');
    PrettyFigureFormat
    
    figure(20) %chi2 plot
    title(['(E0)_{eff} - ',mtd,...
        ' - ',num2str(int32(runtimes/(60*60))),' h - ',...
        ' DT FSD ','Normalization and Shape']);
    nhist(fit_xi2min(fit_xi2min<3*DoF));
    xlabel('$\mathbf{\chi^2}$/14 dof','Interpreter','latex','FontSize',18)
    ylabel('')
    dim = [0.7 0.6 0.3 0.3];
    str = {['N_{fit} = ',num2str(nfit)],['<\chi^2> = ',num2str(round(mean(fit_xi2min(fit_xi2min<3*DoF)),2))],...
        ['Time Dist. = ',mtd],['Runtime = ',num2str(runtimes/(60*60)),'h'],['RMS = ',num2str(rms(fit_xi2min(fit_xi2min<3*DoF)))]};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','red');
    PrettyFigureFormat
    
    
    fprintf(2,'--------------------------------------------------------------\n');
    fprintf(2,'sigma = %g [meV] \n', sigmahat*1e3);
    fprintf(2,'--------------------------------------------------------------\n');
    
 
    fig1str_name = sprintf('./results/E0eff_%s_%gs_%gfits_FSD_fit.eps',mtd,runtimes,nfit);
    fig11str_name = sprintf('./results/E0eff_%s_%gs_%gfits_FSD_fit.fig',mtd,runtimes,nfit);
    fig1str = strcat(savepath, fig1str_name);
    fig11str = strcat(savepath, fig11str_name);
    fprintf('publishing %s ...\n',fig1str);
    publish_figure(10,fig1str);
    export_fig(fig11str);
    fig2str_name = sprintf('./results/E0eff_%s_%gs_%gfits_FSD.eps',mtd,runtimes,nfit);
    fig21str_name = sprintf('./results/E0eff_%s_%gs_%gfits_FSD.fig',mtd,runtimes,nfit);
    fig2str = strcat(savepath,fig2str_name);
    fig21str = strcat(savepath, fig21str_name);
    publish_figure(20,fig2str);
    export_fig(fig21str);
    fprintf('publishing %s ...\n',fig2str);

%% Results
% fprintf(2,'--------------------------------------------------------------\n');
% disp('               TA%                bias                sigma')
% fprintf(2,'--------------------------------------------------------------\n');
% disp([WGTS_CD_MolPerCm2_RelErr,muhat*1e3,sigmahat*1e3]);
% fprintf(2,'--------------------------------------------------------------\n');
% 
% figure(999)
% subplot(2,1,1)
% strtitle = sprintf('(E0)eff - MTD=%s - RunTime=%g - %g fits',mtd,runtimes,nfit);
% errorbar(100*WGTS_CD_MolPerCm2_local/(2.5e17),muhat*1e3,sigmahat*1e3,'s','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','blue','LineWidth',2,'Color','Black');
% grid on
% %xlabel('Column Density Fluctuation (%)')
% xlabel('Nominal Column Density (%)')
% ylabel('(E0)_{eff}   bias/std (meV)')
% grid on
% set(gca,'FontSize',16);
% title(strtitle);
% PrettyFigureFormat
% 
% subplot(2,1,2)
% strtitle = sprintf('(E0)eff - MTD=%s - RunTime=%.1f h - %g fits',mtd,runtimes/(60*60),nfit);
% stairs(100*WGTS_CD_MolPerCm2_local/(2.5e17),sigmahat*1e3,'s','MarkerSize',1,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',2);
% grid on
% %xlabel('Column Density Fluctuation (%) ')
% xlabel('Nominal Column Density (%)')
% ylabel('(E0)_{eff}   uncertainty (meV)')
% grid on
% set(gca,'FontSize',16);
% hold on
% %xx = min(WGTS_CD_MolPerCm2_RelErr*100):0.1:max(WGTS_CD_MolPerCm2_RelErr*100);
% %yy = spline(WGTS_CD_MolPerCm2_RelErr*100,sigmahat*1e3,xx);
% xx = min(100*WGTS_CD_MolPerCm2_local/(2.5e17)):0.1:max(100*WGTS_CD_MolPerCm2_local/(2.5e17));
% yy = spline(100*WGTS_CD_MolPerCm2_local/(2.5e17),sigmahat*1e3,xx);
% plot(100*WGTS_CD_MolPerCm2_local/(2.5e17),sigmahat*1e3,'s',xx,yy,'LineWidth',2,'Color','black')
% hold off
% PrettyFigureFormat
% fig999str = sprintf('./results/E0eff_%s_%gs_%gfits_RFfluctuation_columndensity.eps',mtd,runtimes,nfit);
% fig9991str = sprintf('./results/E0eff_%s_%gs_%gfits_RFfluctuation_columndensity.fig',mtd,runtimes,nfit);
% publish_figure(999,fig999str);
% export_fig(fig9991str);
% save('./results/fite0_RFCovMat_columndensities_Workspace.mat'); %save whole workspace