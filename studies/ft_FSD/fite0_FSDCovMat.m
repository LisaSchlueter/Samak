
%Initialization
addpath(genpath('../../../Samak2.0'));
close all

runtimes = 24*60*60*3;
mtd = 'FT-TL4';
day_str = sprintf('-%udays/',runtimes/(24*60*60));
savefolder = strcat('./results/', mtd, day_str);  

% Fits 
nfit   = 500;            % Per Trials 
npar   = 6;              % Number of Fit Parameters
Fitter = 'Minuit';
clear FitCovErrMatrix ; global FitCovErrMatrix;
    
% FSD Case Studied Case studied
FSDNorm_RelErr = 1e-2.*[0 1 2 5 10];
FSDShape_RelErr = 1e-2.*[0 1 2 5 10];
Ntimes = numel(FSDNorm_RelErr);
% General Parameters
fign = 1;

display = 'OFF';
FitFlag         = 'ON';
savespectrum    = 'ON';
plotresults     = 'ON';
plotdist        = 'OFF';
saveresults     = 'ON';
savedist        = 'ON';

% Initialize output matrices (E0 Fit Results / Column Density Trials)
muhat = zeros(Ntimes,1); sigmahat = zeros(Ntimes,1);
nfitsNotConverged = zeros(Ntimes,1); % Number of fits, which dont converge
nconv = 0; % Number of WGTS_CD_MolPerCm2_local for which all fit dont converge
% Define  Study Object (TBD Object)
CovMat_rho = cell(Ntimes,1);

for ii =1:1:Ntimes 
    A=InitKatrin_ft_RF('TD', mtd, 'TimeSec',runtimes);
% Run Parameters
    A.ComputeTBDDS(); A.ComputeTBDIS('IStype','SIMP');    
    fit_p = zeros(nfit,10);
    fit_xi2min = zeros(nfit,1);
    
    % Compute RF Covariance Matrix
    myEffects = struct('FSD','ON');
    C = CovarianceMatrix('StudyObject',A, 'nTrials',1000,'SysEffect',myEffects,...
       'FSDNorm_RelErr',FSDNorm_RelErr(ii) , 'FSDShape_RelErr', FSDShape_RelErr(ii), ...
       'RecomputeFlag','ON');
    C.ComputeCM_FSD; 
    savename = sprintf('FSDerr-%.2f',FSDNorm_RelErr(ii));
    C.PlotCM('saveplot','ON','savename',savename);
    C.DecomposeCM('saveplot','ON');
    % Define Fit Covariance Matrix : Stat. + Syst. 
    FitCovErrMatrix = diag(A.TBDIS) + ((C.CovMat));
    
    % Define Model, again...
    A=InitKatrin_ft_RF('TD', mtd, 'TimeSec',runtimes);
    A.ComputeTBDDS(); A.ComputeTBDIS('IStype','SIMP'); 

    % Loop of fits
    progressbar('E0,m2,B,N Fits');
    DoF = A.nqU - 4;  
    D = InitKatrin_ft_RF('TD', mtd, 'TimeSec',runtimes);
    for kk = 1:nfit
        % Compute & Stat Fluctuate Data
        D.ComputeTBDDS(); D.ComputeTBDIS('IStype','SIMP');      
        D.AddStatFluctTBDIS;                                                                                            % Fluctuate the Data: stat
        Data = [D.qU,D.TBDIS,sqrt(D.TBDIS)];
       
        % Init Fit Parameters
        i_mnu      = 0; i_Q        = 0;
        i_B        = D.TBDIS(end)./(D.TimeSec*D.qUfrac(end));
        i_N        = 0;
        i_m4       = 0; i_s4       = 0;
        DataTBD = {Data,A};
        ParIni = [i_mnu i_Q i_B i_N i_m4 i_s4];
        switch Fitter
            case 'Minuit'
                % Initializing Fit             
                parnames = ['mSq dQ B dN m4 s4'];
                tmparg = sprintf(['set pri -10 ;'...
                    'fix  5 6 ;'...
                    'set now; minos ; ' ]);
                Args = {ParIni, DataTBD, '-c',tmparg};
                [par, err, chi2min, errmat] = fminuit('NuMassChi2E0_Cov',Args{:});
                
            case 'Matlab'            
                options = optimoptions('fminunc','Algorithm','quasi-newton',...
                    'OptimalityTolerance',1e-7,'StepTolerance',1e-7,'FiniteDifferenceType','central');
                %TBDfun = @(xx) NuMassChi2E0(xx,DataL);
                TBDfun = @(xx) NuMassChi2E0_Cov(xx,DataTBD); 
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
        if chi2min > 3*DoF % if fit didnt converge/ is too bad
        fit_p(kk,:) = ones(1,10)*NaN;  fit_xi2min(kk) = NaN;  
        else
        fit_p(kk,:) = [par(1:4) err(1:4) chi2min DoF];  fit_xi2min(kk) = chi2min;
        progressbar(kk/nfit);        
        end
    end    
    nfitsNotConverged(ii) = sum(isnan(fit_xi2min)); %number of not converged fits
    if nfitsNotConverged(ii) == nfit %when no fit converges
        muhat(ii) = NaN;
        sigmahat(ii) = NaN;
        nconv = nconv + 1;
        fprintf(2,'For %g mol/cm2 no fit converged! :( \n',WGTS_CD_MolPerCm2_local(ii))
        WGTS_CD_MolPerCm2_local(ii) = NaN; 
    else
        %% plot Results
        %if max(fit_xi2min)<20*DoF % plot only when all nfits converged
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
        str = {['N_{fit} = ',num2str(nfit)],['\mu =',num2str(int32(muhat(ii)*1e3)),' meV'],['\sigma = ',num2str(int32(sigmahat(ii)*1e3)),' meV'],...
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
        savename = sprintf('fitfull_%gmolPercm2-%.2ferr',WGTS_CD_MolPerCm2_local(ii),WGTS_CD_MolPerCm2_RelErr);
        save([strcat(savefolder,savename),mtd,'_WGTS.mat']);%,'fit_full','WGTS_CD_MolPerCm2_RelErr','muhat','sigmahat','WGTS_CD_MolPerCm2_local');
        fig1str_name = sprintf('E0eff_%s_%gs_%gfits_rho%g_%.2ferr_fit.eps',mtd,runtimes,nfit,WGTS_CD_MolPerCm2_local(ii),WGTS_CD_MolPerCm2_RelErr);
        fig11str_name = sprintf('E0eff_%s_%gs_%gfits_rho%g_%.2ferr_fit.fig',mtd,runtimes,nfit,WGTS_CD_MolPerCm2_local(ii),WGTS_CD_MolPerCm2_RelErr);
        fig1str = strcat(savefolder,fig1str_name);
        fig11str = strcat(savefolder,fig11str_name);
        fprintf('publishing %s ...\n',fig1str);
        publish_figure(ii+10,fig1str);
        export_fig(fig11str);
        fig2str_name = sprintf('E0eff_%s_%gs_%gfits_rho%g_%.2ferr_chi2.eps',mtd,runtimes,nfit,WGTS_CD_MolPerCm2_local(ii),WGTS_CD_MolPerCm2_RelErr);
        fig21str_name = sprintf('E0eff_%s_%gs_%gfits_rho%g_%.2ferr_chi2.fig',mtd,runtimes,nfit,WGTS_CD_MolPerCm2_local(ii),WGTS_CD_MolPerCm2_RelErr);
        fig2str = strcat(savefolder,fig2str_name);
        fig21str = strcat(savefolder,fig21str_name);
        publish_figure(ii+20,fig2str);
        export_fig(fig21str);
        fprintf('publishing %s ...\n',fig2str);
        %else
        %   fprintf(2,'Fits didnt converge with %g molPercm2',WGTS_CD_MolPerCm2_local(ii))
        %end
        close all %close all figures
    end
  
end
%% Results
% plot only when more than 1 WGTS_CD_MolPerCm2 converged
if nconv >= (Ntimes-1)
    fprintf(2,'Fits didnt converge for all (or all but 1) column densites \n')
else 
figure(999)
subplot(2,1,1)
strtitle = sprintf('(E0)eff - MTD=%s - RunTime=%.1f h - %g fits',mtd,runtimes/(60*60),nfit);
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
xlabel('Nominal Column Density (%)')
ylabel('(E0)_{eff}   uncertainty (meV)')
grid on
set(gca,'FontSize',16);
hold on
xx = min(100*WGTS_CD_MolPerCm2_local/(2.5e17)):0.1:max(100*WGTS_CD_MolPerCm2_local/(2.5e17));
yy = spline(100*WGTS_CD_MolPerCm2_local/(2.5e17),sigmahat*1e3,xx);
plot(100*WGTS_CD_MolPerCm2_local/(2.5e17),sigmahat*1e3,'s',xx,yy,'LineWidth',2,'Color','black')
hold off
PrettyFigureFormat
fig999str_name = sprintf('E0eff_%s_%gs_%gfits_RFfluctuation_columndensity.eps',mtd,runtimes,nfit);
fig9991str_name = sprintf('E0eff_%s_%gs_%gfits_RFfluctuation_columndensity.fig',mtd,runtimes,nfit);
fig999str = strcat(savefolder,fig999str_name);
fig9991str = strcat(savefolder,fig9991str_name);
publish_figure(999,fig999str);
export_fig(fig9991str);
end
