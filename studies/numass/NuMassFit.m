function ul90 = NuMassFit(varargin)
%
%            FIT INTEGRAL SPECTRUM
%           Fit Active Neutrino Mass
%
%          Th. Lasserre - CEA Saclay
%                February 2017
%

    % Initialization
    clear par ; clear parnomix;
    addpath(genpath('../../../Samak2.0'));

    % Parser
    p = inputParser;
    p.addParameter('fign',1,@(x)isfloat(x) && x>0);
    p.addParameter('pub','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('nfit',100,@(x)isfloat(x) && x>=0);
    p.addParameter('mnuSq_t',(0.)^2,@(x)isfloat(x));
    p.addParameter('Fitter','Matlab',@(x)ismember(x,{'Minuit','Matlab'}));
    p.addParameter('CovMat','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('SystCorr','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.parse(varargin{:});
    fign       =    p.Results.fign;
    pub        =    p.Results.pub;
    display    =    p.Results.display;
    nfit       =    p.Results.nfit;
    mnuSq_t    =    p.Results.mnuSq_t;
    Fitter     =    p.Results.Fitter;
    CovMat     =    p.Results.CovMat;
    SystCorr   =    p.Results.SystCorr;

    clear FitCovErrMatrix; global FitCovErrMatrix;

    % Parametrization: True Value
    global A ; A=NuMass_InitKATRIN(); A.mnuSq_i=mnuSq_t;
    A.ComputeTBDDS(); A.ComputeTBDIS(); 
    
    switch CovMat
        case 'ON'
            nCMTrials = 10000;
            fName = sprintf('./data/CovMat/WGTSMACE_CovMat_%gTrials_%s_Year_%g.mat',...
                nCMTrials,A.TD,A.TimeYear);
            if exist(fName, 'file') == 2
                CM=importdata(fName);
            else
                A.ComputeTBDISCovarianceMatrix('nTrials',nCMTrials);
                CM=importdata(fName);
            end
    end
    
    % Loop on fits
    npar       = 6;
    fit_p      = ones(npar,nfit);
    fit_xi2min = ones(1,nfit);

    switch Fitter
        case 'Minuit'
        case 'Matlab'
            tbdisf = @(mb,e0b,bb,nb,qu) A.ComputeTBDISf(qu,mb,e0b,bb,nb,0,0);
            tmp = @(mb,e0b,bb,nb,e) interp1(A.qU,tbdisf(mb,e0b,bb,nb,A.qU),e);
            myf = fittype(@(mb,e0b,bb,nb,qu) tmp(mb,e0b,bb,nb,qu),...
                'independent','qu','coefficients', {'mb','e0b','bb','nb'});
            %,'Robust','Bisquare',
    end
    
    % Number of Fits
    tic
    progressbar('...Neutrino Mass Fitting...');
    for f=1:1:nfit
        progressbar(f/nfit);
        
        A.ComputeTBDDS(); A.ComputeTBDIS(); 
        %AddStatFluctTBDIS(A);
        AddStatSystFluctTBDIS(A,'CovMat',CovMat,'SystCorr',SystCorr);
        Data = [A.qU,A.TBDIS,A.TBDISE];

        % Random Init
        i_mnuSq    = (0.2*randn)^2;     
        i_E0       = 0.05*randn;
        i_B        = 0.1*A.BKG_RateAllFPDSec*randn; 
        i_N        = 1e-3*randn;
        i_m4       = 0; 
        i_s4       = 0;
    
        switch Fitter
            case 'Minuit'
                % Initializing Fit
                ParIni = [i_mnuSq i_E0 i_B i_N i_m4 i_s4];
                parnames = ['mSq dE0 B dN m4 s4'];
                tmparg = sprintf(['set pri -10 ;'...
                    'fix  5 6 ;'...
                    'set now; min ; ' ]);
                Args = {ParIni, Data, '-c',tmparg};
                  
                switch CovMat
                    case 'OFF'
                        fprintf(2,'No Detector Systematics - Stat. + Pulls\n');
                        [par, err, chi2min, errmat] = fminuit('NuMassChi2',Args{:});
                    case 'ON'
                        fprintf(2,'Stat + Detector Systematics (Covariance Matrix)\n');
                        switch SystCorr
                            case 'OFF'
                        FitCovErrMatrix = diag(diag(CM.WGTSMACE_CovMat)) + diag(A.TBDIS); 
                            case 'ON'
                        FitCovErrMatrix = (CM.WGTSMACE_CovMat) + diag(A.TBDIS);
                        end
                        [ par, err, chi2min, errmat ] = fminuit('NuMassChi2CovMat',Args{:});
                end
                
                fit_p(:,f)=par;  fit_xi2min(f) = chi2min;
         
            case 'Matlab'
                opts = fitoptions('StartPoint',[i_mnuSq,i_E0,i_B,i_N],...
                    'Method', 'NonLinearLeastSquares',...
                    'Algorithm','Levenberg-Marquardt','Display','off',...
                    'Weights',(1./A.TBDIS));
                [fit1,gof,fitinfo] = fit(A.qU,A.TBDIS,myf,opts);
                par = [fit1.mb fit1.e0b fit1.bb fit1.nb 0 0];
                fit_p(:,f)=par; err = [0 0 0 0 0 0];
                chi2min = gof.sse;  % Sum of squares due to error - Equivalent to Chi2
                fit_xi2min(f) = chi2min;
        end
        
        switch display
            case 'ON'
                fprintf(2,'-----------------------------------------------------------------------------------\n');
                fprintf('  Processing \t= %g %% \n',f/nfit*100);
                fprintf(2,'-----------------------------------------------------------------------------------\n');
                fprintf(2,'  m^2 \t= %g \t± \t %g \tmeV^2 \t - Init: %g \tmeV^2\n', (A.mnuSq_i+par(1))*1e6,err(1)*1e6,i_mnuSq*1e6);
                fprintf(2,'  dE0 \t= %g \t± \t %g \tmeV \t - Init: %g \tmeV\n', (par(2))*1e3,err(2)*1e3,i_E0*1e3);
                fprintf(2,'  dB \t= %g \t± \t %g \tmcps \t - Init: %g \tmcps\n', (A.BKG_RateSec_i+par(3))*1e3,err(3)*1e3,i_B*1e3);
                fprintf(2,'  dN \t= %g \t± \t %g \t \t - Init: %g \n', par(4),err(4),i_N);
                fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,A.nqU-4);
                fprintf(2,'-----------------------------------------------------------------------------------\n');
                
                figure(10000)
                parnomix = zeros(1,4,1); parnomix = par; parnomix(1) = -A.mnuSq;
                hfit = plot(Data(:,1)-A.Q_i,...
                    NuMassModel4par(par)./NuMassModel4par(parnomix)-1,...
                    'Color','Black','LineWidth',1,'LineStyle','-');
                hold on;
                
                switch CovMat
                    case 'ON'
                e = sqrt((diag(CM.WGTSMACE_CovMat)) + (A.TBDIS));
                    case 'OFF'
                e = Data(:,3);
                end
                
                hdata = errorbar(Data(:,1)-A.Q,...
                    Data(:,2)./NuMassModel4par(parnomix)-1,...
                    e./NuMassModel4par(parnomix),...
                    'ks','MarkerSize',5,'MarkerFaceColor',0.9*[1 1 1],'Color','Black','LineWidth',1);
                hold off;
                grid on
                xlabel('qU-E_0 (eV)','FontSize',14);
                ylabel('Spectral distorsion','FontSize',14);
                set(gca,'FontSize',12);
                axis([min(A.qU-A.Q_i) max(A.qU-A.Q_i)+1 min((Data(:,2)./NuMassModel4par(parnomix)-1)) max(Data(:,3)./NuMassModel4par(parnomix))*3])
                mydata = sprintf('Data: m_{eff}=%.2f eV \n',sqrt(abs(A.mnuSq_i)));
                myfit = sprintf('Fit: m_{eff}=%.2f \\pm %.2f eV (\\chi2 / dof=%.1f/%.0f)',sqrt(abs(A.mnuSq_i+par(1))),err(1),chi2min,A.nqU-4);
                legend([hdata hfit],mydata,myfit,'Location','NorthWest') ;
        end
    end
    toc
    
    %% Plot Results
    if nfit<5
        return;
    end
    
    figure(fign)
    subplot(2,1,1)
    hfit1 = plot(Data(:,1)-A.Q_i,NuMassModel4par(par),...
        'LineWidth',1,'LineStyle','-','Color','Black');
    hold on;
    hfit2 = plot(Data(:,1)-A.Q_i,NuMassModel4par(par),...
        'LineWidth',1,'LineStyle','-','Color','Black');
    hdata = errorbar(Data(:,1)-A.Q_i,Data(:,2),Data(:,3),...
        'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
    %errorbar_tick(hdata,200);
    hold off
    grid on
    xlabel('qU-E_0 (eV)','FontSize',14);
    ylabel('Counts','FontSize',14);
    title(sprintf('KATRIN - %s - %g y - %g mcps',A.TD,A.TimeYear,A.BKG_RateSec_i*1e3));
    set(gca,'FontSize',12);
    set(gca,'yscale','lin');
    mydata = sprintf('Data: m_{eff}=%.2f eV \n',sqrt(abs(A.mnuSq_i)));
    myfit = sprintf('Fit: m_{eff}=%.2f \\pm %.2f eV',sqrt(abs(A.mnuSq_i+par(1))),err(1));
    mychi2 = sprintf('\\chi2 / dof=%.1f/%.0f\n',chi2min,A.nqU-4);
    legend([hdata hfit1 hfit2],mydata,myfit,mychi2,'Location','NorthEast') ;
    axis([min(A.qU-A.Q_i) max(A.qU-A.Q_i)+1 min(Data(:,2)) max(Data(:,2))*1.1])
    
    
    subplot(2,1,2)
    parnomix = zeros(1,4,1); parnomix = par; parnomix(1) = -A.mnuSq;
    hfit = plot(Data(:,1)-A.Q_i,...
        NuMassModel4par(par)./NuMassModel4par(parnomix)-1,...
        'Color','Black','LineWidth',1,'LineStyle','-');
    hold on;
    hdata = errorbar(Data(:,1)-A.Q,...
        Data(:,2)./NuMassModel4par(parnomix)-1,...
        Data(:,3)./NuMassModel4par(parnomix),...
        'ks','MarkerSize',5,'MarkerFaceColor',0.9*[1 1 1],'Color','Black','LineWidth',1);
    hold off;
    grid on
    xlabel('qU-E_0 (eV)','FontSize',14);
    ylabel('Spectral distorsion','FontSize',14);
    set(gca,'FontSize',12);
    a= legend([hdata hfit],'Data/No Mixing - 1','Fit/No Mixing - 1','Location','NorthWest');
    axis([min(A.qU-A.Q_i) max(A.qU-A.Q_i)+1 min((Data(:,2)./NuMassModel4par(parnomix)-1)) max(Data(:,3)./NuMassModel4par(parnomix))*3])
    switch pub
        case 'ON'
           myname = sprintf('./figures/f1-katrin_%g-mcps_%g-numsq.eps',...
                A.BKG_RateAllFPDSec*1e3,mnuSq_t); 
            fprintf(2,'publish: %s\n',myname);publish_figure(fign,myname);
    end
    
    figure(fign+1)
    subplot(2,2,1);
    title('Squared Mass','FontSize',12)
    nhist(fit_p(1,:)*1e6,'text','pdf','color','sequential');xlabel('meV^2','FontSize',8); grid on
    subplot(2,2,2);
    title('Endpoint','FontSize',12)
    nhist((fit_p(2,:))*1e3,'text','pdf','color','sequential');xlabel('meV','FontSize',8); grid on
    subplot(2,2,3);
    title('Background','FontSize',12)
    nhist(fit_p(3,:)*1e3,'text','pdf','color','sequential');xlabel('mcps','FontSize',8); grid on
    subplot(2,2,4);
    title('Normalization','FontSize',12)
    nhist(((fit_p(4,:))),'text','pdf','color','sequential');xlabel('fraction','FontSize',8); grid on
    switch pub
        case 'ON'
           myname = sprintf('./figures/f2-katrin_%g-mcps_%g-numsq.eps',...
                A.BKG_RateAllFPDSec*1e3,mnuSq_t); 
            fprintf(2,'publish: %s\n',myname);publish_figure(fign+1,myname);
    end
    
    figure(fign+2)
    ndhist(fit_p(1,:)*1e6,((fit_p(2,:)))*1e3);
    colorbar
    ylabel('Endpoint (meV)','FontSize',14);
    xlabel('Mass squared (meV^2)','FontSize',14);
    R1 = corrcoef(fit_p(1,:),((fit_p(2,:))));
    switch pub
        case 'ON'
           myname = sprintf('./figures/f3-katrin_%g-mcps_%g-numsq.eps',...
                A.BKG_RateAllFPDSec*1e3,mnuSq_t); 
            fprintf(2,'publish: %s\n',myname);publish_figure(fign+2,myname);
    end
    
    figure(fign+3)
    [t,N,X] = nhist(((fit_p(1,:))*1e6));
    % limit on neutrino mass
    ul90 = sqrt(X(find(cumsum(N)>0.9.*nfit, 1 )));
    fprintf(2,'90 CL Estimate for Active Neutrino Mass Upper Bound: %g meV \n',ul90);
    mystd = std(fit_p(1,:)*1e6); %meV^2
    fprintf(2,'Variance for Active Neutrino Mass Upper Bound: %g meV \n',sqrt(mystd));
    hold on
    x = [ul90^2 ul90^2];
    y = [0 max(t)];
    line(x,y,'Color','red','LineStyle','--')
    hold off
    ylabel('Bin counts','FontSize',14);
    xlabel('Mass squared (meV^2)','FontSize',14);
    title(sprintf('KATRIN - %s - %g y - %g mcps - 90CL = %.2f meV',A.TD,A.TimeYear,A.BKG_RateSec_i*1e3,ul90));
    PrettyFigureFormat;
    
    figure(fign+4)
    title(sprintf('KATRIN - %s - %g y - %g mcps',A.TD,A.TimeYear,A.BKG_RateSec_i*1e3));
    [t,N,X] = nhist(fit_xi2min);
    ylabel('Bin counts','FontSize',14);
    xlabel('\chi^2','FontSize',14);
    PrettyFigureFormat;
    
    % Display Confirguration    
    switch display
        case 'ON'
            A.DisplayTDBInfo();
    end
end
