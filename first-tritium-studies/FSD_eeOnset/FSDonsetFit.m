function FSDonsetFit(varargin)
%
%            FIT INTEGRAL SPECTRUM
%   Fit the onset of electron excited states
%
%          Th. Lasserre - CEA Saclay
%                April 2017
%

    % Initialization
    clear par ; clear parnomix;
    addpath('fminuit'); %digits(6);

    % Parser
    p = inputParser;
    p.addParameter('fign',1,@(x)isfloat(x) && x>0);
    p.addParameter('pub','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('nfit',1,@(x)isfloat(x) && x>=0);
    p.addParameter('mnuSq_t',(0)^2,@(x)isfloat(x));
    p.addParameter('TimeYear',1,@(x)isfloat(x) && x>=0);
    p.addParameter('Fitter','Matlab',@(x)ismember(x,{'Minuit','Matlab'}));
    p.addParameter('CovMat','ON',@(x)ismember(x,{'ON','OFF'}));
    p.parse(varargin{:});
    fign       =    p.Results.fign;
    pub        =    p.Results.pub;
    display    =    p.Results.display;
    nfit       =    p.Results.nfit;
    mnuSq_t    =    p.Results.mnuSq_t;
    TimeYear   =    p.Results.TimeYear;
    Fitter     =    p.Results.Fitter;
    CovMat     =    p.Results.CovMat;

    clear FitCovErrMatrix;
    global FitCovErrMatrix;

    % Parametrization: True Value
    global A ; A=InitKATRIN('TimeYear',TimeYear,'TTFSD','SAENZ','DTFSD','OFF','HTFSD','OFF'); A.mnuSq_i=mnuSq_t;
    A.ComputeTBDDS(); A.ComputeTBDIS(); % A.ComputeTBDISCovarianceMatrix('nTrials',5000);

    global B ; B=InitKATRIN('TimeYear',TimeYear,'TTFSD','SAENZ','DTFSD','OFF','HTFSD','OFF'); B.mnuSq_i=mnuSq_t;
    B.ComputeTBDDS(); B.ComputeTBDIS();
    
    switch CovMat
        case 'ON'
            switch B.TD
                case 'Flat60'
                    CM=importdata('./data/CovMat/WGTSMACE_CovMat_5000Trials_Flat60_Year_0.0273785.mat');
            end
    end
    
    % Loop on fits
    npar       = 8;
    fit_p      = ones(npar,nfit);
    fit_xi2min = ones(1,nfit);

    % Init
    i_mnu      = 0;  i_Q        = 0;
    i_B        = 0;  i_N        = 0;
    i_m4       = 0;  i_s4       = 0;
    i_gs       = 0;  i_es       = 0; % ground/excited states bias 
    
    switch Fitter
        case 'Minuit'
        case 'Matlab'
            tbdisf = @(mb,qb,bb,nb,gsb,esb,qu) A.ComputeTBDISf(qu,mb,qb,bb,nb,gsb,esb);
            tmp = @(mb,qb,bb,nb,gsb,esb,e) interp1(A.qU,tbdisf(mb,qb,bb,nb,gsb,esb,A.qU),e);
            
            myf = fittype(@(mb,qb,bb,nb,gsb,esb,qu) tmp(mb,qb,bb,nb,gsb,esb,qu),...
                'independent','qu','coefficients', {'mb','qb','bb','nb','gsb','esb'});
            opts = fitoptions('StartPoint',[0,0,0,0,0,0],...
                'Method', 'NonlinearLeastSquares','Display','off',...
                'Weights',1./A.TBDIS);
            
%             myf = fittype(@(nb,esb,qu) tmp(0,0,0,nb,esb,qu),...
%                 'independent','qu','coefficients', {'nb','esb'});
%             opts = fitoptions('StartPoint',[0,0],...
%                 'Method', 'NonlinearLeastSquares','Display','iter',...
%                 'Weights',1./A.TBDIS);
    end
    
    % Number of Fits
    tic
    progressbar('FSD Onset Fits');
    for f=1:1:nfit
        progressbar(f/nfit);
        
        B.ComputeTBDDS(); B.ComputeTBDIS(); AddStatFluctTBDIS(B);
        Data = [B.qU,B.TBDIS,B.TBDISE];
        
        switch Fitter
            case 'Minuit'
                % Initializing Fit
                ParIni = [i_mnu i_Q i_B i_N i_m4 i_s4 i_gs i_es];
                parnames = ['mSq dQ B dN m4 s4 gs es'];
                tmparg = sprintf(['set pri -10 ;'...
                    'fix  5 6 ;'...
                    'set now; min ; ' ]);
                Args = {ParIni, Data, '-c',tmparg};
                
                switch CovMat
                    case 'OFF'
                        fprintf(2,'No Detector Systematics - Stat. + Pulls\n');
                        [par, err, chi2min, errmat] = fminuit('FSDonsetChi2',Args{:});
                    case 'ON'
                        fprintf(2,'Stat + Detector Systematics (Covariance Matrix)\n');
                        FitCovErrMatrix = CM.WGTSMACE_CovMat + diag(B.TBDIS);
                        [ par, err, chi2min, errmat ] = fminuit('FSDonsetChi2CovMat',Args{:});
                end
                
                fit_p(:,f)=par;  fit_xi2min(f) = chi2min;
                
            case 'Matlab'
                [fit1,gof,fitinfo] = fit(B.qU,B.TBDIS,myf,opts);
                par = [fit1.mb fit1.qb fit1.bb fit1.nb 0 0 fit1.esb];
                %par = [0 0 0 fit1.nb 0 0 fit1.esb];
                fit_p(:,f)=par;
                err = [0 0 0 0 0 0 0];
                chi2min = gof.sse; fit_xi2min(f) = chi2min;
        end
        switch display
            case 'ON'
                fprintf(2,'--------------------------------------------------------------\n');
                fprintf(2,'  m^2 \t= %g \t± \t %g \tmeV^2 \n', (A.mnuSq_i+par(1))*1e6,err(1)*1e6);
                fprintf(2,'  dQ \t= %g \t± \t %g \tmeV \n', (par(2))*1e3,err(2)*1e3);
                fprintf(2,'  dB \t= %g \t± \t %g \tmcps \n', (A.BKG_RateSec_i+par(3))*1e3,err(3)*1e3);
                fprintf(2,'  dN \t= %g \t± \t %g \n', par(4),err(4));
                fprintf(2,'  dGS \t= %g \t± \t %g \n', par(7),err(7));
                fprintf(2,'  dES \t= %g \t± \t %g \n', par(8),err(8));
                fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,A.nqU-4);
                fprintf(2,'--------------------------------------------------------------\n');
        end
    end
    toc
    
    %% Plot Results
    if nfit<5
        return;
    end
    figure(fign)
    subplot(2,1,1)
    hfit = plot(Data(:,1)-A.Q_i,NuMassModel4parGSESbias(par),...
        'LineWidth',1,'LineStyle','-','Color','Black');
    hold on;
    switch CovMat
        case 'OFF'
            hdata = errorbar(Data(:,1)-A.Q_i,Data(:,2),Data(:,3),...
                'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        case 'ON'
            hdata = errorbar(Data(:,1)-A.Q_i,Data(:,2),Data(:,3),...
                'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
            hdatacm = errorbar(Data(:,1)-A.Q_i,Data(:,2),sqrt(diag(CM.WGTSMACE_CovMat)),...
                'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1,'Color','Red');
    end
    errorbar_tick(hdata,200);    errorbar_tick(hdatacm,200);
    hold off
    grid on
    xlabel('qU-E_0 (eV)','FontSize',14);
    ylabel('Counts','FontSize',14);
    title(sprintf('KATRIN - %s - %g y - %g mcps',A.TD,A.TimeYear,A.BKG_RateSec_i*100));
    set(gca,'FontSize',12);
    set(gca,'yscale','lin');
    mydata = sprintf('Data: m_{eff}=%.2f eV \n',sqrt(abs(A.mnuSq_i)));
    myfit = sprintf('Fit: m_{eff}=%.2f \\pm %.2f eV',sqrt(abs(A.mnuSq_i+par(1))),err(1));
    mychi2 = sprintf('\\chi2 / dof=%.1f/%.0f\n',chi2min,A.nqU-3);
    a = legend([hdata hfit hfit],mydata,myfit,mychi2,'Location','NorthEast') ; legend(a,'boxoff');
    axis([min(A.qU-A.Q_i) max(A.qU-A.Q_i)+1 min(Data(:,2)) max(Data(:,2))*1.1])
    
    
    subplot(2,1,2)
    parnomix = zeros(1,4,1); parnomix = par; parnomix(1) = -A.mnuSq;
    hfit = plot(Data(:,1)-A.Q_i,...
        NuMassModel4parGSESbias(par)./NuMassModel4parGSESbias(parnomix)-1,...
        'Color','Black','LineWidth',1,'LineStyle','-');
    hold on;
    
    switch CovMat
        case 'OFF'
            hdata = errorbar(Data(:,1)-A.Q,...
                Data(:,2)./NuMassModel4parGSESbias(parnomix)-1,...
                Data(:,3)./NuMassModel4parGSESbias(parnomix),...
                'ks','MarkerSize',5,'MarkerFaceColor',0.9*[1 1 1],'Color','Black','LineWidth',1);
            errorbar_tick(hdata,200);
            
        case 'ON'
            hdata = errorbar(Data(:,1)-A.Q,...
                Data(:,2)./NuMassModel4parGSESbias(parnomix)-1,...
                Data(:,3)./NuMassModel4parGSESbias(parnomix),...
                'ks','MarkerSize',5,'MarkerFaceColor',0.9*[1 1 1],'Color','Black','LineWidth',1);
            errorbar_tick(hdata,200);
            hdatacm = errorbar(Data(:,1)-A.Q,...
                Data(:,2)./NuMassModel4parGSESbias(parnomix)-1,...
                sqrt(diag(CM.WGTSMACE_CovMat))./NuMassModel4parGSESbias(parnomix),...
                'ks','MarkerSize',5,'MarkerFaceColor',0.9*[1 1 1],'Color','Black','LineWidth',1,'Color','Red');
            errorbar_tick(hdata,200);  errorbar_tick(hdatacm,200);
            
    end
    
   
    hold off;
    grid on
    xlabel('qU-E_0 (eV)','FontSize',14);
    ylabel('Spectral distorsion','FontSize',14);
    set(gca,'FontSize',12);
    a= legend([hdata hfit],'Data/No Mixing - 1','Fit/No Mixing - 1','Location','SouthWest');
    legend(a,'boxoff');
    axis([min(A.qU-A.Q_i) max(A.qU-A.Q_i)+1 min((Data(:,2)./NuMassModel4parGSESbias(parnomix)-1)) max(Data(:,3)./NuMassModel4parGSESbias(parnomix))*3])
    switch pub
        case 'ON'
           myname = sprintf('./figures/f1-katrinonsetfsd_%g-mcps_%g-numsq-100eV.eps',...
                A.BKG_RateAllFPDSec*1e3,mnuSq_t); 
            fprintf(2,'publish: %s\n',myname);publish_figure(fign,myname);
    end
    
    figure(fign+1)
    subplot(2,2,1);
    title('mass squared','FontSize',12)
    nhist(fit_p(1,:),'text','pdf','color','sequential');xlabel('meV^2','FontSize',8); grid on
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
           myname = sprintf('./figures/f2-katrinonsetfsd_%g-mcps_%g-numsq-100eV.eps',...
                A.BKG_RateAllFPDSec*1e3,mnuSq_t); 
            fprintf(2,'publish: %s\n',myname);publish_figure(fign+1,myname);
    end
    
    figure(fign+2)
    title('Ground/Electronic Excited States Normalization','FontSize',12)
    subplot(1,2,1);
    nhist(fit_p(7,:)*100,'text','pdf','color','sequential');xlabel('GS (Percentage)','FontSize',8); grid on
    subplot(1,2,2);
    nhist(fit_p(8,:)*100,'text','pdf','color','sequential');xlabel('ES (Percentage)','FontSize',8); grid on   
    switch pub
        case 'ON'
            myname = sprintf('./figures/f3-katrinonsetfsd_%g-mcps_%g-numsq-100eV.eps',...
                A.BKG_RateAllFPDSec*1e3,mnuSq_t);
            fprintf(2,'publish: %s\n',myname);publish_figure(fign+2,myname);
    end
    
    switch CovMat
        case 'ON'
            figure(fign+3)
            imagesc(CM.WGTSMACE_CovMat);
            colorbar
            ylabel('qU bins (eV)','FontSize',14);
            xlabel('qU bins (eV)','FontSize',14);
            PrettyFigureFormat
            figure(fign+4)
            imagesc(diag(B.TBDIS));
            colorbar
            ylabel('qU bins (eV)','FontSize',14);
            xlabel('qU bins (eV)','FontSize',14);
            PrettyFigureFormat
            switch pub
                case 'ON'
                    myname = sprintf('./figures/f4-katrinonsetfsd_%g-mcps_%g-numsq-100eV.eps',...
                        A.BKG_RateAllFPDSec*1e3,mnuSq_t);
                    fprintf(2,'publish: %s\n',myname);publish_figure(fign+3,myname);
                        myname = sprintf('./figures/f5-katrinonsetfsd_%g-mcps_%g-numsq-100eV.eps',...
                        A.BKG_RateAllFPDSec*1e3,mnuSq_t);
                    fprintf(2,'publish: %s\n',myname);publish_figure(fign+4,myname);
            end
    end
    
    
     mystdgs = std(fit_p(7,:)*100); %meV^2
    fprintf(2,'Variance for GS: %g %% \n',sqrt(mystdgs));
     mystdes = std(fit_p(8,:)*100); %meV^2
    fprintf(2,'Variance for ES: %g %% \n',sqrt(mystdes));
end
