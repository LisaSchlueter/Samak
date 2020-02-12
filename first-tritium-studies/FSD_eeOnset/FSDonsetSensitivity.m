function ul90 = FSDonsetSensitivity(varargin)
%
%            KATRIN INTEGRAL SPECTRUM
%             FSD Onset Sensitivity
%
%          Th. Lasserre - CEA Saclay
%                April 2017
%

    % Initialization
    clear par ; clear parnomix;
    addpath('fminuit'); addpath('tools'); digits(6);

    % Parser
    p = inputParser;
    p.addParameter('fign',1,@(x)isfloat(x) && x>0);
    p.addParameter('pub','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('eeP_min',0.,@(x)isfloat(x) && x>=0);
    p.addParameter('eeP_max',.5,@(x)isfloat(x) && x>=0);
    p.addParameter('eeP_nTrials',10,@(x)isfloat(x) && x>=0);
    p.addParameter('eeP_t',(0)^2,@(x)isfloat(x));
    p.addParameter('Fitter','Minuit',@(x)ismember(x,{'Minuit','Matlab'}));
    p.addParameter('CovMat','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('TimeYear',10/365.25,@(x)isfloat(x) && x>=0);

    p.parse(varargin{:});
    fign          =    p.Results.fign;
    pub           =    p.Results.pub;
    display       =    p.Results.display;
    eeP_min       =    p.Results.eeP_min;
    eeP_max       =    p.Results.eeP_max;
    eeP_nTrials   =    p.Results.eeP_nTrials;
    eeP_t         =    p.Results.eeP_t;
    Fitter        =    p.Results.Fitter;
    CovMat        =    p.Results.CovMat;
    TimeYear      =    p.Results.TimeYear;

    clear FitCovErrMatrix;
    global FitCovErrMatrix;
    
    % Parametrization: True Value
    global A ; A=InitKATRIN('TimeYear',TimeYear,'TTFSD','SAENZ','DTFSD','OFF','HTFSD','OFF'); A.mnuSq_i=0;
    A.ComputeTBDDS(); A.ComputeTBDIS();

    % Data
    B=InitKATRIN('TimeYear',TimeYear,'TTFSD','SAENZ','DTFSD','OFF','HTFSD','OFF'); B.mnuSq_i=0;
    B.ComputeTBDDS(); B.ComputeTBDIS(); %AddStatFluctTBDIS(B);
    Data = [B.qU,B.TBDIS,B.TBDISE];

    switch CovMat
        case 'ON'
            switch B.TD
                case 'Flat60'
                    CM=importdata('./data/CovMat/WGTSMACE_CovMat_5000Trials_Flat60_Year_0.0273785.mat');
            end
    end
    
    % Loop on fits
    npar       = 8;
    fit_p      = zeros(npar,eeP_nTrials);
    fit_xi2min = zeros(1,eeP_nTrials);

    % Init
    i_mnu      = 0; i_Q        = 0;
    i_B        = 0; i_N        = 0;
    i_m4       = 0; i_s4       = 0;
    i_gs       = 0; i_es       = 0; % ground/excited states bias

    switch Fitter
        case 'Minuit'
        case 'Matlab'
    end
    
    tic
    % Scan grid of active electron neutrino mass %fliplr
    nScanSteps=1:1:eeP_nTrials; counter=0;
    eePPt = (linspace(B.TTNormES_i*0.97,B.TTNormES_i*1.03,numel(nScanSteps)));
    % m^2 > 0
    progressbar('FSD Onset Sensitivity Scan');
    for f=nScanSteps
        counter=counter+1;
        progressbar(f/eeP_nTrials);
        A.TTNormES_i = eePPt(f);
        A.TTNormGS_i = (1 - A.TTNormES_i);
        %A.TTES_bias = -eePPt(f); 
        fprintf(2,'--------------------------------------------------------------\n');
        fprintf(2,'Initializing:  Model(eeP) = %.2f %%  Model(gsP) = %.2f %% (%g)\n - Data(eeP) = %.2f %% \n',...
            A.TTNormES_i*100,A.TTNormGS_i*100,(A.TTNormES_i+A.TTNormGS_i),B.TTNormES_i*100);
        fprintf(2,'--------------------------------------------------------------\n');

        switch Fitter
            case 'Minuit'
                % Initializing Fit
                ParIni = [i_mnu i_Q i_B i_N i_m4 i_s4 i_gs i_es];
                parnames = ['mSq dQ B dN m4 s4 gs es'];
                tmparg = sprintf(['set pri -10 ;'...
                    'fix 5 6  8;'...
                    'set now; min ; ' ]);
                Args = {ParIni, Data, '-c',tmparg};
                
                 switch CovMat
                    case 'OFF'
                        fprintf(2,'No Detector Systematics - Stat. + Pulls\n');
                        [par, err, chi2min, errmat] = fminuit('FSDonsetChi2',Args{:});
                    case 'ON'
                        fprintf(2,'Stat + Detector Systematics (Covariance Matrix)\n');
                        FitCovErrMatrix = CM.WGTSMACE_CovMat + diag(B.TBDIS)/1e9;
                        [ par, err, chi2min, errmat ] = fminuit('FSDonsetChi2CovMat',Args{:});
                 end
                 
                %[par, err, chi2min, errmat] = fminuit('FSDonsetChi2',Args{:});
                fit_p(:,f)=par;  fit_xi2min(f) = chi2min;
                
            case 'Matlab'
                
        end
        
        switch display
            case 'ON'
                fprintf(2,'--------------------------------------------------------------\n');
                fprintf('  Processing \t= %g %% \n',f/numel(nScanSteps)*100);
                fprintf(2,'--------------------------------------------------------------\n');
                fprintf(2,'  m^2 \t= %g \t± \t %g \tmeV^2 \n', (A.mnuSq_i+par(1))*1e6,err(1)*1e6);
                fprintf(2,'  dQ \t= %g \t± \t %g \tmeV \n', (par(2))*1e3,err(2)*1e3);
                fprintf(2,'  B \t= %g \t± \t %g \tmcps \n', (A.BKG_RateSec_i+par(3))*1e3,err(3)*1e3);
                fprintf(2,'  dN \t= %g \t± \t %g \n', par(4),err(4));
                fprintf(2,'  dGS \t= %g \t± \t %g \n', par(7),err(7));
                fprintf(2,'  dES \t= %g \t± \t %g \n', par(8),err(8));
                fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,A.nqU-4);
                fprintf(2,'--------------------------------------------------------------\n');
        end
        
%         CL=2.71;
%         if (chi2min>CL)
%             fprintf(2,'  Chi2 \t= %g / %g dof --> eeP = %g %%\n \n',chi2min,A.nqU-4,A.TTNormES_i*1e2);
%         end
%         if (chi2min>CL*1.1)
%             break;
%         end
    end
        
    %% Plot Results
    B.DisplayTDBInfo();
    figure(fign)
    x = eePPt(1:counter)*1e2;
    y = fit_xi2min(1:counter);
    plot(x,y,'LineWidth',2,'LineStyle','-','Color','Black');
    xlabel('Electronics Excitation Probability'); ylabel('\Delta\chi ^2 ');
    grid on
    func = @(x) interp1( y(1:counter) - min(y(1:counter)) , eePPt(1:counter)*1e2 , x);
    CL=1; eePPtCL = func(CL);
    legIn=sprintf('P(GS) = %.2f +- %.2f %% - P(ES) uncertainty = %.2f %% (1\\sigma)  \n',...
        B.TTNormGS_i*100,2,abs(B.TTNormES_i*100-eePPtCL));
    a = legend({legIn},'Location','northwest','FontSize',10); legend(a,'boxoff');
    %line([eePPtCL eePPtCL],[min(y) max(y)],'LineWidth',2,'Color','Red','LineStyle','--');
    line([min(x) max(x)],[CL CL],'LineWidth',2,'Color','Red','LineStyle','--');
    axis([min(x) max(x) min(y) max(y)]);
    PrettyFigureFormat;

    switch pub
        case 'ON'
            myname = sprintf('./figures/f1-katrineeonset_%g-mcps_%s-%.0f_TimeYear.eps',...
                A.BKG_RateAllFPDSec*1e3,A.TD,A.TimeYear);
            fprintf(2,'publish: %s\n',myname);publish_figure(fign,myname);
            
    end
            
end
