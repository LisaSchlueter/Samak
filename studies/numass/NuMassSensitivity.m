function [ul90m ul90p] = NuMassSensitivity(varargin)
%
%            KATRIN INTEGRAL SPECTRUM
%         Active Neutrino Mass Sensitivity
%
%          Th. Lasserre - CEA Saclay
%                April 2017
%

    % Initialization
    clear par ; clear parnomix;
    addpath(genpath('../../../Samak2.0'));

    % Parser
    p = inputParser;
    p.addParameter('fign',1,@(x)isfloat(x) && x>0);
    p.addParameter('pub','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('Display','ON',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('mnuSq_min',.01,@(x)isfloat(x) && x>=0);
    p.addParameter('mnuSq_max',.4,@(x)isfloat(x) && x>=0);
    p.addParameter('mnuSq_nTrials',100,@(x)isfloat(x) && x>=0);
    p.addParameter('mnuSq_t',(0)^2,@(x)isfloat(x));
    p.addParameter('Fitter','Matlab',@(x)ismember(x,{'Minuit','Matlab'}));
    p.addParameter('PS_Wein93','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('TD','DR30',@(x)ismember(x,{'DR20','DR30','DR40','DR50','Kle15','Kle15Ext','Kle15Dyn','Flat20','Flat30','Flat50','Flat60','Flat100'}));
    p.addParameter('BKG_RateAllFPDSec',10e-3,@(x)isfloat(x) && x>0);
    p.addParameter('CovMat','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('SystCorr','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.parse(varargin{:});
    fign              =    p.Results.fign;
    pub               =    p.Results.pub;
    Display           =    p.Results.Display;
    mnuSq_min         =    p.Results.mnuSq_min;
    mnuSq_max         =    p.Results.mnuSq_max;
    mnuSq_nTrials     =    p.Results.mnuSq_nTrials;
    mnuSq_t           =    p.Results.mnuSq_t;
    Fitter            =    p.Results.Fitter;
    PS_Wein93         =    p.Results.PS_Wein93;
    TD                =    p.Results.TD;
    BKG_RateAllFPDSec =    p.Results.BKG_RateAllFPDSec;
    CovMat            =    p.Results.CovMat;
    SystCorr          =    p.Results.SystCorr;

    % Parametrization: True Value
    global A ; A=InitKATRIN('PS_Wein93',PS_Wein93,'TD',TD,'BKG_RateAllFPDSec',BKG_RateAllFPDSec); A.mnuSq_i=mnuSq_t;
    A.ComputeTBDDS(); A.ComputeTBDIS();

    clear FitCovErrMatrix;
    global FitCovErrMatrix;
    
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
    
    % Data
    B=InitKATRIN('PS_Wein93',PS_Wein93,'TD',TD,'BKG_RateAllFPDSec',BKG_RateAllFPDSec); B.mnuSq_i=0;
    B.ComputeTBDDS(); B.ComputeTBDIS(); %AddStatFluctTBDIS(B);
    Data = [B.qU,B.TBDIS,B.TBDISE];

    % Loop on fits
    npar       = 6;
    fit_p      = zeros(npar,mnuSq_nTrials);
    fit_xi2min = zeros(1,mnuSq_nTrials);

    % Init
    i_mnu      = 0; i_Q        = 0;
    i_B        = 0; i_N        = 0;
    i_m4       = 0; i_s4       = 0;
    
    switch Fitter
        case 'Minuit'
        case 'Matlab'
    end
    
    tic
    % Scan grid of active electron neutrino mass
    nScanSteps=1:1:mnuSq_nTrials; counter=0;
    mact = linspace(mnuSq_min,mnuSq_max,numel(nScanSteps)); mact(1)=0;
    % m^2 > 0
    for f=nScanSteps
        counter=counter+1;
        A.mnuSq_i = mact(f)^2;
        switch Display
            case 'ON'
                fprintf(2,'--------------------------------------------------------------\n');
                fprintf(2,'Initializing m_{eff} = %g meV\n',sqrt(A.mnuSq)*1e3);
                fprintf(2,'--------------------------------------------------------------\n');
        end
        
        switch Fitter
            case 'Minuit'
                % Initializing Fit
                ParIni = [0 i_Q i_B i_N i_m4 i_s4];
                parnames = ['mSq dQ B dN m4 s4'];
                tmparg = sprintf(['set pri -10 ;'...
                    'fix 1 5 6 ;'...
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
                tbdisf = @(qb,bb,nb,qu) A.ComputeTBDISf(qu,0,qb,bb,nb,0,0);
                tmp = @(qb,bb,nb,e) interp1(A.qU,tbdisf(qb,bb,nb,A.qU),e);
                myf = fittype(@(qb,bb,nb,qu) tmp(qb,bb,nb,qu),...
                    'independent','qu','coefficients', {'qb','bb','nb'});
                opts = fitoptions('StartPoint',[0,0,0],...
                    'Method', 'NonlinearLeastSquares','Display','off',...
                    'Weights',1./B.TBDIS);
                [fit1,gof,fitinfo] = fit(B.qU,B.TBDIS,myf,opts);
                par = [0 fit1.qb fit1.bb fit1.nb 0 0];
                fit_p(:,f)=par;
                err = [0 0 0 0 0 0];
                chi2min = gof.sse; fit_xi2min(f) = chi2min;
        end
        switch Display
            case 'ON'
                fprintf(2,'--------------------------------------------------------------\n');
                fprintf('  Processing \t= %g %% \n',f/numel(nScanSteps)*100);
                fprintf(2,'--------------------------------------------------------------\n');
                fprintf('  m^2 \t= %g \t± \t %g \tmeV^2 \n', (A.mnuSq_i+par(1))*1e6,err(1)*1e6);
                fprintf('  dQ \t= %g \t± \t %g \tmeV \n', (par(2))*1e3,err(2)*1e3);
                fprintf('  B \t= %g \t± \t %g \tmcps \n', (A.BKG_RateSec_i+par(3))*1e3,err(3)*1e3);
                fprintf('  dN \t= %g \t± \t %g \n', par(4),err(4));
                fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,A.nqU-4);
                fprintf(2,'--------------------------------------------------------------\n');
        end
        
        CL=2.71;
        if (chi2min>CL)
            switch Display
                case 'ON'
                    fprintf(2,'  Chi2 \t= %g / %g dof --> m_{eff} = %g meV\n \n',chi2min,A.nqU-4,sqrt(A.mnuSq)*1e3);
            end
        end
        if (chi2min>CL*1.1)
            break;
        end
    end

    %% Plot Results
    switch Display
        case 'ON'
            B.DisplayTDBInfo();
    end
    figure(fign)
    subplot(1,2,2)
    x = mact(1:counter)*1e3;
    y = fit_xi2min(1:counter);
    plot(x,y,'LineWidth',2,'LineStyle','-','Color','Black');
    xlabel('m_{eff}^2 > 0 - m_{eff}(meV)'); ylabel('\Delta\chi ^2 ');
    grid on
    func = @(x) interp1( y(1:counter) - min(y(1:counter)) , mact(1:counter)*1e3 , x);
    CL=2.71; mactCL = func(CL); ul90p = mactCL;
    legIn=sprintf('m_{eff} < %.0f meV @90%% CL\n',ul90p);
    a = legend({legIn},'Location','northwest','FontSize',10); legend(a,'boxoff');
    line([mactCL mactCL],[min(y) max(y)],'LineWidth',2,'Color','Red','LineStyle','--');
    line([min(x) max(x)],[CL CL],'LineWidth',2,'Color','Red','LineStyle','--');
    axis([min(x) max(x) min(y) max(y)]);
    %publish_figure(fign,'tmpPos.eps')
    PrettyFigureFormat;

    % Scan grid of active electron neutrino mass - NEGATIVE
    nScanSteps=1:1:mnuSq_nTrials; counter=0;
    mact = linspace(-mnuSq_min,-mnuSq_max,numel(nScanSteps)); mact(1)=0;
    % m^2 < 0
    for f=nScanSteps
        counter=counter+1;
        A.mnuSq_i = -mact(f)^2; 
 switch Display
            case 'ON'
                fprintf(2,'--------------------------------------------------------------\n');
        fprintf(2,'Initializing m_{eff} = %g meV\n',-sqrt(-A.mnuSq)*1e3);
        fprintf(2,'--------------------------------------------------------------\n');
 end
        switch Fitter
            case 'Minuit'
                % Initializing Fit
                ParIni = [0 i_Q i_B i_N i_m4 i_s4];
                parnames = ['mSq dQ B dN m4 s4'];
                tmparg = sprintf(['set pri -10 ;'...
                    'fix 1 5 6 ;'...
                    'set now; min ; ' ]);
                Args = {ParIni, Data, '-c',tmparg};
                [par, err, chi2min, errmat] = fminuit('NuMassChi2',Args{:});
                fit_p(:,f)=par;  fit_xi2min(f) = chi2min;
                
            case 'Matlab'
                tbdisf = @(qb,bb,nb,qu) A.ComputeTBDISf(qu,0,qb,bb,nb,0,0);
                tmp = @(qb,bb,nb,e) interp1(A.qU,tbdisf(qb,bb,nb,A.qU),e);
                myf = fittype(@(qb,bb,nb,qu) tmp(qb,bb,nb,qu),...
                    'independent','qu','coefficients', {'qb','bb','nb'});
                opts = fitoptions('StartPoint',[0,0,0],...
                    'Method', 'NonlinearLeastSquares','Display','off',...
                    'Weights',1./B.TBDIS);
                [fit1,gof,fitinfo] = fit(B.qU,B.TBDIS,myf,opts);
                par = [0 fit1.qb fit1.bb fit1.nb 0 0];
                fit_p(:,f)=par;
                err = [0 0 0 0 0 0];
                chi2min = gof.sse; fit_xi2min(f) = chi2min;
        end
        switch Display
            case 'ON'
                fprintf(2,'--------------------------------------------------------------\n');
                fprintf('  Processing \t= %g %% \n',f/numel(nScanSteps)*100);
                fprintf(2,'--------------------------------------------------------------\n');
                fprintf('  m^2 \t= %g \t± \t %g \tmeV^2 \n', (A.mnuSq_i+par(1))*1e6,err(1)*1e6);
                fprintf('  dQ \t= %g \t± \t %g \tmeV \n', (par(2))*1e3,err(2)*1e3);
                fprintf('  B \t= %g \t± \t %g \tmcps \n', (A.BKG_RateSec_i+par(3))*1e3,err(3)*1e3);
                fprintf('  dN \t= %g \t± \t %g \n', par(4),err(4));
                fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,A.nqU-4);
                fprintf(2,'--------------------------------------------------------------\n');
        end
        
        CL=2.71;
        if (chi2min>CL)
            switch Display
                case 'ON'
                    fprintf(2,'  Chi2 \t= %g / %g dof --> m_{eff} = %g meV\n \n',chi2min,A.nqU-4,-sqrt(-A.mnuSq)*1e3);
            end
        end
        
        if (chi2min>CL*1.1)
            break;
        end
    end
    toc
    
    %% Plot Results
 switch Display
                case 'ON'
                    B.DisplayTDBInfo();
 end
    %figure(fign+1)
    subplot(1,2,1)
    x = mact(1:counter)*1e3;                   
    y = fit_xi2min(1:counter);          
    plot(x,y,'LineWidth',2,'LineStyle','-','Color','Black');
    xlabel('m_{eff}^2 < 0 - m_{eff}(meV)'); ylabel('\Delta\chi ^2 ');
    grid on
    func = @(x) interp1( y(1:counter) - min(y(1:counter)) , mact(1:counter)*1e3 , x);
    CL=2.71; mactCL = func(CL); ul90m = mactCL;
    legIn=sprintf('m_{eff} > %.0f meV @90%% CL\n',ul90m);
    a = legend({legIn},'Location','northeast','FontSize',10); legend(a,'boxoff');
    line([mactCL mactCL],[min(y) max(y)],'LineWidth',2,'Color','Red','LineStyle','--');
    line([min(x) max(x)],[CL CL],'LineWidth',2,'Color','Red','LineStyle','--');
    axis([min(x) max(x) min(y) max(y)]);
    mytitle = sprintf('KATRIN - %g mcps - %s - %g years - [Wein93] Phase Space: %s',...
        A.BKG_RateAllFPDSec*1e3,A.TD,A.TimeYear,A.PS_Wein93);
    mtit(mytitle);
    PrettyFigureFormat;
    
    switch pub
        case 'ON'
            myname = sprintf('./figures/f1-katrinnumass_%g-mcps_%s-%.0f_TimeYear-Wein93_%s.eps',...
                A.BKG_RateAllFPDSec*1e3,A.TD,A.TimeYear,A.PS_Wein93);
            fprintf(2,'publish: %s\n',myname);publish_figure(fign,myname);
            
    end
            
end
