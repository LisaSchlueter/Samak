function [fit_p , errmat, myModel, myData] = FT_FitFSD(varargin)
% Fit First Tritium stacked-pixel single run for
%  - E0 
%  - m
%  - B
%  - N
%  - Pgs --> FSD probability to decay to ground state - new here
%
% Input:
%  - run             run number (first tritium)
%  - Mode            'Data','Sim' 
%  - StatFluc        'ON','OFF' (for simulation mode only)
%  - FSD             'eeON' --> fixed Pgs=0.574 
%                    'eeOFF'--> fixed Pgs=1
%                    'eeFit'--> free Pgs (Pgs+Pes=1)
%  - nfit            number of fits
%  - chi2:           'chi2Stat', 'chi2CM', 'chi2CMFrac','chi2CMShape'...
%  - fitter:         'minuit', 'matlab' 
%  - exclDataStart:  7=400eV, 9=200eV, ... 
%  - displayFit:    'ON','OFF' 
% Output:
% - fit_p: fit parameters
% - errmat: error matrix of fit parameters
% - myModel: Model object (TBD Class)
% - myData: Data structure
% 
% T. Lasserre, Pablo Morales, Lisa Schlueter
% Last Update: June 24 2018
%

close all
addpath(genpath('../../../samak'));


% Parser
p = inputParser;
p.addParameter('run',40667);
p.addParameter('Mode','Data',@(x)ismember(x,{'Data','Sim'}));
p.addParameter('StatFluc','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('FSD','eeFit',@(x)ismember(x,{'eeON','eeOFF','eeFit'}));
p.addParameter('nfit',1,@(x)isfloat(x));
p.addParameter('chi2','chi2Stat',@(x)ismember(x,{'chi2Stat', 'chi2CM', 'chi2CMFrac','chi2CMShape', 'chi2P','chi2Nfix'}));
p.addParameter('fitter','minuit',@(x)ismember(x,{'minuit','matlab'}));
p.addParameter('exclDataStart',7,@(x)isfloat(x));
p.addParameter('displayFit','ON',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});

run            = p.Results.run;
Mode           = p.Results.Mode;
StatFluc       = p.Results.StatFluc;
FSD            = p.Results.FSD;
nfit           = p.Results.nfit;
chi2           = p.Results.chi2;
fitter         = p.Results.fitter;
exclDataStart  = p.Results.exclDataStart;
displayFit     = p.Results.displayFit;


  % Lower  Bin 9  = -200 V
    FitLB = exclDataStart;
    % Higher Bin 26 = + 40 V
    FitHB = 26;
    % Data
    
% Read datafile
datafile = sprintf('../../tritium-data/mat/%gex2.mat',run);
runfile  = load(datafile);

% MTD
reffile = sprintf('Run%gex2.mat',run);
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
fit_p = zeros(nfit,14);
muhat = zeros(Ntimes,1);
sigmahat = zeros(Ntimes,1);
fit_xi2min = zeros(nfit,1);

% Build Model
Sref = ref_runsummaries(run,'ex2','ISCS','Theory','recomputeRF','OFF');
Sref.ComputeTBDDS;Sref.ComputeTBDIS

switch Mode
    case 'Sim'
        switch StatFluc
            case 'ON'
                Sref.AddStatFluctTBDIS;
        end
        Count    = Sref.TBDIS; ErrCount = sqrt(Count);
        DataRef  = [Sref.qU,Count,ErrCount];Data = DataRef;
    case 'Data'
        Count    = (runfile.TBDIS(1:end));
        ErrCount = sqrt(Count);
        DataRef  = [Sref.qU,Count,ErrCount];Data = DataRef;
end


% Define Covariance Matrix - Signal
clear FitCovErrMatrixS ; matrixS = load('CovMat_Run40257_0.1rhod_0.02BField_0.02ISX_0.03FSD_0.023TASR_TC.mat');
covmatfrac = matrixS.CovMatFracCombi;
%FitCovErrMatrixS = ((Count .*  covmatfrac .* Count'))+diag(Count);
FitCovErrMatrixS = diag(Count);

% Define Covariance MAtrix - Background
%clear FitCovErrMatrixB ; %matrixB = load('covmatBfracRun4057.mat');
%FitCovErrMatrixBfrac = matrixB.covmatBfrac;
FitCovErrMatrixBfrac=[];

% Sum
ScaleB      = Sref.BKG_RateSec*Sref.qUfrac.*Sref.TimeSec;
%FitCovMatTB = FitCovErrMatrixS(FitLB:end,FitLB:end) + ( ScaleB(FitLB:end).*FitCovErrMatrixBfrac(FitLB:end,FitLB:end).*ScaleB(FitLB:end)');
%FitCovMatTB = FitCovErrMatrixS  + ( ScaleB .*FitCovErrMatrixBfrac.*ScaleB');
FitCovMatTB = FitCovErrMatrixS ;
%FitCovMatTB = diag(Count);


% Loop
progressbar('E0,m2,B,N,DTGS,DTES Fits');
for kk = 1:nfit
    tic
    
    switch Mode
        case 'Sim'
            SetFitBias(Sref,0); Sref.ComputeTBDDS; Sref.ComputeTBDIS;
            switch StatFluc
                case 'ON'
            Count = mvnrnd(Sref.TBDIS,FitCovErrMatrixS,1)';  ErrCount = sqrt(Count);
                case 'OFF'
            Count = Sref.TBDIS;  ErrCount = sqrt(Count);
            end
            DataRef  = [Sref.qU,Count,ErrCount];Data = DataRef;
    end
    
    % Init Fit Parameters
    i_mnu      = 0; i_Q        = -1;
    i_B        = 0;
    i_N        = 0;
    switch FSD
        case 'eeOFF'
            i_DTGS     = +Sref.DTNormES_i; i_DTES       = -Sref.DTNormES_i;
        case 'eeON'
            i_DTGS     = 0; i_DTES       = 0;
        case 'eeFit'
            i_DTGS     = 0; i_DTES       = 0;
    end
    
    DataTBD = {DataRef,Sref,FitCovErrMatrixS,FitLB,FitHB,FitCovErrMatrixBfrac};
    
    switch Fitter
        case 'Minuit'
            % Initializing Fit
            ParIni = [i_mnu i_Q i_B i_N i_DTGS i_DTES];
            parnames = ['mSq dQ B dN i_DTGS i_DTES'];
            
             switch FSD
                 case {'eeOFF','eeON'}
                     tmparg = sprintf(['set pri -10 ;'...
                         'fix  1 5 6;'...
                         'set now; minos ; imp' ]);
                     FreePar=4;
                 case 'eeFit'
                     tmparg = sprintf(['set pri -10 ;'...
                         'fix 1; set now; minos  ; imp ' ]);
                     FreePar=5;
             end
            DoF = numel(Count(FitLB:FitHB)) - FreePar; 
            
            Args = {ParIni, DataTBD, '-c',tmparg};
            [par, err, chi2min, errmat] = fminuit('NuMassChi2FSD',Args{:});
            
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
            fprintf('  DTGS      = %g +- %g\n', Sref.DTNormGS_i+par(5),err(5));
            fprintf('  DTES      = %g +- %g\n', Sref.DTNormES_i+par(6),err(6));
            fprintf('  Chi2/dof  = %g/%g\n', chi2min, numel(Sref.qU(FitLB:FitHB))-FreePar);
            fprintf('===============================================\n');
            
            switch displayFit 
                case 'ON'
                fig = figure('Name','Fit and Residuals','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
                
                subplot(2,2,[1 2])
                hdata = errorbar(Data(FitLB:FitHB,1),Data(FitLB:FitHB,2)./Sref.qUfrac(FitLB:FitHB)./Sref.TimeSec,Data(FitLB:FitHB,3)./Sref.qUfrac(FitLB:FitHB)./Sref.TimeSec,...
                    'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
                hold on
                hdata2 = boundedline(Data(FitLB:FitHB,1),Data(FitLB:FitHB,2)./Sref.qUfrac(FitLB:FitHB)./Sref.TimeSec,sqrt(diag(FitCovMatTB(FitLB:FitHB,FitLB:FitHB)))./Sref.qUfrac(FitLB:FitHB)./Sref.TimeSec,...
                    'alpha','transparency',0.4,'cmap',rgb('DarkCyan'));
                m = NuMassModelFSD(par,Sref);
                hfit1 = plot(Data(FitLB:FitHB,1),m(FitLB:FitHB)./Sref.qUfrac(FitLB:FitHB)./Sref.TimeSec,...
                    'Color','Black','LineWidth',1,'LineStyle','-');
                hfit3 = line([Data(1,1) Data(end,1)],[(Sref.BKG_RateSec_i+par(3)) (Sref.BKG_RateSec_i+par(3))],'LineStyle','--','Color','Red');
                hold off
                grid on
                xlabel('qU (eV)','FontSize',22);
                ylabel('Counts/qU','FontSize',22);
                set(gca,'FontSize',22);
                set(gca,'yscale','log');
                mydata = sprintf('%s: m=%.2f eV - E0=%.2f eV \n',...
                    Mode,sqrt(Sref.mnuSq_i),18575);
                myfit = sprintf('Fit: \\chi2 / dof=%.2f/%.0f \n m^2=%.3f \t\\pm %.3f eV \n E0=%.3f \t\\pm %.3f eV \n B=%.3f \t\\pm %.3f \n N=%.3f \t\\pm %.3f \n GS=%.3f \t\\pm %.3f',...
                    chi2min,numel(Sref.qU(FitLB:FitHB))-2,((Sref.mnuSq_i+par(1))),err(1),(par(2)),err(2),...
                    (Sref.BKG_RateSec_i+par(3))*1e3,err(3),par(4),err(4),Sref.DTNormGS_i+par(5),err(5));
                mychi2 = sprintf('%s',Mode);
                legend([hdata  hfit1 hfit3],mychi2,myfit,'Offset','Location','SouthWest') ;
                legend('boxoff')
                axis(...
                    [min(Sref.qU(FitLB:FitHB)) ...
                    max(Sref.qU(FitLB:FitHB))+1 ...
                    0.*min(Data(FitLB:FitHB,2)./Sref.qUfrac(FitLB:FitHB)./Sref.TimeSec) ...
                    max(Data(FitLB:FitHB,2)./Sref.qUfrac(FitLB:FitHB)./Sref.TimeSec)*1.2]);
                title(sprintf(' Run %g - %s - Electronic Excitations Onset Study (Samak)',...
                    run,Mode),'FontSize',22);
                PrettyFigureFormat;
                
                subplot(2,2,[3 4])
                hdata = plot(Data(FitLB:FitHB,1),(Data(FitLB:FitHB,2)-m(FitLB:FitHB))./sqrt(diag(FitCovMatTB(FitLB:FitHB,FitLB:FitHB))),...
                    'o','MarkerEdgeColor','k');
                hold on
                line([Data(FitLB,1) Data(FitHB,1)],[0 0 ],'LineStyle','--','Color','Blue');
                %            hdata2 = boundedline(Data(:,1),Data(:,2)-NuMassModelFSD(par,Sref),sqrt(diag(FitCovErrMatrix)),'s','alpha');
                [lstat psys]  = boundedline(Data(FitLB:FitHB,1),Data(FitLB:FitHB,1)-Data(FitLB:FitHB,1),sqrt(diag(FitCovMatTB(FitLB:FitHB,FitLB:FitHB)))./sqrt(diag(FitCovMatTB(FitLB:FitHB,FitLB:FitHB))),...
                    'alpha','transparency',0.4,'cmap',rgb('DarkCyan'));
                lstat.LineStyle= '--'; lsys.LineStyle= '--';
                hdata = plot(Data(FitLB:FitHB,1),(Data(FitLB:FitHB,2)-m(FitLB:FitHB))./sqrt(diag(FitCovMatTB(FitLB:FitHB,FitLB:FitHB))),...
                    'o','MarkerEdgeColor','k');
                hold off
                xlim([min(Sref.qU(FitLB:FitHB)) max(Sref.qU(FitLB:FitHB))+1]);
                grid on
                xlabel('qU (eV)','FontSize',22);
                ylabel('Residuals','FontSize',22);
                set(gca,'FontSize',22);
                set(gca,'yscale','lin');
                PrettyFigureFormat;
                figname=sprintf('figures/FSDonset_run%g_excl%g_%s.pdf',run,exclDataStart,FSD);
                publish_figurePDF(fig,figname);
            end
            
            % Endpoint - Pgs correlation
            switch displayFit
                case 'ON'
                    %fig = figure('Name','Fit and Residuals','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
                    corEOPGS=[errmat(2,2), errmat(2,5);errmat(5,2), errmat(5,5)];
                    clear plt;
                    elpt=ellipsedata(corEOPGS,[par(2)+18575,par(5)+Sref.DTNormGS_i],100,[1,2,3]);
                    plt=Plot(elpt(:,1:2:end),elpt(:,2:2:end));
                    strtitle = sprintf('Run %g',run);
                    plt.Title = strtitle; % plot title
                    plt.XLabel = 'Effective Endpoint (eV)'; % xlabel
                    plt.YLabel = 'P_{gs}'; %ylabel
                    plt.FontSize = 16;
                    plt.Legend = {'1 \sigma','2 \sigma','3 \sigma'}; % legends
                    plt.LegendLoc = 'SouthWest';
                    figname=sprintf('figures/FSDonset_run%g_excl%g_%s_E0PGS.png',run,exclDataStart,FSD);
                    plt.export(figname);
                    PrettyFigureFormat
            end
            
    end
    toc
    myModel = Sref;
    myData  = Count; 
    fit_p(kk,:) = [par(1:6) err(1:6) chi2min DoF];  fit_xi2min(kk) = chi2min;
    disp(kk);
    progressbar(kk/nfit);
end



