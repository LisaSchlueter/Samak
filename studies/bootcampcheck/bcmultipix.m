%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Endpoint sensitivity study
% Detector as 
% Configuration for First Tritium May

% P. Morales 2018
% Last update 23/04/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initialization
clear; close all;
addpath(genpath('../../../Samak2.0'));

load('runmat1_10.mat');

%Parameters
nfit            = 1;

nTeBinFac       = 5;
Fitter          = 'Matlab';
UseCovMat       = 'OFF'; 

Mode            = 'Sim';
TD              = '1_10';

% Computing
savespectrum    = 'OFF';
FitFlag         = 'ON';
% Displaying
display         = 'ON';

% Saving
saveresults     = 'OFF';
savepdfs        = 'OFF';
plotfpdviewer   = 'OFF';

% Parametrization:
mnuSq_t = (0.00)^2;

npixel = 148;
pixels = 1:148;
MACE_Bmax_T = 6;
MACE_Ba_T = 6e-4;

qUdata = TBDISallPixels(:,1:3:(3*npixel-2));
countsdata = TBDISallPixels(:,2:3:(3*npixel-1));
errordata = TBDISallPixels(:,3:3:(3*npixel));

qUdata = qUdata(:,pixels);
countsdata = countsdata(:,pixels);
errordata = errordata(:,pixels);


% Initialization
opt_bin = {...
    'nTeBinningFactor',nTeBinFac};

opt_wgtsmace = {...
    'KTFFlag','Compute',...
    'WGTS_B_T',WGTS_B_T}; %MACE+WGTSIS %SSCW_DCOfficial %SSCW_DCperqU

opt_bkg = {...
    'BKG_Flag','XmasData',...
    'BKG_Type','FLAT'};

opt_fsd= {'DTFSD',DTFSD,...
    'HTFSD',HTFSD,...
    'TTFSD',TTFSD};

opt_katrin = {...
    'TD',TD,...
    'TimeSec',TimeSec,...
    'Mode',Mode,...
    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2};

opt_doppler = {...
    'DopplerEffectFlag',DopplerEffectFlag};

opt_fpd = {'FPD_Eff_pix',FPD_MeanEff*ones(npixel,1),...
    'FPD_Pixel',pixels};

opt_mace = {...
    'MACE_Bmax_T',MACE_Bmax_T,...
    'Pixel_MACE_Ba_TCorr',MACE_Ba_TallPixels/MACE_Ba_T,...
    'MACE_Ba_T',MACE_Ba_T};

fit_p = zeros(nfit,4*length(pixels)+6);

if strcmp(FitFlag,'ON')
    for ii = 1:nfit
        tic

        Data = [qUdata(:),countsdata(:),errordata(:)];
        
%        A = ref_bcmultipix(...
%             opt_katrin{:},...
%             opt_wgtsmace{:},...
%             opt_fsd{:},...
%             opt_bkg{:},...
%             opt_bin{:},...
%             'mnuSq_i',mnuSq_t,...
%             opt_doppler{:},...
%             opt_fpd{:},...
%             opt_mace{:},...
%             'UseParallelRF','ON');
%         
%         A.ComputeTBDDS(); A.ComputeTBDIS();
%         A.AddStatFluctTBDIS();
%         qUdata = A.qU; countsdata = A.TBDIS;
%         errordata = A.TBDISE;
%         
%         Data = [A.qU(:),A.TBDIS(:),A.TBDISE(:)];
        
        Afit = ref_bcmultipix(...
            opt_katrin{:},...
            opt_wgtsmace{:},...
            opt_fsd{:},...
            opt_bkg{:},...
            opt_bin{:},...
            'mnuSq_i',mnuSq_t,...
            opt_doppler{:},...
            opt_fpd{:},...
            opt_mace{:},...
            'UseParallelRF','ON',...
            'BKG_Flag','ON',...
            'BKG_RateAllFPDSec',0);
        
        Afit.ComputeTBDDS();
        Afit.ComputeTBDIS();
        
        % Init
        i_mnu      = 0;
        i_Q        = 0;
        i_B        = countsdata(end-3,:)/(Afit.qUfrac(end)*Afit.TimeSec);
        i_N        = (countsdata(1,:)-countsdata(end,:))./...
            (Afit.TBDIS(1,:)-countsdata(end,:)) - 1;
%         i_N = zeros(1,Afit.nPixels); 
        
        DoF = Afit.nqU*Afit.nPixels - length(i_B) - length(i_N) - 2;
        
        DataTBD = {Data,Afit};
        ParIni = [i_mnu i_Q i_B i_N];
        
        switch Fitter
            case 'Minuit'
                %                 % Initializing Fit
                %
                %                 parnames = ['mSq dQ B dN'];
                %                 % tmparg = sprintf(['set pri 1;'...
                %                 % 'fix  5 6 ; min; minos; imp; scan 1 ']);
                %                 tmparg = sprintf(['set pri -10;'...
                %                     'set now; min;']);
                %                 Args = {ParIni, DataTBD, '-c',tmparg};
                %                 fprintf(2,'No Detector Systematics - Stat.\n');
                %                 [par, err, chi2min, errmat] = fminuit('Chi2E0_first',Args{:});
                
                fprintf('Not working for more than 150 parameters. Use Matlab as Fitter. \n');
            case 'Matlab'
                fprintf('----------BEGIN FIT-------------- \n');
                options = optimoptions('fminunc','Algorithm','quasi-newton',...
                    'OptimalityTolerance',1e-7,'StepTolerance',1e-8,...
                    'FiniteDifferenceType','central',...
                    'MaxFunctionEvaluations',1e6,'UseParallel',true,...
                    'Display','iter');
                switch UseCovMat
                    case 'OFF'
                        TBDfun = @(xx) Chi2BC(xx,DataTBD);
                    case 'ON'
                        TBDfun = @(xx) Chi2E0mpixCM(xx,DataTBD);
                end
                [par,chi2min,exitflag,output,grad,Hessian] = fminunc(TBDfun,ParIni,options);
                if exitflag ~= 1
                    [par,chi2min,exitflag,output,grad,Hessian] = fminunc(TBDfun,par,options);
                end
                errmat = 0.5*Hessian;
                varcov = inv(errmat);
                err = sqrt(diag(varcov))';
                norms_fit = par(3+Afit.nPixels:3+2*Afit.nPixels-1);
                bcks_fit = par(3:2+Afit.nPixels);
        end
        
        switch display
            case 'ON'
                mnuSq_report = Afit.mnuSq_i+par(1);
                mnu_report = sqrt(abs(mnuSq_report));
                err_mnu = sqrt(abs((mnuSq_report + err(1)))) - mnu_report;
                BCK_fit = sum(bcks_fit);
                BCK_err = sqrt(sum(bcks_fit.^2));
                N_ave = mean(abs(norms_fit));
                N_err = sqrt(sum(norms_fit.^2))/Afit.nPixels;
                fprintf(2,'----------------------------------------------\n');
                fprintf('  Processing t= %g %% \n',ii/nfit*100);
                fprintf('===============================================\n');
                fprintf('  m^2       = %g +/- %g eV^2\n', mnuSq_report,err(1));
                fprintf('  m         = %g +/- %g eV\n', mnu_report,err_mnu);
                fprintf('  dQ        = %g +/- %g eV\n', par(2),err(2));
                fprintf('  B total   = %g +/- %g cps\n', BCK_fit,BCK_err);
                fprintf('  dN average= %g +/- %g\n', N_ave,N_err);
                fprintf('  Chi2/dof  = %g/%g\n', chi2min, DoF);
                fprintf('===============================================\n');
        end
        toc
        fit_p(ii,:) = [par err chi2min DoF];
        disp(ii);
    end
end

if strcmp(saveresults,'ON'); save(['fit_bcmultipix_','BC_check','.mat'],'fit_p'); end

%% Plotting in PDF

if strcmp(plotfpdviewer,'ON')
    pixel_norms = NaN(148,1);
    pixel_norms(Afit.FPD_Pixel) = norms_fit; 
    FPDViewer(pixel_norms,'ReDrawSkeleton','ON');
    title(['Normalization bias ',TD,' ',num2str(TimeSec), 's']);
    if ispc
        export_fig('plots/FPDNorm.pdf','-pdf');
    elseif isunix || ismac
        publish_figure(1,'plots/bcmultipix001.eps');
    end
    
    pixel_bcks = NaN(148,1);
    pixel_bcks(Afit.FPD_Pixel) = bcks_fit;
    FPDViewer(pixel_bcks,'ReDrawSkeleton','ON');
    title(['Background ',TD,' ',num2str(TimeSec), 's']);
    
    if ispc
        export_fig('plots/FPDBack.pdf','-pdf');
    elseif isunix || ismac
        publish_figure(2,'plots/bcmultipix002.eps');
    end
end
%% from here

if strcmp(savepdfs,'ON')
    M_IS = Afit.TBDISallPixels;
    M_ISE = Afit.TBDISEallPixels;
    
    D_IS = countsdata;
    D_ISE = errordata;    

    plotName = cell(1,npixel);
    for p = 1:npixel
        figure(p+2);
        set(gcf, 'Visible', 'off');
        hold on
        hline = plot(Afit.qU - Afit.Q, M_IS(:,p),...
            'LineWidth',1,'LineStyle','-','Color','Black');
        hchi2 = plot(Afit.qU - Afit.Q, M_IS(:,p),...
            'LineWidth',1,'LineStyle','-','Color','Black');
        hdata = errorbar(Afit.qU - Afit.Q_i, D_IS(:,p), D_ISE(:,p),...
            'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        hold off
        %mydata = '$$\\mathrm{Data:} E_{0\\mathrm{eff}} = 18575 \\mathrm{eV}$$';
        %mydata = text(0.5,0.5,mydata,'interpreter','latex');
        %myline = sprintf('$$\\mathrm{Fit: }\\delta E_{0\\mathrm{eff}} = %.3f \\pm %.3f \\mathrm{eV}$$',par(2),err(2));
        %myline = text(0.5,0.5,myline,'interpreter','latex');
        mydata = ['Data: E_0eff = 18575 eV'];
        myline = sprintf('Fit: \\delta E_0eff = %.3f \\pm %.3f eV',par(2),err(2));
        mychi2 = sprintf('\\chi^2/dof = %.1f/%.0f\n',chi2min,DoF);
        lgd_fit = legend([hline hchi2 hdata],myline,mychi2,mydata,'Location','Best') ;
        lgd_fit.FontSize = 8;
        xlabel('qU-E_0 (eV)','FontSize',14);
        ylabel('Counts','FontSize',14);
        title(sprintf('KATRIN - %s - pix %g - %0.0f s ',TD,p,TimeSec));
        set(gca,'FontSize',12);
        PrettyFigureFormat;
        plotName(1,p) = {sprintf('plots/Spectra%g.pdf',p)};
        if ispc
            export_fig(plotName{1,p},'-pdf');
        elseif isunix || ismac
            spectra_figname = sprintf('plots/bcmultipix%03d.eps',p);
            publish_figure(p+2,spectra_figname);
        end
        %print('plots/AllSpectra','-dpdf');
        close(gcf);
    end
    if ispc
        delete('plots/SpectraAllPixels.pdf');
        append_pdfs('plots/SpectraAllPixels.pdf',...
            'plots/FPDNorm.pdf','plots/FPDBack.pdf', plotName{:})
    elseif isunix || ismac
        PATH = getenv('PATH');
        setenv('PATH', [PATH,':/usr/local/bin/']);
        cd plots
        command = 'gs -sDEVICE=pdfwrite -sOutputFile="SpectraAllPixelsunix.pdf" -dNOPAUSE -dEPSCrop -c "<</Orientation 0>> setpagedevice" -f bcmultipix*.eps -c quit';
        unix(command);
        unix('rm bcmultipix*.eps');
        cd ..
    end
end


