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

%Parameters
nfit            = 1;
runtime        = (1)*60*60*24;
%runtime        = 210.2189781*18;

nTeBinFac       = 5;
TD              = 'FT-TL4';
Fitter          = 'Matlab';
FPD_Pixel       = 1:148;

FSD             = 'DOSS';
BCK             = 'XmasData';
Mode            = 'Read';

% Computing
savespectrum    = 'OFF';
FitFlag         = 'OFF';
UseCovMat       = 'OFF';
% Displaying
display         = 'OFF';

% Saving
saveresults     = 'OFF';
savepdfs        = 'OFF';
plotfpdviewer   = 'OFF';

% Parametrization: True Value
mnuSq_t = (0.00)^2;

% Initialization
opt_bin = {...
    'nTeBinningFactor',nTeBinFac};

opt_wgtsmace = {...
    'KTFFlag','Compute'}; %MACE+WGTSIS %SSCW_DCOfficial %SSCW_DCperqU

opt_bkg = {...
    'BKG_Flag',BCK};

opt_fsd= {'DTFSD',FSD};

opt_katrin = {...
    'TD',TD,...
    'TimeSec',runtime,...
    'Mode',Mode};

opt_fpd = {...
    'FPD_Pixel',FPD_Pixel};

A = ref_ftmultipix(...
    opt_katrin{:},...
    opt_wgtsmace{:},...
    opt_fsd{:},...
    opt_bkg{:},...
    opt_bin{:},...
    opt_fpd{:},...
    'mnuSq_i',mnuSq_t,...
    'UseParallelRF','ON',...
    'FPD_Segmentation','OFF');

%DATA = HDF5Reader('run',40257,'usetotal','ON');


fit_p = zeros(nfit,4*A.nPixels+6);

if strcmp(UseCovMat,'ON')
    if exist('results/CM.mat','file') == 2
        CM = load('results/CM.mat');
        MegaCovMatFrac = CM.MegaCovMatFrac;
    else
        A_CM = ref_ftmultipix(...
            opt_katrin{:},...
            opt_wgtsmace{:},...
            opt_fsd{:},...
            opt_bkg{:},...
            opt_bin{:},...
            'mnuSq_i',mnuSq_t,...
            'UseParallelRF','OFF',...
            'FPD_Segmentation','OFF');
        
        %2. Define your CM Class with your number of trials, choose SysEffect,...
        SysEffect = struct(...
            'RF_EL','OFF',...  % Response Function(RF) EnergyLoss
            'RF_BF','OFF',...  % R F B-Fields
            'RF_RX','OFF',...  % RF Column Density, Cross Section
            'FSD','ON',...
            'TASR','ON');     % Tritium Activity Stack Runs
        C = CovarianceMatrix('StudyObject',A_CM,'nTrials',1500,...
            'RecomputeFlag','ON',...
            'RunFlag','OFF',...
            'SysEffect',SysEffect,...
            'MACE_Ba_T_RelErr',1e-2,...
            'WGTS_B_T_RelErr',1e-2,...
            'WGTS_CD_MolPerCm2_RelErr',0.1,...
            'ISXsection_RelErr',1e-2);
        %C.WGTS_TASR_RelErr = 0.01;
        
        %%4.Compute Covariance Matrix for your SysEffect; here RF (Response Function)
        %C.ComputeCM_RF();
        C.ComputeCM_FSD;
        % C.PlotCM();        %Plot the CM
        % disp(C.CovMatFile)
        % C.ComputeFracCM('Mode','CM2Frac'); %Compute the fractional CM
        % C.ComputeCMRank(); %Compute the Rank of the CM as a first sanity check
        
        % Build mega Covariance Matrix
        CovMatpixCell = repmat({C.CovMatFrac},1,length(1));
        MegaCovMatFrac = blkdiag(CovMatpixCell{:});
    end
    
    A.ComputeTBDDS();
    A.ComputeTBDIS();
    
    MegaCovMat = A.TBDIS(:).*MegaCovMatFrac.*A.TBDIS(:)';
    
    % Define Covariance Matri for fit
    FitCovErrMatrix = diag(A.TBDIS(:)) + MegaCovMat;
    
end
%%
test = 1;
if test == 1
    A.ComputeTBDDS();
    A.ComputeTBDIS();
    A.AddStatFluctTBDIS();
    Data = [...
        A.qU, ...
        A.TBDIS, ...
        A.TBDISE];
    Afit = ref_ftmultipix(...
        opt_katrin{:},...
        opt_wgtsmace{:},...
        opt_fsd{:},...
        opt_fpd{:},...
        'BKG_Flag','ON',...
        'BKG_Type','FLAT',...
        'BKG_RateAllFPDSec',0,...
        opt_bin{:},...
        'mnuSq_i',mnuSq_t,...
        'UseParallelRF','ON',...
        'FPD_Segmentation','OFF',...
        'FPD_Pixel',[4,6,9]);
    Afit.ComputeTBDDS();
    Afit.ComputeTBDIS();  
    
    F = FITC('SO',Afit,'DATA',Data,'fitter','matlab','pulls',[4 inf inf inf]);
    
    
    %par = F.RESULTS{1}; err = F.RESULTS{2};
    %chi2min = F.RESULTS{3};
    %DoF = Afit.nqU*Afit.nPixels - length(1) - length(1) - 2;
    
end

if strcmp(FitFlag,'ON')
    for ii = 1:nfit
        
        tic
        A.ComputeTBDDS();
        A.ComputeTBDIS();
        A.AddStatFluctTBDIS();
        Data = [...
            A.qU(:) , ...
            A.TBDIS(:), ...
            A.TBDISE(:)];
        % Init
        i_mnu      = 0;
        i_Q        = 0;
        i_B        = A.TBDIS(end,:)/(A.qUfrac(end)*A.TimeSec);
        i_N        = zeros(1,A.nPixels);
        
        DoF = A.nqU*A.nPixels - length(i_B) - length(i_N) - 2;
        
        Afit = ref_ftmultipix(...
            opt_katrin{:},...
            opt_wgtsmace{:},...
            opt_fsd{:},...
            opt_fpd{:},...
            'BKG_Flag','ON',...
            'BKG_Type','FLAT',...
            'BKG_RateAllFPDSec',0,...
            opt_bin{:},...
            'mnuSq_i',mnuSq_t,...
            'UseParallelRF','ON',...
            'FPD_Segmentation','OFF');%'MULTIPIXEL');
        Afit.ComputeTBDDS();
        Afit.ComputeTBDIS();
        
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
                %                 fprintf(2,'No Detector S7ystematics - Stat.\n');
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
                        DataTBD = {Data,Afit};
                        TBDfun = @(xx) Chi2E0_first_allPixels(xx,DataTBD);
                    case 'ON'
                        DataTBD = {Data,Afit,FitCovErrMatrix};
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
                
                forpdf{1} = sprintf('===============================================\n');
                forpdf{2} = sprintf('  m^2       = %g +/- %g eV^2\n', mnuSq_report,err(1));
                forpdf{3} = sprintf('  dQ        = %g +/- %g eV\n', par(2),err(2));
                forpdf{4} = sprintf('  B total   = %g +/- %g cps\n', BCK_fit,BCK_err);
                forpdf{5} = sprintf('  dN average= %g +/- %g\n', N_ave,N_err);
                forpdf{6} = sprintf('  Chi2/dof  = %g/%g\n', chi2min, DoF);
                forpdf{7} = sprintf('===============================================\n');
        end
        toc
        fit_p(ii,:) = [par err chi2min DoF];
        disp(ii);
    end
end

if strcmp(saveresults,'ON'); save(['results/fit_ftmultipix_','CM','.mat'],'fit_p'); end

if false
%% Plotting in PDF
fprintf('Setting figure visibility to OFF. \n')
% Plots from the focal plane detector
fprintf('Plot NORM per pixel in FPD \n')
pixel_norms = NaN(148,1);
pixel_norms(Afit.FPD_Pixel) = norms_fit;
figure(1);
 set(gcf,'Visible','off');      
 set(gcf,'color','white');
FPDViewer(pixel_norms,'ReDrawSkeleton','ON');
title(['Normalization bias ',TD,' ',num2str(runtime), 's']);
if ispc; export_fig('plots/FPDNorm.pdf','-pdf');
elseif isunix || ismac; publish_figure(2,'plots/ftmultipix001.eps'); end
fprintf('Plot NORM done \n')

fprintf('Plot BCK per pixel in FPD \n')
pixel_bcks = NaN(148,1);
pixel_bcks(Afit.FPD_Pixel) = bcks_fit;
figure(2);
 set(gcf,'Visible','off');
 set(gcf,'color','white');
FPDViewer(pixel_bcks,'ReDrawSkeleton','ON');
title(['Background ',TD,' ',num2str(runtime), 's']);

if ispc; export_fig('plots/FPDBCK.pdf','-pdf');
elseif isunix || ismac; publish_figure(2,'plots/ftmultipix002.eps'); end
fprintf('Plot BCK done \n')

% Fit results summarize
figure(3)
 set(gcf,'Visible','off');              
text(0.5,0.5,forpdf,'FontName','FixedWidth','HorizontalAlignment','center')
if ispc; export_fig('plots/FitResults.pdf','-pdf');
elseif isunix || ismac; publish_figure(2,'plots/ftmultipix003.eps'); end
fprintf('Printing fit results done \n')

% Distributions of the parameters per pixel
figure(4);
 set(gcf,'Visible','off'); 
 set(gcf,'color','white');
bcks_residuals = bcks_fit - A.BKG_RateSec;
[mu_b,sigma_b] = normfit(bcks_residuals);
hold on
[~,Nbins,Wbins] = nhist(bcks_residuals);
x_b = linspace(min(bcks_residuals),max(bcks_residuals),length(Nbins)*Afit.nTeBinningFactor);
binWidth = Wbins(2) - Wbins(1);
plot(x_b,Afit.nPixels*binWidth*normpdf(x_b,mu_b,sigma_b));
hold off
annotation('textbox',[0.6 0.5 0.2 0.2],'String',...
    {['\sigma = ',num2str(sigma_b),' cps'],['mean = ',num2str(mu_b),' cps']},'FitBoxToText','on');
title('\delta background [cps]');
xlabel('\delta background [cps]');
if ispc; export_fig('plots/BCKdist.pdf','-pdf');
elseif isunix || ismac; publish_figure(2,'plots/ftmultipix004.eps'); end
fprintf('Plot BCK dist done \n')

figure(5);
 set(gcf,'Visible','off');    
 set(gcf,'color','white');
[mu_N,sigma_N] = normfit(norms_fit);
hold on
[~,Nbins,Wbins] = nhist(norms_fit);
x_N = linspace(min(norms_fit),max(norms_fit),length(Nbins)*Afit.nTeBinningFactor);
binWidth = Wbins(2) - Wbins(1);
plot(x_N,Afit.nPixels*binWidth*normpdf(x_N,mu_N,sigma_N));
hold off
annotation('textbox',[0.6 0.5 0.2 0.2],'String',...
    {['\sigma = ',num2str(sigma_N),''],['mean = ',num2str(mu_N),'']},'FitBoxToText','on');
title('normalization');
xlabel('normalization');
if ispc; export_fig('plots/NORMdist.pdf','-pdf');
elseif isunix || ismac; publish_figure(2,'plots/ftmultipix005.eps'); end
fprintf('Plot norm dist done \n')

% Plots of correlation from background and normalization
figure(6)
 set(gcf,'Visible','off');    
 set(gcf,'color','white');
hc = corrplotm([bcks_fit',norms_fit'],'varNames',{'BCK','NORM'});
    ax = axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
    set(get(ax,'Title'),'Visible','on');
    ax.Title.String = {' ',['Correlation Multipixel Fit']};
if ispc; export_fig('plots/CORRBCKNorm.pdf','-pdf');
elseif isunix || ismac; publish_figure(2,'plots/ftmultipix006.eps'); end
fprintf('Plot correlation done \n')

% Plots from each of the spectra of each pixel

M_IS = Afit.TBDIS;
M_ISE = Afit.TBDISE;
D_IS = A.TBDIS;
D_ISE = A.TBDISE;

plotNameSpectraPixel = cell(1,length(FPD_Pixel));
for p = 1:length(FPD_Pixel)
    fprintf('Plotting data, spectrum and residuals pixel %0.d \n',p)
    figure(p+6);
     set(gcf,'Visible','off'); 
     set(gcf,'color','white');
    % Plot of the data points, spectra and error bars
    subplot(2,1,1)
    hold on
    hline = plot(Afit.qU(:,p) - Afit.Q, M_IS(:,p),...
        'LineWidth',1,'LineStyle','-','Color','Black');
    hdata = errorbar(A.qU(:,p) - A.Q_i, D_IS(:,p), D_ISE(:,p),...
        'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
    hold off
    mydata = ['Data: E_0eff = 18575 eV'];
    myline = sprintf('Fit: \\delta E_0eff = %.3f \\pm %.3f eV',par(2),err(2));
    lgd_fit = legend([hdata hline],'Data','Model',...
        'Location','Best');
    %lgd_fit = legend([hline hchi2 hdata],myline,mychi2,mydata,'Location','Best') ;
    lgd_fit.FontSize = 8;
    xlabel('qU-E_0 (eV)','FontSize',14);
    ylabel('Counts','FontSize',14);
    title(sprintf('KATRIN - %s - pix %g - %0.0f s ',TD,p,runtime));
    set(gca,'FontSize',12);
    PrettyFigureFormat;
    
    % Plot of the residuals according to uncertainty of each data point
    subplot(2,1,2)
    hold on
    hresiduals = scatter(A.qU(:,p)-Afit.Q,...
        ((D_IS(:,p)-M_IS(:,p))./D_ISE(:,p)),...
        'ks','LineWidth',1);
    hresidualserr = errorbar(A.qU(:,p)-Afit.Q,...
        (D_IS(:,p)-M_IS(:,p))./D_ISE(:,p),...
        D_ISE(:,p)./D_ISE(:,p),'ks','LineWidth',1);
    plot(linspace(A.qUmin-Afit.Q,A.qUmax-Afit.Q,A.nqU),zeros(A.nqU,1),'k')
    hold off
    grid on
    xlabel('qU-E_0 (eV)','FontSize',14);
    ylabel('Residuals','FontSize',14);
    set(gca,'FontSize',12);
    lgd_res = legend([hresidualserr],'(Data-Fit)/uncertainty',...
        'Location','Best');
    lgd_res.FontSize = 8;
    PrettyFigureFormat;
        
    % Saving the generated plot
    if ispc
        plotNameSpectraPixel(1,p) = {sprintf('plots/spectrumpix%g.pdf',p)};
        export_fig(plotNameSpectraPixel{1,p},'-pdf');
    elseif isunix || ismac
        spectra_figname = sprintf('plots/ftmultipix%03d.eps',p+6);
        publish_figure(p+6,spectra_figname);
    end
end

% Fit results per pixel

latex_preamble = {'\documentclass[10pt,a4paper]{article}';...
'\usepackage[utf8]{inputenc}';...
'\usepackage{amsmath}';...
'\usepackage{amsfonts}';...
'\usepackage{amssymb}';...
'\begin{document}';...
'\begin{small}'};...

fid = fopen('latex_preamble.tex','w');
fprintf(fid,'%s\n', latex_preamble{:});
fclose(fid);

results_size = length(bcks_fit);
index1 = ceil(results_size/3);
index2 = 2*index1;
index3 = results_size;
space = cell(index1,1);
space(:) = {' '};
extraentries = cell(3*index1-index3,3);
horizontallabels = {'Pix.','Bck.','Norm.',' ','Pix.','Bck.','Norm.',' ','Pix.','Bck.','Norm.';...
    ' ','[mcps]',' ',' ',' ','[mcps]',' ',' ',' ','[mcps]',' '};
numformatlatex = {'%0.d','%0.4g','%0.3g','%0.d','%0.d','%0.4g','%0.3g','%0.d','%0.d','%0.4g','%0.3g','%0.d'};
bcks_fit = bcks_fit*1e3;

celltableforlatex = [num2cell([(1:index1)',bcks_fit(1:index1)',norms_fit(1:index1)']),space,...
    num2cell([(index1+1:index2)',bcks_fit(index1+1:index2)',norms_fit(index1+1:index2)']),space,...
    [num2cell([(index2+1:index3)',bcks_fit(index2+1:index3)',norms_fit(index2+1:index3)']);extraentries]];
latextable(celltableforlatex,'name','fitresultstable.tex','Horiz',horizontallabels,'Hline',[0,2,NaN],...
    'Vline',[1,2,3,5,6,7,9,10],'format',numformatlatex);
bcks_fit = bcks_fit/1e3;

latex_end = {'\end{small}';...
    '\end{document}'};

fid = fopen('latex_end.tex','w');
fprintf(fid,'%s\n', latex_end{:});
fclose(fid);

system('copy latex_preamble.tex+fitresultstable.tex+latex_end.tex fitresultslatex.tex')

system('pdflatex fitresultslatex.tex')

% Save all plots in one PDF and erase the individual PDFs

if ispc
    delete('plots/SpectraAllPixels.pdf');
    append_pdfs('plots/SpectraAllPixels.pdf',...
        'plots/FitResults.pdf','fitresultslatex.pdf','plots/FPDNorm.pdf','plots/FPDBCK.pdf',...
        'plots/BCKdist.pdf','plots/NORMdist.pdf','plots/CORRBCKNorm.pdf',...
        plotNameSpectraPixel{:})
   delete('plots/spectrumpix*.pdf','plots/FPDNorm.pdf','plots/FPDBCK.pdf','plots/FitResults.pdf',...
        'plots/BCKdist.pdf','plots/NORMdist.pdf','plots/CORRBCKNorm.pdf','plots/FitResultsPix.pdf');
elseif isunix || ismac
    PATH = getenv('PATH');
    setenv('PATH', [PATH,':/usr/local/bin/']);
    cd plots
    command = 'gs -sDEVICE=pdfwrite -sOutputFile="SpectraAllPixelsunix.pdf" -dNOPAUSE -dEPSCrop -c "<</Orientation 0>> setpagedevice" -f ftmultipix*.eps -c quit';
    unix(command);
    unix('rm ftmultipix*.eps');
    cd ..
end

fprintf('Saving PDFs done. \n')
close all
end