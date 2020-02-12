function [par err, chi2min, ndof, errmat] = KrL3LineFitData_CovMat(varargin)
%
% FIT INTEGRAL SPECTRUM 83mKr L3-32 Line
% 
% With HV Ripple Systematics - Covariance Matrix Approach
% 
%January 2018, Lisa Schlueter

addpath(genpath('../../../Samak2.0'));
% Initialization
clear par ; clear parnomix; close all

% Parser
p = inputParser;
p.addParameter('fign',1,@(x)isfloat(x) && x>0);
p.addParameter('pub','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('plots','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('webpub','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('CovMat','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('CovMatCompute',1e3,@(x)isfloat(x) && x>0); % Number of Spectra for Covariance Matrix 
p.addParameter('KrObject','',@(x) isa(x,'Kr')); %Object of Krypon class (Kr.m)

p.parse(varargin{:});
fign             =    p.Results.fign;
pub              =    p.Results.pub;
plots            =    p.Results.plots;
webpub           =    p.Results.webpub;
CovMat           =    p.Results.CovMat;
CovMatCompute    =    p.Results.CovMatCompute; % <1000 = 'OFF'
A                =    p.Results.KrObject;

% Compute Spectra
A.ComputeKrDS; A.ComputeKrIS;

        % Read Data
        [qU krdatamap] = readKrData('TD', A.TD);
        Count    = (krdatamap(A.FPD_Pixel,1,:));
        CountErr = (krdatamap(A.FPD_Pixel,2,:));
        Data = [A.qU(:),Count(:),CountErr(:)];

% Initialize Covariance Matrix
StatErrMatrix  = (diag(CountErr(:))); %stat error only
switch CovMat
    case 'OFF'
        FitCovErrMatrix = StatErrMatrix.^2;
    case 'ON'
        if CovMatCompute < 1000
           % CovMatFile = sprintf('CovMat/CovMatrix_KrL3LineFit_Pixel%u.mat', A.FPD_Pixel);
            CovMatFile = sprintf('CovMat_uniform/CovMatrix_KrL3LineFit_Pixel1_uniform2V.mat');
            CM=importdata(CovMatFile);
            fprintf(2,'Reading Covariance Matrix File\n');
            MyCovMat = CM;
        else
            LineE_i_tmp    = A.L3_32_E_i; % Before Building Covariance Matrix
            LineW_i_tmp    = A.L3_32_W_i;

            MyCovMat = BuildCovMat_HVRipple(A, CovMatCompute);
            %mystr = sprintf('../hv2018paper/CovMat/CovMatrix_hv2018paper_Pixel%u.mat', A.FPD_Pixel);
            %save(mystr,'MyCovMat','-mat');      
           
            A.L3_32_E_i = LineE_i_tmp; %set back to initial value from InitKrKATRIN
            A.L3_32_W_i  = LineW_i_tmp;
        end 
       FitCovErrMatrix = MyCovMat + StatErrMatrix.^2;    
end

%% Fit

% Parametrization: Initial Values
A.HVRipplesP2PV = 0.44; %back to initial/central value (CV)
A.ComputeKrIS;

% Line Parameters Initialization
LineE_i    = A.L3_32_E_i;
LineW_i    = A.L3_32_W_i;
% Background / Amplitude Initialization
LinePhi0_i = (Count(1)-Count(end));  A.L3_32_Phi0_i  = LinePhi0_i;
BKG_i      =  Count(end);            A.BKG_RateSec_i = BKG_i;
% Fit Parameters -Bias- Intitialization
npar       = 5;
i_E        = 0;
i_W        = 0;
i_Bkg      = 0;
i_Phi0     = 0;
i_N        = 0;

A.DisplayKrInfo();
fprintf(2,'---------------------------------------------------\n');
fprintf(2,'  Initial Parameters for CovMat Fit of Pixel %u \n', A.FPD_Pixel);
fprintf(2,'---------------------------------------------------\n');
fprintf(2,'  Init E \t= %.3f eV  \n',   A.L3_32_E_i);
fprintf(2,'  Init W \t= %.3f meV \n',   A.L3_32_W_i);
fprintf(2,'  Init Phi0 \t= %.3f  \n',   A.L3_32_Phi0_i);
fprintf(2,'  Init Bkg \t= %.3f cps \n', A.BKG_RateSec_i );
fprintf(2,'---------------------------------------------------\n');


% Fit calling Minuit
ParIni   = [i_E i_W i_Phi0 i_Bkg i_N];
parnames = ['E W Phi0 B N'];
tmparg = sprintf(['set pri 1 ; fix 2 5 ; minos ;']);
DataKr= {Data, A, FitCovErrMatrix};
Args = {ParIni, DataKr, '-c',tmparg};
[par, err, chi2min, errmat] = pminuit('krl332minuit_chi2cov',Args);

% Fit Results (for Display only)
E_Fit             = (LineE_i+par(1));
E_FitError        = (err(1));

W_Fit             = (LineW_i+par(2))*1e3;
W_FitError        = (err(2))*1e3;

Phi0_Fit          = (A.L3_32_Phi0_i+par(3));%*(A.TimeSec.*A.qUfrac(1));
Phi0_FitError     = err(3);%/(A.TimeSec.*A.qUfrac(1));

Offset_Fit        = (A.BKG_RateSec_i+par(4));%*(A.TimeSec.*A.qUfrac(1));
Offset_FitError   = err(4);%/(A.TimeSec.*A.qUfrac(1));

fprintf(2,'----------------------------------------------------\n');
fprintf(2,'Fit Kr83m L3-32 Line \n');
fprintf(2,'----------------------------------------------------\n');
fprintf(2,' Pixel %u \n', A.FPD_Pixel);
fprintf(2,'  E \t= %.3f \t \t %g \t eV \n',E_Fit,E_FitError);
fprintf(2,'  W \t= %.3f \t ± \t %g \t meV \n',W_Fit,W_FitError);
fprintf(2,'  Phi0 \t= %.3f \t ± \t %g \t \n',Phi0_Fit,Phi0_FitError);
fprintf(2,'  Bkg \t= %.3f \t ± \t %g \t cps \n',Offset_Fit,Offset_FitError);
fprintf(2,'  N \t= %g \t ± \t %g \t fixed \n',(par(5)),err(5));
ndof = (A.nqU-4);
fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,ndof);
fprintf(2,'----------------------------------------------------\n');

%% Plot Fit and Residuals
switch plots
    case 'ON'
switch webpub
    case 'OFF'
        
        fig = figure('Name','Fit and Residuals','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
        subplot(2,2,1)
        hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),...
            'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        hold on
        hfit1 = plot(Data(:,1),krl332minuit_modelint(par,A),...
            'Color','Black','LineWidth',1,'LineStyle','-');
        hfit3 = line([Data(1,1) Data(end,1)],[par(4) par(4) ],'LineStyle','--','Color','Red');
        hcov  = errorbar(Data(:,1),Data(:,2),sqrt(diag(FitCovErrMatrix)),...
            'ks','MarkerSize',3,'MarkerFaceColor','RED','Color','RED','LineWidth',1);
        hold off
        grid on
        xlabel('qU (eV)','FontSize',14);
        ylabel('Counts/qU','FontSize',14);
        set(gca,'FontSize',14);
        set(gca,'yscale','lin');
        mydata = sprintf('Data: E=%.2f eV - W=%.2f eV \n',...
            LineE_i,LineW_i);
        myfit = sprintf(' Fit: \\chi2 / dof=%.1f/%.0f \n E=%.3f \t\\pm %.3f eV \n W=%.3f \t\\pm %.3f meV \n A=%.3f \t\\pm %.3f \n O=%.3f \t\\pm %.3f cps',...
            chi2min,A.nqU-4,E_Fit,E_FitError,W_Fit,W_FitError,...
            Phi0_Fit,Phi0_FitError,Offset_Fit,Offset_FitError);
        mychi2 = sprintf('Data');
        legend([hcov hdata  hfit1 hfit3],'Stat.+Syst',mychi2,myfit,'Offset','Location','NorthEast') ;
        axis([min(A.qU) max(A.qU)+1 0.*min(Data(:,2)) max(Data(:,2))*1.2])
        title(sprintf('KATRIN Gaseous Krypton 83m - %s Data and Fit',...
            A.TD),'FontSize',14);
        
        subplot(2,2,3)
        hdata = errorbar(Data(:,1),Data(:,2)-krl332minuit_modelint(par,A),Data(:,3),...
            'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        hold on
        line([Data(1,1) Data(end,1)],[0 0 ],'LineStyle','--','Color','Blue');
        hold off
        grid on
        xlabel('qU (eV)','FontSize',14);
        ylabel('Residuals','FontSize',14);
        set(gca,'FontSize',14);
        set(gca,'yscale','lin');
        axis([min(A.qU) max(A.qU)+1 min(Data(:,2)-krl332minuit_modelint(par,A))*2 max(Data(:,2)-krl332minuit_modelint(par,A))*2]);
        
        subplot(2,2,2)
        hfit2 = plot(A.Te,(krl332minuit_modeldiff(par,A)),...
            'Color','Black','LineWidth',1,'LineStyle','-');
        grid on
        xlabel('Te (eV)','FontSize',14);
        ylabel('Counts/Te','FontSize',14);
        set(gca,'FontSize',14);
        set(gca,'yscale','lin');
        title(sprintf('Differential Spectrum - Doppler %s (fixed \\sigma=%.1f eV)',...
            A.DopplerEffectFlag),'FontSize',14);
        
%         subplot(2,2,4)
%         %figure('Name','Fit Correlation Matrix','NumberTitle','off','rend','painters','pos',[10 10 500 500])
%         %title(sprintf('KATRIN Krypton 83m Fit - Correlation Matrix - %s - Doppler %s\n',A.TD,A.DopplerEffectFlag));
%         [h cor] = corrplot(errmat(1:4,1:4));
%         xticklabels({'Position','Width','Amplitude','Offset'})
%         yticklabels({'Position','Width','Amplitude','Offset'})
%         disp(cor);
        
        switch pub
            case 'ON'
                myname = sprintf('./figures/krl332minuit_ref.eps');
                fprintf(2,'publish: %s\n',myname);
                publish_figure(999+fign,myname);
        end
        
    case 'ON'
        
        %% Plot Data
        clear x y e data1
        x = Data(:,1)';
        y = Data(:,2)';
        e = Data(:,3)';
        data1 = {...
            struct(...
            'x', x, ...
            'y', y, ...
            'mode', 'markers', ...
            'name', 'measured', ...
            'error_y', struct(...
            'type', 'data', ...
            'array', e, ...
            'visible', true,...
            'color', '#85144B', ...
            'thickness', 1.5, ...
            'width', 3, ...
            'opacity', 1), ...
            'marker', struct(...
            'color', '#85144B', ...
            'size', 8),...
            'type', 'scatter')...
            };
        response = plotly(data1, struct('filename', 'Integral Spectrum and Fit', 'fileopt', 'overwrite','offline', true,'open',true));
        plot_url3 = response.url;
        disp(plot_url3);

        %% Plot Residuals
        clear x y e data2
        x = Data(:,1)';
        y = (Data(:,2)-krl332minuit_modelint(par,A))';
        e = Data(:,3)';
        data2 = {...
            struct(...
            'x', x, ...
            'y', y, ...
            'mode', 'markers', ...
            'name', 'measured', ...
            'error_y', struct(...
            'type', 'data', ...
            'array', e, ...
            'visible', true,...
            'color', '#85144B', ...
            'thickness', 1.5, ...
            'width', 3, ...
            'opacity', 1), ...
            'marker', struct(...
            'color', '#85144B', ...
            'size', 8),...
            'type', 'scatter')...
            };
        response = plotly(data2, struct('filename', 'Residuals', 'fileopt', 'overwrite','offline', true,'open',true));
        plot_url1 = response.url;
        disp(plot_url1);
        
        %% Plot Diff. Spectrum
        figwww1 = figure();
        plot(A.Te,(krl332minuit_modeldiff(par,A)));
        grid on
        xlabel('Te (eV)','FontSize',14);
        ylabel('Counts/Te','FontSize',14);
        set(gca,'FontSize',14);
        set(gca,'yscale','lin');
        title(sprintf('Differential Spectrum - Doppler %s (fixed \\sigma=%.1f eV)',...
            A.DopplerEffectFlag),'FontSize',14);
        response = fig2plotly(figwww1, 'filename', 'Differential Spectrum');
        plotly_url2 = response.url;        
        disp(plotly_url2);

end
end
end
