function [par err, chi2min, ndof, fig] = krl332minuit_fitcovm(varargin)
%
% FIT INTEGRAL SPECTRUM
% 83mKr L3-32 Line
%
% With HV Ripple Systematics
% Covariance Matrix Approach
%
% Th. Lasserre - CEA Saclay
% December 2017
%

% Initialization
clear par ; clear parnomix; close all

% Parser
p = inputParser;
p.addParameter('fign',1,@(x)isfloat(x) && x>0);
p.addParameter('pub','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('webpub','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('CovMat','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('CovMatCompute',1e3-1,@(x)isfloat(x) && x>0);
p.parse(varargin{:});
fign             =    p.Results.fign;
pub              =    p.Results.pub;
webpub           =    p.Results.webpub;
CovMat           =    p.Results.CovMat;
CovMatCompute    =    p.Results.CovMatCompute; % <1000 = 'OFF'

% Parametrization: True Value
global A ; A =krl332minuit_init(); %Object of Krypon class (Kr.m)

%Switch on HV Ripple 
A.HVRipples     = 'ON';
A.HVRipplesP2PV = 1;

% Set Time to 1 s
A.TimeYear = 1/365.25/86400 ; A.TimeSec = 1;
fprintf(2,'MoSFitterJuly2017: Time=%.2f sec\n',A.TimeSec);

% Set qUfrac to 1 per bin
A.qUfrac  = ones(numel(A.qUfrac),1);
fprintf(2,'MoSFitterJuly2017: Fraction of time in qU bins = %s \n',...
    num2str(A.qUfrac'));

% Compute Spectra
A.ComputeKrDS(); A.ComputeKrIS();

% Fit Parameters -Bias- Intitialization
npar       = 5;
i_E        = 0;
i_W        = 0;
i_Bkg      = 0;
i_Phi0     = 0;
i_N        = 0;

% Read Data
[qU krdatamap] = readKrData('TD','MoSFitterJuly2017');
Count    = (krdatamap(1,1,:));
CountErr = (krdatamap(1,2,:));
Data = [A.qU(:),Count(:),CountErr(:)];
fprintf(2,'MoSFitterJuly2017: Data - Min=%.2f Max=%.2f Ndata=%.0f\n',...
    min(Count(:)), max(Count(:)),numel(Count(:)));

% Initialize Covariance Matrix
clear FitCovErrMatrix; global FitCovErrMatrix;
StatErrMatrix  = (diag(CountErr(:)));
switch CovMat
    case 'OFF'
        FitCovErrMatrix = StatErrMatrix.^2;
    case 'ON'
        if CovMatCompute < 1000
            CM=importdata('./krl332minuit_fitcovm.mat');
            fprintf(2,'Reading Covariance Matrix File\n');
            MyCovMat = CM;
        else
            A.L3_32_E_i     = 30472.3  ;
            A.L3_32_W_i     = 0.996    ;
            A.L3_32_Phi0_i  = 726.34   ;
            A.BKG_RateSec_i = 276.41   ;
            A.qUfrac  = ones(numel(A.qUfrac),1);
            A.ComputeKrDS();
            
            HVRipplesValue = ...
                A.HVRipplesP2PV ...
                + A.HVRipplesP2PV/10*randn(CovMatCompute,1); % Scan HVRipple Values
            myKrIS      = zeros(CovMatCompute,A.nqU);        % for storing IS
            
            tic
            progressbar('Generate KrIS Covariance Matrix');
            for i=1:CovMatCompute
                progressbar(i/CovMatCompute);
                % change hv ripple and compute integral spectrum
                A.HVRipplesP2PV = HVRipplesValue(i);
                A.ComputeKrIS; myKrIS(i,:)  = A.KrIS;
            end
            
            MyCovMat = cov(myKrIS);
            mystr = sprintf('./krl332minuit_fitcovm.mat');
            save(mystr,'MyCovMat','-mat');
            
            figS = figure('Name','Covariance Matrix Building - Spectra','NumberTitle','off','rend','painters','pos',[300 300 400 400]);
            hdataC = errorbar(Data(:,1),Data(:,2),Data(:,3),...
                'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
            hold on
            for i=1:CovMatCompute
                hcovC  = plot(Data(:,1),myKrIS(i,:),'Color',rgb('orange'),'LineWidth',1);
            end
            hold off
            grid on
            xlabel('qU (eV)','FontSize',14);
            ylabel('Counts','FontSize',14);
            set(gca,'FontSize',14);
            set(gca,'yscale','lin');
            
        end
        
        FitCovErrMatrix = MyCovMat + StatErrMatrix.^2;
        figC = figure('Name','Covariance Matrix','NumberTitle','off','rend','painters','pos',[100 100 600 200]);
        subplot(1,3,1)
        imagesc(StatErrMatrix.^2);
        colormap(copper);
        colorbar;
        title('Stat. Fluct.');
        pbaspect([1 1 1])
        
        subplot(1,3,2)
        imagesc(MyCovMat);
        colormap(copper);
        colorbar;
        title('HV Ripple CM - No Stat. Fluct.');
        pbaspect([1 1 1])
        
        subplot(1,3,3)
        imagesc(FitCovErrMatrix);
        colormap(copper);
        colorbar;
        title('HV Ripple CM - With Stat. Fluct.');
        pbaspect([1 1 1])
end

% Parametrization: True Value
A =krl332minuit_init();
A.HVRipples     = 'ON'; A.HVRipplesP2PV = 1;
A.TimeYear = 1/365.25/86400 ; A.TimeSec = 1;
A.qUfrac  = ones(numel(A.qUfrac),1);
A.ComputeKrDS(); A.ComputeKrIS();

% Line Parameters Initialization
LineE_i    = A.L3_32_E_i;
LineW_i    = A.L3_32_W_i;
LinePhi0_i = A.L3_32_Phi0_i;
Pixel      = 0;

% Background / Amplitude Initialization
LinePhi0_i = (Count(1)-Count(end));  A.L3_32_Phi0_i  = LinePhi0_i;
BKG_i      =  Count(end);            A.BKG_RateSec_i = BKG_i;
fprintf(2,'MoSFitterJuly2017: L3_32_Phi0_i=%.2f \t BKG_RateSec_i=%.2f \n',...
    A.L3_32_Phi0_i,A.BKG_RateSec_i);
A.DisplayKrInfo();

fprintf(2,'---------------------------------------------------\n');
fprintf(2,'  Initial Parameters \n');
fprintf(2,'---------------------------------------------------\n');
fprintf(2,'  Init E \t= %.3f eV  \n',   LineE_i);
fprintf(2,'  Init W \t= %.3f meV \n',   LineW_i);
fprintf(2,'  Init Phi0 \t= %.3f  \n',   LinePhi0_i);
fprintf(2,'  Init Bkg \t= %.3f cps \n', BKG_i );
fprintf(2,'---------------------------------------------------\n');


% Fit calling Minuit
ParIni   = [i_E i_W i_Phi0 i_Bkg i_N];
parnames = ['E W Phi0 B N'];
tmparg = sprintf(['set pri -100 ; fix 5 ; set now; min ; imp']);
Args = {ParIni, Data, '-c',tmparg};
[par, err, chi2min, errmat] = fminuit('krl332minuit_chi2cov',Args{:});

% Definition of Amplitude / Background
E_Fit             = (LineE_i+par(1));
E_FitError        = (err(1));

W_Fit             = (LineW_i+par(2))*1e3;
W_FitError        = (par(2))*1e3;

Phi0_Fit          = (A.L3_32_Phi0_i+par(3))*(A.TimeSec.*A.qUfrac(1));
Phi0_FitError     = err(3)/(A.TimeSec.*A.qUfrac(1));

Offset_Fit        = (A.BKG_RateSec_i+par(4))*(A.TimeSec.*A.qUfrac(1));
Offset_FitError   = err(4)/(A.TimeSec.*A.qUfrac(1));

fprintf(2,'----------------------------------------------------\n');
fprintf(2,'Fit Kr83m L3-32 Line \n');
fprintf(2,'----------------------------------------------------\n');
fprintf(2,'  E \t= %.3f \t� \t %g \t eV \n',E_Fit,E_FitError);
fprintf(2,'  W \t= %.3f \t� \t %g \t meV \n',W_Fit,W_FitError);
fprintf(2,'  Phi0 \t= %.3f \t� \t %g \t \n',Phi0_Fit,Phi0_FitError);
fprintf(2,'  Bkg \t= %.3f \t� \t %g \t cps \n',Offset_Fit,Offset_FitError);
fprintf(2,'  N \t= %g \t� \t %g \t fixed \n',(par(5)),err(5));
%             fprintf(2,'  E \t= %.3f \t� \t %g \t eV \n',   (par(1)),err(1));
%             fprintf(2,'  W \t= %.3f \t� \t %g \t meV \n',  (par(2))*1e3,err(2)*1e3);
%             fprintf(2,'  Phi0 \t= %.3f \t� \t %g \t    \n',(par(3)),err(3));
%             fprintf(2,'  Bkg \t= %.3f \t� \t %g \t cps \n',(par(4)),err(4));
%             fprintf(2,'  N \t= %g \t� \t %g \t    \n',     (par(5)),err(5));
ndof = A.nqU-4;fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,ndof);
fprintf(2,'----------------------------------------------------\n');

%% Plot Fit and Residuals

switch webpub
    case 'OFF'
        
        fig = figure('Name','Fit and Residuals','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
        subplot(2,2,1)
        hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),...
            'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        hold on
        hfit1 = plot(Data(:,1),krl332minuit_modelint(par),...
            'Color','Black','LineWidth',1,'LineStyle','-');
        hfit3 = line([Data(1,1) Data(end,1)],[par(4) par(4) ],'LineStyle','--','Color','Red');
        hcov  = errorbar(Data(:,1),Data(:,2),sqrt(diag(FitCovErrMatrix)),...
            'ks','MarkerSize',3,'MarkerFaceColor',rgb('orange'),'Color',rgb('orange'),'LineWidth',1);
        hold off
        grid on
        xlabel('qU (eV)','FontSize',14);
        ylabel('Counts/qU','FontSize',14);
        set(gca,'FontSize',14);
        set(gca,'yscale','lin');
        mydata = sprintf('Data: E=%.2f eV - W=%.2f eV \n',...
            LineE_i,LineW_i);
        myfit = sprintf(' Fit: \\chi2 / dof=%.1f/%.0f \n E=%.3f \t\\pm %.3f eV \n W=%.3f \t\\pm %.3f meV \n A=%.3f \t\\pm %.3f \n O=%.3f \t\\pm %.3f',...
            chi2min,A.nqU-4,E_Fit,E_FitError,W_Fit,W_FitError,...
            Phi0_Fit,Phi0_FitError,Offset_Fit,Offset_FitError);
        mychi2 = sprintf('Data');
        legend([hcov hdata  hfit1 hfit3],'Stat.+Syst',mychi2,myfit,'Offset','Location','NorthEast') ;
        axis([min(A.qU) max(A.qU)+1 0.*min(Data(:,2)) max(Data(:,2))*1.2])
        title(sprintf('KATRIN Gaseous Krypton 83m - %s Data and Fit',...
            A.TD),'FontSize',14);
        
        subplot(2,2,3)
        hdata = errorbar(Data(:,1),Data(:,2)-krl332minuit_modelint(par),Data(:,3),...
            'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        hold on
        line([Data(1,1) Data(end,1)],[0 0 ],'LineStyle','--','Color','Blue');
        hold off
        grid on
        xlabel('qU (eV)','FontSize',14);
        ylabel('Residuals','FontSize',14);
        set(gca,'FontSize',14);
        set(gca,'yscale','lin');
        axis([min(A.qU) max(A.qU)+1 min(Data(:,2)-krl332minuit_modelint(par))*2 max(Data(:,2)-krl332minuit_modelint(par))*2]);
        
        subplot(2,2,2)
        hfit2 = plot(A.Te,(krl332minuit_modeldiff(par)),...
            'Color','Black','LineWidth',1,'LineStyle','-');
        grid on
        xlabel('Te (eV)','FontSize',14);
        ylabel('Counts/Te','FontSize',14);
        set(gca,'FontSize',14);
        set(gca,'yscale','lin');
        title(sprintf('Differential Spectrum - Doppler %s (fixed \\sigma=%.1f eV)',...
            A.DopplerEffectFlag,A.StandDev),'FontSize',14);
        
        subplot(2,2,4)
        %figure('Name','Fit Correlation Matrix','NumberTitle','off','rend','painters','pos',[10 10 500 500])
        %title(sprintf('KATRIN Krypton 83m Fit - Correlation Matrix - %s - Doppler %s\n',A.TD,A.DopplerEffectFlag));
        [h cor] = corrplot(errmat(1:4,1:4));
        xticklabels({'Position','Width','Amplitude','Offset'})
        yticklabels({'Position','Width','Amplitude','Offset'})
        disp(cor);
        
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
        y = (Data(:,2)-krl332minuit_modelint(par))';
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
        plot(A.Te,(krl332minuit_modeldiff(par)));
        grid on
        xlabel('Te (eV)','FontSize',14);
        ylabel('Counts/Te','FontSize',14);
        set(gca,'FontSize',14);
        set(gca,'yscale','lin');
        title(sprintf('Differential Spectrum - Doppler %s (fixed \\sigma=%.1f eV)',...
            A.DopplerEffectFlag,A.StandDev),'FontSize',14);
        response = fig2plotly(figwww1, 'filename', 'Differential Spectrum');
        plotly_url2 = response.url;        
        disp(plotly_url2);

end

end
