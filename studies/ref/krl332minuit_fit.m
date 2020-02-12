function [par , err , chi2min, ndof , fig] = krl332minuit_fit(varargin)
%
% FIT INTEGRAL SPECTRUM
% 83mKr L3-32 Line
%
% Th. Lasserre - CEA Saclay
% December 2017
%

% Initialization
clear par ; clear parnomix; close all
%addpath('../../MATLAB-Online-master/')
%addpath(genpath('../../../Samak2.0'));

% Parser
p = inputParser;
p.addParameter('fign',1,@(x)isfloat(x) && x>0);
p.addParameter('pub','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('webpub','OFF',@(x)ismember(x,{'ON','OFF'}));
p.parse(varargin{:});
fign             =    p.Results.fign;
pub              =    p.Results.pub;
webpub           =    p.Results.webpub;

% Parametrization: True Value
%global A ; 
A = krl332minuit_init();
A.ComputeKrDS(); A.ComputeKrIS();

% Line Parameters Initialization
LineE_i    = A.L3_32_E_i;
LineW_i    = A.L3_32_W_i;
%LinePhi0_i = A.L3_32_Phi0_i;
%Pixel      = 0;

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

% Set Time to 1 s
A.TimeYear = 1/365.25/86400 ; A.TimeSec = 1;
fprintf(2,'MoSFitterJuly2017: Time=%.2f sec\n',A.TimeSec);

% Set qUfrac to 1 per bin
A.qUfrac  = ones(numel(A.qUfrac),1);
fprintf(2,'MoSFitterJuly2017: Fraction of time in qU bins = %s \n',...
    num2str(A.qUfrac'));

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
%tmparg = sprintf(['set pri -10 ; fix 5 ; min ; imp ']);
%tmparg = sprintf(['fix 5 ; minos ; imp ; SCAn 1 ; SCAn 2 ; SCAn 3 ; SCAn 4 ; CONtour 1 2 ; CONtour 1 3 ; CONtour 1 4 ; CONtour 2 3 ; CONtour 2 4 , CONtour 3 4']);
tmparg = sprintf([' set pri -10 ; fix 5 ; min ; imp   ']);
DataL = {Data,A};
Args = {ParIni, DataL,'-n',parnames, '-c',tmparg};
[par, err, chi2min, errmat] = pminuit('krl332minuit_chi2',Args);
drawnow('update')

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
fprintf(2,'  E \t= %.3f \t± \t %g \t eV \n',E_Fit,E_FitError);
fprintf(2,'  W \t= %.3f \t± \t %g \t meV \n',W_Fit,W_FitError);
fprintf(2,'  Phi0 \t= %.3f \t± \t %g \t \n',Phi0_Fit,Phi0_FitError);
fprintf(2,'  Bkg \t= %.3f \t± \t %g \t cps \n',Offset_Fit,Offset_FitError);
fprintf(2,'  N \t= %g \t± \t %g \t fixed \n',(par(5)),err(5));
%             fprintf(2,'  E \t= %.3f \t± \t %g \t eV \n',   (par(1)),err(1));
%             fprintf(2,'  W \t= %.3f \t± \t %g \t meV \n',  (par(2))*1e3,err(2)*1e3);
%             fprintf(2,'  Phi0 \t= %.3f \t± \t %g \t    \n',(par(3)),err(3));
%             fprintf(2,'  Bkg \t= %.3f \t± \t %g \t cps \n',(par(4)),err(4));
%             fprintf(2,'  N \t= %g \t± \t %g \t    \n',     (par(5)),err(5));
ndof = A.nqU-4;fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,ndof);
fprintf(2,'----------------------------------------------------\n');
drawnow('update');

%% Plot Fit and Residuals
close all

switch webpub
    case 'OFF'
        
        fig = figure('Name','Fit and Residuals','NumberTitle','off','rend','painters','pos',[15 15 1260 760]);
        
        subplot(2,2,1)
        hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),...
            'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        hold on
        hfit1 = plot(Data(:,1),krl332minuit_modelint(par,A),...
            'Color','Black','LineWidth',1,'LineStyle','-');
        hfit3 = line([Data(1,1) Data(end,1)],[Offset_Fit Offset_Fit ],'LineStyle','--','Color','Red');
        hold off
        grid on
        xlabel('qU (eV)','FontSize',14);
        ylabel('Counts/qU','FontSize',14);
        set(gca,'FontSize',14);
        set(gca,'yscale','lin');
        mydata = sprintf('Data: E=%.2f eV - W=%.2f eV \n',...
            LineE_i,LineW_i);
        myfit = sprintf(' Fit: \\chi2 / dof=%.1f/%.0f \n E=%.3f \t\\pm %.3f eV \n W=%.3f \t\\pm %.3f eV \n A=%.3f \t\\pm %.3f \n O=%.3f \t\\pm %.3f',...
            chi2min,A.nqU-4,E_Fit,E_FitError,W_Fit,W_FitError,...
            Phi0_Fit,Phi0_FitError,Offset_Fit,Offset_FitError);
        mychi2 = sprintf('Data');
        legend([hdata  hfit1 hfit3],mychi2,myfit,'Offset','Location','NorthEast') ;
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
            A.DopplerEffectFlag,A.DE_sigma),'FontSize',14);
        
        subplot(2,2,4)
        %figure('Name','Fit Correlation Matrix','NumberTitle','off','rend','painters','pos',[10 10 500 500])
        %title(sprintf('KATRIN Krypton 83m Fit - Correlation Matrix - %s - Doppler %s\n',A.TD,A.DopplerEffectFlag));
        [h cor] = corplot(errmat(1:4,1:4));
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
            A.DopplerEffectFlag,A.DE_sigma),'FontSize',14);
        response = fig2plotly(figwww1, 'filename', 'Differential Spectrum');
        plotly_url2 = response.url;        
        disp(plotly_url2);

end

end
