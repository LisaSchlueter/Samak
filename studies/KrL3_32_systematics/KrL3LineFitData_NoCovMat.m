function [par, err, fit_xi2min] = KrL3LineFitData_NoCovMat(varargin)
%
% FIT INTEGRAL SPECTRUM DATA
% 83mKr L3-32 Line

% Written for a quick Fit without Systematic Effects

% Th. Lasserre - CEA Saclay
% Last Updated: January 2018
%
orange=1/255*[255,140,0];

% Initialization
clear par ; clear parnomix;
% Parser
p = inputParser;

p.addParameter('plots', 'OFF', @(x)ismember(x,{'ON', 'OFF'}));
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('Fitter','Minuit',@(x)ismember(x,{'Minuit','Matlab'}));
p.addParameter('KrObject','',@(x) isa(x,'Kr'));
p.parse(varargin{:});

display           =    p.Results.display;
Fitter            =    p.Results.Fitter;
plots             =    p.Results.plots;
A                 =    p.Results.KrObject;
% Parametrization: True Value

LineE_i    = A.L3_32_E_i;
LineW_i    = A.L3_32_W_i;

% Loop on fits
npar       = 5; fit_p      = ones(npar,1); fit_xi2min = ones(1,1);

% Init
i_E        = 0; i_W        = 0; i_Bkg      = 0; i_Phi0     = 0;i_N        = 0;

        % Read Data File
        [~, krdatamap] = readKrData('TD',A.TD);
        Count   = (krdatamap(A.FPD_Pixel,1,:));
        CountErr = (krdatamap(A.FPD_Pixel,2,:));
        Data = [A.qU(:), Count(:), CountErr(:)];

% Background / Amplitude Initialization / Fill Data
  LinePhi0_i = (Count(1)-Count(end));  A.L3_32_Phi0_i  = LinePhi0_i;
  BKG_i      =  Count(end);            A.BKG_RateSec_i = BKG_i;  

%disp(Data);

switch display
    case 'ON'    
        A.DisplayKrInfo;
        fprintf(2,'--------------------------------------------------------------\n');
        fprintf(2,'  Initial Parameters \n');
        fprintf(2,'--------------------------------------------------------------\n');
        fprintf(2,'  Init E \t= %.3f eV  \n',   LineE_i);
        fprintf(2,'  Init W \t= %.3f meV \n',   LineW_i);
        fprintf(2,'  Init Phi0 \t= %.3f  \n',   LinePhi0_i);
        fprintf(2,'  Init Bkg \t= %.3f cps \n', BKG_i );
        fprintf(2,'--------------------------------------------------------------\n');      
end

switch Fitter
    case 'Matlab'
        krisf = @(eb,wb,pb,bb,nb,qu) A.ComputeKrISf(qu,eb,wb,pb,bb,nb);
        tmpf = @(eb,wb,pb,bb,nb,e) interp1(A.qU,krisf(eb,wb,pb,bb,nb,A.qU),e);
        myf = fittype(@(eb,wb,pb,bb,nb,qu) tmpf(eb,wb,pb,bb,nb,qu),...
            'independent','qu','coefficients',{'eb','wb','pb','bb'},'problem','nb');
        opts = fitoptions('TolFun',1e-6,'StartPoint',[0 0 0 0],...
            'Method', 'NonLinearLeastSquares',...
            'Display','off',...
            'Weights',(1./A.KrIS));
end

switch Fitter
    case 'Minuit'
        % Initializing Fit
        ParIni = [i_E i_W i_Phi0 i_Bkg i_N];
        parnames = ['E W Phi0 B N'];
       %tmparg = sprintf(['set pri -10 ; fix 2 5; set now; min ; imp']);
        tmparg = sprintf(['set pri 1 ; fix 5 ; minos ;']);
        DataKr = {Data, A};  %Data + KrObject 
        Args = {ParIni, DataKr, '-c',tmparg};
        fprintf(2,'No Systematics - Stat. + Pulls\n');
        %[par, err, chi2min, errmat] = fminuit('KrLineChi2', Args{:});
        [par, err, chi2min, errmat] = pminuit('KrLineChi2', Args);
        fit_xi2min = chi2min; 
        fit_p =par;  
       
       
    case 'Matlab'
        [fit1,gof,~] = fit(A.qU,A.KrIS,myf,opts,'problem',0);
        par = [fit1.eb fit1.wb fit1.pb fit1.bb fit1.nb];
        fit_p(:)=par; err = [0 0 0 0 0];
        chi2min = gof.sse; fit_xi2min(f) = chi2min;
end


% Definition of Amplitude / Background: For Display
E_Fit             = (LineE_i+par(1));
E_FitError        = (err(1));

W_Fit             = (LineW_i+par(2))*1e3;
W_FitError        = (err(2))*1e3;

Phi0_Fit          = (A.L3_32_Phi0_i+par(3));%*(A.TimeSec.*A.qUfrac(1));
Phi0_FitError     = err(3);%/(A.TimeSec.*A.qUfrac(1));

Offset_Fit        = (A.BKG_RateSec_i+par(4));%*(A.TimeSec.*A.qUfrac(1));
Offset_FitError   = err(4);%/(A.TimeSec.*A.qUfrac(1));

ndof = A.nqU-4;

switch display
    case 'ON'
        fprintf(2,'----------------------------------------------------\n');
        fprintf(2,'Fit Kr83m L3-32 Line Pixel %u \n', A.FPD_Pixel);
        fprintf(2,'----------------------------------------------------\n');
        fprintf(2,'  E \t= %.3f \t  +-\t %g \t eV \n',E_Fit,E_FitError);
        fprintf(2,'  W \t= %.3f \t  +- \t %g \t meV \n',W_Fit,W_FitError);
        fprintf(2,'  Phi0 \t= %.3f \t  +- \t %g \t \n',Phi0_Fit,Phi0_FitError);
        fprintf(2,'  Bkg \t= %.3f \t  +- \t %g \t cps \n',Offset_Fit,Offset_FitError);
        fprintf(2,'  N \t= %g \t  +- \t %g \t fixed \n',par(5),err(5));
        fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,ndof);
        fprintf(2,'----------------------------------------------------\n');
end

%% Plot Results 
switch plots
    case 'ON'
        
fig = figure('NumberTitle','off','rend','painters' ,'pos',[10 10 1300 1000]);
set(gcf,'Color','w') %background white
subplot(2,1,1)
 hdata = errorbar(Data(:,1)*1e-03,Data(:,2),Data(:,3),...
     's','MarkerSize',3,'Color', orange );%,'LineWidth',1);
 hold on
 hfit1 = plot(Data(:,1)*1e-03,KrLineModel4par(par,A),...
     'Color','Red','LineStyle','-', 'LineWidth',1.5);
 hfit3 = line([Data(1,1)*1e-03 Data(end,1)*1e-03],[Offset_Fit Offset_Fit],'LineStyle','--','LineWidth',1.5);
 hold off
 grid on
 xlabel('qU (keV)','FontSize',14);
 ylabel('$\dot{\textrm{\textbf{N}}} (\textrm{\textbf{cps}})$','FontSize',14, 'Interpreter','latex');
 set(gca,'FontSize',14);
 set(gca,'yscale','lin');
 %mydata = sprintf('Data: E=%.2f eV - W=%.2f eV \n',...
  %   LineE_i,LineW_i);
 myfit = sprintf('Fit: \\chi2 / dof=%.1f/%.0f\n E= %.3f ± %.3f eV \n W=%.0f  ± %.0f meV \n A=%.3f  ± %.3f \n O=%.3f  ± %.3f cps',...
     chi2min,A.nqU-4,E_Fit,E_FitError,W_Fit,W_FitError,...
     Phi0_Fit,Phi0_FitError,Offset_Fit,Offset_FitError);
 mychi2 = sprintf('Data');
 legend([hdata  hfit1 hfit3],mychi2,myfit,'Offset','Location','NorthEast') ;
 xlim([min(A.qU)*1e-03 max(A.qU)*1e-03]);% 0.*min(Data(:,2)) max(Data(:,2))*1.2])
 title(sprintf('KATRIN Gaseous Krypton 83m -{%s} Pixel %u Data and Fit (No Cov. Matrix)',...
     'KrL3-32', A.FPD_Pixel),'FontSize',14);
 
 subplot(2,1,2)
 plot(A.qU*1e-03, (Count(:)-KrLineModel4par(par,A))./CountErr(:), 'x', 'Color', orange );
 hold on
 refline = line([A.qU(1)*1e-03, A.qU(end)*1e-03], [0, 0], 'LineStyle', '--', 'LineWidth', 1.5);
 hold off
 grid on
 xlabel('qU (keV)','FontSize',14);
 ylabel('\textrm{\textbf{Residuals}}/$\mathbf{\sigma}$','FontSize',14, 'Interpreter','latex');
 set(gca,'FontSize',14);
 set(gca,'yscale','lin');
 xlim([min(A.qU)*1e-03 max(A.qU)*1e-03]);%*1e-03 min(Data(:,2)-KrLineModel4par(par))*2 max(Data(:,2)-KrLineModel4par(par))*2]);
end
end
